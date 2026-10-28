"""Read-only reference checks shared by marker queries and liftover preparation."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
import gzip
import hashlib
from pathlib import Path
import re
import struct


class ReferenceError(ValueError):
    pass


@dataclass(frozen=True)
class FaiEntry:
    length: int
    offset: int
    bases: int
    width: int


def nonempty(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ReferenceError(f"required reference file is missing or empty: {path}")


def dictionary_paths(path: Path) -> tuple[Path, ...]:
    name = str(path)
    for suffix in (".fasta.gz", ".fna.gz", ".fas.gz", ".fa.gz", ".fasta", ".fna", ".fas", ".fa"):
        if name.endswith(suffix):
            canonical = Path(name[:-len(suffix)] + ".dict")
            break
    else:
        canonical = path.with_suffix(".dict")
    return tuple(dict.fromkeys((canonical, path.with_suffix(".dict"))))


def fasta_format(path: Path) -> str:
    nonempty(path)
    with path.open("rb") as fh:
        header = fh.read(12)
        if header[:2] != b"\x1f\x8b":
            if not header.startswith(b">"):
                raise ReferenceError(f"not a FASTA file: {path}")
            return "plain"
        if len(header) == 12 and header[2] == 8 and header[3] & 4:
            extra = fh.read(struct.unpack("<H", header[10:12])[0])
            offset = 0
            while offset + 4 <= len(extra):
                tag, length = extra[offset:offset+2], struct.unpack("<H", extra[offset+2:offset+4])[0]
                stop = offset + 4 + length
                if stop > len(extra):
                    break
                if tag == b"BC" and length == 2:
                    return "bgzf"
                offset = stop
    raise ReferenceError(f"gzip FASTA is not BGZF: {path}; stage uncompressed or bgzip and rebuild indexes")


def read_fai(path: Path) -> dict[str, FaiEntry]:
    nonempty(path)
    entries = {}
    previous_offset = -1
    with path.open(encoding="utf-8") as fh:
        for line, raw in enumerate(fh, 1):
            fields = raw.rstrip("\r\n").split("\t")
            try:
                if len(fields) != 5 or not fields[0] or fields[0] in entries:
                    raise ValueError
                length, offset, bases, width = map(int, fields[1:])
                if length <= 0 or offset <= previous_offset or bases <= 0 or width < bases or bases > length:
                    raise ValueError
            except ValueError:
                raise ReferenceError(f"invalid or duplicate FASTA index entry at {path}:{line}") from None
            entries[fields[0]] = FaiEntry(length, offset, bases, width)
            previous_offset = offset
    return entries


def check_gzi(path: Path) -> None:
    nonempty(path)
    with path.open("rb") as fh:
        first = fh.read(8)
        if len(first) != 8:
            raise ReferenceError(f"invalid BGZF .gzi index: {path}")
        count, = struct.unpack("<Q", first)
        if path.stat().st_size != 8 + 16 * count:
            raise ReferenceError(f"invalid BGZF .gzi index size: {path}")
        previous = (-1, -1)
        for _ in range(count):
            current = struct.unpack("<QQ", fh.read(16))
            if current[0] <= previous[0] or current[1] <= previous[1]:
                raise ReferenceError(f"non-monotonic BGZF .gzi index: {path}")
            previous = current


def indexed_metadata(path: Path) -> dict[str, FaiEntry]:
    fmt = fasta_format(path)
    entries = read_fai(Path(f"{path}.fai"))
    if fmt == "bgzf":
        check_gzi(Path(f"{path}.gzi"))
    else:
        with path.open("rb") as fh:
            for name, entry in entries.items():
                fh.seek(max(0, entry.offset - 1))
                if entry.offset == 0 or fh.read(1) not in {b"\n", b"\r"}:
                    raise ReferenceError(f"FASTA index offset does not start a sequence line for {name!r}: {path}")
                if entry.width > entry.bases:
                    fh.seek(entry.offset + entry.bases)
                    ending = fh.read(entry.width - entry.bases)
                    # faidx may retain a nominal EOL width on a final one-line
                    # sequence even though the file ends immediately after it.
                    at_eof = (not ending and entry.length == entry.bases and
                              entry.offset + entry.length == path.stat().st_size and
                              entry.width - entry.bases <= 2)
                    if ending not in {b"\n", b"\r\n"} and not at_eof:
                        raise ReferenceError(f"FASTA index line layout is stale for {name!r}: {path}")
    return entries


def read_dictionary(path: Path) -> dict[str, tuple[int, str]]:
    nonempty(path)
    result = {}
    with path.open(encoding="utf-8") as fh:
        for line, raw in enumerate(fh, 1):
            if not raw.startswith("@SQ\t"):
                continue
            fields = {}
            try:
                for field in raw.rstrip("\r\n").split("\t")[1:]:
                    key, value = field.split(":", 1)
                    if key in fields:
                        raise ValueError
                    fields[key] = value
                name, length = fields["SN"], int(fields["LN"])
                checksum = fields.get("M5", "").lower()
                if not name or length <= 0 or name in result or (checksum and not re.fullmatch(r"[0-9a-f]{32}", checksum)):
                    raise ValueError
            except (ValueError, KeyError):
                raise ReferenceError(f"invalid or duplicate @SQ at {path}:{line}") from None
            result[name] = (length, checksum)
    if not result:
        raise ReferenceError(f"dictionary has no @SQ entries: {path}")
    return result


def check_dictionary(path: Path, entries: dict[str, FaiEntry]) -> dict[str, tuple[int, str]]:
    dictionary = read_dictionary(path)
    if list(dictionary) != list(entries) or any(dictionary[n][0] != entries[n].length for n in entries):
        raise ReferenceError(f"dictionary sequence names, lengths or order differ from FASTA index: {path}")
    return dictionary


@contextmanager
def indexed_reference(path: Path):
    entries = indexed_metadata(path)
    try:
        import pysam
    except ImportError:
        raise ReferenceError("indexed reference validation requires pysam (python3 -m pip install pysam)") from None
    try:
        fasta = pysam.FastaFile(str(path))
    except (OSError, ValueError) as exc:
        raise ReferenceError(f"cannot open indexed reference {path}: {exc}") from None
    with fasta:
        # Small boundary/line-wrap probes catch broken layouts without scanning a genome.
        for name, entry in entries.items():
            for pos in {0, entry.length - 1, min(entry.bases, entry.length - 1)}:
                try:
                    base = fasta.fetch(name, pos, pos + 1).upper()
                except (KeyError, OSError, ValueError) as exc:
                    raise ReferenceError(f"unusable FASTA index for {name!r}: {exc}") from None
                if len(base) != 1 or base not in "ACGTNRYKMSWBDHV":
                    raise ReferenceError(f"FASTA index/layout mismatch for {name!r} in {path}")
        yield lru_cache(maxsize=4096)(fasta.fetch)


def check_sequence_md5(path: Path, entries: dict[str, FaiEntry], dictionary: dict[str, tuple[int, str]]) -> None:
    """Optional full sequential scan; metadata/probes alone do not establish sequence identity."""
    if any(not digest for _, digest in dictionary.values()):
        raise ReferenceError("full sequence validation requires M5 for every dictionary sequence")
    actual, actual_index = {}, {}
    unterminated_single_lines = set()
    name, digest, length = None, None, 0
    offset, line_bases, width, position, previous_line = 0, 0, 0, 0, None
    opener = gzip.open if fasta_format(path) == "bgzf" else open
    with opener(path, "rb") as fh:
        for raw in fh:
            start = position
            position += len(raw)
            if raw.startswith(b">"):
                if name is not None:
                    actual[name] = (length, digest.hexdigest())
                    actual_index[name] = FaiEntry(length, offset, line_bases, width)
                words = raw[1:].split()
                if not words:
                    raise ReferenceError(f"empty FASTA header: {path}")
                name = words[0].decode("utf-8")
                if name in actual:
                    raise ReferenceError(f"duplicate FASTA sequence name: {name!r}")
                digest, length = hashlib.md5(), 0
                previous_line = None
                offset, line_bases, width = 0, 0, 0
            else:
                if digest is None:
                    raise ReferenceError(f"invalid FASTA content: {path}")
                bases = b"".join(raw.split()).upper()
                if not bases:
                    raise ReferenceError(f"blank sequence line prevents unambiguous FASTA indexing: {path}")
                if previous_line is not None and previous_line != (line_bases, width):
                    raise ReferenceError(f"irregular FASTA line wrapping: {path}")
                if length == 0:
                    offset, line_bases, width = start, len(bases), len(raw)
                elif len(bases) > line_bases:
                    raise ReferenceError(f"irregular FASTA line wrapping: {path}")
                previous_line = (len(bases), len(raw))
                digest.update(bases)
                length += len(bases)
    if name is not None:
        actual[name] = (length, digest.hexdigest())
        actual_index[name] = FaiEntry(length, offset, line_bases, width)
        if length == line_bases and width == line_bases:
            unterminated_single_lines.add(name)
    layout_ok = list(actual_index) == list(entries)
    if layout_ok:
        for name, declared in entries.items():
            observed = actual_index[name]
            if observed == declared:
                continue
            if (name not in unterminated_single_lines or observed.length != declared.length or
                    observed.offset != declared.offset or observed.bases != declared.bases or
                    not observed.width <= declared.width <= observed.width + 2):
                layout_ok = False
                break
    if list(actual) != list(entries) or actual != dictionary or not layout_ok:
        raise ReferenceError(f"FASTA content differs from index/dictionary M5 values: {path}")
