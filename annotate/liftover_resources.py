"""Shared liftover reference discovery and trusted local shell configuration."""

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_NAMES = {
    "hs1": ("chm13v2.0_maskedY_rCRS.fa", "chm13v2.0_maskedY_rCRS.fa.gz"),
    "hg38": ("GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna",
             "GCA_000001405.15_GRCh38_no_alt_analysis_set_masked.fna.gz",
             "hg38.fa", "hg38.fa.gz"),
}
EXTENDED_NAMES = {
    "hs1": ("chm13v2.0_maskedY_rCRS.H3630.fa", "chm13v2.0_maskedY_rCRS.H3630.fa.gz",
            "chm13v2.0_maskedY_rCRSDH3630O1102.fa"),
    "hg38": ("hg38p14.H3630.fa", "hg38p14.H3630.fa.gz", "hg38p14DH3630O1102.fa"),
}
EXTENDED_PATTERN = re.compile(r"H3630|O1102|hs38d1|oral[_-]?microbiome|metagenom|IMGT", re.I)


class ResourceError(ValueError):
    pass


def config_environment(environ=None, cwd=None):
    """Run the same POSIX-shell loader as the shell drivers, without parsing shell in Python.

    Configuration is executable, trusted local code. Its diagnostics go to stderr;
    NUL-delimited environment output preserves spaces, newlines and '=' in values.
    Do not print this environment: it can contain private paths or credentials.
    """
    result = subprocess.run(
        ["sh", "-eu", "-c", '. "$1/lib_liftover.sh"; load_liftover_config >&2; env -0',
         "liftover-config", str(SCRIPT_DIR)], env=environ, cwd=cwd, capture_output=True,
    )
    if result.stderr:
        print(os.fsdecode(result.stderr), end="", file=sys.stderr)
    if result.returncode:
        raise ResourceError("could not load trusted local liftover configuration")
    return dict(os.fsdecode(value).split("=", 1)
                for value in result.stdout.split(b"\0") if value)


def reference_roots(root=None, *, environ=None, cwd=None, repo_root=None):
    env = os.environ if environ is None else environ
    cwd = Path.cwd() if cwd is None else Path(cwd)
    repo = SCRIPT_DIR.parent if repo_root is None else Path(repo_root)
    configured = root or env.get("REF_ROOT") or env.get("BIO_TOOLS_INDEX_DIR") or env.get("INDEX_DIR")
    candidates = ([Path(configured).expanduser()] if configured else [])
    if env.get("Y_RESOURCES_DIR"):
        candidates.append(Path(env["Y_RESOURCES_DIR"]).expanduser() / "ref")
    candidates += [cwd / "resources/ref", repo / "mapping/index"]
    roots = []
    for path in candidates:
        if not path.is_absolute():
            path = cwd / path
        path = path.absolute()
        if path not in roots:
            roots.append(path)
    return roots


def checked_reference(path, allow_extended=False):
    # Check both names: symlinking a known extended reference does not opt it in.
    if not allow_extended and any(EXTENDED_PATTERN.search(p.name) for p in (path, path.resolve())):
        raise ResourceError(f"reference appears extended/decoy-enriched: {path}; "
                            "use a non-extended reference or explicitly pass --allow-extended-ref")
    if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
        raise ResourceError(f"not a readable, nonempty reference: {path}")
    return path


def find_reference(assembly, explicit=None, root=None, *, allow_extended=False,
                   environ=None, cwd=None, repo_root=None):
    """Choose a path only; format, sidecars and sequence identity are separate checks."""
    if assembly not in REFERENCE_NAMES:
        raise ResourceError(f"unsupported assembly: {assembly}")
    env = os.environ if environ is None else environ
    cwd = Path.cwd() if cwd is None else Path(cwd)
    explicit = explicit or env.get(f"REF_{assembly.upper()}")
    if explicit:
        path = Path(explicit).expanduser()
        return checked_reference((path if path.is_absolute() else cwd / path).absolute(), allow_extended)
    roots = reference_roots(root, environ=env, cwd=cwd, repo_root=repo_root)
    # Prefer a base analysis reference anywhere in the search roots over an
    # extended fallback, even when the latter has been explicitly permitted.
    groups = [REFERENCE_NAMES[assembly]]
    if allow_extended:
        groups.append(EXTENDED_NAMES[assembly])
    for names in groups:
        for directory in roots:
            for name in names:
                path = directory / name
                if path.exists() or path.is_symlink():
                    return checked_reference(path, allow_extended)
    raise ResourceError(f"--ref-{assembly} not set and auto-detect failed; searched: "
                        + ", ".join(map(str, roots))
                        + ". Supply a reference explicitly; extended names require --allow-extended-ref.")


def main():
    parser = argparse.ArgumentParser(description="Select a liftover reference without modifying it.")
    parser.add_argument("--assembly", required=True, choices=tuple(REFERENCE_NAMES))
    parser.add_argument("--reference")
    parser.add_argument("--ref-root")
    parser.add_argument("--allow-extended-ref", action="store_true")
    args = parser.parse_args()
    try:
        print(find_reference(args.assembly, args.reference, args.ref_root,
                             allow_extended=args.allow_extended_ref))
    except (ResourceError, OSError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
