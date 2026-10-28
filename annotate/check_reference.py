#!/usr/bin/env python3
"""Validate reference metadata and indexed access without creating or replacing files."""
from __future__ import annotations

import argparse
from pathlib import Path

from reference_utils import (ReferenceError, check_dictionary, check_sequence_md5,
                             dictionary_paths, fasta_format, indexed_metadata, indexed_reference)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", type=Path)
    modes = ap.add_mutually_exclusive_group()
    modes.add_argument("--format-only", action="store_true")
    modes.add_argument("--index-only", action="store_true")
    modes.add_argument("--print-dict-path", action="store_true")
    ap.add_argument("--check-md5", action="store_true", help="also scan all sequence bytes and compare dictionary M5 values")
    args = ap.parse_args()
    if args.check_md5 and (args.format_only or args.index_only or args.print_dict_path):
        ap.error("--check-md5 requires full dictionary validation")
    try:
        if args.print_dict_path:
            print(dictionary_paths(args.reference)[0])
            return 0
        if args.format_only:
            print(fasta_format(args.reference))
            return 0
        entries = indexed_metadata(args.reference)
        dictionaries = []
        if not args.index_only:
            dictionaries = [p for p in dictionary_paths(args.reference) if p.exists() or p.is_symlink()]
            if not dictionaries:
                raise ReferenceError(f"missing .dict file for {args.reference}")
            variants = [check_dictionary(dictionary, entries) for dictionary in dictionaries]
            for current in variants[1:]:
                for name in entries:
                    first, other = variants[0][name][1], current[name][1]
                    if first and other and first != other:
                        raise ReferenceError(f"dictionary aliases have conflicting M5 for {name!r}")
            if args.check_md5:
                if any(not value[1] for variant in variants for value in variant.values()):
                    raise ReferenceError("full sequence validation requires M5 in all dictionary aliases")
                check_sequence_md5(args.reference, entries, variants[0])
        with indexed_reference(args.reference):
            pass
        print(f"Validated {len(entries)} reference sequences" + (" and dictionary" if dictionaries else ""))
    except (ReferenceError, OSError, UnicodeError) as exc:
        raise SystemExit(f"ERROR: {exc}") from None
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
