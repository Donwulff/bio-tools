#!/usr/bin/env python3
"""
Repair EAGER/AdapterRemoval-prefixed BAM read names and pairing.

This tool expects queryname-sorted BAM input where read names are grouped
by base name (prefix-stripped) when pairing. If your BAM still has
M_/F_/R_ prefixes, run this tool in --strip-prefix mode, then
queryname-sort, then run again in pairing mode.

Recommended workflow:
  ./util/eager_repair_bam.py --strip-prefix input.bam -o stripped.bam
  samtools sort -n -o stripped.qname.bam stripped.bam
  ./util/eager_repair_bam.py --pair stripped.qname.bam -o repaired.bam
  samtools fixmate -m repaired.bam repaired.fixmate.bam
  samtools sort -o repaired.coord.bam repaired.fixmate.bam

Requires: pysam
"""

import argparse
import sys
import re

__version__ = "0.3.0"

try:
    import pysam
except ImportError as exc:
    sys.stderr.write("error: pysam is required (pip install pysam)\n")
    raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Repair EAGER/AdapterRemoval-prefixed BAM read names and pairing."
    )
    parser.add_argument(
        "input_bam",
        help="BAM input path or '-' for stdin. For pairing, input must be name-sorted.",
    )
    parser.add_argument(
        "-o",
        "--output-bam",
        help="Output BAM path or '-' for stdout (default: <input>.repaired.bam).",
    )
    parser.add_argument(
        "--output-uncompressed",
        action="store_true",
        help="Write uncompressed BAM (fast for pipes).",
    )
    parser.add_argument(
        "--prefix-regex",
        default=r"^(?P<pfx>[MFR])_",
        help="Regex to capture prefix letter as group 'pfx'.",
    )
    parser.add_argument(
        "--prefix-tag",
        default="XQ",
        help="SAM tag to store origin (default: XQ). Uses M/F/R/U/D.",
    )
    parser.add_argument(
        "--duplicate-merged",
        action="store_true",
        help="Duplicate M_ reads into pseudo-pairs. Duplicate gets XQ:D.",
    )
    parser.add_argument(
        "--dup-orientation",
        choices=("same", "flip"),
        default="same",
        help="Orientation for duplicate mate: same as original or flipped (default: same).",
    )
    parser.add_argument(
        "--add-legacy-suffix",
        action="store_true",
        help="Append /1 or /2 to query names for paired reads.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--strip-prefix",
        action="store_true",
        dest="strip_prefix",
        help="Strip M_/F_/R_ prefix and add prefix tag; no pairing changes.",
    )
    mode.add_argument(
        "--pair",
        action="store_true",
        dest="pair_only",
        help="Pair reads using existing prefix tags; do not strip prefix.",
    )
    mode.add_argument(
        "--strip-prefix-only",
        action="store_true",
        dest="strip_prefix",
        help=argparse.SUPPRESS,
    )
    mode.add_argument(
        "--pair-only",
        action="store_true",
        dest="pair_only",
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def clone_rec(rec, header):
    return pysam.AlignedSegment.fromstring(rec.to_string(), header)


def group_records(bam, prefix_re, prefix_tag, use_qname_prefix, strip_prefix):
    current_base = None
    bucket = []
    for rec in bam:
        pfx = None
        if prefix_tag and rec.has_tag(prefix_tag):
            tag = rec.get_tag(prefix_tag)
            if tag in ("M", "F", "R"):
                pfx = tag

        base = rec.query_name
        if use_qname_prefix:
            m = prefix_re.match(rec.query_name)
            if m:
                pfx = pfx or m.group("pfx")
                if strip_prefix:
                    base = rec.query_name[m.end():]
                else:
                    base = rec.query_name

        if current_base is None:
            current_base = base
        if base != current_base:
            yield current_base, bucket
            bucket = []
            current_base = base

        bucket.append((pfx, rec))

    if bucket:
        yield current_base, bucket


def add_suffix(name: str, suffix: str) -> str:
    if name.endswith("/1") or name.endswith("/2"):
        return name
    return f"{name}{suffix}"


def main() -> int:
    args = parse_args()
    prefix_re = re.compile(args.prefix_regex)
    if args.strip_prefix:
        use_qname_prefix = True
        strip_prefix = True
    elif args.pair_only:
        use_qname_prefix = False
        strip_prefix = False
    else:
        use_qname_prefix = True
        strip_prefix = True
    do_pairing = not args.strip_prefix

    output_bam = args.output_bam
    if not output_bam:
        if args.input_bam.endswith(".bam"):
            output_bam = args.input_bam[:-4] + ".repaired.bam"
        else:
            output_bam = args.input_bam + ".repaired.bam"

    in_mode = "rb"
    out_mode = "wb0" if args.output_uncompressed else "wb"

    with pysam.AlignmentFile(args.input_bam, in_mode) as bam_in:
        header_dict = bam_in.header.to_dict()
        pg_id = "eager_repair_bam"
        existing_pg = {pg.get("ID") for pg in header_dict.get("PG", []) if "ID" in pg}
        if pg_id in existing_pg:
            suffix = 1
            while f"{pg_id}.{suffix}" in existing_pg:
                suffix += 1
            pg_id = f"{pg_id}.{suffix}"
        header_dict.setdefault("PG", []).append(
            {
                "ID": pg_id,
                "PN": "eager_repair_bam",
                "VN": __version__,
                "CL": " ".join(sys.argv),
            }
        )
        out_header = pysam.AlignmentHeader.from_dict(header_dict)
        with pysam.AlignmentFile(output_bam, out_mode, header=out_header) as bam_out:
            for base, records in group_records(
                bam_in, prefix_re, args.prefix_tag, use_qname_prefix, strip_prefix
            ):
                if not do_pairing:
                    # Strip/retag only. Optionally duplicate merged reads.
                    for pfx, rec in records:
                        if not args.pair_only:
                            rec.query_name = base
                            if pfx in ("M", "F", "R"):
                                rec.set_tag(args.prefix_tag, pfx, value_type="Z")
                            else:
                                rec.set_tag(args.prefix_tag, "U", value_type="Z")

                        if pfx == "M" and args.duplicate_merged:
                            rec1 = rec
                            rec2 = clone_rec(rec, bam_in.header)
                            if args.dup_orientation == "flip":
                                rec2.is_reverse = not rec1.is_reverse
                            for idx, r in enumerate((rec1, rec2), start=0):
                                r.is_paired = True
                                r.is_read1 = (idx == 0)
                                r.is_read2 = (idx == 1)
                                r.is_proper_pair = False
                                r.set_tag(
                                    args.prefix_tag,
                                    "M" if idx == 0 else "D",
                                    value_type="Z",
                                )
                                if args.add_legacy_suffix:
                                    r.query_name = add_suffix(
                                        r.query_name, "/1" if idx == 0 else "/2"
                                    )
                            # Set mate fields after orientation is finalized.
                            rec1.mate_is_unmapped = rec2.is_unmapped
                            rec2.mate_is_unmapped = rec1.is_unmapped
                            rec1.mate_is_reverse = rec2.is_reverse
                            rec2.mate_is_reverse = rec1.is_reverse
                            rec1.next_reference_id = rec2.reference_id
                            rec1.next_reference_start = rec2.reference_start
                            rec2.next_reference_id = rec1.reference_id
                            rec2.next_reference_start = rec1.reference_start
                            if args.dup_orientation == "flip":
                                tlen = rec1.reference_length or rec1.query_length or 0
                                rec1.template_length = tlen
                                rec2.template_length = -tlen
                            else:
                                rec1.template_length = 0
                                rec2.template_length = 0
                            bam_out.write(rec1)
                            bam_out.write(rec2)
                        else:
                            bam_out.write(rec)
                    continue

                by_prefix = {"M": [], "F": [], "R": [], None: []}
                for pfx, rec in records:
                    by_prefix[pfx].append(rec)

                has_f = len(by_prefix["F"]) > 0
                has_r = len(by_prefix["R"]) > 0
                mate_present = has_f and has_r

                # Process F/R reads (paired)
                for rec in by_prefix["F"]:
                    if not args.pair_only:
                        rec.query_name = base
                        rec.set_tag(args.prefix_tag, "F", value_type="Z")
                    rec.is_paired = True
                    rec.is_read1 = True
                    rec.is_read2 = False
                    rec.mate_is_unmapped = not mate_present
                    if args.add_legacy_suffix:
                        rec.query_name = add_suffix(rec.query_name, "/1")
                    bam_out.write(rec)

                for rec in by_prefix["R"]:
                    if not args.pair_only:
                        rec.query_name = base
                        rec.set_tag(args.prefix_tag, "R", value_type="Z")
                    rec.is_paired = True
                    rec.is_read1 = False
                    rec.is_read2 = True
                    rec.mate_is_unmapped = not mate_present
                    if args.add_legacy_suffix:
                        rec.query_name = add_suffix(rec.query_name, "/2")
                    bam_out.write(rec)

                # Process merged reads
                for rec in by_prefix["M"]:
                    if not args.pair_only:
                        rec.query_name = base
                        rec.set_tag(args.prefix_tag, "M", value_type="Z")
                    if not args.duplicate_merged:
                        rec.is_paired = False
                        rec.is_read1 = False
                        rec.is_read2 = False
                        rec.mate_is_unmapped = False
                        rec.mate_is_reverse = False
                        rec.next_reference_id = -1
                        rec.next_reference_start = -1
                        rec.template_length = 0
                        bam_out.write(rec)
                        continue

                    # Duplicate into pseudo-pairs (tagged)
                    rec1 = rec
                    rec2 = clone_rec(rec, bam_in.header)
                    if args.dup_orientation == "flip":
                        rec2.is_reverse = not rec1.is_reverse
                    for idx, r in enumerate((rec1, rec2), start=0):
                        r.is_paired = True
                        r.is_read1 = (idx == 0)
                        r.is_read2 = (idx == 1)
                        r.is_proper_pair = False
                        r.set_tag(
                            args.prefix_tag,
                            "M" if idx == 0 else "D",
                            value_type="Z",
                        )
                        if args.add_legacy_suffix:
                            r.query_name = add_suffix(r.query_name, "/1" if idx == 0 else "/2")
                    # Set mate fields after orientation is finalized.
                    rec1.mate_is_unmapped = rec2.is_unmapped
                    rec2.mate_is_unmapped = rec1.is_unmapped
                    rec1.mate_is_reverse = rec2.is_reverse
                    rec2.mate_is_reverse = rec1.is_reverse
                    rec1.next_reference_id = rec2.reference_id
                    rec1.next_reference_start = rec2.reference_start
                    rec2.next_reference_id = rec1.reference_id
                    rec2.next_reference_start = rec1.reference_start
                    if args.dup_orientation == "flip":
                        tlen = rec1.reference_length or rec1.query_length or 0
                        rec1.template_length = tlen
                        rec2.template_length = -tlen
                    else:
                        rec1.template_length = 0
                        rec2.template_length = 0
                    bam_out.write(rec1)
                    bam_out.write(rec2)

                # Unprefixed reads: pass through with optional tag
                for rec in by_prefix[None]:
                    if not args.pair_only:
                        rec.set_tag(args.prefix_tag, "U", value_type="Z")
                    bam_out.write(rec)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
