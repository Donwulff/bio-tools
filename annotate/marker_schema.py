"""Shared marker-status table contract and ancestral/derived allele resolution."""
from __future__ import annotations


OUTPUT_COLUMNS = ("chrom", "pos", "id", "ref", "alt", "aa", "hg", "isogg",
                  "gt", "gq", "dp", "source", "status", "reason")
STATUSES = frozenset({"derived", "ancestral", "nocall", "ambiguous", "filtered", "other", "unknown"})
CALLABLE = frozenset({"derived", "ancestral"})


def derived_allele(ref: str, alt: str, ancestral: str, derived: str = ".") -> str:
    """Use the explicit allele when available, otherwise the resolver's first-ALT rule."""
    if derived and derived != ".":
        return derived
    first_alt = alt.split(",")[0]
    if ancestral == ref:
        return first_alt
    if ancestral == first_alt:
        return ref
    return "."
