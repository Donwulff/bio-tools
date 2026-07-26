#!/usr/bin/env python3
"""Self-contained checks for ytree.py and y_path_rank.py's TreeScorer.

Run: python3 annotate/test_ytree.py    (no pytest, no fixtures on disk)

The placement cases below encode calls this repository has actually made from
reads -- the Iceman terminal at G-L166 while derived at the refuted G-Z6208, the
Swiss Neolithic men stopping at G-Z6219, the Sardinians stopping at G-PF3239.
If the scorer ever stops reproducing those, it is the scorer that changed, and
this file says so before any result table does.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ytree

REPO = Path(__file__).resolve().parent.parent
TREE = REPO / "markers" / "tree_local.tsv"

HEADER = "\t".join(ytree.COLUMNS)

failures: list[str] = []


def check(cond: bool, what: str) -> None:
    if cond:
        print(f"  ok   {what}")
    else:
        print(f"  FAIL {what}")
        failures.append(what)


def check_raises(text: str, what: str, fragment: str = "") -> None:
    """The tree file is load-time validated; a malformed one must never load."""
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        ytree.load_tree(path)
    except ytree.TreeError as exc:
        if fragment and fragment not in str(exc):
            check(False, f"{what} (raised, but not about {fragment!r}: {exc})")
        else:
            check(True, what)
    else:
        check(False, f"{what} (loaded without error)")
    finally:
        Path(path).unlink()


def row(marker: str, status: str, ref: str = "C", alt: str = "A") -> list[str]:
    """A marker_status.tsv row: 13 columns, ID at 2, REF/ALT at 3/4, STATUS at 12."""
    return ["chrY", "1", marker, ref, alt, ".", ".", ".", ".", ".", ".", ".", status]


def score(calls: dict[str, str], tree: ytree.Tree) -> ytree.TreeScorer:
    s = ytree.TreeScorer(tree)
    for marker, status in calls.items():
        r = row(marker, status)
        s.add(r[2], ytree.status_from_call(r[12]))
    return s


def test_real_tree_loads() -> ytree.Tree:
    print("markers/tree_local.tsv")
    tree = ytree.load_tree(TREE)
    check(tree.root == "G-P287", f"root is G-P287 (got {tree.root})")
    check(tree.nodes["G-Z6219"].parent == "G-PF3239", "G-Z6219 hangs below G-PF3239")
    check(tree.nodes["G-L166"].parent == "G-Z6219", "G-L166 hangs below G-Z6219")
    check(tree.nodes["G-Z6499"].parent == "G-L166", "G-Z6499 hangs below G-L166")
    check(tree.nodes["G-Z6208"].status == "refuted", "G-Z6208 is marked refuted")
    check(tree.lookup("L166") == ("G-L166", ytree.ROLE_DEFINING),
          "L166 is a defining marker of G-L166")
    check(tree.lookup("Z6134") == ("G-L166", ytree.ROLE_EQUIVALENT),
          "Z6134 is an equivalent of G-L166, not decisive")
    check(tree.lookup("Z6135") == ("G-Z6219", ytree.ROLE_EQUIVALENT),
          "Z6135 is an unregistered equivalent of G-Z6219")
    check(tree.lookup("M201") is None, "an unlisted marker resolves to nothing")
    check("G-Z6208" not in tree.path_nodes("G-Z6208"),
          "a refuted node never appears in a path")
    return tree


def test_validation() -> None:
    print("validation refuses malformed trees")
    base = f"{HEADER}\nG-A\t-\tpublished\tA1\t-\troot\n"
    check_raises(base + "G-B\tG-A\tpublished\tA1\t-\tdup\n",
                 "a marker on two nodes is rejected", "already assigned")
    check_raises(base + "G-B\tG-NOPE\tpublished\tB1\t-\tx\n",
                 "an unknown parent is rejected", "unknown parent")
    check_raises(base + "G-B\t-\tpublished\tB1\t-\tx\n",
                 "two roots are rejected", "exactly one root")
    check_raises(f"{HEADER}\nG-A\tG-B\tpublished\tA1\t-\tx\nG-B\tG-A\tpublished\tB1\t-\tx\n",
                 "a cycle is rejected")
    check_raises(base + "G-B\tG-A\tmaybe\tB1\t-\tx\n",
                 "an invented status is rejected", "not one of")
    check_raises(base + "G-B\tG-A\tpublished\t-\tB1\tx\n",
                 "a node with no defining marker is rejected", "no defining")
    check_raises(base + "G-A\tG-A\tpublished\tA9\t-\tx\n",
                 "a duplicate node name is rejected", "listed twice")
    check_raises("G-A\t-\tpublished\tA1\t-\tx\n",
                 "a missing header is rejected", "header")


def test_placements(tree: ytree.Tree) -> None:
    print("placement from read-level calls")

    # The Iceman: derived through L166 including the refuted Z6208, ancestral at
    # both branches below L166.
    iceman = score({
        "P287": "derived", "P15": "derived", "PF3147": "derived", "L91": "derived",
        "PF3239": "derived", "Z6219": "derived", "FGC5671": "derived",
        "L166": "derived", "L167": "derived",
        "Z6208": "derived",
        "Z6499": "ancestral", "Z6494": "ancestral", "FGC5687": "ancestral",
    }, tree)
    check(iceman.terminal() == ["G-L166"],
          f"Iceman terminates at G-L166 (got {iceman.terminal()})")
    check(iceman.call("G-Z6208") == "DERIVED",
          "derived at Z6208 is still reported as derived")
    check("G-Z6208" not in iceman.terminal(),
          "a refuted node cannot be a terminal placement even when derived")
    check(iceman.excluded() == ["G-Z6494", "G-Z6499"],
          f"both branches below L166 are excluded (got {iceman.excluded()})")

    # Oberbipp / Aesch / Muttenz: derived at Z6219, ancestral at the two L166
    # transversions. This is the split that put Z6219 above L166 in the first place.
    swiss = score({
        "PF3239": "derived", "Z6219": "derived",
        "L166": "ancestral", "L167": "ancestral",
    }, tree)
    check(swiss.terminal() == ["G-Z6219"],
          f"Swiss Neolithic men terminate at G-Z6219 (got {swiss.terminal()})")
    check(swiss.excluded() == ["G-L166"], "G-L166 is excluded for them")

    # The three Sardinians: their own evidence stops at PF3239.
    sardinia = score({
        "PF3239": "derived", "FGC5671": "ancestral", "L166": "ancestral",
    }, tree)
    check(sardinia.terminal() == ["G-PF3239"],
          f"Sardinians terminate at G-PF3239 (got {sardinia.terminal()})")
    check(sardinia.excluded() == ["G-Z6219", "G-L166"],
          f"both G-Z6219 and G-L166 are excluded (got {sardinia.excluded()})")

    # Coverage gaps above the call must not sink it -- most of the backbone is
    # uncovered in most of these libraries.
    gappy = score({"Z6219": "derived", "L166": "derived"}, tree)
    check(gappy.terminal() == ["G-L166"],
          "an uncovered backbone does not block a deeper placement")

    # An ancestral read above a derived one does.
    contradicted = score({"PF3239": "ancestral", "L166": "derived"}, tree)
    check(contradicted.terminal() == [],
          f"an ancestral ancestor blocks placement (got {contradicted.terminal()})")

    # Both states at one node's defining markers.
    conflict = score({
        "PF3239": "derived", "Z6219": "derived",
        "L166": "derived", "L167": "ancestral",
    }, tree)
    check(conflict.call("G-L166") == "CONFLICT", "L166 derived + L167 ancestral is CONFLICT")
    check(conflict.terminal() == ["G-Z6219"],
          f"a conflicted node does not terminate (got {conflict.terminal()})")

    # Equivalents are counted but never decisive.
    equiv_only = score({"Z6219": "derived", "Z6134": "derived"}, tree)
    check(equiv_only.eq_derived["G-L166"] == 1, "an equivalent marker is counted")
    check(equiv_only.terminal() == ["G-Z6219"],
          f"an equivalent alone cannot deepen a call (got {equiv_only.terminal()})")

    # Two sibling branches supported at once is a result, not a tie to break.
    both = score({"Z6219": "derived", "L166": "derived",
                  "Z6494": "derived", "Z6499": "derived"}, tree)
    check(sorted(both.terminal()) == ["G-Z6494", "G-Z6499"],
          f"two supported siblings are both reported (got {both.terminal()})")

    empty = score({}, tree)
    check(empty.terminal() == [], "no coverage yields no placement")


def test_newick(tree: ytree.Tree) -> None:
    print("newick export")
    nwk = ytree.to_newick(tree)
    check(nwk.endswith(";"), "newick is terminated")
    check(nwk.count("(") == nwk.count(")"), "newick parentheses balance")
    check("G-Z6208" not in nwk, "a refuted node is not exported by default")
    check("G-Z6208" in ytree.to_newick(tree, include_refuted=True),
          "include_refuted=True exports it")
    check(":" not in ytree.to_newick(tree, lengths=False), "lengths can be omitted")
    check("G-Z6219" in nwk and "G-Z6499" in nwk, "provisional nodes are exported")


def main() -> int:
    tree = test_real_tree_loads()
    test_validation()
    test_placements(tree)
    test_newick(tree)
    print()
    if failures:
        print(f"{len(failures)} FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
