"""Load, validate and export the local Y tree (markers/tree_local.tsv).

Marker infrastructure in this directory has always handled two things: named
markers resolved through the YBrowse index, and bare coordinates. Neither
expresses *topology*. That G-Z6219 sits above G-L166 and G-Z6499 below it was
stated only in prose, which means no script could act on it and no script could
be wrong about it either. This module is the fix: one file states the tree, one
loader validates it, and everything downstream reads the same statement.

The status column is what makes that safe. Two of these nodes are not in any
published tree and one published label is refuted here, so a consumer that
ignores status will happily report a provisional node as settled. Callers get
`status` on every node and `Tree.path_nodes()` refuses to route through refuted
ones; deciding what to do with `provisional` and `putative` is the caller's job,
but it cannot claim not to have been told.

Imported by scripts in this directory; Python puts the script's own directory on
sys.path, so a plain `import ytree` works regardless of the caller's cwd. Also
runnable directly to check the file or emit Newick.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Sequence

STATUSES = ("published", "provisional", "putative", "refuted")

# Mapping from ylib.site_call()'s vocabulary onto the three states a node can be
# scored from. `low_power_1read_ancestral` deliberately lands in "nocall": one
# ancestral read is one molecule, and treating it as ancestral would convert "we
# could not tell" into "we tested and it was negative" -- the exact confusion
# ylib.site_call() was written to prevent. Anything not named here is a nocall.
DERIVED_CALLS = frozenset({"DERIVED", "DERIVED_1read_transversion", "derived"})
ANCESTRAL_CALLS = frozenset({"ancestral"})
MIXED_CALLS = frozenset({"mixed"})

# Nocalls that nonetheless lean one way. These do not change any call -- the rule
# stays exactly as ylib.site_call() fixed it -- but they are counted and reported,
# because a node called DERIVED off one transversion read while another of its
# defining markers has a lone ancestral read is a placement worth a second look,
# and silently reporting it as a clean derived call would be the same failure this
# project keeps documenting in other people's tables.
WEAK_ANCESTRAL_CALLS = frozenset({"low_power_1read_ancestral"})
WEAK_DERIVED_CALLS = frozenset({"nocall_damage_prone_1read", "nocall_1read_transition"})


def status_from_call(call: str) -> str:
    """One of derived / ancestral / mixed / nocall."""
    call = (call or "").strip()
    if call in DERIVED_CALLS:
        return "derived"
    if call in ANCESTRAL_CALLS:
        return "ancestral"
    if call in MIXED_CALLS:
        return "mixed"
    return "nocall"

# Statuses a terminal placement may rest on. `putative` is excluded because it
# means "not registered and/or not tested" -- see the header of tree_local.tsv.
DECISIVE_STATUSES = ("published", "provisional")

COLUMNS = ("node", "parent", "status", "defining", "equivalent", "evidence")

ROLE_DEFINING = "defining"
ROLE_EQUIVALENT = "equivalent"


class TreeError(ValueError):
    """The tree file is malformed. Always fatal -- never fall back to a guess."""


class Node:
    __slots__ = ("name", "parent", "status", "defining", "equivalent",
                 "evidence", "children", "depth")

    def __init__(self, name: str, parent: str | None, status: str,
                 defining: List[str], equivalent: List[str], evidence: str):
        self.name = name
        self.parent = parent
        self.status = status
        self.defining = defining
        self.equivalent = equivalent
        self.evidence = evidence
        self.children: List[str] = []
        self.depth = 0

    @property
    def markers(self) -> List[str]:
        return self.defining + self.equivalent

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<Node {self.name} status={self.status} depth={self.depth}>"


class Tree:
    def __init__(self, nodes: Dict[str, Node], root: str, path: Path):
        self.nodes = nodes
        self.root = root
        self.path = path
        # marker name -> (node name, role). One marker, one node, enforced at load.
        self.markers: Dict[str, tuple[str, str]] = {}
        for n in nodes.values():
            for m in n.defining:
                self.markers[m] = (n.name, ROLE_DEFINING)
            for m in n.equivalent:
                self.markers[m] = (n.name, ROLE_EQUIVALENT)

    def __len__(self) -> int:
        return len(self.nodes)

    def lookup(self, marker: str) -> tuple[str, str] | None:
        return self.markers.get(marker)

    def ancestors(self, name: str) -> List[str]:
        """Node names from the parent of `name` up to the root, nearest first."""
        out: List[str] = []
        cur = self.nodes[name].parent
        while cur is not None:
            out.append(cur)
            cur = self.nodes[cur].parent
        return out

    def descendants(self, name: str) -> Iterator[str]:
        for child in self.nodes[name].children:
            yield child
            yield from self.descendants(child)

    def preorder(self, start: str | None = None) -> Iterator[Node]:
        """Depth-first, children in file order -- so output matches the file."""
        stack = [start or self.root]
        while stack:
            node = self.nodes[stack.pop()]
            yield node
            stack.extend(reversed(node.children))

    def path_nodes(self, name: str) -> List[str]:
        """Root-to-node path, refuted nodes excluded.

        A refuted node is a label in circulation, not a place on the tree, so it
        may appear in a report but never in a path.
        """
        chain = [name] + self.ancestors(name)
        chain.reverse()
        return [n for n in chain if self.nodes[n].status != "refuted"]


def _split_markers(field: str, where: str) -> List[str]:
    field = (field or "").strip()
    if field in ("", "-", "."):
        return []
    out: List[str] = []
    for part in field.split(","):
        m = part.strip()
        if not m:
            continue
        if m in out:
            raise TreeError(f"{where}: marker {m!r} listed twice")
        out.append(m)
    return out


def load_tree(path: str | Path) -> Tree:
    path = Path(path)
    if not path.exists():
        raise TreeError(f"tree file not found: {path}")

    nodes: Dict[str, Node] = {}
    order: List[str] = []
    seen_markers: Dict[str, str] = {}
    header_seen = False

    with path.open("r", encoding="utf-8", newline="") as fh:
        for lineno, raw in enumerate(fh, 1):
            if raw.startswith("#") or not raw.strip():
                continue
            row = raw.rstrip("\n").split("\t")
            if not header_seen:
                if tuple(c.strip() for c in row[:len(COLUMNS)]) != COLUMNS:
                    raise TreeError(
                        f"{path}:{lineno}: header must be {'/'.join(COLUMNS)}, got {row}")
                header_seen = True
                continue
            if len(row) < len(COLUMNS):
                raise TreeError(
                    f"{path}:{lineno}: expected {len(COLUMNS)} columns, got {len(row)}")

            name, parent, status = (row[0].strip(), row[1].strip(), row[2].strip())
            where = f"{path}:{lineno} ({name})"
            if not name:
                raise TreeError(f"{where}: empty node name")
            if name in nodes:
                raise TreeError(f"{where}: node listed twice")
            if status not in STATUSES:
                raise TreeError(
                    f"{where}: status {status!r} not one of {'/'.join(STATUSES)}")

            defining = _split_markers(row[3], where)
            equivalent = _split_markers(row[4], where)
            if not defining:
                raise TreeError(f"{where}: no defining markers -- a node with no "
                                "defining marker cannot be called")
            for m in defining + equivalent:
                if m in seen_markers:
                    raise TreeError(
                        f"{where}: marker {m!r} already assigned to {seen_markers[m]}")
                seen_markers[m] = name

            nodes[name] = Node(name, parent if parent not in ("-", "") else None,
                               status, defining, equivalent, row[5].strip())
            order.append(name)

    if not header_seen:
        raise TreeError(f"{path}: no header row found")
    if not nodes:
        raise TreeError(f"{path}: no nodes")

    roots = [n for n in order if nodes[n].parent is None]
    if len(roots) != 1:
        raise TreeError(f"{path}: expected exactly one root, found {roots or 'none'}")
    root = roots[0]

    for name in order:
        parent = nodes[name].parent
        if parent is None:
            continue
        if parent not in nodes:
            raise TreeError(f"{path}: node {name!r} has unknown parent {parent!r}")
        nodes[parent].children.append(name)

    # Cycle check by walking every node to the root under a bounded step count.
    for name in order:
        cur, steps = nodes[name].parent, 0
        while cur is not None:
            steps += 1
            if steps > len(nodes):
                raise TreeError(f"{path}: cycle in parent chain at {name!r}")
            cur = nodes[cur].parent
        nodes[name].depth = steps

    return Tree(nodes, root, path)


def to_newick(tree: Tree, lengths: bool = True, include_refuted: bool = False) -> str:
    """Newick with internal node labels; branch length = count of defining markers.

    Branch length here is a marker count, not a coalescent time or a rate-scaled
    distance. It is meaningful to a placement tool that walks branches counting
    diagnostic sites and meaningless to anything that reads it as time. Pass
    lengths=False if the consumer would misread it.

    Refuted nodes are dropped by default, for the same reason path_nodes() skips
    them: exporting one would assert a branch this project says does not exist.
    """
    def render(name: str) -> str:
        node = tree.nodes[name]
        kids = [c for c in node.children
                if include_refuted or tree.nodes[c].status != "refuted"]
        inner = f"({','.join(render(c) for c in kids)})" if kids else ""
        label = f"{inner}{name}"
        if lengths:
            label += f":{len(node.defining)}"
        return label

    return render(tree.root) + ";"


def marker_names(field: str) -> List[str]:
    """Marker names in an ID field. Compound IDs are common: the pileup tables
    carry `L91,S285,PF3246` and `Z6516,FGC5675` as single cells."""
    field = (field or "").strip().strip('"')
    if not field or field == ".":
        return []
    return [t.strip() for t in re.split(r"[;,]", field) if t.strip()]


class TreeScorer:
    """Place a sample on the local tree from its per-marker calls.

    Deliberately blind to label columns. The ranking in y_path_rank.py works from
    HG/ISOGG labels, which is exactly the layer this project keeps finding to be
    wrong -- a label one node deeper than its own evidence. This scorer maps
    marker names to nodes and calls each node from its own *defining* markers
    only, so the two can disagree; when they do, the disagreement is the result.

    Shared by y_path_rank.py (VCF-derived marker_status tables) and
    y_tree_place.py (read-level pileup tables) for the reason ylib.py exists: two
    consumers of one rule drift apart if each keeps its own copy.
    """

    def __init__(self, tree: Tree, name: str = ""):
        self.tree = tree
        self.name = name
        self.def_derived: Dict[str, int] = defaultdict(int)
        self.def_derived_tv: Dict[str, int] = defaultdict(int)
        self.def_ancestral: Dict[str, int] = defaultdict(int)
        self.def_nocall: Dict[str, int] = defaultdict(int)
        self.eq_derived: Dict[str, int] = defaultdict(int)
        self.eq_ancestral: Dict[str, int] = defaultdict(int)
        self.eq_nocall: Dict[str, int] = defaultdict(int)
        self.def_weak_ancestral: Dict[str, int] = defaultdict(int)
        self.def_weak_derived: Dict[str, int] = defaultdict(int)
        self.unmapped: set[str] = set()

    def add(self, id_field: str, status: str, transversion: bool = False,
            weak: str = "") -> None:
        """Record one marker observation.

        `status` is derived/ancestral/mixed/nocall. `weak` is "ancestral" or
        "derived" when the observation is a nocall that leans that way -- counted
        for reporting, never for the call.
        """
        for name in marker_names(id_field):
            hit = self.tree.lookup(name)
            if hit is None:
                self.unmapped.add(name)
                continue
            node, role = hit
            defining = role == ROLE_DEFINING
            if status in ("derived", "mixed"):
                (self.def_derived if defining else self.eq_derived)[node] += 1
                if defining and transversion:
                    self.def_derived_tv[node] += 1
            if status in ("ancestral", "mixed"):
                (self.def_ancestral if defining else self.eq_ancestral)[node] += 1
            if status == "nocall":
                (self.def_nocall if defining else self.eq_nocall)[node] += 1
                if defining and weak == "ancestral":
                    self.def_weak_ancestral[node] += 1
                elif defining and weak == "derived":
                    self.def_weak_derived[node] += 1
            # A marker belongs to one node, so the first hit is the only hit.
            return

    def add_call(self, id_field: str, call: str, mutation_class: str = "") -> None:
        """Record one row of a pileup table, translating ylib.site_call()'s vocabulary."""
        call = (call or "").strip()
        weak = ("ancestral" if call in WEAK_ANCESTRAL_CALLS else
                "derived" if call in WEAK_DERIVED_CALLS else "")
        self.add(id_field, status_from_call(call),
                 transversion=mutation_class.strip() == "transversion", weak=weak)

    def call(self, node: str) -> str:
        d, a = self.def_derived[node], self.def_ancestral[node]
        if d and a:
            return "CONFLICT"
        if d:
            return "DERIVED"
        if a:
            return "ANCESTRAL"
        return "NOCALL"

    def supported(self, node: str) -> bool:
        """Derived here, and nothing above contradicts it.

        NOCALL ancestors are permitted: absence of coverage is not evidence of
        absence, and in these libraries most of the backbone is uncovered in most
        samples. A single ANCESTRAL or CONFLICT ancestor is disqualifying.
        """
        if self.tree.nodes[node].status == "refuted":
            return False
        if self.call(node) != "DERIVED":
            return False
        return all(self.call(a) in ("DERIVED", "NOCALL")
                   for a in self.tree.ancestors(node))

    def terminal(self) -> List[str]:
        """Deepest supported nodes with no supported node below them.

        Returns a list: more than one means the reads place the sample on two
        branches at once, which is a result about the data, not something to
        resolve by picking the higher score.
        """
        supported = [n for n in self.tree.nodes if self.supported(n)]
        return [n for n in supported
                if not any(self.supported(d) for d in self.tree.descendants(n))]

    def excluded(self) -> List[str]:
        """Nodes the sample is demonstrably NOT at, nearest the root first."""
        out = [n for n in self.tree.nodes
               if self.call(n) == "ANCESTRAL"
               and self.tree.nodes[n].status != "refuted"]
        return sorted(out, key=lambda n: (self.tree.nodes[n].depth, n))

    def conflicts(self) -> List[str]:
        return sorted((n for n in self.tree.nodes if self.call(n) == "CONFLICT"),
                      key=lambda n: (self.tree.nodes[n].depth, n))

    def ladder(self, node: str) -> str:
        return " > ".join(
            f"{n}{format_status(self.tree.nodes[n].status)}({self.call(n)[0]})"
            for n in self.tree.path_nodes(node))

    def on_path(self) -> set[str]:
        out: set[str] = set()
        for t in self.terminal():
            out.update(self.tree.path_nodes(t))
        return out

    def rows(self) -> Iterator[List]:
        """Per-node rows in file order, matching NODE_COLUMNS."""
        terminal, on_path = set(self.terminal()), self.on_path()
        for node in self.tree.preorder():
            n = node.name
            yield [
                self.name, node.depth, n, node.parent or "-", node.status,
                self.call(n),
                self.def_derived[n], self.def_derived_tv[n], self.def_ancestral[n],
                self.def_nocall[n], self.eq_derived[n], self.eq_ancestral[n],
                self.eq_nocall[n],
                "yes" if n in on_path else "no",
                "yes" if n in terminal else "no",
                self.def_weak_ancestral[n], self.def_weak_derived[n],
                ";".join(self.caveats(n)) or "-",
                node.evidence,
            ]

    def caveats(self, node: str) -> List[str]:
        """Reasons to distrust a call that the rules nonetheless allow."""
        out: List[str] = []
        if self.call(node) == "DERIVED" and self.def_weak_ancestral[node]:
            out.append(f"weak_ancestral={self.def_weak_ancestral[node]}")
        if self.call(node) == "ANCESTRAL" and self.def_weak_derived[node]:
            out.append(f"weak_derived={self.def_weak_derived[node]}")
        if self.call(node) == "DERIVED" and self.def_derived[node] == 1:
            out.append("single_marker")
        return out

    def summary(self) -> str:
        """One line: what this sample is, in the form a report can quote."""
        terminal = self.terminal()
        if not terminal:
            excluded = self.excluded()
            why = f"; excluded at {', '.join(excluded)}" if excluded else ""
            return f"no placement{why}"
        if len(terminal) > 1:
            return f"AMBIGUOUS: {', '.join(terminal)} supported at once"
        t = terminal[0]
        node = self.tree.nodes[t]
        out = (f"{t} [{node.status}] "
               f"({self.def_derived[t]} derived, {self.def_derived_tv[t]} transversion)")
        caveats = self.caveats(t)
        if caveats:
            out += f" [{'; '.join(caveats)}]"
        return out


NODE_COLUMNS = (
    "sample", "depth", "node", "parent", "status", "call",
    "def_derived", "def_derived_transversion", "def_ancestral", "def_nocall",
    "eq_derived", "eq_ancestral", "eq_nocall", "on_path", "terminal",
    "def_weak_ancestral", "def_weak_derived", "caveats", "evidence",
)


def format_status(status: str) -> str:
    """One-character suffix so a node's standing survives being pasted anywhere."""
    return {"published": "", "provisional": "*", "putative": "?", "refuted": "!"}[status]


STATUS_LEGEND = "status suffixes: * provisional, ? putative, ! refuted (not on any path)"


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Check or export the local Y tree.")
    ap.add_argument("--tree", default="markers/tree_local.tsv", help="tree TSV")
    ap.add_argument("--newick", help="write Newick to this path (- for stdout)")
    ap.add_argument("--no-lengths", action="store_true",
                    help="omit branch lengths from the Newick")
    ap.add_argument("--markers-out", help="write a marker/node/role TSV here")
    args = ap.parse_args(argv)

    try:
        tree = load_tree(args.tree)
    except TreeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    counts: Dict[str, int] = {s: 0 for s in STATUSES}
    for node in tree.nodes.values():
        counts[node.status] += 1
    print(f"{tree.path}: {len(tree)} nodes, root {tree.root}, "
          + ", ".join(f"{counts[s]} {s}" for s in STATUSES if counts[s]))
    print(f"{len(tree.markers)} markers, no duplicates, no cycles")

    for node in tree.preorder():
        indent = "  " * node.depth
        marks = f"{len(node.defining)} defining"
        if node.equivalent:
            marks += f", {len(node.equivalent)} equivalent"
        print(f"  {indent}{node.name}{format_status(node.status)}  [{marks}]")
    print(f"  ({STATUS_LEGEND})")

    if args.markers_out:
        with open(args.markers_out, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["marker", "node", "role", "node_status"])
            for node in tree.preorder():
                for m in node.defining:
                    w.writerow([m, node.name, ROLE_DEFINING, node.status])
                for m in node.equivalent:
                    w.writerow([m, node.name, ROLE_EQUIVALENT, node.status])
        print(f"wrote {args.markers_out}")

    if args.newick:
        nwk = to_newick(tree, lengths=not args.no_lengths)
        if args.newick == "-":
            print(nwk)
        else:
            Path(args.newick).write_text(nwk + "\n", encoding="utf-8")
            print(f"wrote {args.newick}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
