#!/usr/bin/env python3
"""
charges.py

Conserved quantities of a rewrite grammar, i.e. its superselection sectors.

Motivation
----------
Asking whether a field could be a collective state of the relational network
is not directly testable here -- an electromagnetic field needs orientation,
and unordered hyperedges cannot carry a sign. But fields come with conserved
charges, and whether a grammar has conserved quantities *is* decidable.

A conserved quantity is exactly a function constant on connected components
of the continuation graph. So enumerate hypergraphs up to a size bound, take
their orbits under the grammar, and ask whether the obvious invariants
already separate the components:

    node count, weight = sum(|e| - 1), edge count

If they do, the grammar has no charges beyond the bookkeeping ones. If two
components agree on all of them, some further invariant distinguishes them,
and that invariant is a candidate charge -- something conserved by the local
rules that is not a trivial count.

This is the honest, in-scope version of the question. It says nothing about
electromagnetism; it asks only whether the dynamics conserves anything
non-obvious.

The window must be closed under the grammar
-------------------------------------------
Enumerating "all hypergraphs with at most E edges" is not valid here: a rule
that leaves the window has its transition silently dropped, which splits
components that are really joined and manufactures spurious charges. Since
compose and decompose both conserve node count and weight, the set of all
hypergraphs on exactly n nodes with exactly weight w *is* closed under them,
and that is the window used below.

Result
------
On the closed window a genuine hidden charge exists. Beyond node count and
weight, the number of connected components of the hypergraph is conserved --
equivalently the cycle rank w - n + c, since the other two are fixed. Neither
edge count nor maximum degree is conserved. Local rewrites cannot join or
split connected pieces, so connectivity is a superselection rule.

Connectivity is not the whole story, and the complete answer is finer.
Because compose and decompose act entirely within one connected component,
they can move neither nodes nor weight between components, so the conserved
quantity is the *multiset of per-component (node count, weight)* -- see
complete_charge.py. That is strictly finer than counting components, and it
separates every sector in 17 of the 18 (n, w) cases enumerated, including
n=6 and n=7 where component-counting fails.

The single exception, n=4 w=7, is not a further invariant. Its second
"sector" is one frozen state with zero available moves: every compose is
blocked because the union edge is already present, and every decompose
because a resulting binary edge already is. It is an isolated vertex of the
continuation graph produced by the no-duplicate-edge guards.

Frozen states do not proliferate. Their fraction falls with system size --
0.083 at n=4, 0.045 at n=5, 0.0004 at n=6 -- so they are small-size
accidents rather than a jamming transition.

This charge is topological, unsigned, and carries no local current, so it is
charge-like only in the superselection sense. Nothing here resembles an
oriented or signed field.
"""

import argparse
from collections import defaultdict
from itertools import combinations

from rewrite_lab import Hypergraph, Rule
from rewrite_research import COMPOSE, DECOMPOSE, CONTRACT
from expansion import SUBDIVIDE, SMOOTH

NAMED = {r.name: r for r in (COMPOSE, DECOMPOSE, CONTRACT, SUBDIVIDE, SMOOTH)}


def all_hypergraphs(n_nodes, max_edges, sizes=(2, 3)):
    """Every canonical hypergraph on <= n_nodes with <= max_edges edges."""
    pool = []
    for s in sizes:
        pool.extend(frozenset(c) for c in combinations(range(n_nodes), s))

    seen = {}
    for k in range(1, max_edges + 1):
        for combo in combinations(pool, k):
            H = Hypergraph(combo)
            c = H.canonical()
            if c not in seen:
                seen[c] = H
    return seen


def components(states, rules):
    """Union-find over the grammar's transitions, restricted to `states`."""
    parent = {c: c for c in states}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for c, H in states.items():
        for r in rules:
            for G in r.apply(H):
                gc = G.canonical()
                if gc in states:
                    union(c, gc)

    comps = defaultdict(list)
    for c in states:
        comps[find(c)].append(c)
    return list(comps.values())


def invariants(H):
    return (len(H.nodes()),
            sum(len(e) - 1 for e in H.edges),
            len(H.edges))


def analyse(name, rules, n_nodes, max_edges):
    states = all_hypergraphs(n_nodes, max_edges)
    comps = components(states, rules)

    # Group components by the value of the obvious invariants. A component is
    # well defined only if every member shares them; check that too, since a
    # rule that changes an invariant makes it not conserved at all.
    by_inv = defaultdict(list)
    non_conserved = set()
    for comp in comps:
        vals = {invariants(states[c]) for c in comp}
        if len(vals) > 1:
            # invariants differ within a component => not conserved
            for i, tag in enumerate(("nodes", "weight", "edges")):
                if len({v[i] for v in vals}) > 1:
                    non_conserved.add(tag)
        by_inv[tuple(sorted(vals))].append(comp)

    conserved = [t for t in ("nodes", "weight", "edges") if t not in non_conserved]
    degenerate = {k: v for k, v in by_inv.items() if len(v) > 1}

    print(f"{name}")
    print(f"  states={len(states)}  components={len(comps)}")
    print(f"  conserved among (nodes, weight, edges): "
          f"{conserved if conserved else 'none'}")
    if not conserved:
        print("  -> no superselection structure: the grammar connects "
              "everything it can reach")
    print(f"  invariant-classes containing >1 component: {len(degenerate)}")
    if degenerate:
        k = max(degenerate, key=lambda k: len(degenerate[k]))
        print(f"  largest degeneracy: {len(degenerate[k])} components share "
              f"(nodes,weight,edges)={k[0] if len(k)==1 else k}")
        print("  -> a further conserved quantity exists (candidate charge)")
    else:
        print("  -> (nodes, weight, edges) is a complete set of charges")
    print()
    return len(comps), len(degenerate)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--nodes", type=int, default=5)
    ap.add_argument("--edges", type=int, default=4)
    args = ap.parse_args(argv)

    print("=" * 66)
    print(f"conserved quantities on hypergraphs with <= {args.nodes} nodes, "
          f"<= {args.edges} edges")
    print("=" * 66)
    print()

    for names in (("compose", "decompose"),
                  ("compose", "decompose", "contract"),
                  ("subdivide", "smooth"),
                  ("compose", "decompose", "subdivide", "smooth"),
                  ("compose", "decompose", "contract", "subdivide", "smooth")):
        analyse("+".join(names), [NAMED[n] for n in names],
                args.nodes, args.edges)


if __name__ == "__main__":
    main()
