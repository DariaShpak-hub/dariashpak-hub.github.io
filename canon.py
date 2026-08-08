#!/usr/bin/env python3
"""
canon.py

Canonical labelling of finite unlabelled hypergraphs by individualization-
refinement, replacing brute-force minimisation over all n! relabellings.

Method
------
1. Colour nodes by an isomorphism-invariant local signature (degree together
   with the multiset of incident edge sizes).
2. Refine to an equitable partition by Weisfeiler-Leman: a node's new colour
   is its old colour plus the multiset, over incident edges, of the edge size
   and the colours of that edge's other members. Iterate until the number of
   colour classes stops growing.
3. If the partition is discrete, every node has a unique rank and relabelling
   by that rank yields a canonical form directly.
4. Otherwise pick a target cell -- the non-singleton cell of lowest colour,
   a choice that depends only on structure -- individualize each of its nodes
   in turn, refine again, recurse, and take the lexicographic minimum over
   the resulting leaves.

Refinement prunes almost all of the n! search in practice while step 4
preserves exactness on symmetric graphs, where refinement alone would stall.

A note on the renumbering inside _refine: colours are compressed to small
integers by sorted signature order, which is safe *here* because the result
of this module is a relabelled edge set and all cross-graph comparison
happens on that edge set. Compressing colours per-graph and then comparing
the colours themselves across graphs would not be safe -- it destroys the
correspondence between colour values in different graphs.
"""

from collections import defaultdict
from itertools import permutations


def _incidence(edges):
    inc = defaultdict(list)
    for e in edges:
        for v in e:
            inc[v].append(e)
    return inc


def _compress(sig):
    order = {s: i for i, s in enumerate(sorted(set(sig.values())))}
    return {v: order[s] for v, s in sig.items()}


def _refine(colour, inc):
    """Iterate WL refinement to an equitable partition."""
    ncls = len(set(colour.values()))
    while True:
        sig = {}
        for v, es in inc.items():
            s = sorted((len(e), tuple(sorted(colour[u] for u in e if u != v)))
                       for e in es)
            sig[v] = (colour[v], tuple(s))
        colour = _compress(sig)
        n2 = len(set(colour.values()))
        if n2 == ncls:
            return colour
        ncls = n2


def _individualize(colour, v):
    """Split v out of its cell, ordering it just after its former cell-mates."""
    raw = {u: (2 * c + (1 if u == v else 0)) for u, c in colour.items()}
    return _compress({u: (c,) for u, c in raw.items()})


def _leaf(colour, edges):
    return tuple(sorted(tuple(sorted(colour[v] for v in e)) for e in edges))


def _search(colour, inc, edges):
    colour = _refine(colour, inc)

    cells = defaultdict(list)
    for v, c in colour.items():
        cells[c].append(v)

    target = None
    for c in sorted(cells):
        if len(cells[c]) > 1:
            target = c
            break

    if target is None:
        return _leaf(colour, edges)

    best = None
    for v in sorted(cells[target]):
        cand = _search(_individualize(colour, v), inc, edges)
        if best is None or cand < best:
            best = cand
    return best


def canonical_form(edges):
    """
    Canonical representative of a hypergraph under node relabelling.

    ``edges`` is any iterable of iterables of hashable node labels. Returns a
    sorted tuple of sorted node-index tuples, identical for isomorphic inputs.
    """
    edges = [frozenset(e) for e in edges]
    inc = _incidence(edges)
    if not inc:
        return ()

    init = {v: (len(es), tuple(sorted(len(e) for e in es)))
            for v, es in inc.items()}
    return _search(_compress(init), inc, edges)


def canonical_bruteforce(edges):
    """Reference implementation: minimise over all n! relabellings."""
    edges = [frozenset(e) for e in edges]
    nodes = sorted({v for e in edges for v in e})
    if not nodes:
        return ()
    best = None
    for p in permutations(range(len(nodes))):
        m = dict(zip(nodes, p))
        img = tuple(sorted(tuple(sorted(m[v] for v in e)) for e in edges))
        if best is None or img < best:
            best = img
    return best
