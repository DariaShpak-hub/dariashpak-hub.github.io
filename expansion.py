#!/usr/bin/env python3
"""
expansion.py

Growth of accessible state space in a rewrite universe, and whether
"accelerated expansion" is even well defined without a derived clock.

The proposal
------------
Define an accessible-volume entropy on a depth-k ball around the seed,

    S_acc(k) = ln |B_k|,      B_k = states reachable within k rewrites,

read an effective scale factor off it, and ask whether that scale factor
accelerates.

Two obstructions, both computable rather than philosophical
-----------------------------------------------------------
1. Finiteness. The continuation graph of a weight-conserving rule set is a
   fixed finite set, so |B_k| rises and then saturates. Saturating growth
   decelerates by construction: no such rule set can produce sustained
   expansion, whatever observable is read off it. Expansion requires rules
   that create nodes.

2. Reparameterization. "Acceleration" is the sign of a second derivative
   with respect to a time coordinate, and that sign is not invariant. For
   any monotone reclocking k -> tau(k) the sign of d^2 a / d tau^2 can be
   changed at will. Rewrite depth is a choice of clock, not a fact about the
   rule set, so a positive result under it means nothing until the clock is
   itself derived from the system. This module measures the same expansion
   under two clocks -- rewrite depth, and an entropic clock tau = ln|B_k| --
   to show the sign of the acceleration is a property of the choice.
"""

import argparse
import math

from rewrite_lab import Hypergraph, Rule, continuation_graph, path_seed
from rewrite_research import COMPOSE, DECOMPOSE, RULE_SETS


def subdivide_edge(H):
    """
    {x,y} -> {x,z},{z,y} with z a fresh node.

    A node-creating rule: the only kind that can make the accessible state
    space grow without bound. Weight increases by one per application, so it
    has no inverse in the bundled rule set.
    """
    nodes = H.nodes()
    fresh = (max(nodes) + 1) if nodes else 0
    seen = set()

    for e in H.edges:
        if len(e) != 2:
            continue
        x, y = sorted(e)
        new_edges = set(H.edges)
        new_edges.remove(e)
        new_edges.add(frozenset((x, fresh)))
        new_edges.add(frozenset((fresh, y)))
        G = Hypergraph(new_edges)
        c = G.canonical()
        if c not in seen:
            seen.add(c)
            yield G


SUBDIVIDE = Rule("subdivide", subdivide_edge)


def smooth_node(H):
    """
    {x,z},{z,y} -> {x,y} where z has degree exactly 2: the inverse of
    subdivide.

    Adding this makes the expanding sector reversible, which is what allows
    unbounded accessible-volume growth to coexist with zero entropy
    production -- an infinite symmetric graph has detailed balance and still
    has |B_k| growing without bound.
    """
    inc = {}
    for e in H.edges:
        for v in e:
            inc.setdefault(v, []).append(e)

    seen = set()
    for z, es in inc.items():
        if len(es) != 2 or any(len(e) != 2 for e in es):
            continue
        a, b = es
        x = next(iter(a - {z}), None)
        y = next(iter(b - {z}), None)
        if x is None or y is None or x == y:
            continue
        new_edges = set(H.edges)
        new_edges.discard(a)
        new_edges.discard(b)
        if frozenset((x, y)) in new_edges:
            continue
        new_edges.add(frozenset((x, y)))
        G = Hypergraph(new_edges)
        c = G.canonical()
        if c not in seen:
            seen.add(c)
            yield G


SMOOTH = Rule("smooth", smooth_node)


def ball_volumes(seed, rules, max_depth, max_states=400000):
    """|B_k| for k = 0..max_depth: states within k rewrites of the seed."""
    seen = {seed.canonical()}
    frontier = [seed]
    vols = [1]

    for _ in range(max_depth):
        nxt = []
        for H in frontier:
            for r in rules:
                for G in r.apply(H):
                    c = G.canonical()
                    if c not in seen:
                        seen.add(c)
                        nxt.append(G)
        vols.append(len(seen))
        frontier = nxt
        if not frontier or len(seen) >= max_states:
            break
    return vols


def second_difference(xs, ys):
    """Discrete d^2 y / dx^2 on a possibly non-uniform grid."""
    out = []
    for i in range(1, len(xs) - 1):
        h1 = xs[i] - xs[i - 1]
        h2 = xs[i + 1] - xs[i]
        if h1 <= 0 or h2 <= 0:
            out.append(float("nan"))
            continue
        d1 = (ys[i] - ys[i - 1]) / h1
        d2 = (ys[i + 1] - ys[i]) / h2
        out.append(2.0 * (d2 - d1) / (h1 + h2))
    return out


def fit_growth(vols):
    """Return (poly_exponent, exp_rate) fitted over the non-saturated range."""
    ks = list(range(len(vols)))
    use = [(k, v) for k, v in zip(ks, vols) if k >= 1 and v > 1]
    if len(use) < 3:
        return float("nan"), float("nan")

    def slope(xs, ys):
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        return (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
                if sxx > 0 else float("nan"))

    lk = [math.log(k) for k, _ in use]
    lv = [math.log(v) for _, v in use]
    k = [float(k) for k, _ in use]
    return slope(lk, lv), slope(k, lv)


def report(name, seed, rules, max_depth):
    vols = ball_volumes(seed, rules, max_depth)
    ks = list(range(len(vols)))
    poly, expo = fit_growth(vols)

    print(f"{name}")
    print(f"  |B_k| = {vols}")
    print(f"  growth fit: |B_k| ~ k^{poly:.3f}   or   e^({expo:.3f} k)")

    # Effective scale factor. With no derived notion of spatial dimension the
    # honest choice is a = |B_k|^(1/d) for a fitted d; d only rescales a, it
    # cannot change the sign of the acceleration, which is the point at issue.
    d = poly if poly and poly == poly and poly > 0.2 else 1.0
    a = [v ** (1.0 / d) for v in vols]

    # Clock 1: rewrite depth.
    acc_k = second_difference([float(x) for x in ks], a)
    # Clock 2: entropic time tau = ln|B_k|.
    tau = [math.log(v) for v in vols]
    acc_t = second_difference(tau, a)

    def sgn(seq):
        s = [x for x in seq if x == x]
        if not s:
            return "n/a"
        p = sum(1 for x in s if x > 1e-12)
        n = sum(1 for x in s if x < -1e-12)
        return f"{p} positive, {n} negative"

    print(f"  a''  under depth clock  : {sgn(acc_k)}")
    print(f"  a''  under entropic clock: {sgn(acc_t)}")
    print(f"  saturated: {vols[-1] == vols[-2] if len(vols) > 1 else False}")
    print()
    return vols, acc_k, acc_t


def star(m):
    """Star with m branches: unlike a path, its edges are not interchangeable,
    so subdivision produces a multiset of branch lengths rather than a single
    number, and the orbit is m-1 dimensional instead of 1 dimensional."""
    return Hypergraph([{0, i} for i in range(1, m + 1)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--depth", type=int, default=14)
    ap.add_argument("--path", type=int, default=6)
    args = ap.parse_args(argv)

    print("=" * 72)
    print("accessible-volume growth |B_k| and the clock-dependence of a''")
    print("=" * 72)
    print()

    seed = path_seed(args.path)
    report("G1 (compose+decompose): weight-conserving, finite",
           seed, RULE_SETS["G1-reversible"]["rules"], args.depth)
    report("G1 + subdivide: node-creating, unbounded",
           path_seed(2), [COMPOSE, DECOMPOSE, SUBDIVIDE], args.depth)
    report("subdivide only: node-creating, no reversible sector",
           path_seed(2), [SUBDIVIDE], args.depth)
    report("subdivide+smooth: node-creating AND reversible",
           star(3), [SUBDIVIDE, SMOOTH], args.depth)


if __name__ == "__main__":
    main()

