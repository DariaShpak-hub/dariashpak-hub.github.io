#!/usr/bin/env python3
"""
boundary_entropy.py

Does conditioning on the conserved sector suppress bulk degrees of freedom?

The question is whether the connectivity / cycle-rank charge found in
charges.py actually eliminates interior freedom, rather than whether S ~ A
holds at scales we cannot reach. For a region R with interior edges X and
crossing (boundary) edges Y, measure

    S(X | Y, C)              conditional interior entropy in sector C
    s_int = S(X | Y, C) / V  interior entropy density, V = interior slots
    I(X : Y) = S(X) + S(Y) - S(X, Y)

If s_int tends to a positive constant the system has ordinary bulk
information and the holographic branch fails. If s_int falls toward zero
while the boundary entropy stays finite, the constraint is genuinely
determining the interior from the boundary.

Labeled, not canonical
----------------------
A region is a set of node labels, which is not well defined across
isomorphism classes -- node 0 of one canonical form need not correspond to
node 0 of another. So this module works with *labeled* hypergraphs on a fixed
site set, the ordinary setup for entanglement-style region analysis, and the
ensemble is uniform over labeled graphs at fixed node count, weight, and
connectivity.

The volume proxy V is the number of possible interior edges, C(r,2) + C(r,3),
i.e. the count of bulk degrees of freedom available -- not the number
realised, which would build the answer into the denominator.

Result: volume law. The holographic branch fails.
-------------------------------------------------
s_int does not fall to zero in any sector. At n=7, w=5 it plateaus at roughly
0.186 (sector 2), 0.177 (sector 3) and 0.127 (sector 4) as r grows, and in
the most constrained sector it turns back upward rather than continuing down.
S(X|Y,C) grows linearly in the number of interior slots -- 0.34, 1.05, 2.17,
3.88, 6.50 against slots 1, 4, 10, 20, 35 -- which is a volume law with a
reduced coefficient.

Conditioning on the conserved sector does remove information: at r=6 the
conditional interior entropy falls from 6.98 pooled to 6.50, 6.18 and 4.46 in
sectors 2, 3 and 4, a reduction of 7% to 36%. So the connectivity charge is
doing real work. But it rescales the entropy density without changing the
scaling class.

The mutual information diagnostic agrees. Boundary determination would need
I(X:Y) ~ S(X); measured, I(X:Y) is 1.64 to 2.36 against S(X) of 6.82 to 8.15,
so the boundary fixes roughly a fifth of the interior information. The bulk
carries ordinary independent information.
"""

import argparse
import math
from collections import defaultdict
from itertools import combinations


def enumerate_labeled(n, w):
    """All labeled edge-sets on n sites with weight exactly w."""
    binary = [frozenset(c) for c in combinations(range(n), 2)]
    ternary = [frozenset(c) for c in combinations(range(n), 3)]
    for t in range(w // 2 + 1):
        b = w - 2 * t
        if b < 0:
            continue
        for tt in combinations(ternary, t):
            for bb in combinations(binary, b):
                yield frozenset(tt) | frozenset(bb)


def connectivity(edges, n):
    """Number of connected components over all n sites (isolated count)."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in edges:
        vs = sorted(e)
        for v in vs[1:]:
            a, b = find(vs[0]), find(v)
            if a != b:
                parent[a] = b
    return len({find(v) for v in range(n)})


def entropy(counts):
    tot = sum(counts)
    if tot <= 0:
        return 0.0
    return -sum((c / tot) * math.log(c / tot) for c in counts if c > 0)


def analyse(n, w, r):
    """Return per-sector statistics for region R = {0..r-1}."""
    R = frozenset(range(r))
    # joint[(sector)][(boundary)][interior] = count
    joint = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for edges in enumerate_labeled(n, w):
        c = connectivity(edges, n)
        interior = frozenset(e for e in edges if e <= R)
        bound = frozenset(e for e in edges if (e & R) and not (e <= R))
        joint[c][bound][interior] += 1

    slots = len(list(combinations(range(r), 2))) + \
        len(list(combinations(range(r), 3)))

    out = []
    for c, by_bound in sorted(joint.items()):
        total = sum(sum(d.values()) for d in by_bound.values())
        if total < 50:
            continue
        # S(X|Y) = sum_y p(y) H(X|Y=y)
        s_cond = 0.0
        x_counts = defaultdict(int)
        y_counts = defaultdict(int)
        xy_counts = []
        for y, d in by_bound.items():
            ny = sum(d.values())
            s_cond += (ny / total) * entropy(list(d.values()))
            y_counts[y] += ny
            for x, k in d.items():
                x_counts[x] += k
                xy_counts.append(k)
        sx = entropy(list(x_counts.values()))
        sy = entropy(list(y_counts.values()))
        sxy = entropy(xy_counts)
        out.append({"sector": c, "states": total, "slots": slots,
                    "S_cond": s_cond, "s_int": s_cond / slots if slots else 0.0,
                    "S_X": sx, "S_Y": sy, "I": sx + sy - sxy})

    # pooled over all sectors: conditioning on C removed
    pooled = defaultdict(lambda: defaultdict(int))
    for c, by_bound in joint.items():
        for y, d in by_bound.items():
            for x, k in d.items():
                pooled[y][x] += k
    tot = sum(sum(d.values()) for d in pooled.values())
    s_pool = sum((sum(d.values()) / tot) * entropy(list(d.values()))
                 for d in pooled.values()) if tot else 0.0
    return out, s_pool, slots


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--n", type=int, nargs="+", default=[6, 7])
    ap.add_argument("--w", type=int, default=5)
    args = ap.parse_args(argv)

    for n in args.n:
        print("=" * 74)
        print(f"n = {n} sites, weight = {args.w}")
        print("=" * 74)
        print("%3s %6s %7s %8s %9s %9s %9s %9s %9s" % (
            "r", "slots", "sector", "states", "S(X|Y,C)", "s_int",
            "S(X|Y)pool", "I(X:Y)", "S(X)"))
        for r in range(2, n):
            rows, s_pool, slots = analyse(n, args.w, r)
            for row in rows:
                print("%3d %6d %7d %8d %9.4f %9.4f %9.4f %9.4f %9.4f" % (
                    r, slots, row["sector"], row["states"], row["S_cond"],
                    row["s_int"], s_pool, row["I"], row["S_X"]))
        print()


if __name__ == "__main__":
    main()
