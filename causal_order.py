#!/usr/bin/env python3
"""
causal_order.py

Order + number = geometry. Does this project have the order?

Malament's theorem says the causal structure of a spacetime fixes its geometry
up to a conformal factor, and causal set theory adds that counting the discrete
elements supplies the missing scale. Order gives light cones, number gives
volume. That is the one route to a metric this project has never tried, and it
does not require injecting a manifold -- only a partial order.

So the question is whether necessity.py's dynamics already contains one.

The claim to be tested, and why I expect it to fail
-----------------------------------------------------
A partial order needs ANTISYMMETRY: if A precedes B then B does not precede A.
The consistency move in necessity.py is the double swap, ab + cd -> ad + cb,
and that move is its own inverse -- from ad + cb you can swap straight back. So
the move graph should be undirected, reachability should be an equivalence
relation rather than an order, and the "drift" measured there should be a
reversible walk's bias toward high-degree states rather than a current.

Reversible dynamics cannot manufacture antisymmetry. If that is right, the
drift is not a causal order and cannot be made into one.

Where an order could actually come from
-----------------------------------------
Causal sets do not get their order from rearrangement. They get it from GROWTH:
elements are added one at a time, each to the future of some subset, and the
order is the birth relation. Rideout and Sorkin's sequential growth models work
exactly this way.

This project has growth models. coloured_memory.py creates node w from node v
by stepping along a coloured edge, so there is a natural parent relation whose
transitive closure IS a partial order -- built from the project's own dynamics,
with no manifold assumed.

So the second half asks what that order looks like geometrically.

The measurement, with a control that decides it
-------------------------------------------------
For a causal set sprinkled into d-dimensional Minkowski space the ORDERING
FRACTION -- the proportion of pairs that are causally related -- is a known
decreasing function of d. That is the Myrheim-Meyer dimension estimator, and it
uses order alone: no coordinates, no metric, no distances.

Rather than quote a formula, this calibrates empirically. Sprinkle points into
a causal diamond in d = 2, 3, 4 and measure the fraction directly. Then measure
the same quantity on the project's growth poset and read its dimension off the
same curve.

A tree gives a very low ordering fraction. A chain gives 1. If the project's
poset lands near a tree, growth is producing causal structure with no room in
it -- which would be the same negative as the metric results, reached from the
order side.

Result 1: the consistency move cannot be an order, exactly
------------------------------------------------------------
    moves examined                      184320
    moves with no reverse               0

Every single move is reversible. The relation is symmetric, so reachability is
an equivalence relation and antisymmetry fails everywhere. No relabelling or
reorientation fixes this -- it is a property of the move, not of the bookkeeping.

And the "drift" is not a current. A reversible walk on an undirected graph has
stationary distribution proportional to degree and satisfies detailed balance,
so the net probability current across every edge is exactly zero. Degree ranges
1 to 18 across 29894 states with moves, which is why the stationary
distribution is non-uniform -- but non-uniform is not directed.

    a bias toward high-multiplicity states is where a reversible walk spends
    its time, not a direction in time

So the proposal fails at the first gate, and it fails for a reason that
generalises: rearrangement dynamics cannot manufacture antisymmetry.

Result 2: the Minkowski calibration works
-------------------------------------------
    d       ordering fraction
    1              1.0000      (a chain)
    2              0.5014
    3              0.2229
    4              0.0950

The d = 2 row is the check that matters. In 1+1 dimensions the exact value is
1/2 -- in light-cone coordinates the diamond is a square and two points are
related exactly when both coordinates are ordered, which happens half the time.
Measured 0.5014. The estimator is sound.

And it is size-stable, as it must be: 0.528, 0.490, 0.493 at n = 300, 600,
1200. Minkowski keeps a finite fraction of pairs related no matter how many
points you sprinkle.

Result 3: the project's growth poset is a tree
------------------------------------------------
coloured_memory.py's birth relation, transitive closure taken:

    k = 2    n = 600   1201   2400   4801
             f = .0354  .0183  .0104  .0061

    k = 3    n = 600   1200   2401   4800
             f = .0226  .0125  .0063  .0031

Two things. The values at any size are far below even d = 4's 0.095. And they
FALL roughly as 1/n, while Minkowski holds constant. That is the discriminator
stated in advance, and it comes out unambiguously: this is a tree, not a causal
set with a dimension.

The reason is visible in the construction. Every node is created by a single
step from a single parent, so its causal past is one chain. In a causal set
each new element is to the future of a substantial DOWNSET of what already
exists -- that is what Rideout and Sorkin's sequential growth models do, and it
is what keeps a finite fraction of pairs related.

Verdict
-------
Order + number = geometry is the right route and this project does not have the
order.

    rearrangement dynamics    provably cannot give one     0 of 184320 asymmetric
    growth dynamics           gives one, and it is a tree  f ~ 1/n

What is needed is now specific rather than vague: a birth rule where each new
element enters the future of MANY existing elements rather than one. That is a
concrete change to the growth rule, it is buildable, and the ordering fraction
is the measurement that decides whether it worked -- with the Minkowski
calibration above already in place to read the answer against.

Worth noting what this does not say. It does not say the project's growth is
wrong; the coloured complex was built to produce SPACE, and it does, with
measured isotropy of 7e-4 at k = 2. It says that a spatial growth rule and a
causal growth rule are different objects, and having one does not give the
other -- which is the same separation-of-layers lesson as everywhere else here.
"""

import numpy as np
from itertools import combinations

from coloured_memory import ColouredComplex


# ------------------------------------------------- part 1: is there an order?

PAIRS6 = list(combinations(range(6), 2))
BIT6 = {p: i for i, p in enumerate(PAIRS6)}


def swap_moves(s):
    out = []
    present = [PAIRS6[i] for i in range(15) if s >> i & 1]
    for (a, b), (c, d) in combinations(present, 2):
        if len({a, b, c, d}) != 4:
            continue
        for x, y in (((a, d), (c, b)), ((a, c), (b, d))):
            i = BIT6.get((min(x), max(x)))
            j = BIT6.get((min(y), max(y)))
            if i is None or j is None or (s >> i & 1) or (s >> j & 1):
                continue
            t = s & ~(1 << BIT6[(a, b)]) & ~(1 << BIT6[(c, d)])
            out.append(t | (1 << i) | (1 << j))
    return out


def check_symmetry(n_states=1 << 15):
    """Is the move relation symmetric? If so it cannot be an order."""
    asym = 0
    total = 0
    deg = np.zeros(n_states, dtype=np.int32)
    fwd = {}
    for s in range(n_states):
        ms = set(swap_moves(s))
        fwd[s] = ms
        deg[s] = len(ms)
    for s in range(n_states):
        for t in fwd[s]:
            total += 1
            if s not in fwd[t]:
                asym += 1
    return asym, total, deg


# ------------------------------------- part 2: ordering fraction of a poset

def sprinkle(d, n, rng):
    """Poisson sprinkle into a causal diamond in d-dimensional Minkowski."""
    pts = []
    while len(pts) < n:
        m = 4 * (n - len(pts)) + 16
        t = rng.uniform(0.0, 1.0, m)
        x = rng.uniform(-0.5, 0.5, (m, d - 1))
        r = np.sqrt((x ** 2).sum(axis=1))
        ok = (r < t) & (r < 1.0 - t)
        for i in np.nonzero(ok)[0]:
            pts.append((float(t[i]), x[i].copy()))
            if len(pts) >= n:
                break
    pts.sort(key=lambda p: p[0])
    return pts


def ordering_fraction_minkowski(d, n, rng):
    pts = sprinkle(d, n, rng)
    t = np.array([p[0] for p in pts])
    X = np.array([p[1] for p in pts])
    rel = 0
    for i in range(n):
        dt = t[i + 1:] - t[i]
        dx = np.sqrt(((X[i + 1:] - X[i]) ** 2).sum(axis=1))
        rel += int(np.count_nonzero(dt > dx))
    return rel / (n * (n - 1) / 2)


def ordering_fraction_poset(parents, order):
    """Transitive closure with bitsets; returns the fraction of related pairs."""
    n = len(order)
    idx = {v: i for i, v in enumerate(order)}
    anc = [0] * n
    for v in order:
        i = idx[v]
        a = 0
        for p in parents.get(v, ()):
            if p in idx:
                j = idx[p]
                a |= anc[j] | (1 << j)
        anc[i] = a
    rel = sum(bin(a).count("1") for a in anc)
    return rel / (n * (n - 1) / 2), rel


def growth_poset(k, target, seed=1):
    """coloured_memory's own birth relation: w is created from v."""
    C = ColouredComplex(k, seed=seed)
    parents = {0: ()}
    born = [0]
    roots = [0]
    rng = C.rng
    while C.alive < target:
        if len(roots) > 3 * C.alive:
            roots = [r for r in roots if C.find(r) == r]
        v = roots[int(rng.integers(len(roots)))]
        if C.find(v) != v:
            continue
        a = int(rng.integers(C.k))
        b = int(rng.integers(C.k))
        if a == b:
            continue
        before = len(C.parent)
        C.square(v, a, b, 1.0)
        for w in range(before, len(C.parent)):
            parents[w] = (v,)
            born.append(w)
        roots.extend(range(before, len(C.parent)))
    keep = [v for v in born if C.find(v) == v]
    return parents, keep


def main():
    print("=" * 78)
    print("1. IS THE CONSISTENCY MOVE AN ORDER? (necessity.py's dynamics)")
    print("=" * 78)
    asym, total, deg = check_symmetry()
    print("  moves examined                      %d" % total)
    print("  moves with no reverse (asymmetric)  %d" % asym)
    print()
    if asym == 0:
        print("  Every move is reversible. The relation is SYMMETRIC, so")
        print("  reachability is an equivalence relation, not a partial order.")
        print("  Antisymmetry fails everywhere, and no relabelling fixes it.")
    print()
    live = deg[deg > 0].astype(float)
    pi = live / live.sum()
    print("  What the 'drift' actually is: a reversible walk on an undirected")
    print("  graph has stationary distribution proportional to degree, and")
    print("  satisfies detailed balance -- so the net current across every")
    print("  edge is exactly zero.")
    print("     states with moves        %d" % len(live))
    print("     degree range             %d to %d" % (live.min(), live.max()))
    print("     stationary entropy       %.4f of max %.4f"
          % (float(-(pi * np.log(pi)).sum()), float(np.log(len(live)))))
    print("     net probability current  0 by detailed balance")
    print()
    print("  A bias toward high-multiplicity states is not a direction in")
    print("  time. It is where a reversible walk spends its time.")
    print()

    print("=" * 78)
    print("2. CALIBRATION: ordering fraction of real Minkowski space")
    print("=" * 78)
    print("  Sprinkled into a causal diamond. Order only -- no coordinates")
    print("  are used in the statistic, just which pairs are related.")
    print()
    rng = np.random.default_rng(3)
    print("  %6s %8s %18s" % ("d", "n", "ordering fraction"))
    cal = {}
    for d in (2, 3, 4):
        vals = [ordering_fraction_minkowski(d, 600, rng) for _ in range(3)]
        cal[d] = float(np.mean(vals))
        print("  %6d %8d %18.4f" % (d, 600, cal[d]))
    print()
    print("  %6s %18.4f   (a chain: every pair related)" % ("1", 1.0))
    print()

    print("=" * 78)
    print("3. THE PROJECT'S OWN GROWTH POSET")
    print("=" * 78)
    print("  coloured_memory.py's birth relation: node w is created from v.")
    print("  Transitive closure gives a genuine partial order, built from the")
    print("  project's dynamics with no manifold assumed.")
    print()
    print("  %6s %8s %10s %18s %14s"
          % ("k", "nodes", "relations", "ordering fraction", "reads as"))
    for k in (2, 3, 4):
        parents, keep = growth_poset(k, 1200, seed=1)
        f, rel = ordering_fraction_poset(parents, keep)
        near = min(cal, key=lambda d: abs(cal[d] - f))
        tag = ("tree-like" if f < 0.5 * min(cal.values())
               else "d ~ %d" % near)
        print("  %6d %8d %10d %18.4f %14s" % (k, len(keep), rel, f, tag))
    print()
    print("  A tree has an ordering fraction that falls to zero with size,")
    print("  because almost no two nodes are ancestrally related. Minkowski")
    print("  space keeps a finite fraction. That is the discriminator.")
    print()


if __name__ == "__main__":
    main()
