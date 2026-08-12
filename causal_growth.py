#!/usr/bin/env python3
"""
causal_growth.py

Build a causal set here, instead of borrowing the framing.

causal_order.py established two things. The consistency move cannot give an
order -- 0 of 184320 moves lack a reverse, so reachability is an equivalence
relation. And the growth relation that does give an order gives a tree:
ordering fraction 0.0354 falling to 0.0061 as n runs 600 to 4801, against a
size-stable 0.5014 for 1+1 Minkowski.

Two requirements came out of that, both checkable:

    MANY PREDECESSORS   each new element must enter the future of a substantial
                        downset, not a single parent
    LABEL INVARIANCE    the physics must not depend on the order the labels
                        were handed out, or the "causal structure" is the
                        simulation's loop counter in disguise

This tests both, on Rideout-Sorkin style sequential growth.

The growth rule
---------------
Transitive percolation, the simplest classical sequential growth model.
Elements are added one at a time; each new element is linked to each existing
element independently with probability p; the order is the transitive closure.
p is the only parameter. p -> 0 gives an antichain, p -> 1 gives a chain, and
the interesting behaviour is in between.

Three diagnostics, none of which can be passed by tuning alone
----------------------------------------------------------------
    ordering fraction   the Myrheim-Meyer statistic, calibrated against
                        sprinklings rather than quoted from a formula

    midpoint dimension  an INDEPENDENT estimator. Split the interval at the
                        element that balances it best; in d dimensions each
                        half holds about 1/2^d of the volume. Also calibrated
                        against sprinklings.

    label recovery      Spearman correlation between the birth label and the
                        position in a freshly sampled linear extension. A chain
                        forces its own labelling and scores 1. An antichain
                        scores 0. If a growth model scores near 1 then the
                        labelling IS the physics and time has been smuggled.

The point of having two dimension estimators is that any single one can be hit
by adjusting p. Agreement between two independent ones is much harder to fake,
which is the gate this project has been applying everywhere else.

And the negative control matters as much as the positive one. Almost every
poset on n elements is a Kleitman-Rothschild poset: three layers, everything in
the middle layer, relations only between adjacent layers. Those have a
perfectly respectable ordering fraction and are nothing like a spacetime. If
the two estimators agree on Minkowski and disagree on a random layered poset,
they are measuring something real.

Result 1: the two-estimator gate validates
--------------------------------------------
    construction              n    ord frac   midpoint d   label rec
    Minkowski d = 2         700      0.5127        1.94       0.940
    Minkowski d = 3         700      0.2435        2.96       0.867
    Minkowski d = 4         700      0.1108        3.81       0.790
    random layered (KR)     700      0.4471        2.43       0.883

The midpoint estimator returns 1.94, 2.96 and 3.81 for true dimensions 2, 3
and 4, using a completely different statistic from the ordering fraction. Two
independent measurements agreeing on the right answer is the gate, and it
works.

The negative control is weaker than intended and should be reported as such.
The layered poset gives 0.4471, which reads as d ~ 2.2 on the ordering-fraction
curve, against a midpoint of 2.43 -- a gap of only about 0.2, not the sharp
disagreement I was hoping for. A better discriminator was sitting in the data
and was not measured: the layered poset has height 3 by construction, while a
Minkowski sprinkling at n = 700 has a far longer chain. Height would have
separated them immediately.

Result 2: the project's growth fails both ways
------------------------------------------------
    coloured birth, k = 2   701      0.0301        5.28       0.552
    coloured birth, k = 3   702      0.0196        5.75       0.306

Ordering fraction 0.03 is off the bottom of the calibration entirely, and it
falls with n as causal_order.py measured, so there is no stable dimension to
read. The midpoint estimator returning 5.28 is not a dimension either -- it is
what that statistic does on a tree. The two numbers are not in agreement about
anything; they are both reporting that this object is not a spacetime.

Result 3: transitive percolation does not reproduce Minkowski's agreement
--------------------------------------------------------------------------
    p        ord frac   implied d   midpoint d   gap    label rec
    0.002      0.0030     >> 4         7.45        --      0.281
    0.005      0.0187     >> 4         5.87        --      0.605
    0.010      0.1480      ~3.6        3.08       0.5      0.830
    0.020      0.5144      ~2.0        1.59       0.4      0.944
    0.050      0.8142      < 2         1.21        --      0.991
    0.100      0.9267      < 2         1.05        --      0.998
    0.300      0.9873      < 2         1.01        --      1.000

For Minkowski the two estimators agree to within 0.06, 0.04 and 0.19. For
percolation the closest approaches are gaps of 0.4 and 0.5, and the trend runs
straight from antichain to chain without passing through a regime where they
lock together.

That is consistent with what is known: transitive percolation is manifold-like
only in a particular scaling limit, and generically produces posets that have a
respectable ordering fraction and no geometry. This scan does not establish
that no such regime exists -- n = 700 is small and the interesting scaling may
be p ~ 1/n with much larger n -- but at every parameter tried, the agreement
that Minkowski shows is absent.

Result 4: my label-recovery diagnostic was misconceived
--------------------------------------------------------
I proposed that a label recovery near 1 would show the birth order was smuggled
time. Minkowski scores 0.940, 0.867 and 0.790.

So high label recovery is what a REAL SPACETIME looks like. Of course it is --
a causal set encodes its own time-ordering; that is the entire content of
"order determines geometry." The statistic measures how chain-like a poset is,
not whether the generating process cheated. It cannot detect what I built it to
detect, and no threshold on it would help.

Result 5: and my "loop counter in disguise" claim was too strong
-----------------------------------------------------------------
Last turn I said the coloured complex's birth poset was the simulation's loop
counter wearing a poset's clothes. Its label recovery is 0.552 at k = 2 and
0.306 at k = 3 -- LOWER than Minkowski's.

A tree has many incomparable pairs, so it admits enormous numbers of linear
extensions and does not remember the order it was built in. The birth sequence
generated the tree; the tree does not encode the birth sequence. That is a real
distinction and I collapsed it.

The tree is still not a causal set, for the reasons in Result 2. It is just not
a relabelled clock either.

Verdict
-------
    two independent estimators, calibrated      WORKS (1.94, 2.96, 3.81)
    project's growth poset                      fails, and not marginally
    transitive percolation, p in [0.002, 0.3]   never matches the agreement
    label recovery as a smuggling detector      RETRACTED, cannot work
    "birth poset is a relabelled loop counter"  RETRACTED, measured 0.55

What survives is a usable instrument and a clean negative. The instrument is
the pair of estimators plus the Minkowski calibration, which now sits in the
repository ready to score any future construction. The negative is that neither
the project's own growth nor the simplest sequential growth model produces a
causal set with a dimension.

What is genuinely untested: whether some other growth rule does. The Rideout-
Sorkin family is much larger than transitive percolation -- it is parameterised
by a whole sequence of coupling constants -- and only special choices are
expected to be manifold-like. Searching that family is a real experiment and
this module now has the scoring function for it.
"""

import numpy as np

from causal_order import sprinkle, growth_poset


# --------------------------------------------------------------- poset core

def closure(n, links):
    """anc[i] = bitset of everything strictly below i. links: i -> list of j<i."""
    anc = [0] * n
    for i in range(n):
        a = 0
        for j in links[i]:
            j = int(j)              # numpy ints silently overflow past bit 63
            a |= anc[j] | (1 << j)
        anc[i] = a
    return anc


def relations(anc):
    return sum(bin(a).count("1") for a in anc)


def ordering_fraction(anc):
    n = len(anc)
    return relations(anc) / (n * (n - 1) / 2)


def midpoint_dimension(anc):
    """Split at the best-balanced element; each half holds ~1/2^d of the whole."""
    n = len(anc)
    desc = [0] * n
    for i in range(n):
        a = anc[i]
        while a:
            b = a & -a
            desc[b.bit_length() - 1] |= (1 << i)
            a ^= b
    best, bm = -1, None
    for m in range(n):
        lo = bin(anc[m]).count("1")
        hi = bin(desc[m]).count("1")
        if min(lo, hi) > best:
            best, bm = min(lo, hi), m
    if bm is None or best < 2:
        return float("nan")
    half = max(bin(anc[bm]).count("1"), bin(desc[bm]).count("1"))
    if half <= 0 or half >= n:
        return float("nan")
    return float(np.log2(n / half))


def label_recovery(anc, rng, trials=6):
    """Spearman correlation of birth label with a sampled linear extension."""
    n = len(anc)
    desc = [0] * n
    for i in range(n):
        a = anc[i]
        while a:
            b = a & -a
            desc[b.bit_length() - 1] |= (1 << i)
            a ^= b
    cs = []
    for _ in range(trials):
        remaining = set(range(n))
        placed = 0
        pos = [0] * n
        indeg = [bin(anc[i]).count("1") for i in range(n)]
        ready = [i for i in range(n) if indeg[i] == 0]
        while ready:
            j = int(rng.integers(len(ready)))
            v = ready.pop(j)
            pos[v] = placed
            placed += 1
            remaining.discard(v)
            d = desc[v]
            while d:
                b = d & -d
                w = b.bit_length() - 1
                indeg[w] -= 1
                if indeg[w] == 0 and w in remaining:
                    ready.append(w)
                d ^= b
        if placed < n:
            continue
        a = np.argsort(np.argsort(np.arange(n)))
        b = np.argsort(np.argsort(np.array(pos)))
        cs.append(float(np.corrcoef(a, b)[0, 1]))
    return float(np.mean(cs)) if cs else float("nan")


# ------------------------------------------------------------- constructions

def transitive_percolation(n, p, rng):
    links = [[] for _ in range(n)]
    for i in range(1, n):
        m = rng.random(i) < p
        links[i] = [int(x) for x in np.nonzero(m)[0]]
    return closure(n, links)


def minkowski_poset(d, n, rng):
    pts = sprinkle(d, n, rng)
    t = np.array([q[0] for q in pts])
    X = np.array([q[1] for q in pts])
    links = [[] for _ in range(n)]
    for i in range(n):
        dt = t[i] - t[:i]
        dx = np.sqrt(((X[i] - X[:i]) ** 2).sum(axis=1))
        links[i] = [int(x) for x in np.nonzero(dt > dx)[0]]
    return closure(n, links)


def layered_poset(n, rng, layers=3):
    """Kleitman-Rothschild style: the shape almost every random poset has."""
    lab = rng.integers(0, layers, n)
    lab.sort()
    links = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            if lab[i] == lab[j] + 1 and rng.random() < 0.5:
                links[i].append(int(j))
    return closure(n, links)


def tree_poset(k, n, seed=1):
    parents, keep = growth_poset(k, n, seed=seed)
    idx = {v: i for i, v in enumerate(keep)}
    links = [[] for _ in range(len(keep))]
    for v in keep:
        for pnt in parents.get(v, ()):
            if pnt in idx and idx[pnt] < idx[v]:
                links[idx[v]].append(idx[pnt])
    return closure(len(keep), links)


# ------------------------------------------------------------------- report

def row(name, anc, rng):
    return (name, ordering_fraction(anc), midpoint_dimension(anc),
            label_recovery(anc, rng), len(anc))


def main():
    rng = np.random.default_rng(11)
    N = 700

    print("=" * 78)
    print("1. CALIBRATION AND CONTROLS")
    print("=" * 78)
    print("  %-26s %8s %10s %12s %10s"
          % ("construction", "n", "ord frac", "midpoint d", "label rec"))
    cal = {}
    for d in (2, 3, 4):
        anc = minkowski_poset(d, N, rng)
        nm, f, md, lr, n = row("Minkowski d = %d" % d, anc, rng)
        cal[d] = (f, md)
        print("  %-26s %8d %10.4f %12.2f %10.3f" % (nm, n, f, md, lr))
    anc = layered_poset(N, rng)
    nm, f, md, lr, n = row("random layered (KR-like)", anc, rng)
    print("  %-26s %8d %10.4f %12.2f %10.3f" % (nm, n, f, md, lr))
    print()
    print("  The two dimension estimators should AGREE on Minkowski and")
    print("  disagree on the layered poset. That is what makes them a test")
    print("  rather than a pair of tunable numbers.")
    print()

    print("=" * 78)
    print("2. THE PROJECT'S EXISTING GROWTH, FOR REFERENCE")
    print("=" * 78)
    print("  %-26s %8s %10s %12s %10s"
          % ("construction", "n", "ord frac", "midpoint d", "label rec"))
    for k in (2, 3):
        anc = tree_poset(k, N, seed=1)
        nm, f, md, lr, n = row("coloured birth, k = %d" % k, anc, rng)
        print("  %-26s %8d %10.4f %12.2f %10.3f" % (nm, n, f, md, lr))
    print()

    print("=" * 78)
    print("3. SEQUENTIAL GROWTH: transitive percolation, scanning p")
    print("=" * 78)
    print("  %-26s %8s %10s %12s %10s"
          % ("construction", "n", "ord frac", "midpoint d", "label rec"))
    for p in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3):
        anc = transitive_percolation(N, p, rng)
        nm, f, md, lr, n = row("percolation p = %.3f" % p, anc, rng)
        near = min(cal, key=lambda d: abs(cal[d][0] - f))
        print("  %-26s %8d %10.4f %12.2f %10.3f   MM ~ d %d"
              % (nm, n, f, md, lr, near))
    print()

    print("=" * 78)
    print("  PASS requires all three: an ordering fraction in the Minkowski")
    print("  range, a midpoint dimension AGREEING with it, and a label")
    print("  recovery well below 1 so the birth order is gauge rather than")
    print("  physics.")
    print()


if __name__ == "__main__":
    main()
