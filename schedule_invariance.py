#!/usr/bin/env python3
"""
schedule_invariance.py

Which observables survive a change of update schedule?

The rewrite system has no time. Step counts are path lengths along a
traversal, not a physical coordinate, so any quantity defined per step --
H = dS/dt, w_eff, "the peak is at t=3696" -- may be an artifact of how the
simulation was scheduled rather than a fact about the system.

The test: keep the transition structure X_i -> X_j fixed and vary only the
schedule, then see what is unchanged. Two kinds of change are applied.

    laziness      p <- (1-a) p + a P p, for several a.  A pure slowing-down:
                  the same kernel, iterated differently.
    rate bias     compose, decompose and escape fire at different relative
                  rates.  This reweights which transitions are taken, not
                  just how fast.

Observables are compared *not* against the step counter but against an
intrinsic quantity, the mean weight <W>, which is a function of state. If the
S-versus-<W> curves from different schedules collapse onto one another, the
relation is structural. If they separate, it belonged to the schedule.

Expectation, stated in advance: laziness should collapse, since it only
rescales how far the distribution moves per step. Rate bias changes which
paths carry the mass and need not collapse. Purely structural quantities --
the reachable set, the sectors, the set of weights reached -- must be
identical under every schedule, and serve as the control.

Result: the expectation holds, and it separates the observables cleanly.

Laziness collapses. lazy=0.5 and lazy=0.8 agree on S at matched <W> to about
0.005 nats across the whole grid, so the S-versus-<W> relation does not care
how far the distribution is advanced per step.

lazy=0.0 sits below both by 0.691 to 0.698 nats, which is ln 2 = 0.69315.
That is not a schedule dependence but the bipartite-periodicity artifact:
compose and decompose each flip the parity of the edge count, so the
non-lazy walk occupies one parity class at a time, halving its effective
support and costing exactly ln 2. The same artifact appeared earlier in this
project as a non-converging power iteration.

Rate bias does not collapse. Changing the relative firing rates of compose,
decompose and escape moves S at matched <W> by up to 2.245 nats, far outside
the laziness spread. Notably the disagreement shrinks as <W> grows -- 2.25 at
<W>=5 down to 0.18 at <W>=10 -- so the schedule dependence is concentrated in
the early cascade and washes out as the system approaches its ceiling.

The consequence for the cosmological reading is direct. Step counts and
anything defined per step were already suspect; this shows the S-versus-<W>
relation is robust to *reparameterisation* of the schedule but not to
*reweighting* of it. An emergent-time construction must therefore fix the
relative rates of the rules, not merely the step size, before any expansion
history read off it can be called a property of the system.
"""

import numpy as np
from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = 10
N = 1 << 20
BIN_MASK = (1 << NB) - 1

COMPOSE = [(i, j, IDX[EDGES[i] | EDGES[j]])
           for i in range(NB) for j in range(i + 1, NB)
           if len(EDGES[i] & EDGES[j]) == 1]
DECOMPOSE = []
for k in range(NB, 20):
    e = EDGES[k]
    for mid in sorted(e):
        o = sorted(e - {mid})
        DECOMPOSE.append((k, IDX[frozenset((mid, o[0]))],
                          IDX[frozenset((mid, o[1]))]))

A = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int16); T = np.zeros(N, np.int16)
for i in range(NB):
    B += ((A >> i) & 1).astype(np.int16)
for i in range(NB, 20):
    T += ((A >> i) & 1).astype(np.int16)
W = (B + 2 * T).astype(np.float64)


def channels():
    """Separate (src,dst) arrays per rule type, so rates can be varied."""
    cs, cd = [], []
    for i, j, k in COMPOSE:
        m = (((A >> i) & 1) == 1) & (((A >> j) & 1) == 1) & (((A >> k) & 1) == 0)
        s = A[m]; cs.append(s); cd.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
    ds, dd = [], []
    for k, i, j in DECOMPOSE:
        m = (((A >> k) & 1) == 1) & (((A >> i) & 1) == 0) & (((A >> j) & 1) == 0)
        s = A[m]; ds.append(s); dd.append((s & ~(1 << k)) | (1 << i) | (1 << j))
    gate = (A & BIN_MASK) == 0
    s0 = A[gate]
    es, ed = [], []
    for i in range(NB, 20):
        free = ((s0 >> i) & 1) == 0
        s = s0[free]; es.append(s); ed.append(s | (1 << i))
    return ((np.concatenate(cs), np.concatenate(cd)),
            (np.concatenate(ds), np.concatenate(dd)),
            (np.concatenate(es), np.concatenate(ed)))


CH = channels()


def kernel(rates):
    """Row-normalised weighted kernel: rate r_c, r_d, r_e per channel."""
    src = np.concatenate([CH[i][0] for i in range(3)])
    dst = np.concatenate([CH[i][1] for i in range(3)])
    wt = np.concatenate([np.full(len(CH[i][0]), rates[i], dtype=np.float64)
                         for i in range(3)])
    tot = np.bincount(src, weights=wt, minlength=N)
    live = tot > 0
    return src, dst, wt, tot, live


def shannon(p):
    q = p[p > 0]
    return float(-(q * np.log(q)).sum())


def trajectory(rates, lazy, steps):
    src, dst, wt, tot, live = kernel(rates)
    seed = 0
    for e in ((0, 1), (1, 2), (2, 3), (3, 4)):
        seed |= 1 << IDX[frozenset(e)]
    p = np.zeros(N); p[seed] = 1.0
    out = []
    share = np.zeros(N)
    for _ in range(steps):
        np.divide(p, tot, out=share, where=live)
        moved = np.bincount(dst, weights=wt * share[src], minlength=N)
        p = (1 - lazy) * moved + lazy * p + (1 - lazy) * np.where(live, 0.0, p)
        p /= p.sum()
        out.append((float((p * W).sum()), shannon(p)))
    return np.array(out)


def collapse(curves, labels, grid):
    """Interpolate S onto a common <W> grid and report the spread."""
    print("%9s" % "<W>", *["%11s" % l for l in labels], "%11s" % "spread")
    rows = []
    for g in grid:
        vals = []
        for c in curves:
            w, s = c[:, 0], c[:, 1]
            if w[0] <= g <= w[-1]:
                vals.append(float(np.interp(g, w, s)))
            else:
                vals.append(float("nan"))
        good = [v for v in vals if v == v]
        spread = (max(good) - min(good)) if len(good) > 1 else float("nan")
        rows.append(spread)
        print("%9.2f" % g, *["%11.5f" % v for v in vals], "%11.5f" % spread)
    ok = [r for r in rows if r == r]
    print("  max spread over grid: %.5f nats" % (max(ok) if ok else float("nan")))
    print()


if __name__ == "__main__":
    grid = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]

    print("=" * 78)
    print("A. laziness only  (same kernel, different step size)")
    print("=" * 78)
    lazies = [0.0, 0.5, 0.8]
    curves = [trajectory((1, 1, 1), a, 2500) for a in lazies]
    collapse(curves, ["lazy=%.1f" % a for a in lazies], grid)

    print("=" * 78)
    print("B. rate bias  (same transitions, different relative rates)")
    print("=" * 78)
    rate_sets = [(1, 1, 1), (3, 1, 1), (1, 3, 1), (1, 1, 5)]
    curves = [trajectory(r, 0.5, 2500) for r in rate_sets]
    collapse(curves, ["c%d d%d e%d" % r for r in rate_sets], grid)
