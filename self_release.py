#!/usr/bin/env python3
"""
self_release.py

A self-triggering constraint release: the system unlocks its own sector.

constraint_release.py opened the sector by switching the rule set by hand at
a chosen step. That is an external intervention, and a universe has no
exterior to perform it. Here the rule set is fixed for the entire run and the
escape rule carries a guard that is false in the prepared state:

    unlock:  if the state has NO binary edges, any binary edge may appear

compose and decompose conserve the weight b + 2t exactly, so the seed is
locked in its sector. The seed is all binary (b = 4, t = 0), so the guard
fails and the unlock is unavailable. To satisfy it the system must first
consume its own structure, composing every binary edge into ternaries until
b = 0. Only then does adding a binary edge become possible, which changes the
weight and breaks the conservation law that was confining it.

Nothing is switched. The ignition condition is reached, or not, by the
dynamics alone.

Because this is run on distributions rather than single trajectories, the
release is not a step but a leak: whatever probability sits on the b = 0
states escapes each step, so the escape rate is set by the equilibrium
occupation of the bottleneck. The expected signature is delayed ignition --
dS/dt decaying as the sector equilibrates, then turning back up as mass finds
the bottleneck, rather than the instantaneous jump a hand-thrown switch
gives.

Results
-------
Self-triggering works: escape begins at t=3 with no rule set ever switched,
against a control that stays locked in 290 states with P(escaped) = 0 for the
whole run.

Two things did not go as predicted.

There is no delayed ignition. The gate is not a rare fluctuation -- b = 0
states hold 1/12 of the equilibrium mass, and the seed reaches them in two
composes -- so escape starts almost immediately instead of after a waiting
time. Building a long fuse would need a gate condition of genuinely small
equilibrium weight.

The cascade does not continue; it terminates itself after exactly one escape,
by parity. The gate needs b = 0, which forces w = 2t and hence even weight,
while the escape adds a binary edge and raises the weight by one. So a system
starting at weight 4 escapes to weight 5, and weight 5 contains zero gate
states (against 45 at weight 4 and 120 at weight 6). The door opens once and
then locks behind it, permanently. Support settles at 1562 = 290 states of
weight 4 plus 1272 of weight 5, with 99.996% of the mass in the latter.

The entropy is not monotone. It rises to 7.283 while spread across both
sectors, then falls back through ln(1272) = 7.148 toward 7.08 as the mass
concentrates into the final sector alone -- an entropy pulse rather than a
step up. There is no contradiction with the H-theorem: nothing removes a
binary edge once b = 1, so the chain is irreversible, its stationary
distribution goes as degree rather than uniform, and the equilibrium entropy
of the final sector is below ln of its size.
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


def base_moves():
    a = np.arange(N, dtype=np.int64)
    S, D = [], []
    for i, j, k in COMPOSE:
        m = (((a >> i) & 1) == 1) & (((a >> j) & 1) == 1) & (((a >> k) & 1) == 0)
        s = a[m]
        S.append(s); D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
    for k, i, j in DECOMPOSE:
        m = (((a >> k) & 1) == 1) & (((a >> i) & 1) == 0) & (((a >> j) & 1) == 0)
        s = a[m]
        S.append(s); D.append((s & ~(1 << k)) | (1 << i) | (1 << j))
    return S, D


def unlock_moves():
    """Guarded escape: fires only on states with no binary edge."""
    a = np.arange(N, dtype=np.int64)
    gate = (a & BIN_MASK) == 0
    s0 = a[gate]
    S, D = [], []
    for i in range(NB):
        S.append(s0); D.append(s0 | (1 << i))
    return S, D


def build(with_unlock):
    S, D = base_moves()
    if with_unlock:
        s2, d2 = unlock_moves()
        S += s2; D += d2
    src, dst = np.concatenate(S), np.concatenate(D)
    outdeg = np.bincount(src, minlength=N).astype(np.float64)
    return src, dst, outdeg, outdeg > 0


idx = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int16)
T = np.zeros(N, np.int16)
for i in range(NB):
    B += ((idx >> i) & 1).astype(np.int16)
for i in range(NB, 20):
    T += ((idx >> i) & 1).astype(np.int16)
W = B + 2 * T
GATE = (idx & BIN_MASK) == 0


def shannon(p):
    q = p[p > 0]
    return float(-(q * np.log(q)).sum())


def run(label, with_unlock, steps=140):
    src, dst, od, live = build(with_unlock)

    seed = 0
    for e in ((0, 1), (1, 2), (2, 3), (3, 4)):
        seed |= 1 << IDX[frozenset(e)]
    p = np.zeros(N); p[seed] = 1.0
    w0 = W[seed]

    print(f"{label}   seed weight={w0}, b={B[seed]}, t={T[seed]}")
    print("%5s %10s %10s %10s %11s %11s" % (
        "t", "S(p)", "dS/dt", "|supp|", "P(gate)", "P(escaped)"))
    prev = shannon(p)
    for t in range(1, steps + 1):
        w = np.zeros(N)
        np.divide(p, od, out=w, where=live)
        nxt = np.bincount(dst, weights=0.5 * w[src], minlength=N)
        nxt += 0.5 * p
        nxt[~live] += 0.5 * p[~live]
        p = nxt
        s = shannon(p)
        if t <= 8 or t % 10 == 0:
            print("%5d %10.5f %10.6f %10d %11.3e %11.3e" % (
                t, s, s - prev, int((p > 1e-15).sum()),
                p[GATE].sum(), p[W != w0].sum()))
        prev = s
    print()


if __name__ == "__main__":
    run("CONTROL  compose+decompose only", with_unlock=False)
    run("SELF-RELEASE  fixed rule set, guarded unlock", with_unlock=True)
