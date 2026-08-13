#!/usr/bin/env python3
"""
constraint_release.py

A constraint-release transition: does opening a conserved sector produce a
burst of entropy production?

The setup is the thermodynamic content of a backdraft, with no external
reservoir needed. compose and decompose conserve the weight b + 2t exactly,
so a prepared state is locked inside a small sector of the 2^20 configuration
space and equilibrates there. Switching on a rule that changes the weight --
here a single-edge toggle -- removes the conservation law and makes the rest
of the space reachable. Nothing is added to the system: the same states
existed all along, they were simply unreachable.

Measured across the switch:

    S(p_t)      Shannon entropy of the distribution
    |supp|      number of states carrying weight
    dS/dt       entropy production rate
    S(q)        entropy of the charge distribution q = weight mod 3

The control never opens the constraint, which separates a genuine release
burst from ordinary continued relaxation.

Result
------
The transition is sharp and the control isolates it. Both runs are identical
up to t=30: support locks at 290 states by t=12, entropy saturates near 5.38,
and dS/dt decays smoothly to 0.024. The control continues that decay to
0.00016 by t=90 with the support never leaving 290 states.

Releasing the constraint at t=30 makes dS/dt jump from 0.024 to 1.188 in a
single step, a factor of 50, with the support going 290 -> 1622 immediately
and reaching all 1048576 states by t=84. The charge entropy rises from 0 to
exactly ln 3 = 1.09861, so the weight-mod-3 charge fully randomises across
its three sectors. Sixty steps after release the production rate is still
0.029 against the control's 0.00016, a factor of 184.

Caveats
-------
That entropy rises when a larger region becomes reachable is not itself a
discovery; the content is the sharpness, and the control showing it is the
constraint change rather than continued relaxation. The size of the jump is
also not intrinsic -- releasing later means a lower pre-release rate and so a
larger ratio.

The important limitation is that the door is opened from outside. The rule
set is switched by hand at t=30, which is an external intervention, and a
universe has no exterior to open it. A self-triggering version would need the
rule's applicability to depend on the state, so the system unlocks its own
sector. That is the difference between a backdraft and a genuine
cosmological transition, and it is a well-posed next construction rather
than a fixable detail of this one.
"""

import numpy as np
from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = 10
N = 1 << 20

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


def move_arrays(toggles=()):
    src_all = np.arange(N, dtype=np.int64)
    S, D = [], []
    for i, j, k in COMPOSE:
        m = (((src_all >> i) & 1) == 1) & (((src_all >> j) & 1) == 1) & \
            (((src_all >> k) & 1) == 0)
        s = src_all[m]
        S.append(s); D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
    for k, i, j in DECOMPOSE:
        m = (((src_all >> k) & 1) == 1) & (((src_all >> i) & 1) == 0) & \
            (((src_all >> j) & 1) == 0)
        s = src_all[m]
        S.append(s); D.append((s & ~(1 << k)) | (1 << i) | (1 << j))
    for i in toggles:
        S.append(src_all); D.append(src_all ^ (1 << i))
    return np.concatenate(S), np.concatenate(D)


def weights():
    idx = np.arange(N, dtype=np.int64)
    b = np.zeros(N, np.int16); t = np.zeros(N, np.int16)
    for i in range(NB):
        b += ((idx >> i) & 1).astype(np.int16)
    for i in range(NB, 20):
        t += ((idx >> i) & 1).astype(np.int16)
    return b + 2 * t


W = weights()
Q = (W % 3).astype(np.int8)


def step(p, src, dst, outdeg, live):
    """One lazy step: stay with probability 1/2, else uniform successor."""
    w = np.zeros(N)
    np.divide(p, outdeg, out=w, where=live)
    nxt = np.bincount(dst, weights=0.5 * w[src], minlength=N)
    nxt += 0.5 * p
    nxt[~live] += 0.5 * p[~live]
    return nxt


def shannon(p):
    q = p[p > 0]
    return float(-(q * np.log(q)).sum())


def charge_entropy(p):
    c = np.array([p[Q == k].sum() for k in range(3)])
    c = c[c > 0]
    return float(-(c * np.log(c)).sum())


def run(label, release_at, steps=90):
    closed = move_arrays()
    opened = move_arrays(tuple(range(NB)))
    od_c = np.bincount(closed[0], minlength=N).astype(np.float64)
    od_o = np.bincount(opened[0], minlength=N).astype(np.float64)

    seed = 0
    for e in ((0, 1), (1, 2), (2, 3), (3, 4)):
        seed |= 1 << IDX[frozenset(e)]
    p = np.zeros(N); p[seed] = 1.0

    print(f"{label}   (seed weight = {W[seed]}, charge q = {Q[seed]})")
    print("%5s %10s %10s %10s %10s %8s" % (
        "t", "S(p)", "dS/dt", "|supp|", "S(q)", "phase"))
    prev = shannon(p)
    for t in range(1, steps + 1):
        if release_at is not None and t > release_at:
            src, dst, od = opened[0], opened[1], od_o
            phase = "open"
        else:
            src, dst, od = closed[0], closed[1], od_c
            phase = "closed"
        p = step(p, src, dst, od, od > 0)
        s = shannon(p)
        if t <= 6 or t % 6 == 0 or (release_at and release_at <= t <= release_at + 6):
            print("%5d %10.5f %10.5f %10d %10.5f %8s" % (
                t, s, s - prev, int((p > 1e-15).sum()), charge_entropy(p), phase))
        prev = s
    print()


if __name__ == "__main__":
    run("CONTROL: constraint never released", release_at=None)
    run("RELEASE at t=30", release_at=30)
