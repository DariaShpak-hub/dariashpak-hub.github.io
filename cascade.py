#!/usr/bin/env python3
"""
cascade.py

Removing the parity lock: does the constraint release become a staircase?

self_release.py found that a guarded escape fires exactly once and then locks
permanently. The gate needs b = 0, which forces even weight, while the escape
adds a binary edge and makes the weight odd, so no gate state is ever
reachable again.

The obvious repair is to have the escape add a *ternary* edge instead, raising
the weight by two and preserving parity. That is variant A below. It has a
flaw worth testing rather than assuming: adding a ternary edge to a state with
b = 0 leaves b = 0, so the gate does not close behind the escape and the
system can climb repeatedly without ever doing work in between. That is a
runaway, not a staircase.

A genuine staircase needs the escape to close its own gate, so that reopening
it costs work. Variant B gates on t = 0 instead: the escape adds a ternary,
which immediately falsifies t = 0, and the system must decompose that ternary
back into binary edges to reopen. Each step then costs a decomposition, and
the climb should terminate when the binary slots run out rather than by parity.

    variant A   gate b == 0, escape adds a ternary   (parity preserved)
    variant B   gate t == 0, escape adds a ternary   (gate self-closing)
    variant V1  gate b == 0, escape adds a binary    (the parity-locked original)

Results
-------
The parity lock is removed and the cascade runs the whole ladder. Variant A
reaches every even weight from 4 to 20 -- the ceiling, all ten ternary slots
filled -- against V1's {4, 5}. Support grows from 1562 to 498586, mean weight
climbs to 9.97 and is still rising at t=300, and the termination is by
exhaustion of edge slots rather than by an arithmetic obstruction.

Two predictions in this file's own preamble were wrong.

A is not a runaway. The worry was that leaving b = 0 satisfied would let the
system climb instantly; instead the mean weight rises slowly and is nowhere
near the ceiling after 300 steps, because decompose keeps pulling mass back
off the gate.

Neither is it a staircase. All nine weight sectors are populated by t=25 and
the entropy rises smoothly with no steps. The system spreads across the ladder
rather than climbing it rung by rung, so "cascade" is right but "staircase"
is not.

Variant B, built specifically to produce steps by closing its own gate,
cascades *less* far rather than more: weights {6, 8, 10, 12} only. Gating on
t = 0 requires every ternary absent, which needs w binary edges, and with ten
binary slots that is impossible beyond weight 10 -- so it halts at 12. Making
the gate self-closing restricts the climb instead of structuring it. B also
drains its own seed weight entirely, losing weight 4 from the support by
t=175.

The entropy pulse of V1 does not recur in A or B within 300 steps. V1 peaked
and fell because it finished draining into a single terminal sector; A and B
keep opening new territory, so their entropy is still monotone rising when
the run ends.
"""

import numpy as np
from collections import Counter
from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = 10
N = 1 << 20
BIN_MASK = (1 << NB) - 1
TER_MASK = ((1 << 20) - 1) ^ BIN_MASK

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

IDXA = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int16)
T = np.zeros(N, np.int16)
for i in range(NB):
    B += ((IDXA >> i) & 1).astype(np.int16)
for i in range(NB, 20):
    T += ((IDXA >> i) & 1).astype(np.int16)
W = B + 2 * T


def base_moves():
    S, D = [], []
    for i, j, k in COMPOSE:
        m = (((IDXA >> i) & 1) == 1) & (((IDXA >> j) & 1) == 1) & \
            (((IDXA >> k) & 1) == 0)
        s = IDXA[m]
        S.append(s); D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
    for k, i, j in DECOMPOSE:
        m = (((IDXA >> k) & 1) == 1) & (((IDXA >> i) & 1) == 0) & \
            (((IDXA >> j) & 1) == 0)
        s = IDXA[m]
        S.append(s); D.append((s & ~(1 << k)) | (1 << i) | (1 << j))
    return S, D


def escape(gate_mask, add_range):
    """Guarded escape: on gated states, any absent edge in add_range appears."""
    s0 = IDXA[gate_mask]
    S, D = [], []
    for i in add_range:
        free = (s0 >> i) & 1 == 0
        s = s0[free]
        S.append(s); D.append(s | (1 << i))
    return S, D


def build(gate_mask, add_range):
    S, D = base_moves()
    s2, d2 = escape(gate_mask, add_range)
    S += s2; D += d2
    src, dst = np.concatenate(S), np.concatenate(D)
    od = np.bincount(src, minlength=N).astype(np.float64)
    return src, dst, od, od > 0


def shannon(p):
    q = p[p > 0]
    return float(-(q * np.log(q)).sum())


def run(label, gate_mask, add_range, steps=300):
    src, dst, od, live = build(gate_mask, add_range)
    seed = 0
    for e in ((0, 1), (1, 2), (2, 3), (3, 4)):
        seed |= 1 << IDX[frozenset(e)]
    p = np.zeros(N); p[seed] = 1.0

    print(f"{label}    seed w={W[seed]}, b={B[seed]}, t={T[seed]}")
    print("%5s %10s %11s %10s %10s %9s" % (
        "t", "S(p)", "dS/dt", "|supp|", "<weight>", "#weights"))
    prev = shannon(p)
    for k in range(1, steps + 1):
        w = np.zeros(N)
        np.divide(p, od, out=w, where=live)
        nxt = np.bincount(dst, weights=0.5 * w[src], minlength=N)
        nxt += 0.5 * p
        nxt[~live] += 0.5 * p[~live]
        p = nxt
        s = shannon(p)
        if k <= 4 or k % 25 == 0:
            supp = p > 1e-15
            print("%5d %10.5f %11.6f %10d %10.4f %9d" % (
                k, s, s - prev, int(supp.sum()),
                float((p * W).sum()), len(set(W[supp].tolist()))))
        prev = s
    supp = p > 1e-15
    wc = Counter(W[supp].tolist())
    print("  weights reached:", sorted(wc))
    print("  mass by weight :", {int(k): round(float(p[(W == k) & supp].sum()), 4)
                                 for k in sorted(wc)})
    print()


if __name__ == "__main__":
    run("V1  gate b==0, escape adds BINARY  (parity-locked)",
        (IDXA & BIN_MASK) == 0, range(0, NB))
    run("A   gate b==0, escape adds TERNARY (parity preserved)",
        (IDXA & BIN_MASK) == 0, range(NB, 20))
    run("B   gate t==0, escape adds TERNARY (gate self-closing)",
        (IDXA & TER_MASK) == 0, range(NB, 20))
