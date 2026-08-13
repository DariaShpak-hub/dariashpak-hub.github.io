#!/usr/bin/env python3
"""
memory.py

How much does a state remember about how it began?

Quantified as the mutual information between "which initial condition was
prepared" and "where the system is now":

    I(seed ; X_t) = H( mean_i p_i ) - mean_i H( p_i )

with a uniform prior over seeds. I = ln k means the k seeds stay perfectly
distinguishable forever; I = 0 means the origin has been completely erased.

The prediction from the charge results is sharp: memory should equal exactly
the conserved quantities and nothing else. compose and decompose conserve the
per-component (nodes, weight) multiset, so seeds sharing that invariant must
become indistinguishable, while seeds differing in it can never mix. With
seeds spread over c distinct charge classes the prediction is therefore

    I(seed ; X_t)  ->  ln c        not ln k, and not 0.

Four seeds are used, two in each of two weight classes, so k = 4 and c = 2,
giving predicted I -> ln 2 = 0.6931.

Result, and its precise scope
-----------------------------
I decays from ln 4 to ln 2 = 0.693147, matching to 1.1e-16, with within-class
mutual information reaching exactly 0.000. Note this is computed over the
full 2^20 microstate distribution, not over chosen observables: the two
same-sector seeds converge to the *same distribution*, so no function of the
current state can separate them.

But that is a statement about ensembles, not trajectories, and in a finite
recurrent chain the two differ sharply. The seed here has pi = 2.7778e-3, so
its mean Poincare recurrence time is 1/pi = 360 steps, against an erasure
time of roughly 100 steps. The trajectory returns to its exact initial state
more often than the ensemble takes to forget it. Memory is erased from the
distribution while being revisited by the path.

This matters for cyclic arguments. A conserved charge confines a cycle to one
sector; it does not forbid the cycle. Recurrence guarantees the low-entropy
initial state recurs, and entropy genuinely falls when it does. The real
obstructions to a cyclic cosmology are the recurrence timescale, of order
e^S and therefore e^(10^104) for the observed universe, and the Boltzmann
brain problem: small fluctuations are exponentially more likely than large
ones, so a randomly chosen low-entropy observer should find a minimal
fluctuation rather than an ordered cosmos.
"""

import numpy as np
from itertools import combinations

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

A = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int16); T = np.zeros(N, np.int16)
for i in range(NB):
    B += ((A >> i) & 1).astype(np.int16)
for i in range(NB, 20):
    T += ((A >> i) & 1).astype(np.int16)
W = B + 2 * T


def moves():
    S, D = [], []
    for i, j, k in COMPOSE:
        m = (((A >> i) & 1) == 1) & (((A >> j) & 1) == 1) & (((A >> k) & 1) == 0)
        s = A[m]; S.append(s); D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
    for k, i, j in DECOMPOSE:
        m = (((A >> k) & 1) == 1) & (((A >> i) & 1) == 0) & (((A >> j) & 1) == 0)
        s = A[m]; S.append(s); D.append((s & ~(1 << k)) | (1 << i) | (1 << j))
    src, dst = np.concatenate(S), np.concatenate(D)
    od = np.bincount(src, minlength=N).astype(np.float64)
    return src, dst, od, od > 0


def bits(edges):
    m = 0
    for e in edges:
        m |= 1 << IDX[frozenset(e)]
    return m


def shannon(p):
    q = p[p > 0]
    return float(-(q * np.log(q)).sum())


def main():
    src, dst, od, live = moves()

    seeds = {
        "path  w=4": bits([(0, 1), (1, 2), (2, 3), (3, 4)]),
        "star  w=4": bits([(0, 1), (0, 2), (0, 3), (0, 4)]),
        "cycle w=5": bits([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]),
        "tri+e w=5": bits([(0, 1), (1, 2), (0, 2), (2, 3), (3, 4)]),
    }
    names = list(seeds)
    print("seeds and their conserved weight:")
    for n_ in names:
        print("   %-10s weight = %d" % (n_, W[seeds[n_]]))
    print()

    ps = []
    for n_ in names:
        p = np.zeros(N); p[seeds[n_]] = 1.0
        ps.append(p)

    print("%6s %10s %12s %12s" % ("t", "I(seed;X)", "ln k = ln4", "ln c = ln2"))
    for t in range(1, 601):
        for i in range(len(ps)):
            p = ps[i]
            w = np.zeros(N)
            np.divide(p, od, out=w, where=live)
            nxt = np.bincount(dst, weights=0.5 * w[src], minlength=N)
            nxt += 0.5 * p
            nxt[~live] += 0.5 * p[~live]
            ps[i] = nxt
        if t in (1, 2, 5, 10, 25, 50, 100, 200, 400, 600):
            mix = sum(ps) / len(ps)
            I = shannon(mix) - sum(shannon(p) for p in ps) / len(ps)
            print("%6d %10.6f %12.4f %12.4f" % (t, I, np.log(4), np.log(2)))

    mix = sum(ps) / len(ps)
    I = shannon(mix) - sum(shannon(p) for p in ps) / len(ps)
    print()
    print("final I = %.6f   ln 2 = %.6f   difference = %.2e"
          % (I, np.log(2), abs(I - np.log(2))))
    # within-class distinguishability
    for a, b in ((0, 1), (2, 3)):
        m2 = (ps[a] + ps[b]) / 2
        Iab = shannon(m2) - (shannon(ps[a]) + shannon(ps[b])) / 2
        print("   within-class I(%s vs %s) = %.3e"
              % (names[a], names[b], Iab))


if __name__ == "__main__":
    main()
