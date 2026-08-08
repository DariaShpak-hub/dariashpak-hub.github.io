#!/usr/bin/env python3
"""
null_model.py

Are the structural results special to compose/decompose, or generic?

Every structural claim in this project was measured on one hand-picked rule
set. Without a matched random ensemble there is no way to tell "compose and
decompose conserve a Z_3 charge" from "any rule set of this shape conserves
something". This builds the null.

The null ensemble
-----------------
compose consumes two binary edges sharing exactly one node and produces the
ternary edge equal to their union, giving a table of triples (i, j, k) with
i, j binary and k ternary. decompose is the same table run backwards.

A null rule set keeps the *arity signature* -- same number of triples, same
two-binary-to-one-ternary shape, forwards and backwards -- but randomises the
*geometry*, dropping the requirement that k be the union of i and j. Anything
the real rule set shares with the null is forced by arity; anything it does
not share is contributed by the geometric structure.

Prediction, stated before running: the Z_3 charge should be generic. Any rule
of this arity changes (b, t) by (-2, +1), so b - t changes by exactly 3
whatever edges are involved, and (b-t) mod 3 is conserved for structural
reasons that have nothing to do with compose. Hypergraph connectivity should
be special, since only the union condition keeps a rewrite inside a connected
component.

Result: the prediction holds, and it splits the charge result in two
------------------------------------------------------------------
    weight          real: yes    null: 40/40    generic
    weight mod 3    real: yes    null: 40/40    generic
    edge count      real: no     null:  0/40    (sanity check)
    connectivity    real: yes    null:  0/40    SPECIAL

Weight conservation is pure arity: two binary edges carry weight 2 and one
ternary edge carries weight 2, so any rule of this shape conserves weight
whatever edges it touches. The Z_3 charge is then not an independent finding
at all for this rule set -- it follows from exact weight conservation, since
q = w mod 3. Both were reported earlier as though they were structural
properties of compose and decompose. They are not.

(Z_3 was a genuine independent charge in sector_breaking.py, where the
sigma-symmetrised rule set breaks weight conservation and only the residue
survives. That reading stands. What does not stand is treating it as special
to compose/decompose in the unsymmetrised case.)

Connectivity is the real structural result: no random rule set of matched
arity conserves it, because only the union condition -- k equal to i or j --
keeps a rewrite inside a single connected component. The per-component
(nodes, weight) charge from charges.py therefore survives the null in its
per-component half and dissolves into arity bookkeeping in its weight half.
"""

import numpy as np
from itertools import combinations

N_NODES = 5
EDGES = [frozenset(c) for c in combinations(range(N_NODES), 2)] + \
        [frozenset(c) for c in combinations(range(N_NODES), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = sum(1 for e in EDGES if len(e) == 2)
NE = len(EDGES)
N = 1 << NE

REAL = [(i, j, IDX[EDGES[i] | EDGES[j]])
        for i in range(NB) for j in range(i + 1, NB)
        if len(EDGES[i] & EDGES[j]) == 1]

A = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int8); T = np.zeros(N, np.int8)
for i in range(NB):
    B += ((A >> i) & 1).astype(np.int8)
for i in range(NB, NE):
    T += ((A >> i) & 1).astype(np.int8)
W = (B.astype(np.int16) + 2 * T.astype(np.int16))


def components_all():
    """Number of connected components of the hypergraph, for every state."""
    node_mask = [0] * NE
    for i, e in enumerate(EDGES):
        m = 0
        for v in e:
            m |= 1 << v
        node_mask[i] = m

    reach = [np.full(N, 1 << v, dtype=np.int8) for v in range(N_NODES)]
    for _ in range(N_NODES):
        for v in range(N_NODES):
            r = reach[v]
            for i in range(NE):
                present = ((A >> i) & 1).astype(bool)
                touches = (r & node_mask[i]) != 0
                r = np.where(present & touches, r | node_mask[i], r)
            reach[v] = r

    comp = np.zeros(N, np.int8)
    for v in range(N_NODES):
        r = reach[v].astype(np.int32)
        lowest = r & (-r)                       # lowest set bit
        comp += (lowest == (1 << v)).astype(np.int8)
    return comp


COMP = components_all()


def transitions(triples):
    """Forward (2 binary -> 1 ternary) and backward moves for a triple table."""
    S, D = [], []
    for i, j, k in triples:
        m = (((A >> i) & 1) == 1) & (((A >> j) & 1) == 1) & (((A >> k) & 1) == 0)
        s = A[m]
        S.append(s); D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))
        m2 = (((A >> k) & 1) == 1) & (((A >> i) & 1) == 0) & (((A >> j) & 1) == 0)
        s2 = A[m2]
        S.append(s2); D.append((s2 & ~(1 << k)) | (1 << i) | (1 << j))
    return np.concatenate(S), np.concatenate(D)


def random_triples(rng, count):
    out = set()
    while len(out) < count:
        i, j = rng.choice(NB, 2, replace=False)
        k = int(rng.integers(NB, NE))
        out.add((int(min(i, j)), int(max(i, j)), k))
    return sorted(out)


def invariants(src, dst):
    """Which candidate quantities are conserved by every transition?"""
    return {
        "weight": bool((W[src] == W[dst]).all()),
        "weight mod 3": bool((W[src] % 3 == W[dst] % 3).all()),
        "edge count": bool(((B[src] + T[src]) == (B[dst] + T[dst])).all()),
        "connectivity": bool((COMP[src] == COMP[dst]).all()),
        "node coverage": bool(((COMP[src] == COMP[dst])).all()),
    }


def main():
    rng = np.random.default_rng(20260808)
    keys = ["weight", "weight mod 3", "edge count", "connectivity"]

    src, dst = transitions(REAL)
    real = invariants(src, dst)
    print("real rule set (compose + decompose), %d triples, %d moves"
          % (len(REAL), len(src)))
    for k in keys:
        print("   %-14s conserved: %s" % (k, real[k]))
    print()

    K = 40
    print("null ensemble: %d random rule sets, %d triples each" % (K, len(REAL)))
    tally = {k: 0 for k in keys}
    for _ in range(K):
        tr = random_triples(rng, len(REAL))
        s, d = transitions(tr)
        inv = invariants(s, d)
        for k in keys:
            tally[k] += int(inv[k])
    print("%-14s %10s %14s %s" % ("invariant", "real", "null (of %d)" % K,
                                  "verdict"))
    for k in keys:
        frac = tally[k] / K
        if real[k] and frac > 0.9:
            verdict = "GENERIC - forced by arity"
        elif real[k] and frac < 0.1:
            verdict = "SPECIAL - geometry"
        elif real[k]:
            verdict = "partly generic"
        else:
            verdict = "not conserved even for real"
        print("%-14s %10s %14d %s" % (k, real[k], tally[k], verdict))


if __name__ == "__main__":
    main()
