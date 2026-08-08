#!/usr/bin/env python3
"""
conjugation.py

Can a charge-symmetric rule set produce a motif/anti-motif population
asymmetry on its own?

The answer is fixed before any measurement, and the argument is short. Let
sigma be an involution on states (the conjugation) and suppose the rule set is
sigma-equivariant, meaning H -> G exactly when sigma H -> sigma G. Then sigma
is an automorphism of the continuation graph, any equivariant transition
kernel has pi(sigma H) = pi(H), and for every motif class M

    A_M = [P(M) - P(sigma M)] / [P(M) + P(sigma M)] = 0

identically. This is the same statement that forced the probability flux to
vanish on the reversible rule set: a symmetry of the structure is inherited by
every population built from it.

So the possibilities are exhaustive:

    rule set respects sigma   ->  A_M = 0 exactly, for all motifs
    rule set violates sigma   ->  A_M can be nonzero, but the asymmetry was
                                 put into the rules

Asking for A_M != 0 "without asymmetry in the rules or the initial
conditions" rules out both sources, so the answer is no. Nothing here is
spontaneous symmetry breaking; that requires degenerate vacua and a
thermodynamic limit, neither of which a finite symmetric continuation graph
provides.

The conjugation used
--------------------
On a fixed five-node vertex set the edge complement sigma(e) = V \\ e maps
2-edges to 3-edges and back, and sigma^2 = identity. It is the natural
candidate conjugation for this rule set, exchanging exactly the two edge
species the rule set acts on.

Both branches are verified below: the bundled rule set is shown to violate
sigma, its symmetrisation is shown to respect it, and A_M is measured in
each case.

Result
------
Verified by exact enumeration with bitmask states over the 20 possible edges
on five nodes, from a conjugation-symmetric seed.

    symmetrised rule set : 349412 states, sigma-closed, 0 equivariance
                          violations, N(M) = N(Mbar) = 82385 exactly,
                          A_M = +0.0000000000
    bundled rule set     : 4 states, not sigma-closed, 4/4 violations,
                          A_M = +1.0

The populations are exactly equal in the symmetric case -- 82385 against
82385, not approximately -- and maximally unequal when the rules break sigma.
No spontaneous asymmetry exists in between.

An earlier run of this file reported A_M = -0.13 for the symmetrised rule set.
That was a truncation artifact: exploration stopped at a 200000-state cap
mid-flight, so the partial state set was not sigma-closed. The tell was that
the reported count equalled the cap exactly.

This is not a failure of the model. Sakharov's conditions require C and CP
violation for a matter/antimatter asymmetry, and the experiment reproduces
precisely that requirement: conjugation asymmetry must be put into the
dynamics, it cannot arise from symmetric rules and a symmetric preparation.
"""

from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
SIG = [IDX[frozenset(V - e)] for e in EDGES]          # sigma as index permutation

def sigma(mask):
    out = 0
    for i in range(20):
        if mask >> i & 1:
            out |= 1 << SIG[i]
    return out

# precompute moves
COMPOSE, DECOMPOSE = [], []
for i in range(10):
    for j in range(i + 1, 10):
        a, b = EDGES[i], EDGES[j]
        if len(a & b) == 1:
            COMPOSE.append((i, j, IDX[a | b]))
for k in range(10, 20):
    e = EDGES[k]
    for mid in sorted(e):
        o = sorted(e - {mid})
        DECOMPOSE.append((k, IDX[frozenset((mid, o[0]))], IDX[frozenset((mid, o[1]))]))

def succ_base(m):
    for i, j, k in COMPOSE:
        if (m >> i & 1) and (m >> j & 1) and not (m >> k & 1):
            yield m & ~(1 << i) & ~(1 << j) | (1 << k)
    for k, i, j in DECOMPOSE:
        if (m >> k & 1):
            rest = m & ~(1 << k)
            if not (rest >> i & 1) and not (rest >> j & 1):
                yield rest | (1 << i) | (1 << j)

def succ_sym(m):
    yield from succ_base(m)
    for g in succ_base(sigma(m)):
        yield sigma(g)

def reachable(seed, succ):
    seen = {seed}; stack = [seed]
    while stack:
        m = stack.pop()
        for g in succ(m):
            if g not in seen:
                seen.add(g); stack.append(g)
    return seen

def stats(name, states, succ):
    closed = all(sigma(m) in states for m in states)
    bad = sum(1 for m in states
              if {sigma(g) for g in succ(m)} != set(succ(sigma(m))))
    nM = nMb = 0
    for m in states:
        b = bin(m & 0x3FF).count("1"); t = bin(m >> 10).count("1")
        if b > t: nM += 1
        elif t > b: nMb += 1
    A = (nM - nMb) / (nM + nMb) if nM + nMb else 0.0
    print("  %s" % name)
    print("    states=%d  sigma-closed=%s  equivariance violations=%d"
          % (len(states), closed, bad))
    print("    N(M)=%d N(Mbar)=%d   A_M = %+.10f" % (nM, nMb, A))
    print()

seed = (1 << IDX[frozenset((0, 1))]) | (1 << IDX[frozenset((2, 3, 4))])
assert sigma(seed) == seed, "seed not conjugation-symmetric"

print("exact enumeration, conjugation-symmetric seed")
stats("bundled (compose+decompose)", reachable(seed, succ_base), succ_base)
stats("symmetrised (+ sigma-conjugates)", reachable(seed, succ_sym), succ_sym)
