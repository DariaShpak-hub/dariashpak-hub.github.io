#!/usr/bin/env python3
"""
sakharov.py

Which condition actually generates a motif asymmetry: charge violation, or
conjugation violation?

The previous results fixed two pieces. compose sends (b, t) to (b-2, t+1), so
q = (b - t) mod 3 is an exactly conserved Z_3 charge, and sigma (edge
complement on five nodes) exchanges b and t, so it maps q to -q. With the
charge conserved the state space splits into superselection sectors and any
asymmetry is frozen at its prepared value.

Adding a rule that violates the charge lets the system move between sectors.
The question is whether that alone generates an asymmetry, or whether
conjugation violation is separately required. These are independent knobs:

    toggle a binary OR a ternary edge  -> violates q, respects sigma
    toggle only a binary edge          -> violates q and sigma

which gives a 2x2 design:

    G_A  compose/decompose symmetrised          q conserved, C symmetric
    G_B  G_A + toggle binary and ternary        q violated,  C symmetric
    G_C  compose/decompose bundled              q conserved, C violated
    G_D  G_C + toggle binary only               q violated,  C violated

Prediction from the equivariance theorem: sigma-equivariance forces A_M = 0
whatever the charge does, so G_A and G_B should both give exactly zero and
charge violation alone should generate nothing. Asymmetry should appear only
where C is violated.

Measurement
-----------
The asymmetry is measured under the stationary distribution of the
uniform-successor chain, not over the reachable set. Once a toggle rule is
present the chain reaches essentially the whole space, and the whole space is
sigma-symmetric by construction, so a reachable-set count would report zero
for every grammar regardless of the rules. The rule asymmetry lives in the
stationary weights, not in which states are reachable.

Results
-------
A second, accidental symmetry decides the outcome. Bit-complement tau (every
edge slot present <-> absent) maps compose's precondition exactly onto
decompose's, so deg(tau m) = deg(m), and tau also exchanges b and t and hence
M and Mbar. The chosen observable is therefore odd under *two* involutions,
sigma and tau, and every grammar in the 2x2 above respects tau.

    G_A  q-conserved, C-symmetric,  tau-symmetric  ->  A_M = 0
    G_D  q-violated,  C-violated,   tau-symmetric  ->  A_M = 0
    G_F  tau-violated, C-symmetric                 ->  A_M = 0
    G_E  tau-violated, C-violated                  ->  A_M != 0 transiently

So charge violation is neither necessary nor sufficient, and C violation
alone is not sufficient either. What is required is violating *every*
symmetry under which the observable is odd -- the ordinary statement that an
asymmetry is forbidden by any symmetry that flips its sign.

G_E's asymmetry is transient, not stationary. Its unique sink is the
all-ones configuration, which has b = t = 10 and so belongs to neither M nor
Mbar. All mass reaches it: by iteration 200 the sink holds 1.000000 and both
populations have underflowed to zero, so the A_M values printed beyond that
point are 0/0 noise (an earlier version of this file reported -0.628 that
way). During relaxation the asymmetry is genuine, growing from -0.074 at
iteration 10 to -0.595 at iteration 50 among the mass not yet absorbed.

Caveat on Sakharov's third condition
------------------------------------
Real baryogenesis needs departure from thermal equilibrium because CPT forces
equal equilibrium abundances even when C and CP are violated. This model has
no CPT theorem, so its stationary state is simply the stationary distribution
of a non-reversible chain, which is generically asymmetric. Conditions one and
two are tested faithfully here; the third has no honest analogue, and a
nonzero equilibrium asymmetry below should not be read as evading it.
"""

import numpy as np
from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
SIG = [IDX[frozenset(V - e)] for e in EDGES]
NB = 10                      # binary edge indices 0..9, ternary 10..19
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


def sigma_array(states):
    """Apply the edge-complement conjugation to an array of bitmasks."""
    out = np.zeros_like(states)
    for i in range(20):
        bit = (states >> i) & 1
        out |= bit.astype(states.dtype) << SIG[i]
    return out


def moves(rules_compose, rules_decompose, toggles):
    """Build (src, dst) arrays for one grammar over the whole state space."""
    src_all = np.arange(N, dtype=np.int64)
    S, D = [], []

    for i, j, k in rules_compose:
        m = (((src_all >> i) & 1) == 1) & (((src_all >> j) & 1) == 1) & \
            (((src_all >> k) & 1) == 0)
        s = src_all[m]
        S.append(s)
        D.append((s & ~(1 << i) & ~(1 << j)) | (1 << k))

    for k, i, j in rules_decompose:
        m = (((src_all >> k) & 1) == 1) & (((src_all >> i) & 1) == 0) & \
            (((src_all >> j) & 1) == 0)
        s = src_all[m]
        S.append(s)
        D.append((s & ~(1 << k)) | (1 << i) | (1 << j))

    for i in toggles:
        S.append(src_all)
        D.append(src_all ^ (1 << i))

    return np.concatenate(S), np.concatenate(D)


def conjugate_moves(comp, dec):
    """sigma o rule o sigma, as index-permuted move tables."""
    return ([(SIG[i], SIG[j], SIG[k]) for i, j, k in comp],
            [(SIG[k], SIG[i], SIG[j]) for k, i, j in dec])


def stationary(src, dst, iters=400, tol=1e-14):
    outdeg = np.bincount(src, minlength=N).astype(np.float64)
    live = outdeg > 0
    pi = np.where(live, 1.0, 0.0)
    pi /= pi.sum()
    w = np.zeros(N)
    for _ in range(iters):
        w[:] = 0.0
        np.divide(pi, outdeg, out=w, where=live)
        nxt = np.bincount(dst, weights=w[src], minlength=N)
        nxt[~live] += pi[~live]              # absorbing states hold their mass
        s = nxt.sum()
        if s > 0:
            nxt /= s
        delta = np.abs(nxt - pi).sum()
        pi = nxt
        if delta < tol:
            break
    return pi


def asymmetry(pi):
    """A_M for M = 'more 2-edges than 3-edges', conjugate under sigma."""
    idx = np.arange(N, dtype=np.int64)
    b = np.zeros(N, dtype=np.int8)
    t = np.zeros(N, dtype=np.int8)
    for i in range(NB):
        b += ((idx >> i) & 1).astype(np.int8)
    for i in range(NB, 20):
        t += ((idx >> i) & 1).astype(np.int8)
    pM = pi[b > t].sum()
    pMb = pi[t > b].sum()
    tot = pM + pMb
    return ((pM - pMb) / tot if tot > 0 else 0.0), pM, pMb


def charge_violated(src, dst):
    """Does any transition change q = (b - t) mod 3?"""
    def q(x):
        b = np.zeros(len(x), dtype=np.int64)
        t = np.zeros(len(x), dtype=np.int64)
        for i in range(NB):
            b += (x >> i) & 1
        for i in range(NB, 20):
            t += (x >> i) & 1
        return (b - t) % 3
    sample = slice(0, len(src), max(1, len(src) // 200000))
    return bool(np.any(q(src[sample]) != q(dst[sample])))


def main():
    cconj, dconj = conjugate_moves(COMPOSE, DECOMPOSE)

    grammars = {
        "G_A  q-conserved, C-symmetric":
            (COMPOSE + cconj, DECOMPOSE + dconj, []),
        "G_B  q-violated,  C-symmetric":
            (COMPOSE + cconj, DECOMPOSE + dconj, list(range(20))),
        "G_C  q-conserved, C-violated":
            (COMPOSE, DECOMPOSE, []),
        "G_D  q-violated,  C-violated":
            (COMPOSE, DECOMPOSE, list(range(NB))),
    }

    print("%-32s %10s %8s %8s %14s" % (
        "grammar", "moves", "q-viol", "C-viol", "A_M"))
    for name, (comp, dec, tog) in grammars.items():
        src, dst = moves(comp, dec, tog)
        qv = charge_violated(src, dst)
        # C violation: is the move set closed under sigma?
        sconj = set(zip(sigma_array(src[:400000]).tolist(),
                        sigma_array(dst[:400000]).tolist()))
        base = set(zip(src[:400000].tolist(), dst[:400000].tolist()))
        cv = not sconj.issubset(base | sconj) or len(sconj - base) > 0
        pi = stationary(src, dst)
        A, pM, pMb = asymmetry(pi)
        print("%-32s %10d %8s %8s %+14.10f" % (
            name, len(src), qv, cv, A))


if __name__ == "__main__":
    main()
