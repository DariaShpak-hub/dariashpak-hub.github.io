#!/usr/bin/env python3
"""
sector_breaking.py

Can a conjugation-symmetric rule set break the symmetry spontaneously?

conjugation.py showed that a sigma-equivariant rule set has A_M = 0 exactly.
That argument is about the *ensemble*, and it is exactly the argument that
spontaneous symmetry breaking evades in ordinary physics: a symmetric
Hamiltonian has vanishing average magnetisation while ergodicity breaks and
any single realisation sits in one asymmetric sector.

The sharper question is therefore about sectors rather than averages. Take the
full state space, decompose it into connected components of the continuation
graph under a sigma-equivariant rule set, and ask how sigma acts on the set of
components.

    sigma fixes every component
        every sector is self-conjugate, A_M = 0 within each one as well as
        overall, and no preparation can retain an asymmetry.

    sigma maps some component to a different one
        those two components are mirror images. A system prepared in one can
        never reach the other, so it carries N_M != N_Mbar forever, from
        symmetric rules and a symmetric rule set. That is spontaneous
        symmetry breaking in the ergodicity-breaking sense, and it is not
        asymmetry inserted into the dynamics.

States are bitmasks over the 20 possible edges on five nodes, so the entire
space of 2^20 configurations is enumerated exactly rather than sampled.

Result: superselection, not symmetry breaking
---------------------------------------------
sigma-paired components do exist, which looks at first like the mechanism
being sought. The full space splits into three large components of 349520,
349520 and 349412 states plus 124 frozen singletons, and the two equal-sized
ones are sigma-mirrors carrying A_M = +0.1894 and -0.1894.

But they are distinguished by a conserved charge, so this is not spontaneous.
The three components are exactly the residue classes of (b - t) mod 3, where
b and t count 2-edges and 3-edges: compose sends (b, t) to (b-2, t+1), so
b - t falls by exactly 3 and the residue is invariant. sigma exchanges b and
t, hence maps the charge q to -q, which is why q=1 and q=2 are mirrors while
q=0 is self-conjugate and has A_M = 0.

Genuine spontaneous symmetry breaking requires the mirror sectors to be
degenerate in every conserved quantity, with ergodicity broken by suppressed
tunnelling rather than by a conservation law. Here a Z_3 charge forbids the
transition outright. Choosing a sector is choosing the charge, and the charge
is fixed by the initial condition -- exactly the source the question rules
out.

The structure mirrors baryon number correctly. A conserved charge exists,
conjugation flips its sign, and no asymmetry can be generated from a
symmetric start without violating it. That is Sakharov's first condition,
reproduced rather than evaded.
"""

from collections import defaultdict
from itertools import combinations

V = frozenset(range(5))
EDGES = [frozenset(c) for c in combinations(range(5), 2)] + \
        [frozenset(c) for c in combinations(range(5), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
SIG = [IDX[frozenset(V - e)] for e in EDGES]
N = 1 << 20

SIGMA = [0] * N
for m in range(N):
    o = 0
    mm = m
    while mm:
        i = (mm & -mm).bit_length() - 1
        o |= 1 << SIG[i]
        mm &= mm - 1
    SIGMA[m] = o

COMPOSE, DECOMPOSE = [], []
for i in range(10):
    for j in range(i + 1, 10):
        if len(EDGES[i] & EDGES[j]) == 1:
            COMPOSE.append((i, j, IDX[EDGES[i] | EDGES[j]]))
for k in range(10, 20):
    e = EDGES[k]
    for mid in sorted(e):
        o = sorted(e - {mid})
        DECOMPOSE.append((k, IDX[frozenset((mid, o[0]))],
                          IDX[frozenset((mid, o[1]))]))


def succ_base(m):
    for i, j, k in COMPOSE:
        if (m >> i & 1) and (m >> j & 1) and not (m >> k & 1):
            yield m & ~(1 << i) & ~(1 << j) | (1 << k)
    for k, i, j in DECOMPOSE:
        if m >> k & 1:
            rest = m & ~(1 << k)
            if not (rest >> i & 1) and not (rest >> j & 1):
                yield rest | (1 << i) | (1 << j)


def succ_sym(m):
    yield from succ_base(m)
    s = SIGMA[m]
    for g in succ_base(s):
        yield SIGMA[g]


def main():
    parent = list(range(N))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for m in range(N):
        rm = find(m)
        for g in succ_sym(m):
            rg = find(g)
            if rm != rg:
                parent[rg] = rm
                rm = find(m)

    comps = defaultdict(list)
    for m in range(N):
        comps[find(m)].append(m)

    roots = {find(m): find(m) for m in range(N)}
    self_conj = paired = 0
    pair_examples = []
    for r, members in comps.items():
        img = find(SIGMA[members[0]])
        if img == r:
            self_conj += 1
        else:
            paired += 1
            if len(pair_examples) < 3:
                pair_examples.append((r, img, len(members)))

    print("full space: %d configurations, %d connected components"
          % (N, len(comps)))
    print("  components fixed by sigma        : %d" % self_conj)
    print("  components mapped to a different : %d" % paired)
    print()

    def asym(members):
        nM = nMb = 0
        for m in members:
            b = bin(m & 0x3FF).count("1")
            t = bin(m >> 10).count("1")
            if b > t:
                nM += 1
            elif t > b:
                nMb += 1
        tot = nM + nMb
        return ((nM - nMb) / tot if tot else 0.0), nM, nMb

    sizes = sorted((len(v) for v in comps.values()), reverse=True)
    print("  component sizes (largest 6): %s" % sizes[:6])
    big = max(comps.values(), key=len)
    A, nM, nMb = asym(big)
    print("  largest component: %d states, N(M)=%d N(Mbar)=%d, A_M=%+.10f"
          % (len(big), nM, nMb, A))
    print("  largest component sigma-closed: %s"
          % (find(SIGMA[big[0]]) == find(big[0])))

    if pair_examples:
        print()
        print("  sigma-paired components exist -> ergodicity breaking")
        for r, img, size in pair_examples:
            A1, a, b = asym(comps[r])
            A2, c, d = asym(comps[img])
            print("    pair sizes %d/%d   A_M = %+.6f vs %+.6f"
                  % (size, len(comps[img]), A1, A2))
    else:
        print()
        print("  every component is self-conjugate -> no ergodicity breaking,")
        print("  so no preparation can retain a motif asymmetry")


if __name__ == "__main__":
    main()
