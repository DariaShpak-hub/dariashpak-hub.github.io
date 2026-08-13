#!/usr/bin/env python3
"""
cosmology_audit.py

The cosmological model tested with time strictly excluded.

Admissibility rule: an observable counts only if it is a function of the
static transition graph. Graph distance qualifies -- it is a metric on the
state space, invariant under both reparameterising and rate-reweighting the
schedule. Step counts, rates, dS/dt, H, q, w_eff, burst factors and
recurrence times do not qualify and appear nowhere below.

The model makes four claims that survive translation into time-free form:

    1  the initial state is special: few accessible configurations
    2  the constraint releases itself: the gate is reachable from within
    3  release opens a much larger space
    4  the arrow follows from the preparation, not from the rules

Each is tested twice: against the rest of the sector, to ask whether the
chosen seed is unusual rather than arbitrary, and against random rule sets of
matched arity, to ask whether the behaviour belongs to compose/decompose or
to any rule of that shape.

The seed's percentile is the number that matters. A model resting on a
"highly constrained initial state" needs its seed to sit in the tail. If the
seed is typical of its sector, claim 1 is not a property of the preparation
at all, and claims 2 to 4 inherit that.

Results: every claim is either generic or done better by random rule sets
------------------------------------------------------------------------
Claim 1 fails outright. The path seed sits at the 51.7th percentile of its
own sector, with dS_B = -0.0154, marginally *above* the equilibrium mean
rather than below it. Its macrostate holds 125 of the sector's 290 states
while the smallest holds 15. Only 5.2% of states in the sector have dS_B > 0
at all. The seed is median, not constrained.

Note this contradicts nothing measured earlier, and that is the point: under
the degree-sequence coarse-graining the same seed occupies a singleton
macrostate, Omega = 1, and is maximally special. Under the arity profile it
occupies the largest block. Specialness is a property of the (state,
coarse-graining) pair, exactly as the arrow is, and cannot be attributed to
the preparation alone.

Claim 2 holds but is generic. The gate lies at graph distance 2 from the
seed, closer than 79% of the sector, so the constraint does release itself.
But 100% of random matched-arity rule sets also have a reachable gate, so
self-release distinguishes nothing about compose/decompose.

Claims 3 and 4 are worse than random. The real rule set opens
ln(sector/Omega) = 0.84 against a null median of 1.18, and gives a static
gradient g(seed) = 0.182 against a null median of 0.705, roughly four times
weaker. Random rule sets of the same arity expand more and have a stronger
entropic gradient than the one the model is built on.

The single respect in which the real rule set stands out is that it is more
confining: sector size 290 against a null median of 662. That follows from
conserving connectivity, which the null model identified as its one genuinely
structural property. It is also precisely what makes it worse on claims 3
and 4. The geometric structure that is real works against the cosmological
reading rather than supporting it.
"""

import numpy as np
from collections import deque, defaultdict
from itertools import combinations

N_NODES = 5
EDGES = [frozenset(c) for c in combinations(range(N_NODES), 2)] + \
        [frozenset(c) for c in combinations(range(N_NODES), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = sum(1 for e in EDGES if len(e) == 2)
NE = len(EDGES)
N = 1 << NE
BIN_MASK = (1 << NB) - 1

REAL = [(i, j, IDX[EDGES[i] | EDGES[j]])
        for i in range(NB) for j in range(i + 1, NB)
        if len(EDGES[i] & EDGES[j]) == 1]

A = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int8); T = np.zeros(N, np.int8)
for i in range(NB):
    B += ((A >> i) & 1).astype(np.int8)
for i in range(NB, NE):
    T += ((A >> i) & 1).astype(np.int8)


def adjacency(triples):
    """Undirected neighbour lists, as a dict of state -> set of states."""
    adj = defaultdict(set)
    for i, j, k in triples:
        m = (((A >> i) & 1) == 1) & (((A >> j) & 1) == 1) & (((A >> k) & 1) == 0)
        s = A[m]
        d = (s & ~(1 << i) & ~(1 << j)) | (1 << k)
        for u, v in zip(s.tolist(), d.tolist()):
            adj[u].add(v); adj[v].add(u)
    return adj


def sector_of(adj, seed):
    seen = {seed}; q = deque([seed])
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):
            if v not in seen:
                seen.add(v); q.append(v)
    return seen


def distances_to(adj, sector, targets):
    """Graph distance from every sector state to the nearest target."""
    dist = {t: 0 for t in targets if t in sector}
    q = deque(dist)
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):
            if v in sector and v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def macro(state):
    """Arity profile: the coarse-graining with real compression."""
    return (int(B[state]), int(T[state]))


def audit(name, triples, seed, verbose=True):
    adj = adjacency(triples)
    sector = sector_of(adj, seed)
    if len(sector) < 8:
        return None

    blocks = defaultdict(list)
    for x in sector:
        blocks[macro(x)].append(x)
    omega = {m: len(v) for m, v in blocks.items()}
    n = len(sector)
    S_eq = sum((omega[m] / n) * np.log(omega[m]) for m in omega)

    # claim 1: is the seed's macrostate unusually small?
    gaps = np.array([S_eq - np.log(omega[macro(x)]) for x in sorted(sector)])
    seed_gap = S_eq - np.log(omega[macro(seed)])
    pct_gap = float((gaps < seed_gap).mean())

    # claim 2: gate reachability, purely structural
    gates = [x for x in sector if (x & BIN_MASK) == 0]
    if gates:
        dist = distances_to(adj, sector, gates)
        seed_d = dist.get(seed, None)
        ds = np.array([dist[x] for x in sorted(sector) if x in dist])
        pct_d = float((ds > seed_d).mean()) if seed_d is not None else float("nan")
    else:
        seed_d, pct_d, ds = None, float("nan"), np.array([0])

    # claim 3: how much larger is the space beyond the gate?
    #          measured structurally as sector size vs the seed's macrostate
    opening = np.log(n) - np.log(omega[macro(seed)])

    # claim 4: static gradient at the seed's macrostate
    acc = 0.0; cnt = 0
    su = np.log(omega[macro(seed)])
    for x in blocks[macro(seed)]:
        for y in adj.get(x, ()):
            if y in sector:
                acc += np.log(omega[macro(y)]) - su; cnt += 1
    g_seed = acc / cnt if cnt else 0.0

    if verbose:
        print("  sector size            : %d   macrostates %d" % (n, len(omega)))
        print("  seed Omega             : %d   (smallest in sector: %d)"
              % (omega[macro(seed)], min(omega.values())))
        print("  claim 1  dS_B(seed)    : %+.4f   percentile %.1f%%"
              % (seed_gap, 100 * pct_gap))
        print("  claim 2  dist to gate  : %s   percentile %.1f%%  (gates: %d)"
              % (seed_d, 100 * pct_d, len(gates)))
        print("  claim 3  ln(sector/Om) : %+.4f" % opening)
        print("  claim 4  g(seed)       : %+.4f" % g_seed)
    return dict(n=n, gap=seed_gap, pct_gap=pct_gap, dist=seed_d,
                pct_d=pct_d, opening=opening, g=g_seed, gates=len(gates))


def main():
    seed = 0
    for e in ((0, 1), (1, 2), (2, 3), (3, 4)):
        seed |= 1 << IDX[frozenset(e)]

    print("=" * 72)
    print("REAL rule set, path seed  (b=%d, t=%d)" % (B[seed], T[seed]))
    print("=" * 72)
    real = audit("real", REAL, seed)
    print()

    print("=" * 72)
    print("Is the SEED special? Same rule set, every state in the sector")
    print("=" * 72)
    adj = adjacency(REAL)
    sector = sorted(sector_of(adj, seed))
    blocks = defaultdict(list)
    for x in sector:
        blocks[macro(x)].append(x)
    omega = {m: len(v) for m, v in blocks.items()}
    n = len(sector)
    S_eq = sum((omega[m] / n) * np.log(omega[m]) for m in omega)
    gaps = np.array([S_eq - np.log(omega[macro(x)]) for x in sector])
    print("  dS_B over the sector: min %+.4f  median %+.4f  max %+.4f"
          % (gaps.min(), np.median(gaps), gaps.max()))
    print("  seed                : %+.4f   -> percentile %.1f%%"
          % (real["gap"], 100 * real["pct_gap"]))
    frac_pos = float((gaps > 0).mean())
    print("  fraction of ALL states with dS_B > 0 : %.3f" % frac_pos)
    print()

    print("=" * 72)
    print("NULL rule sets, matched arity, same seed")
    print("=" * 72)
    rng = np.random.default_rng(11)
    rows = []
    for _ in range(12):
        tr = set()
        while len(tr) < len(REAL):
            i, j = rng.choice(NB, 2, replace=False)
            k = int(rng.integers(NB, NE))
            tr.add((int(min(i, j)), int(max(i, j)), k))
        r = audit("null", sorted(tr), seed, verbose=False)
        if r:
            rows.append(r)
    for key, lab in (("n", "sector size"), ("gap", "dS_B(seed)"),
                     ("pct_gap", "gap percentile"), ("opening", "ln(sector/Om)"),
                     ("g", "g(seed)")):
        vals = [r[key] for r in rows]
        print("  %-16s real %10.4f    null median %10.4f   [%.4f, %.4f]"
              % (lab, real[key], float(np.median(vals)),
                 float(np.min(vals)), float(np.max(vals))))
    gate_frac = float(np.mean([r["gates"] > 0 for r in rows]))
    print("  %-16s real %10s    null %10.0f%% have a reachable gate"
          % ("self-release", "yes" if real["gates"] else "no", 100 * gate_frac))


if __name__ == "__main__":
    main()
