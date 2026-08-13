#!/usr/bin/env python3
"""
gravity_arrow.py

Can the arrow matter observes come from geometry rather than from matter?

The proposal under test: matter does not begin out of equilibrium. Geometry
does. Matter starts already at its own entropy ceiling, and acquires an arrow
only because geometry's relaxation moves the ceiling. This is the Weyl
curvature hypothesis in a form this model can actually evaluate.

To test it the state has to be split, so matter is added to the substrate as
the most passive thing possible: one bit on each present edge. It has no
dynamics of its own in part 1 and only a slow randomising dynamics in part 2.
It is carried by geometry, never driving it. Then the total entropy splits

    S_total = S(geometry) + S(matter | geometry)

and with bits conditionally uniform the second term is exactly

    S(matter | geometry) = ln 2 * E[number of edges]

so matter's arrow is nothing but the expected edge count of the geometry. The
question "does gravitational relaxation give matter an arrow" becomes "does
relaxation raise or lower the edge count", which is decidable in closed form.

Two substrates are compared, and they answer it differently.

Part 1, closed substrate: the arrow runs backwards
--------------------------------------------------
Under compose/decompose the weight w = b + 2t is exactly conserved, so

    edges = b + t = w - t

and the matter ceiling is a strictly decreasing function of the ternary count.
Matter therefore gains entropy if and only if the geometry starts with MORE
ternary edges than it will have at equilibrium -- if and only if geometry
starts clumped and relaxes by dispersing.

That is the opposite of what the proposal needs. A smooth low-entropy start --
the path, the star, the cycle, all binary, t = 0 -- is the most dispersed
geometry in its sector, so relaxation can only raise t, and

    dS(matter) = ln 2 * (t_seed - E_eq[t]) < 0   strictly, for every t=0 seed.

Measured, in nats:

    seed              t   dS(matter)   dS(geometry)   dS(total)
    path              0     -0.4302        +5.6699      +5.2397
    star              0     -0.4302        +5.6699      +5.2397
    cycle             0     -0.7193        +7.1483      +6.4290
    triangle+edge     0     -0.7193        +7.1483      +6.4290
    one clump         1     +0.2629        +5.6699      +5.9328
    two clumps        2     +0.9561        +5.6699      +6.6259
    heavy clump       3     +1.0898        +8.2079      +9.2977

The total rises everywhere, so no law is broken -- gravity's entropy increase
is partly paid for out of matter's. But the sign of matter's share tracks
t_seed against the sector mean, and it is negative for exactly the smooth
seeds the hypothesis is about. The seeds that do give matter a positive arrow
are the concentrated ones, which start clumped and relax by spreading. That is
inverted relative to the hypothesis.

So in a closed substrate the answer is no, and the reason is not an accident
of these rules. It is conservation. If the substrate is fixed, matter's phase
space is bounded once and for all, and every nat matter gains is a nat some
other sector state has to lose. Relaxation redistributes the ceiling; it
cannot raise it.

Part 2, open substrate: the arrow needs nodes to be born
--------------------------------------------------------
Subdivide replaces edge {x,y} by a new node z and edges {x,z}, {z,y}. It
conserves nothing -- nodes and edges both grow -- so the ceiling is not fixed:

    S_max(t) = ln 2 * (E_0 + t)

grows without bound. Now matter is given a slow dynamics, one random bit flip
per step, and its marginals are tracked exactly (each p_e evolves as
p_e <- p_e + (1 - 2 p_e)/E, which is exact for the marginal).

Matter entropy now rises for the whole run, from a start where it was already
saturated. But it never catches the ceiling. At 400 steps from a 4-edge start
the ceiling is 280.031 nats, matter holds 233.612, and the gap is 46.419 and
still widening at a constant 0.1145 nats per step -- the rate has not turned
over anywhere in the run. Each new edge is born ordered while ever more edges
compete for the same one flip per step, so per-edge equilibration slows as 1/E
against a constant birth rate. The gap settles to a fixed FRACTION of the
ceiling, 17.5% at step 25 and 16.6% from step 200 on, so matter stays a
constant proportion behind however long the run goes.

That is the mechanism the proposal is reaching for, and it is not relaxation.
Nothing here relaxes. The arrow is the ceiling receding faster than the
contents can follow, which is Layzer's entropy gap, and it requires the
substrate to grow. Slowing the growth kills it -- at 400 steps:

    subdivide every    edges   ceiling   S_matter    gap
    1 step               404   280.031    233.612   46.419
    2 steps              204   141.402    128.263   13.139
    4 steps              104    72.087     68.307    3.780
    8 steps               54    37.430     36.124    1.306

At one birth per eight steps matter tracks the ceiling to within 3.5%.

Verdict
-------
The hypothesis is half right, and the half it gets wrong is the important one.
Matter's arrow can indeed be sourced entirely by geometry with matter starting
in equilibrium -- part 2 exhibits that. But "gravitational relaxation" is the
wrong mechanism for it. Relaxation toward equilibrium in a closed substrate
gives matter a NEGATIVE arrow whenever geometry starts smooth, which is the
case the hypothesis is about. What produces the positive arrow is expansion of
the substrate itself: nodes being born, phase space opening, and matter
falling behind. Growth, not relaxation.
"""

import numpy as np
from collections import deque, defaultdict
from itertools import combinations
from math import log

N_NODES = 5
EDGES = [frozenset(c) for c in combinations(range(N_NODES), 2)] + \
        [frozenset(c) for c in combinations(range(N_NODES), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = sum(1 for e in EDGES if len(e) == 2)
NE = len(EDGES)
N = 1 << NE
LN2 = log(2.0)

COMPOSE = [(i, j, IDX[EDGES[i] | EDGES[j]])
           for i in range(NB) for j in range(i + 1, NB)
           if len(EDGES[i] & EDGES[j]) == 1]

A = np.arange(N, dtype=np.int64)
B = np.zeros(N, np.int8)
T = np.zeros(N, np.int8)
for i in range(NB):
    B += ((A >> i) & 1).astype(np.int8)
for i in range(NB, NE):
    T += ((A >> i) & 1).astype(np.int8)


def adjacency():
    adj = defaultdict(set)
    for i, j, k in COMPOSE:
        m = (((A >> i) & 1) == 1) & (((A >> j) & 1) == 1) & (((A >> k) & 1) == 0)
        s = A[m]
        d = (s & ~(1 << i) & ~(1 << j)) | (1 << k)
        for u, v in zip(s.tolist(), d.tolist()):
            adj[u].add(v)
            adj[v].add(u)
    return adj


ADJ = adjacency()


def sector_of(seed):
    seen = {seed}
    q = deque([seed])
    while q:
        u = q.popleft()
        for v in ADJ.get(u, ()):
            if v not in seen:
                seen.add(v)
                q.append(v)
    return sorted(seen)


def bits(edges):
    m = 0
    for e in edges:
        m |= 1 << IDX[frozenset(e)]
    return m


SEEDS = {
    "path (smooth)":      bits([(0, 1), (1, 2), (2, 3), (3, 4)]),
    "star (smooth)":      bits([(0, 1), (0, 2), (0, 3), (0, 4)]),
    "cycle (smooth)":     bits([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]),
    "tri+edge (smooth)":  bits([(0, 1), (1, 2), (0, 2), (2, 3), (3, 4)]),
    "one clump":          bits([(0, 1, 2), (2, 3), (3, 4)]),
    "two clumps":         bits([(0, 1, 2), (2, 3, 4)]),
    "heavy clump":        bits([(0, 1, 2), (1, 2, 3), (2, 3, 4)]),
}


def part1():
    print("=" * 74)
    print("PART 1  closed substrate: does relaxing geometry give matter an arrow?")
    print("=" * 74)
    print("  matter = one bit per present edge, conditionally uniform, no dynamics")
    print("  so  S(matter | geometry) = ln2 * E[edges]  exactly.")
    print("  weight w = b + 2t is conserved, hence edges = w - t.")
    print()
    print("  %-18s %5s %5s %8s %9s %10s %10s %10s"
          % ("seed", "t", "edg", "sector", "E_eq[t]", "dS_matt", "dS_geom", "dS_tot"))
    for name, seed in SEEDS.items():
        sec = sector_of(seed)
        n = len(sec)
        t_eq = float(np.mean([T[x] for x in sec]))
        e_eq = float(np.mean([B[x] + T[x] for x in sec]))
        e_seed = int(B[seed] + T[seed])
        dS_m = LN2 * (e_eq - e_seed)
        dS_g = log(n)
        print("  %-18s %5d %5d %8d %9.4f %+10.4f %+10.4f %+10.4f"
              % (name, T[seed], e_seed, n, t_eq, dS_m, dS_g, dS_m + dS_g))
    print()
    print("  every t=0 seed loses matter entropy; only the clumped seeds gain it.")
    print("  the total always rises, so gravity's gain is partly matter's loss.")
    print()


def part2(steps=400, births=(1, 2, 4, 8), e0=4):
    print("=" * 74)
    print("PART 2  open substrate: subdivide spawns nodes, ceiling is not fixed")
    print("=" * 74)
    print("  subdivide: {x,y} -> {x,z},{z,y}, a new node each time.")
    print("  matter starts SATURATED (every bit uniform) and is given one random")
    print("  bit flip per step. Bits on newly born edges start ordered, p = 0.")
    print("  per-edge marginals are exact: p <- p + (1 - 2p)/E.")
    print()
    print("  %-14s %10s %10s %10s %10s"
          % ("birth period", "edges", "ceiling", "S_matter", "gap"))
    for period in births:
        p = np.full(e0, 0.5)                # start at the ceiling
        for s in range(1, steps + 1):
            E = len(p)
            p = p + (1.0 - 2.0 * p) / E     # one flip somewhere, exact marginal
            if s % period == 0:             # a node is born, an edge splits
                k = (s // period) % E
                p = np.concatenate([np.delete(p, k), [p[k], 0.0]])
        h = np.where((p > 0) & (p < 1),
                     -(p * np.log(np.clip(p, 1e-300, 1)) +
                       (1 - p) * np.log(np.clip(1 - p, 1e-300, 1))), 0.0)
        S = float(h.sum())
        ceil = len(p) * LN2
        print("  every %-8d %10d %10.3f %10.3f %10.3f"
              % (period, len(p), ceil, S, ceil - S))
    print()

    print("  the gap over time, at one subdivide per step:")
    print("  %8s %10s %10s %10s %10s"
          % ("step", "edges", "ceiling", "S_matter", "gap"))
    p = np.full(e0, 0.5)
    for s in range(1, steps + 1):
        E = len(p)
        p = p + (1.0 - 2.0 * p) / E
        k = s % E
        p = np.concatenate([np.delete(p, k), [p[k], 0.0]])
        if s in (1, 5, 25, 50, 100, 200, 300, 400):
            h = np.where((p > 0) & (p < 1),
                         -(p * np.log(np.clip(p, 1e-300, 1)) +
                           (1 - p) * np.log(np.clip(1 - p, 1e-300, 1))), 0.0)
            print("  %8d %10d %10.3f %10.3f %10.3f"
                  % (s, len(p), len(p) * LN2, float(h.sum()),
                     len(p) * LN2 - float(h.sum())))
    print()
    print("  matter entropy rises the whole way from an already-saturated start,")
    print("  and never catches the ceiling. Nothing relaxed; the ceiling ran.")
    print()


if __name__ == "__main__":
    part1()
    part2()
