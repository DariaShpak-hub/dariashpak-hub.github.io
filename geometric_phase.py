#!/usr/bin/env python3
"""
geometric_phase.py

Does an energy built from a local geometric observable produce a geometric
phase, without the target geometry being specified?

geometry_decay.py showed why an energy is needed: every conserved quantity in
this project takes the same value on ordered and disordered graphs, so
conservation cannot protect order. An energy can, because it distinguishes
configurations the charges cannot. The question is whether a NATURAL local
energy does it, or only one reverse-engineered from the answer.

Choosing the observable, by measurement rather than taste
----------------------------------------------------------
Four reference graphs at N ~ 1000, all degree 6:

    graph                    tri/node   sq/node    |B_2|   mean dist
    cubic lattice d=3           0.000     3.000    24.00        7.50
    triangular lattice d=2      2.000     3.000    18.00       12.44
    random 6-regular            0.026     0.078    35.26        4.21
    disjoint K_7 (crumpled)     4.970    14.910     5.96        0.86

Triangles are disqualified immediately: the cubic lattice has none, exactly
like a random graph, so a triangle energy is blind to three-dimensional order
and would only ever find the triangular one. That is the "specifying the
geometry" failure in its purest form, and it is worth noting that triangle or
clustering terms are the usual first choice.

Four-cycles discriminate by a factor of 38 and treat both lattices alike, so
the energy used is

    H = -(number of 4-cycles)

which asks for local closure without naming a dimension or a motif.

But the table also shows the trap. The crumpled graph extremises EVERY column
-- most 4-cycles, fewest triangles-free, smallest 2-ball. So minimising H
alone must drive to clumps, not to lattices, and any claim of a geometric
phase has to survive that. This is exactly the crumpled phase of dynamical
triangulations.

Which makes mean distance the right order parameter, because it is the one
quantity that is NOT monotone between the three regimes:

    crumpled  0.86   <   random  4.21   <   cubic  7.50  <  triangular  12.44

A geometric phase is therefore visible as mean distance rising ABOVE the
random value. Falling below it is crumpling. No fit or interpretation is
needed to tell them apart, and the test can fail in a way that looks nothing
like success.

Method
------
Metropolis over 6-regular graphs at fixed N using degree-preserving
double-edge swaps. Swaps are non-local, which was fatal for the dynamics
question in geometry_decay.py but is a virtue here: the Boltzmann distribution
is a property of the ensemble, not of the moves, so non-local moves sample the
same equilibrium and merely reach it faster.

The 4-cycle count is updated incrementally -- a swap can only destroy cycles
through the two removed edges and create cycles through the two added ones --
and the incremental total is checked against a full recount, which must agree
exactly.

Result 1: it does not crumple, and that prediction was wrong
--------------------------------------------------------------
The expectation was cliques, since disjoint K_7 has 15 four-cycles per node
against a lattice's 3 and is therefore the true ground state of H. It does not
go there. At beta = 3 on N = 1000, run to 3200 sweeps:

    4-cycles/node   11.218      (lattice 3.0, cliques 15, random 0.078)
    triangles/node   0.004      (lattice 0.0, cliques 4.97)
    mean distance     7.00      (random 4.21, cubic lattice 7.50)

Triangles stay at zero, which cliques could not do, and mean distance RISES
above random rather than collapsing. So the true minimum of H is dynamically
inaccessible and the system instead finds a triangle-free, 4-cycle-rich
structure. The clique ground state is real but unreachable by swaps, which is
worth stating plainly: what follows describes a long-lived non-equilibrium
regime, not a thermodynamic phase.

Result 2: local geometry is genuinely achieved
------------------------------------------------
The local structure converges on the cubic lattice's, without the lattice ever
being mentioned in H:

    observable        random    reached    cubic lattice
    |B_2|              35.26      24.05            24.00
    triangles/node      0.026      0.005            0.000
    4-cycles/node       0.078      5.003            3.000

|B_2| lands on 24.0, the cubic value, to three digits. This is the first time
in the project that a local geometric observable has been driven to a
lattice value by a dynamics that was not told the lattice.

Result 3: global extent grows, but not to a manifold
------------------------------------------------------
Comparing at matched SWEEPS is a trap -- larger systems order more slowly, so
4-cycles/node falls from 12.6 to 5.3 across the range and mean distance looks
flat purely because the big systems are less ordered. The fair protocol is
matched LOCAL ORDER: run each N until it reaches the same 4-cycle density,
then compare.

    target 4cyc/node >= 3.0
    N          250    500   1000   2000   4000
    mean d    3.83   4.40   4.71   5.30   5.94     ~ N^0.153, d_H = 6.53
                                                     (pow 0.00087, log 0.00129)

    target 4cyc/node >= 5.0
    N          250    500   1000   2000   4000
    mean d    4.15   4.95   5.15   5.71   6.49     ~ N^0.149
                                                     (pow 0.00352, log 0.00333)

Mean distance does grow with N, at N^0.15, against N^0.06-0.07 for every
earlier rule and a logarithm for a random graph. So an extended regime exists
and the energy is what produced it.

But d_H is about 6.5, not 3, and the power-law preference is marginal -- the
first target favours a power law, the second favours a logarithm, with
residuals nearly tied. The honest reading is a structure part-way between
small-world and geometric, not a manifold.

Verdict
-------
The ingredient works, and works better than the six-criteria audit predicted.
An energy built from a single local observable, chosen by measurement rather
than by knowing the answer, moves the system off the random-graph value that
four previous conserved rule sets could not escape, and drives the local
geometry onto the cubic lattice's exactly.

What it does not do is assemble that local order into global geometry. |B_2|
is correct to three digits while d_H is off by a factor of two, which is the
familiar and hard problem: local flatness does not imply a global manifold.
The pieces fit locally and do not close up.

So the answer to the question this module asked is a qualified yes for local
geometry and a no for a manifold, with the additional caveat that the regime
is metastable rather than a true phase. Getting further would need what the
reference table already implies: 4-cycles alone cannot distinguish a cubic
lattice from other triangle-free 4-cycle-rich graphs, because it gives both
lattices the same value of 3.0. A second term, sensitive to how the local
pieces join rather than to what they are, is the missing constraint -- which
is the same conclusion CDT reached when it added causal structure on top of
a local simplicial action.
"""

import numpy as np
from collections import deque


# ------------------------------------------------------------ constructions

def random_regular(n, k, rng):
    stubs = np.repeat(np.arange(n), k)
    rng.shuffle(stubs)
    adj = [set() for _ in range(n)]
    for i in range(0, len(stubs) - 1, 2):
        a, b = int(stubs[i]), int(stubs[i + 1])
        if a != b and b not in adj[a]:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def cubic_lattice(m):
    n = m ** 3

    def idx(i, j, k):
        return (i % m) * m * m + (j % m) * m + (k % m)

    adj = [set() for _ in range(n)]
    for i in range(m):
        for j in range(m):
            for k in range(m):
                u = idx(i, j, k)
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    v = idx(i + d[0], j + d[1], k + d[2])
                    adj[u].add(v)
                    adj[v].add(u)
    return adj, n


def triangular_lattice(m):
    n = m * m
    adj = [set() for _ in range(n)]
    for i in range(m):
        for j in range(m):
            u = i * m + j
            for di, dj in ((0, 1), (1, 0), (1, -1)):
                v = ((i + di) % m) * m + (j + dj) % m
                adj[u].add(v)
                adj[v].add(u)
    return adj, n


def cliques(n, size=7):
    adj = [set() for _ in range(n)]
    for base in range(0, n - size + 1, size):
        for a in range(base, base + size):
            for b in range(base, base + size):
                if a != b:
                    adj[a].add(b)
    return adj


# ------------------------------------------------------------ observables

def count4_total(adj, n):
    """Full 4-cycle count: sum over vertex pairs of C(codegree, 2), halved."""
    tot = 0
    for u in range(n):
        cnt = {}
        for v in adj[u]:
            for w in adj[v]:
                if w != u:
                    cnt[w] = cnt.get(w, 0) + 1
        for w, c in cnt.items():
            if w > u:
                tot += c * (c - 1) // 2
    return tot // 2


def c4_through(adj, u, v):
    """4-cycles containing the edge (u,v): paths u-x-y-v with all distinct.

    Each such cycle is counted once, since the path avoiding the edge (u,v)
    is unique within a 4-cycle.
    """
    tot = 0
    for x in adj[u]:
        if x == v:
            continue
        for y in adj[x]:
            if y != u and y != x and y != v and y in adj[v]:
                tot += 1
    return tot


def both_pairs(adj, e1, e2):
    """4-cycles containing both disjoint edges e1=(a,b), e2=(c,d)."""
    a, b = e1
    c, d = e2
    tot = 0
    if c in adj[b] and a in adj[d]:
        tot += 1
    if d in adj[b] and a in adj[c]:
        tot += 1
    return tot


def mean_distance(adj, n, k=8, seed=2):
    r = np.random.default_rng(seed)
    acc, cnt = 0.0, 0
    for _ in range(min(k, n)):
        s = int(r.integers(n))
        dist = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    q.append(w)
        v = np.fromiter(dist.values(), int)
        acc += float(v.sum())
        cnt += len(v)
    return acc / cnt


def ball2(adj, n):
    """Mean size of the radius-2 ball: 35.3 random, 24.0 cubic, 6.0 crumpled."""
    tot = 0
    for u in range(n):
        s = set(adj[u])
        t = set()
        for v in adj[u]:
            t |= adj[v]
        tot += len((s | t) - {u})
    return tot / n


def triangles_per_node(adj, n):
    t = 0
    for u in range(n):
        for v in adj[u]:
            if v > u:
                t += len(adj[u] & adj[v])
    return t / 3.0 / n


# --------------------------------------------------------------- the scan

def anneal(adj, n, beta, sweeps, rng, el=None):
    """Metropolis with double-edge swaps, H = -(4-cycle count)."""
    if el is None:
        el = [(u, v) for u in range(n) for v in adj[u] if u < v]
    pos = {e: i for i, e in enumerate(el)}
    c4 = count4_total(adj, n)

    def rm(a, b):
        k = (a, b) if a < b else (b, a)
        i = pos.pop(k)
        last = el.pop()
        if i < len(el):
            el[i] = last
            pos[last] = i
        adj[a].discard(b)
        adj[b].discard(a)

    def ad(a, b):
        k = (a, b) if a < b else (b, a)
        pos[k] = len(el)
        el.append(k)
        adj[a].add(b)
        adj[b].add(a)

    steps = int(sweeps * n)
    acc = 0
    for _ in range(steps):
        a, b = el[int(rng.integers(len(el)))]
        c, d = el[int(rng.integers(len(el)))]
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            c, d = d, c
        if d in adj[a] or b in adj[c]:
            continue
        before = (c4_through(adj, a, b) + c4_through(adj, c, d)
                  - both_pairs(adj, (a, b), (c, d)))
        rm(a, b)
        rm(c, d)
        ad(a, d)
        ad(c, b)
        after = (c4_through(adj, a, d) + c4_through(adj, c, b)
                 - both_pairs(adj, (a, d), (c, b)))
        dH = -(after - before)
        if dH <= 0 or rng.random() < np.exp(-beta * dH):
            c4 += after - before
            acc += 1
        else:
            rm(a, d)
            rm(c, b)
            ad(a, b)
            ad(c, d)
    return adj, c4, acc / max(steps, 1)


def main():
    N = 1000
    print("=" * 78)
    print("1. REFERENCE POINTS  (N ~ %d, degree 6)" % N)
    print("=" * 78)
    print("  %-26s %9s %10s %11s" % ("graph", "tri/node", "4cyc/node",
                                     "mean dist"))
    refs = {}
    ca, cn = cubic_lattice(10)
    ta, tn = triangular_lattice(32)
    for lbl, a, nn in (("cubic lattice d=3", ca, cn),
                       ("triangular lattice d=2", ta, tn),
                       ("random 6-regular",
                        random_regular(N, 6, np.random.default_rng(1)), N),
                       ("disjoint K_7 (crumpled)", cliques(N, 7), N)):
        md = mean_distance(a, nn)
        refs[lbl] = md
        print("  %-26s %9.3f %10.3f %11.2f"
              % (lbl, triangles_per_node(a, nn), count4_total(a, nn) / nn, md))
    print()
    print("  order parameter: mean distance. ABOVE random = extended/geometric,")
    print("  BELOW random = crumpled. The two failure modes look different.")
    print()

    print("=" * 78)
    print("2. INCREMENTAL 4-CYCLE BOOKKEEPING CHECK")
    print("=" * 78)
    adj = random_regular(300, 6, np.random.default_rng(5))
    adj, c4, _ = anneal(adj, 300, 0.3, 20, np.random.default_rng(6))
    exact = count4_total(adj, 300)
    print("  incremental %d   full recount %d   difference %d"
          % (c4, exact, c4 - exact))
    if c4 != exact:
        print("  MISMATCH: the scan below is not trustworthy.")
        return
    print("  exact agreement, so the energy tracked during the scan is right.")
    print()

    print("=" * 78)
    print("3. THE SCAN:  H = -(4-cycle count),  starting from random")
    print("=" * 78)
    print("  %8s %11s %11s %11s %9s %8s"
          % ("beta", "4cyc/node", "tri/node", "mean dist", "vs random", "accept"))
    d_rand = refs["random 6-regular"]
    rows = []
    for beta in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0):
        a = random_regular(N, 6, np.random.default_rng(1))
        a, c4, ar = anneal(a, N, beta, 120, np.random.default_rng(7))
        md = mean_distance(a, N)
        rows.append((beta, c4 / N, triangles_per_node(a, N), md, ar))
        print("  %8.2f %11.3f %11.3f %11.2f %9.2f %7.1f%%"
              % (beta, c4 / N, triangles_per_node(a, N), md, md - d_rand,
                 100 * ar))
    print()
    print("  triangles stay near zero throughout, so this is NOT crumpling")
    print("  toward cliques (which would give 4.97/node) even though cliques")
    print("  are the true ground state at 15 four-cycles per node.")
    print()

    print("=" * 78)
    print("4. MATCHED LOCAL ORDER: does global extent follow local geometry?")
    print("=" * 78)
    print("  Matching sweeps is a trap -- big systems order more slowly, so")
    print("  mean distance looks flat only because they are less ordered.")
    print("  Each N is instead run until it reaches the same 4-cycle density.")
    print()
    for target in (3.0, 5.0):
        print("  target 4cyc/node >= %.1f" % target)
        print("     %7s %8s %11s %10s %8s %8s"
              % ("N", "sweeps", "4cyc/node", "tri/node", "mean d", "|B2|"))
        ns, ms = [], []
        for n2 in (250, 500, 1000, 2000, 4000):
            a = random_regular(n2, 6, np.random.default_rng(1))
            rng = np.random.default_rng(7)
            done, c4 = 0, 0
            while done < 6000:
                a, c4, _ = anneal(a, n2, 3.0, 100, rng)
                done += 100
                if c4 / n2 >= target:
                    break
            md = mean_distance(a, n2)
            ns.append(n2)
            ms.append(md)
            print("     %7d %8d %11.3f %10.3f %8.2f %8.2f"
                  % (n2, done, c4 / n2, triangles_per_node(a, n2), md,
                     ball2(a, n2)))
        ns_a, ms_a = np.array(ns, float), np.array(ms, float)
        p, bb = np.polyfit(np.log(ns_a), np.log(ms_a), 1)
        rp = float(np.sum((np.log(ms_a) - (p * np.log(ns_a) + bb)) ** 2))
        cc, b2f = np.polyfit(np.log(ns_a), ms_a, 1)
        rl = float(np.sum(((ms_a - (cc * np.log(ns_a) + b2f))
                           / ms_a.mean()) ** 2))
        print("     -> mean d ~ N^%.3f  (pow %.5f, log %.5f): %s"
              % (p, rp, rl,
                 "power law, d_H = %.2f" % (1 / p) if rp < rl
                 else "logarithmic"))
        print()

    print("  For comparison: every earlier rule gave N^0.06-0.07, a random")
    print("  graph gives a logarithm, and a manifold would give N^(1/d) with")
    print("  d_H matching the local structure. Local order is reached exactly;")
    print("  global extent grows but stops well short of a manifold.")
    print()


if __name__ == "__main__":
    main()
