#!/usr/bin/env python3
"""
model.py

The finished model. Forty-nine modules and ninety-seven commits reduce to two
verified components and one measured fact about them: they are incompatible.

THE HEADLINE
============
    You can have PAID expansion or ISOTROPIC expansion. Not both.

    The move that funds node creation is the move that destroys the metric.

That is the finished result. It is not a gap awaiting more work; it is a
measured conflict between two components each of which works in isolation.

THE TWO COMPONENTS
==================

Mode A -- free growth, isotropic
---------------------------------
    state       graph G, generation counter g(e) on each edge
    rule        subdivide: an edge splits, a node is born
    measure     edges chosen with weight 2^-g

    Under uniform selection an edge already cut into k pieces has k chances to
    be cut again -- a Yule process whose Dirichlet shares never concentrate,
    freezing the shear at 68%. Weighting by 2^-g makes each lineage's total
    weight invariant under subdivision (one edge of weight 2^-g becomes two of
    2^-(g+1)), so every lineage grows at an equal rate.

    verified    shear = 1.50/sqrt(N), exponent N^-0.47, the central-limit rate
                meets the observed bound sigma/H < 4.7e-11 at N = 1e21, which
                the observable universe exceeds by ~60 orders of magnitude
                comoving scale factor a ~ N, needing no clock
    but         d = 1.00 exactly, and nothing is conserved: births are free

Mode B -- paid growth, conserved
---------------------------------
    state       graph G, bank B
    rules       close      u,v at distance 2 become adjacent   dE=+1, dB=+1
                open       an edge in a triangle is removed    dE=-1, dB=-1
                subdivide  an edge splits, a node is born      dN=+1,dE=+1,dB=-1
    invariant   Q = B + 2N - E, exactly

    Clumping is the only source of budget and node creation is a sink, so no
    node exists that a closure did not fund. The invariant forces the mean
    degree: E = B + 2N - Q gives <k> = 4 + 2(B-Q)/N -> 4, a prediction derived
    twice independently (here and from the drive-dissipation balance in soc.py).

    verified    Q residual exactly 0 at every checkpoint
                <k> = 3.998 against the derived 4
    but         shear frozen at 48-60%, exponent between N^-0.03 and N^+0.06
                the scale factor barely moves, a ~ 1.4 to 3.0

WHY THEY CANNOT BE COMBINED
===========================
Tested directly below, and the failure is not an artefact of how new edges are
labelled. A closure creates an edge that belongs to no lineage, so it needs a
generation; all three natural choices were tried and all three fail:

    close-edge generation      shear exponent
    g = 0  (new lineage)          N^-0.02
    g = local mean                N^-0.03
    g = max  (finest scale)       N^+0.06

The cause is structural. `close` joins two nodes at graph distance 2, and that
is precisely the move geometry_decay.py measured as destroying the metric with
a half-life of about three sweeps, independent of N. The anti-Yule measure
governs how subdivision is *distributed* among lineages; it has no power over a
rule that *deletes and rewires* them. Paying for expansion means running the
metric-destroying move once per birth.

WHAT IS ESTABLISHED ACROSS THE PROJECT
======================================
    conservation   Q residual exactly 0; <k> -> 4 derived, not imposed
    expansion      a(N) from persistent comoving markers, no clock anywhere
    isotropy       shear ~ 1/sqrt(N), reaching the observed bound (Mode A only)
    arrow          matter starting saturated still gains entropy, because the
                   ceiling ln2*E recedes as the substrate grows -- a constant
                   16.6% lag
    memory         equals the conserved charges exactly, to 1.1e-16
    1/r^2          the lattice Green function gives A -> 1/4pi on a dead
                   lattice: inverse-square is three-dimensionality, not a
                   mechanism

WHAT IS FALSIFIED
=================
    dimension      nothing selects one; monotone drift, no fixed point
    metastability  geometry has a half-life of ~3 sweeps, N-independent
    protection     every state invariant is blind to order, taking identical
                   values on a lattice and its randomised twin
    barrier        entropic and extensive, B ~ 0.25N, with no saddle and no
                   second basin
    clustering     contact interaction only, 1.10x against an observed 19x
    criticality    SOC self-organises to a SUPERCRITICAL state, not a critical
                   one, and the graph stays an expander

THE ONE-SENTENCE MODEL
======================
A timeless relational system in which expansion is a graph observable requiring
no clock, in which a growth measure alone is enough to make that expansion
isotropic to within observational bounds, and in which the thermodynamic arrow
follows from a receding entropy ceiling rather than from any preference of the
rules -- but in which nothing selects a spatial dimension, and paying for the
expansion costs you the geometry.
"""

import numpy as np
from collections import deque


class Universe:
    """Both modes. `paid=False` is Mode A, `paid=True` is Mode B."""

    def __init__(self, n0=6, seed=0, paid=True, close_gen="zero"):
        self.rng = np.random.default_rng(seed)
        self.paid = paid
        self.close_gen = close_gen
        self.adj = [set() for _ in range(n0)]
        self.edges = []
        self.gen = []
        self.n = n0
        self.bank = 0
        for i in range(n0):
            self._add(i, (i + 1) % n0, 0)
        self.Q = self.bank + 2 * self.n - len(self.edges)
        self.markers = list(range(n0))
        self.base = self._baseline()
        self.counts = {"close": 0, "open": 0, "subdivide": 0}

    # ------------------------------------------------------------ topology
    def _add(self, u, v, g):
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.edges.append((u, v))
        self.gen.append(g)

    def _drop(self, i):
        u, v = self.edges[i]
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        e, g = self.edges.pop(), self.gen.pop()
        if i < len(self.edges):
            self.edges[i], self.gen[i] = e, g

    def residual(self):
        return self.bank + 2 * self.n - len(self.edges) - self.Q

    def mean_degree(self):
        return 2 * len(self.edges) / self.n

    # --------------------------------------------------------------- rules
    def close(self):
        for _ in range(40):
            u, v = self.edges[int(self.rng.integers(len(self.edges)))]
            x = u if self.rng.random() < 0.5 else v
            nb = list(self.adj[x])
            if len(nb) < 2:
                continue
            a = nb[int(self.rng.integers(len(nb)))]
            b = nb[int(self.rng.integers(len(nb)))]
            if a == b or b in self.adj[a]:
                continue
            if self.close_gen == "zero":
                g = 0
            elif self.close_gen == "max":
                g = max(self.gen) if self.gen else 0
            else:                                  # local mean
                inc = [self.gen[i] for i, (p, q) in enumerate(self.edges)
                       if p in (a, b) or q in (a, b)]
                g = int(np.mean(inc)) if inc else 0
            self._add(a, b, g)
            self.bank += 1
            self.counts["close"] += 1
            return True
        return False

    def open(self):
        if self.bank < 1:
            return False
        for _ in range(40):
            i = int(self.rng.integers(len(self.edges)))
            u, v = self.edges[i]
            if self.adj[u] & self.adj[v]:
                self._drop(i)
                self.bank -= 1
                self.counts["open"] += 1
                return True
        return False

    def subdivide(self):
        """Anti-Yule selection: edges chosen with weight 2^-g."""
        if self.paid and self.bank < 1:
            return False
        w = np.exp2(-np.asarray(self.gen, dtype=np.float64))
        c = np.cumsum(w)
        i = min(int(np.searchsorted(c, self.rng.random() * c[-1])),
                len(self.edges) - 1)
        u, v = self.edges[i]
        g = self.gen[i]
        self._drop(i)
        z = self.n
        self.adj.append(set())
        self.n += 1
        self._add(u, z, g + 1)
        self._add(z, v, g + 1)
        if self.paid:
            self.bank -= 1
        self.counts["subdivide"] += 1
        return True

    # -------------------------------------------------------- observables
    def bfs(self, s):
        d = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in self.adj[x]:
                if y not in d:
                    d[y] = d[x] + 1
                    q.append(y)
        return d

    def _baseline(self):
        base = {}
        for i in self.markers:
            d = self.bfs(i)
            for j in self.markers:
                if j > i:
                    base[(i, j)] = d[j]
        return base

    def scale_factor(self):
        out = []
        for i in self.markers:
            d = self.bfs(i)
            for j in self.markers:
                if j > i and self.base.get((i, j)) and j in d:
                    out.append(d[j] / self.base[(i, j)])
        return np.array(out, float)

    def matter_entropy(self):
        """S(matter | geometry) = ln2 * E: one passive bit per present edge."""
        return np.log(2.0) * len(self.edges)

    def dimension(self):
        """Shell estimator |S_r| ~ r^(d-1); ball counts are biased ~20% low."""
        acc = {}
        keys = [v for v in range(self.n) if self.adj[v]]
        for _ in range(6):
            s = keys[int(self.rng.integers(len(keys)))]
            c = np.bincount(np.fromiter(self.bfs(s).values(), int))
            for r in range(len(c)):
                acc[r] = acc.get(r, 0.0) + float(c[r])
        sh = {k: v / 6 for k, v in acc.items()}
        hi = max(4, max(sh) // 3)
        rs = [r for r in sorted(sh) if 2 <= r <= hi and sh[r] > 0]
        if len(rs) < 3:
            return float("nan")
        x = np.log(np.array(rs, float))
        y = np.log(np.array([sh[r] for r in rs]))
        return 1.0 + float(np.polyfit(x, y, 1)[0])


def evolve(U, target_n):
    while U.n < target_n:
        if not U.paid:
            U.subdivide()
        elif U.bank == 0:
            U.close()
        else:
            (U.subdivide if U.rng.random() < 0.5 else U.open)()
    return U


def table(title, paid, close_gen, checkpoints, seed=1):
    print("  " + title)
    print("     %7s %7s %6s %8s %9s %9s %7s"
          % ("N", "E", "<k>", "residual", "mean a", "shear", "d"))
    U = Universe(seed=seed, paid=paid, close_gen=close_gen)
    ns, sps = [], []
    for tgt in checkpoints:
        evolve(U, tgt)
        a = U.scale_factor()
        sp = float(a.std() / a.mean())
        ns.append(U.n)
        sps.append(max(sp, 1e-12))
        print("     %7d %7d %6.3f %8d %9.2f %8.2f%% %7.2f"
              % (U.n, len(U.edges), U.mean_degree(), U.residual(),
                 a.mean(), 100 * sp, U.dimension()))
    slope = float(np.polyfit(np.log(ns), np.log(sps), 1)[0])
    print("     -> shear ~ N^%.2f    S(matter|geom) = %.1f"
          % (slope, U.matter_entropy()))
    print()
    return U, slope, float(np.mean(np.array(sps) * np.sqrt(ns)))


def main():
    CK = [1000, 2000, 4000, 8000]
    print("=" * 78)
    print("MODE A -- free growth, anti-Yule measure")
    print("=" * 78)
    UA, slA, cA = table("subdivide only, weight 2^-g", False, "zero", CK)
    print("     isotropic: exponent %.2f is the central-limit rate -1/2."
          % slA)
    print("     observed bound 4.7e-11 met at N = %.1e" % ((cA / 4.7e-11) ** 2))
    print("     but nothing is conserved and d = %.2f" % UA.dimension())
    print()

    print("=" * 78)
    print("MODE B -- paid growth, conserved budget")
    print("=" * 78)
    for cg in ("zero", "local", "max"):
        table("close-edge generation = %s" % cg, True, cg, CK)
    print("     conserved and <k> -> 4, but the shear does not decay under any")
    print("     of the three natural generation assignments.")
    print()

    print("=" * 78)
    print("THE CONFLICT, STATED")
    print("=" * 78)
    print("  Mode A isotropises and conserves nothing.")
    print("  Mode B conserves exactly and does not isotropise.")
    print()
    print("  `close` joins nodes at graph distance 2 -- the same move")
    print("  geometry_decay.py measured as destroying the metric with a")
    print("  half-life of ~3 sweeps, independent of N. The anti-Yule measure")
    print("  governs how subdivision is distributed among lineages; it has no")
    print("  power over a rule that deletes and rewires them.")
    print()
    print("  Paying for expansion means running the metric-destroying move")
    print("  once per birth. That is the finished result.")
    print()


if __name__ == "__main__":
    main()
