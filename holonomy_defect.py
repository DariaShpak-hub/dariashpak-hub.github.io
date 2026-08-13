#!/usr/bin/env python3
"""
holonomy_defect.py

A locally measurable relational defect in an otherwise intact vacuum.

Two independent routes have now killed the same idea. The p_comm scan in
coloured_metric.py found zero holonomy defects at every commutation probability
including p = 0, because a square that fails to close does not contradict
anything -- the two paths simply reach different nodes, and the loop that would
have detected the disagreement is never formed. Address-frustrated fission
fails the same way and worse: cutting a periodic seam removes exactly L edges,
so Delta beta_1 = -L, which is a global topology change masquerading as a local
defect, and not a fixed-size rewrite vector at all.

    missing loop  !=  nontrivial holonomy

The fix is to separate POSITION from TRANSPORTED FRAME. Give each edge a
transport variable as well as an adjacency, and a loop can close in position
while the frame comes back rotated. The smallest group that does this is Z_2:

    U_e in {+1, -1},    H(gamma) = product of U_e around gamma
    vacuum   H = +1 on every face
    defect   H = -1, with the loop completely intact

This is Z_2 lattice gauge theory, and two of its standard properties matter
here. The flux is gauge invariant -- flipping every edge at a node leaves all
faces unchanged -- so the defect is not an artifact of labelling. And flipping
a chain of edges flips only the faces at its two ends, so defects are created
in PAIRS joined by an invisible string, which is the Z_2 version of a
vortex-antivortex pair.

Why a gauge field alone is not enough
--------------------------------------
U_e is a label. It does not touch adjacency, so graph distances are literally
unchanged and epsilon = 0 identically. A passive decoration passes the holonomy
gate and cannot reach the metric gate. The transport has to couple to the
rewrite rule, and the coupling chosen here is the simplest one available:

    subdivision is allowed only on the boundary of a frustrated face

So a defect is a place where the complex keeps refining. That gives, by
construction and checkable:

    dN = +1, dE = +1 per move       a local fixed-size rewrite vector
    beta_1 unchanged                 subdivision preserves cycle rank exactly
    epsilon >= 0                     paths lengthen, nothing is ever cut
    holonomy preserved               the two new edges carry (U_e, +1)

Every one of those is a gate the seam fission failed.

What is genuinely open
-----------------------
Isotropy, and the prediction is specific. ONE defect must fail: it stretches
the paths that pass through it and leaves the rest alone, which is directional
tearing and gives A of order 1. The claim worth testing is that isotropy is a
property of a dilute GAS of defects, by the central limit theorem:

    A ~ n_d^(-1/2)      at fixed size, more defects means rounder
    A ~ N^(-1/2)        at fixed density, since n_d grows with N

If that holds, isotropy is emergent and statistical rather than built in, and
it arrives at the same central-limit rate that shear_damping.py (N^-0.47) and
coloured_metric.py (N^-0.47 at k = 2) both measured by unrelated mechanisms.
If A stalls at a constant instead, the defect picture fails the isotropy gate
exactly as fission and chords did, and the transport variable buys nothing.

The control is the point of comparison: the same number of subdivisions placed
on uniformly random edges. That is isotropic by construction, so it is the
floor A has to approach if the defect gas is working.

Result 1: the vacuum, the defect, and the loop survives
--------------------------------------------------------
L = 16 torus, N = 256, E = 512, beta_1 = 257, zero frustrated faces.
Flipping a chain of 6 edges gives exactly two frustrated faces, at (4,8) and
(10,8) -- the two ends of the chain, with the interior flipped twice and back.
beta_1 stays 257: nothing was cut.

    Wilson loop                        H
    1x1 far from the defects          +1
    1x1 on a defect                   -1
    2x2 enclosing one                 -1
    8x2 enclosing both                +1

A random gauge transformation changes 0 of 256 face fluxes. So the defect is
gauge invariant, locally detectable by a single plaquette, and -- unlike every
earlier candidate in this project -- it lives in a loop that is completely
intact. The holonomy gate passes, and so do the three gates seam fission
failed: dN = +1 and dE = +1 per move, beta_1 unchanged, nothing ever cut.

Result 2: my own isotropy statistic was measuring sparsity, and the
control caught it
--------------------------------------------------------------------
The first run gave A ~ 6-12 for the defect gas and reported it as anisotropy.
It is not. If only a fraction f of marker pairs move at all, then

    A = sd/mean ~ 1/sqrt(f)

whatever the directional structure is, and the uniform control -- isotropic by
construction -- returned A = 6.7 to 9.6 at the same f. Identical numbers from a
process with no preferred direction at all. The statistic was reading how few
pairs had been touched.

The fix is to keep raising the subdivision budget until essentially every pair
has moved, and only then read A. At L = 24, budget as a multiple of E:

    pairs  subdiv     mode   mean eps      A    moved   1/sqrt(f)
        8     288   defect     0.0118   5.10     5.7%       4.20
        8     288  uniform     0.0403   2.40    26.2%       1.95
        8    1152   defect     0.0150   4.66     6.6%       3.88
        8    1152  uniform     0.2996   0.76    94.2%       1.03
        8    4608   defect     0.0158   4.56     6.6%       3.88
        8    4608  uniform     1.5802   0.37    99.8%       1.00
        8   18432   defect     0.0165   4.21     7.4%       3.68
        8   18432  uniform     6.9036   0.31   100.0%       1.00

       32     288   defect     0.0340   2.71    19.6%       2.26
       32    1152   defect     0.0941   1.85    34.2%       1.71
       32    4608   defect     0.1660   1.91    38.8%       1.61
       32   18432   defect     0.3013   1.98    42.4%       1.54

In the sparse regime A tracks 1/sqrt(f) for BOTH processes, which is the
artifact. Once the control saturates at f = 1 the tracking stops and A finally
means something: 0.31, and still falling with budget.

Result 3: the defect gas fails, and the reason is geometric
-------------------------------------------------------------
The two columns behave completely differently as the budget grows 64-fold:

    uniform    mean eps  0.040 -> 6.90    unbounded, A -> 0.31, f -> 100%
    defect     mean eps  0.012 -> 0.017   SATURATES, A ~ 4.4, f stuck at 7%

At 8 pairs, a 64x increase in refinement moves the mean dilation by 40% and the
moved fraction not at all. At 32 pairs it is better but still saturating: f
levels off near 42% and A near 1.9, against the control's 0.31.

The mechanism is not subtle once seen. Subdividing the boundary of a face makes
that face's edges long, so shortest paths route AROUND the refined region
rather than through it. Localized refinement is invisible to the metric,
because geodesics avoid it. Pouring more subdivisions into the same faces
changes nothing, which is exactly what the saturation shows.

Verdict: transport passes, the coupling fails, and the negative is general
---------------------------------------------------------------------------
The Z_2 transport variable does what it was introduced to do. A closed loop
with nontrivial holonomy exists, it is gauge invariant, it is locally
detectable, and creating it costs a fixed-size local rewrite. That is the first
time in this project that a defect has existed at all -- both incomplete
commutation and seam fission produced missing loops instead.

What fails is the coupling from holonomy to metric, and the failure generalises
beyond this particular choice:

    LOCALIZED REFINEMENT CANNOT PRODUCE EXPANSION.

Shortest paths route around any bounded region that has been made expensive, so
the dilation saturates no matter how much refinement is added. Expansion has to
be spatially extensive -- and the control shows that when it is, it works
properly: unbounded mean dilation and A falling to 0.31 and still dropping.

Which sharpens the standing hole rather than filling it. The project now has a
defect sector and an expansion mechanism, and they are measurably incompatible
in the obvious coupling: the thing that makes a defect local is exactly what
makes it metrically invisible. Any successful version needs the holonomy to
modify transport EVERYWHERE, not to mark a place where refinement happens.

Caveats. Two dimensions only. Z_2 flux defects are points in 2D but loops in
3D, so the defect sector does not carry over unchanged. The marker sample is 24
nodes and A in the saturated regime is measured at one system size, so the
control's 0.31 is a value, not a scaling law.
"""

import numpy as np
from collections import deque


class Z2Lattice:
    """Periodic L x L square lattice with a Z_2 transport variable per edge."""

    def __init__(self, L, seed=0):
        self.L = L
        self.rng = np.random.default_rng(seed)
        self.n = L * L
        self.adj = [set() for _ in range(self.n)]
        self.eu, self.ev, self.U = [], [], []
        self.edge_faces = []
        self.faces = []
        self.H = {}
        self.V = {}
        for y in range(L):
            for x in range(L):
                self.H[(x, y)] = self._add(self.node(x, y),
                                           self.node((x + 1) % L, y))
                self.V[(x, y)] = self._add(self.node(x, y),
                                           self.node(x, (y + 1) % L))
        for y in range(L):
            for x in range(L):
                f = [self.H[(x, y)], self.V[((x + 1) % L, y)],
                     self.H[(x, (y + 1) % L)], self.V[(x, y)]]
                fi = len(self.faces)
                self.faces.append(list(f))
                for e in f:
                    self.edge_faces[e].append(fi)

    def node(self, x, y):
        return y * self.L + x

    def _add(self, u, v, U=1):
        i = len(self.eu)
        self.eu.append(u)
        self.ev.append(v)
        self.U.append(U)
        self.edge_faces.append([])
        self.adj[u].add(v)
        self.adj[v].add(u)
        return i

    # ------------------------------------------------------------- topology
    def n_edges(self):
        return len(self.eu)

    def betti1(self):
        return self.n_edges() - self.n + 1

    def flux(self, fi):
        p = 1
        for e in self.faces[fi]:
            p *= self.U[e]
        return p

    def frustrated(self):
        return [f for f in range(len(self.faces)) if self.flux(f) < 0]

    # --------------------------------------------------------------- defects
    def insert_pair(self, x0, y, m):
        """Flip V(x,y) for x in (x0, x0+m]: flips only the two end faces."""
        for j in range(1, m + 1):
            e = self.V[((x0 + j) % self.L, y)]
            self.U[e] = -self.U[e]
        return (self.node(x0, y), self.node((x0 + m) % self.L, y))

    def gauge_transform(self, rng):
        """U_e -> g_u U_e g_v. Physical content must be untouched."""
        g = rng.choice([-1, 1], size=self.n)
        for i in range(len(self.eu)):
            self.U[i] *= int(g[self.eu[i]]) * int(g[self.ev[i]])

    def wilson(self, x0, y0, a, b):
        """Rectangular loop on the pristine lattice, corner (x0,y0), size axb."""
        L, p = self.L, 1
        for i in range(a):
            p *= self.U[self.H[((x0 + i) % L, y0 % L)]]
            p *= self.U[self.H[((x0 + i) % L, (y0 + b) % L)]]
        for j in range(b):
            p *= self.U[self.V[(x0 % L, (y0 + j) % L)]]
            p *= self.U[self.V[((x0 + a) % L, (y0 + j) % L)]]
        return p

    # -------------------------------------------------------------- dynamics
    def _subdivide_clean(self, e, u, v, Ue, fs):
        """u-v becomes u-z-v carrying (U_e, +1). dN=+1, dE=+1, flux fixed."""
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        z = self.n
        self.adj.append(set())
        self.n += 1
        self.eu[e], self.ev[e], self.U[e] = u, z, Ue
        self.adj[u].add(z)
        self.adj[z].add(u)
        nb = self._add(z, v, 1)
        self.edge_faces[nb] = list(fs)
        for fi in fs:
            self.faces[fi].append(nb)
        return z

    def grow_on_frustration(self, steps, rng):
        fr = self.frustrated()
        if not fr:
            return 0
        done = 0
        for _ in range(steps):
            fi = fr[int(rng.integers(len(fr)))]
            es = self.faces[fi]
            e = es[int(rng.integers(len(es)))]
            self._subdivide_clean(e, self.eu[e], self.ev[e], self.U[e],
                                  list(self.edge_faces[e]))
            done += 1
        return done

    def grow_uniform(self, steps, rng):
        for _ in range(steps):
            e = int(rng.integers(len(self.eu)))
            self._subdivide_clean(e, self.eu[e], self.ev[e], self.U[e],
                                  list(self.edge_faces[e]))
        return steps

    # ---------------------------------------------------------- measurement
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


def marker_distances(G, markers):
    out = {}
    for i in markers:
        d = G.bfs(i)
        for j in markers:
            if j > i and j in d:
                out[(i, j)] = d[j]
    return out


def epsilon_stats(base, now):
    e = np.array([(now[p] - base[p]) / base[p] for p in base
                  if p in now and base[p] > 0], float)
    if e.size < 3:
        return float("nan"), float("nan"), float("nan"), 0.0
    mu, sd = float(e.mean()), float(e.std())
    return mu, sd, (sd / mu if mu > 0 else float("nan")), float((e > 0).mean())


def experiment(L, n_pairs, steps_per_pair=None, seed=0, uniform=False,
               n_markers=24):
    rng = np.random.default_rng(seed)
    G = Z2Lattice(L, seed=seed)
    b1_before = G.betti1()
    markers = sorted(rng.choice(L * L, size=min(n_markers, L * L),
                                replace=False).tolist())
    base = marker_distances(G, markers)
    sep = max(2, L // 4)
    if not uniform:
        for _ in range(n_pairs):
            x0 = int(rng.integers(L))
            y = int(rng.integers(L))
            G.insert_pair(x0, y, sep)
    nfr = len(G.frustrated())
    total = (steps_per_pair or (L * L) // 8) * max(n_pairs, 1)
    if uniform:
        G.grow_uniform(total, rng)
    else:
        G.grow_on_frustration(total, rng)
    now = marker_distances(G, markers)
    mu, sd, A, frac = epsilon_stats(base, now)
    return dict(L=L, N=L * L, pairs=n_pairs, defects=nfr, steps=total,
                b1_before=b1_before, b1_after=G.betti1(),
                mean_eps=mu, sd_eps=sd, A=A, frac_moved=frac)


def main():
    print("=" * 78)
    print("1. THE VACUUM, THE DEFECT, AND THAT THE LOOP SURVIVES")
    print("=" * 78)
    L = 16
    G = Z2Lattice(L, seed=1)
    print("  pristine L=%d torus   N=%d  E=%d  beta_1=%d  frustrated faces=%d"
          % (L, G.n, G.n_edges(), G.betti1(), len(G.frustrated())))
    G.insert_pair(4, 8, 6)
    fr = G.frustrated()
    print("  after one flipped chain of 6 edges: frustrated faces = %d at %s"
          % (len(fr), [(f % L, f // L) for f in fr]))
    print("  beta_1 = %d   (unchanged: nothing was cut)" % G.betti1())
    print()
    print("  Wilson loops on the pristine lattice:")
    for name, (x0, y0, a, b) in (
            ("1x1 far from defects", (12, 2, 1, 1)),
            ("1x1 on a defect     ", (4, 8, 1, 1)),
            ("2x2 enclosing one   ", (3, 7, 2, 2)),
            ("8x2 enclosing both  ", (3, 7, 8, 2))):
        print("     %s   H = %+d" % (name, G.wilson(x0, y0, a, b)))
    print()
    rng = np.random.default_rng(7)
    before = [G.flux(f) for f in range(len(G.faces))]
    G.gauge_transform(rng)
    after = [G.flux(f) for f in range(len(G.faces))]
    print("  random gauge transformation: faces changed = %d of %d"
          % (sum(1 for a, b in zip(before, after) if a != b), len(before)))
    print("  So the defect is gauge invariant -- it is not a labelling")
    print("  artifact, and it is locally detectable by one plaquette.")
    print()

    print("=" * 78)
    print("2. IS A MEASURING ANISOTROPY, OR JUST SPARSITY?")
    print("=" * 78)
    print("  If only a fraction f of marker pairs move at all, then")
    print("  A = sd/mean ~ 1/sqrt(f) whatever the direction structure is.")
    print("  The uniform control is isotropic BY CONSTRUCTION, so if it shows")
    print("  the same A at the same f, then A is reading sparsity, not shape.")
    print()
    L, E0 = 24, 2 * 24 * 24
    print("  %7s %9s %9s %10s %8s %9s %11s"
          % ("pairs", "subdiv", "mode", "mean eps", "A", "moved", "1/sqrt(f)"))
    for npair in (8, 32):
        for mult in (0.25, 1, 4, 16):
            for uni in (False, True):
                tot = int(mult * E0)
                rs = []
                for s in range(3):
                    rng = np.random.default_rng(s)
                    G = Z2Lattice(L, seed=s)
                    mk = sorted(rng.choice(L * L, size=24,
                                           replace=False).tolist())
                    base = marker_distances(G, mk)
                    if uni:
                        G.grow_uniform(tot, rng)
                    else:
                        for _ in range(npair):
                            G.insert_pair(int(rng.integers(L)),
                                          int(rng.integers(L)), max(2, L // 4))
                        G.grow_on_frustration(tot, rng)
                    rs.append(epsilon_stats(base, marker_distances(G, mk)))
                mu = np.nanmean([r[0] for r in rs])
                A = np.nanmean([r[2] for r in rs])
                f = np.nanmean([r[3] for r in rs])
                print("  %7d %9d %9s %10.4f %8.3f %8.1f%% %11.2f"
                      % (npair, tot, "uniform" if uni else "defect", mu, A,
                         100 * f, 1.0 / np.sqrt(max(f, 1e-9))))
        print()

    print("=" * 78)
    print("3. WHY THE DEFECT GAS SATURATES")
    print("=" * 78)
    print("  Refinement confined to a defect boundary makes those edges long,")
    print("  so shortest paths route AROUND the refined region instead of")
    print("  through it. Adding more subdivisions to the same faces then buys")
    print("  nothing, and both mean eps and the moved fraction level off.")
    print()


if __name__ == "__main__":
    main()
