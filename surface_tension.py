#!/usr/bin/env python3
"""
surface_tension.py

A nucleation barrier for geometry: does an interface cost make a lattice
metastable?

geometry_decay.py found the lattice dies in O(N) moves with no metastability,
and diagnosed why: every conserved quantity takes the SAME value on ordered
and disordered graphs, so nothing forbids or even slows the transition.
coherence_field.py then showed an Ising field bolted on top cannot supply the
missing barrier, because in its ordered phase H = -J E depends only on the
bond count and is blind to shape.

Classical nucleation theory says exactly what is missing. A superheated liquid
is metastable because a vapour bubble of radius r costs

    dG(r) = -(4/3) pi r^3 dg  +  4 pi r^2 sigma

bulk gain against SURFACE cost, which creates a critical radius. Below it,
bubbles shrink. The model built here has bulk terms and no surface term, so
ordered regions dissolve from their edges at no charge.

The fix has to avoid the trap that killed the coherence field: the order
parameter must be a function of the GRAPH, not a field living on it. So it is
read off the graph's own 4-cycle structure. For node i,

    q_i = sum over pairs of neighbours (u,w) of |N(u) & N(w)| - 1

counts 4-cycles through i. Measured, it separates perfectly:

    2-torus (degree 4)     q = 4 at every node
    3-torus (degree 6)     q = 12 at every node
    random degree 4        q mean 0.010,  fraction with q >= 3: 0.0000
    random degree 6        q mean 0.152,  fraction with q >= 3: 0.0000

so s_i = [q_i >= 3] is an exact indicator of "locally lattice-like", with no
threshold tuning. The energy is then

    H = -eps * sum_i s_i  +  gamma * #{edges (i,j) with s_i != s_j}

a bulk reward for order and a cost per unit interface. Both terms are
functionals of the graph. Metropolis over the same conserved close/open moves,
N fixed.

A note on what success looks like, because the obvious criterion is wrong.
Classical nucleation gives a lifetime proportional to exp(barrier)/N -- larger
systems decay SOONER in absolute time, because there are more places to
nucleate. So "half-life grows with N" is not the signature. The signature is
exponential growth of the half-life with the surface coupling gamma, and that
is what is measured below.

Result: the surface coupling does nothing, for an instructive reason
----------------------------------------------------------------------
    gamma   half-life   final order   final phi
      0.0        1.21        0.9975       0.154
      0.5        1.00        1.0000       0.168
      1.0        0.98        1.0000       0.183
      2.0        1.12        1.0000       0.188
      4.0        1.12        1.0000       0.188

log(half-life) against gamma has slope 0.007, so the half-life multiplies by
1.01 per unit of surface cost. Flat. There is no barrier, no exponential, and
if anything the mildest gammas decay marginally faster. (Size dependence is
equally flat: half-lives 0.68, 1.21, 1.47 at N = 256, 400, 784 for gamma = 0
and 0.70, 1.12, 1.34 for gamma = 2.)

The "final order" column is the diagnostic. It sits at 1.0000 while phi falls
to 0.15, meaning every node is still locally lattice-like after the metric has
gone 85% of the way to random. No ordered/disordered interface ever forms, so
the surface term has nothing to charge for.

Why no threshold can fix it
-----------------------------
The natural response is that q* = 3 is too lenient against a lattice value of
4. It is not a threshold problem. Tracking the raw distribution as the metric
decays:

    sweeps    phi   q mean   f(q>=3)   f(q>=4)   f(q>=6)
      0.25  0.736     4.05    0.9600    0.8850    0.0550
      1.00  0.545     4.78    0.9150    0.7150    0.3000
      4.00  0.224     8.58    0.9900    0.9575    0.8400
     16.00  0.154     9.36    0.9975    0.9950    0.9075

q RISES from 4.05 to 9.36 -- more than doubling -- while phi falls from 0.74
to 0.15. Local 4-cycle richness INCREASES as the geometry is destroyed,
because `close` manufactures local closure while the accumulated rewiring
destroys long-range structure. Every threshold therefore moves the wrong way,
and an energy rewarding local order would actively drive the destruction
rather than resist it.

What this establishes
-----------------------
The failure is not that the wrong local observable was chosen. It is that the
quantity needing protection has no local signature at all. Long-range metric
structure is destroyed by moves that are locally invisible -- worse than
invisible, locally they look like increasing order.

That generalises the sharpest earlier finding rather than adding to it. In
geometric_phase.py the 4-cycle energy drove |B_2| onto the cubic lattice's
24.00 to three digits while d_H stayed wrong by a factor of two. Local order
and global geometry were decoupled there; here they are measured to be
ANTI-correlated under this dynamics.

So no local energy of this family can supply a nucleation barrier for
geometry. Nucleation works in a superheated liquid because a bubble has a
boundary that is locally detectable -- you can stand at the interface and see
liquid on one side and vapour on the other. The destruction of a graph's
metric has no such boundary. It is delocalised, and a local term cannot price
it.

Which leaves the barrier question where the analogy could not take it. Getting
metastability would need a term sensitive to something non-local -- the
spectral gap, a cycle-space invariant, a global embedding cost -- and at that
point it is no longer a local rewriting model, which was the premise of the
whole programme.
"""

import numpy as np
from collections import deque


def square_torus(L):
    n = L * L
    adj = [set() for _ in range(n)]
    for i in range(L):
        for j in range(L):
            u = i * L + j
            for di, dj in ((1, 0), (0, 1)):
                v = ((i + di) % L) * L + (j + dj) % L
                adj[u].add(v)
                adj[v].add(u)
    return adj, n


def random_regular(n, k, rng):
    stubs = np.repeat(np.arange(n), k)
    rng.shuffle(stubs)
    adj = [set() for _ in range(n)]
    for i in range(0, len(stubs) - 1, 2):
        a, b = int(stubs[i]), int(stubs[i + 1])
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def mean_distance(adj, n, srcs=6, seed=2):
    rng = np.random.default_rng(seed)
    acc, cnt = 0.0, 0
    for _ in range(srcs):
        s = int(rng.integers(n))
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


class Nucleation:
    """close/open dynamics under H = -eps*bulk + gamma*interface."""

    def __init__(self, adj, n, eps, gamma, cap=6, qstar=3):
        self.adj = [set(a) for a in adj]
        self.n = n
        self.eps = eps
        self.gamma = gamma
        self.cap = cap
        self.qstar = qstar
        self.el = []
        self.pos = {}
        for u in range(n):
            for v in self.adj[u]:
                if u < v:
                    self.pos[(u, v)] = len(self.el)
                    self.el.append((u, v))
        self.bank = 0

    # ---------------------------------------------------------- structure
    def _add(self, a, b):
        k = (a, b) if a < b else (b, a)
        if k in self.pos:
            return False
        self.adj[a].add(b)
        self.adj[b].add(a)
        self.pos[k] = len(self.el)
        self.el.append(k)
        return True

    def _del(self, a, b):
        k = (a, b) if a < b else (b, a)
        i = self.pos.pop(k)
        last = self.el.pop()
        if i < len(self.el):
            self.el[i] = last
            self.pos[last] = i
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    # ------------------------------------------------- the order parameter
    def s(self, i):
        nb = list(self.adj[i])
        t = 0
        for a in range(len(nb)):
            na = self.adj[nb[a]]
            for b in range(a + 1, len(nb)):
                t += len(na & self.adj[nb[b]]) - 1
                if t >= self.qstar:
                    return 1
        return 1 if t >= self.qstar else 0

    def patch(self, x, y):
        """Nodes whose s can change: within distance 2 of x or y."""
        p = {x, y}
        for z in (x, y):
            p |= self.adj[z]
            for w in self.adj[z]:
                p |= self.adj[w]
        return p

    def patch_energy(self, P):
        sv = {i: self.s(i) for i in P}
        for i in list(P):
            for j in self.adj[i]:
                if j not in sv:
                    sv[j] = self.s(j)
        e = -self.eps * sum(sv[i] for i in P)
        seen = set()
        for i in P:
            for j in self.adj[i]:
                k = (i, j) if i < j else (j, i)
                if k in seen:
                    continue
                seen.add(k)
                if sv[i] != sv[j]:
                    e += self.gamma
        return e

    # --------------------------------------------------------------- moves
    def close(self, rng):
        for _ in range(30):
            a, b = self.el[int(rng.integers(len(self.el)))]
            v = a if rng.random() < 0.5 else b
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            x = nb[int(rng.integers(len(nb)))]
            y = nb[int(rng.integers(len(nb)))]
            if x == y or y in self.adj[x]:
                continue
            if len(self.adj[x]) >= self.cap or len(self.adj[y]) >= self.cap:
                continue
            self._add(x, y)
            P = self.patch(x, y)          # computed with the edge present
            e_after = self.patch_energy(P)
            self._del(x, y)
            e_before = self.patch_energy(P)
            dH = e_after - e_before
            if dH <= 0 or rng.random() < np.exp(-dH):
                self._add(x, y)
                self.bank += 1
                return True
            return False
        return False

    def open(self, rng):
        if self.bank < 1:
            return False
        for _ in range(30):
            a, b = self.el[int(rng.integers(len(self.el)))]
            if not (self.adj[a] & self.adj[b]):
                continue
            P = self.patch(a, b)          # computed with the edge present
            e_before = self.patch_energy(P)
            self._del(a, b)
            e_after = self.patch_energy(P)
            dH = e_after - e_before
            if dH <= 0 or rng.random() < np.exp(-dH):
                self.bank -= 1
                return True
            self._add(a, b)
            return False
        return False

    def order_fraction(self):
        return sum(self.s(i) for i in range(self.n)) / self.n


def half_life(sw, phi):
    for i in range(1, len(phi)):
        if phi[i] <= 0.5 <= phi[i - 1]:
            lo, hi = phi[i - 1], phi[i]
            if abs(hi - lo) < 1e-12:
                return sw[i]
            f = (lo - 0.5) / (lo - hi)
            return float(np.exp(np.log(sw[i - 1])
                                + f * (np.log(sw[i]) - np.log(sw[i - 1]))))
    return float("nan")


MARKS = [0.25, 0.5, 1, 2, 4, 8, 16, 32]


def decay(L, eps, gamma, seed=3):
    adj, n = square_torus(L)
    U = Nucleation(adj, n, eps, gamma)
    rng = np.random.default_rng(seed)
    d0 = mean_distance(U.adj, n)
    d_rand = mean_distance(random_regular(n, 4, np.random.default_rng(9)), n)
    marks = sorted({max(1, int(s * n)) for s in MARKS})
    rows = []
    turn = 0
    for t in range(1, marks[-1] + 1):
        for _ in range(6):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
        if t in marks:
            d = mean_distance(U.adj, n)
            rows.append(dict(sweeps=t / n, d=d,
                             phi=(d - d_rand) / (d0 - d_rand),
                             order=U.order_fraction()))
    return d0, d_rand, rows


def q_values(U):
    """The raw order-parameter distribution, for diagnosing thresholds."""
    out = []
    for i in range(U.n):
        nb = list(U.adj[i])
        t = 0
        for a in range(len(nb)):
            na = U.adj[nb[a]]
            for b in range(a + 1, len(nb)):
                t += len(na & U.adj[nb[b]]) - 1
        out.append(t)
    return np.array(out)


def q_trajectory(L, seed=3):
    """Track q against the metric, to see whether they move together."""
    adj, n = square_torus(L)
    U = Nucleation(adj, n, 1.0, 0.0)
    rng = np.random.default_rng(seed)
    d0 = mean_distance(U.adj, n)
    dr = mean_distance(random_regular(n, 4, np.random.default_rng(9)), n)
    marks = {int(s * n) for s in (0.25, 1, 4, 16, 32)}
    rows = []
    turn = 0
    for t in range(1, 32 * n + 1):
        for _ in range(6):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
        if t in marks:
            q = q_values(U)
            d = mean_distance(U.adj, n)
            rows.append((t / n, (d - dr) / (d0 - dr), float(q.mean()),
                         float((q >= 3).mean()), float((q >= 4).mean()),
                         float((q >= 6).mean())))
    return rows


def main():
    print("=" * 78)
    print("1. THE ORDER PARAMETER IS READ OFF THE GRAPH, NOT ADDED TO IT")
    print("=" * 78)
    adj, n = square_torus(24)
    U = Nucleation(adj, n, 1.0, 0.0)
    print("  square torus      fraction with s=1 : %.4f" % U.order_fraction())
    R = Nucleation(random_regular(n, 4, np.random.default_rng(1)), n, 1.0, 0.0)
    print("  random 4-regular  fraction with s=1 : %.4f" % R.order_fraction())
    print()
    print("  s_i = [q_i >= 3] where q_i counts 4-cycles through node i.")
    print("  No tuning: the two populations do not overlap at all.")
    print()

    print("=" * 78)
    print("2. HALF-LIFE vs SURFACE COUPLING   (L=20, N=400, eps=1)")
    print("=" * 78)
    print("  Nucleation predicts lifetime ~ exp(barrier)/N, so the signature")
    print("  is exponential growth in gamma, NOT growth with N.")
    print()
    print("  %8s %12s %14s %14s"
          % ("gamma", "half-life", "final order", "final phi"))
    hs = []
    for g in (0.0, 0.5, 1.0, 2.0, 4.0):
        d0, dr, rows = decay(20, 1.0, g)
        h = half_life([r["sweeps"] for r in rows], [r["phi"] for r in rows])
        hs.append((g, h))
        print("  %8.1f %12.2f %14.4f %14.3f"
              % (g, h, rows[-1]["order"], rows[-1]["phi"]))
    print()
    fin = [(g, h) for g, h in hs if not np.isnan(h)]
    if len(fin) >= 3:
        gg = np.array([f[0] for f in fin], float)
        hh = np.array([f[1] for f in fin], float)
        sl = np.polyfit(gg, np.log(hh), 1)[0]
        print("  log(half-life) vs gamma: slope %.3f per unit gamma" % sl)
        print("  -> half-life multiplies by %.2f per unit of surface cost"
              % np.exp(sl))
    print()

    print("=" * 78)
    print("3. WHY: DOES LOCAL ORDER TRACK THE METRIC AT ALL?")
    print("=" * 78)
    print("  If the surface term is to bite, the order parameter must FALL")
    print("  as the metric decays, so that an interface exists.")
    print()
    print("  %8s %8s %9s %9s %9s %9s"
          % ("sweeps", "phi", "q mean", "f(q>=3)", "f(q>=4)", "f(q>=6)"))
    for r in q_trajectory(20):
        print("  %8.2f %8.3f %9.2f %9.4f %9.4f %9.4f" % r)
    print()
    print("  q RISES from 4.05 to 9.36 while phi falls from 0.74 to 0.15.")
    print("  Local 4-cycle richness more than doubles as the metric is")
    print("  destroyed, so every threshold moves the WRONG way and no")
    print("  interface ever forms. The surface term has nothing to act on.")
    print()


if __name__ == "__main__":
    main()
