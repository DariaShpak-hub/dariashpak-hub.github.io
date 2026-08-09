#!/usr/bin/env python3
"""
dimension_attractor.py

Is three a selected dimension, or a free parameter?

Everything downstream of a stable 3D substrate now has a measured answer in
this project: 1/r^2 comes from Gauss's law (self_amplification.py), clumping
comes from attraction, extended structure needs the long range that 1/r^2
supplies (matter_clumping.py). What none of it explains is the 3.

Nobody else derives it either, without either string theory or an anthropic
step. Ehrenfest's 1917 argument -- stable orbits and stable atoms need d <= 3,
sharp wave propagation needs odd d >= 3 -- selects 3 but does not produce it.
Brandenberger-Vafa string gas is the one genuine dynamical derivation, and it
needs string theory. CDT puts 4 in through its simplices, causal sets sprinkle
into a manifold, quantum graphity sets d by its Hamiltonian.

And this project's own measurements point the same way: in six_criteria.py the
degree cap behaved as a DIAL on dimension, cap 6 giving d_H = 4.03 and cap 12
giving 6.34. Dimension came out tunable, with nothing selecting a value.

So the question is whether the dynamics has any preference at all. Tori of
dimension 2, 3 and 4 are all allowed states of the same conserved local
dynamics -- close/open, N fixed, as in geometry_decay.py -- and they match
exactly at N = 4096, since 64^2 = 16^3 = 8^4. Run each and ask:

    if all dimensions decay at the same rate to the same place, dimension is
    a free parameter and no rule in this family will ever select 3

    if one decays more slowly, that is a selection mechanism and its value is
    the prediction

    if low-d seeds drift up and high-d seeds drift down to a common value,
    that value is an attractor

There is a confound to control. A d-torus has degree 2d, so dimension and
degree move together and a rate difference could be either. The control is a
2-torus with diagonals -- the king's graph, dimension 2 but degree 8 -- which
matches the 4-torus in degree while differing in dimension.

Honest limit on the measurement
---------------------------------
The shell estimator needs a scaling range, and at fixed N the side length
falls as N^(1/d), so high dimensions run out of room. Measured on PRISTINE
tori at N ~ 16000:

    d      2       3       4       5       6
    est  2.000   2.946   3.673   4.173   4.432

Exact at d = 2, good at d = 3, already 8% low at d = 4, and useless beyond.
So the scan stops at d = 4, and everything below is reported as a change from
each lattice's own pristine value rather than as an absolute dimension.

Result: dimension is not free, but nothing selects three
----------------------------------------------------------
Pristine calibration at N = 4096, then decay under identical dynamics:

    seed            degree   pristine d   half-life   per edge
    2-torus              4        2.000        2.35      1.173
    3-torus              6        2.921        1.76      0.587
    4-torus              8        3.548        1.30      0.324
    2-torus + diagonals  8        2.000        3.12      0.781

The per-edge column matters: a sweep is N moves but a degree-2d graph has dN
edges, so denser graphs get less attention per edge. The 4-torus dies fastest
despite each of its edges being touched half as often as the 2-torus's.

Dimension is NOT a free parameter with respect to survival. The half-life
falls monotonically with d -- 1.173, 0.587, 0.324 per edge for d = 2, 3, 4 --
and the control separates the causes cleanly:

    same dimension, degree 4 vs 8:      1.173 vs 0.781    ratio 1.5
    same degree, dimension 2 vs 4:      0.781 vs 0.324    ratio 2.4

so both matter and dimension matters more. But the preference runs the wrong
way for cosmology: LOWER dimension survives longer. Extrapolated, this rule
selects d = 1, which is exactly where scale_factor.py's pure subdivision
ended up. Three is nowhere on the curve.

And there is no attractor. The measured dimension drifts UP from every seed
before the metric fails -- 2.00 to 2.29, 2.92 to 3.93, 3.55 to 4.94 -- and
none of them converge. Low-d seeds do not rise toward 3 and high-d seeds do
not fall toward it. Everything inflates toward dimensionlessness at a rate set
by where it started. (Late-time values like 0.107 and -1.420 are artifacts:
once the metric is gone the diameter collapses and there is no scaling window
left to fit.)

What this closes
------------------
This answers the question it was built for, in the negative and for the whole
family. Local conserved rewriting does have a dimensional preference, so the
weaker worry -- that dimension is simply invisible to the dynamics -- is
wrong. But the preference is monotone toward low d, with no stationary point,
so no rule of this kind will select three by dynamics alone. Something outside
this family has to fix it: a compactification mechanism as in string gas
cosmology, a causal structure as in CDT, or an anthropic step as in Ehrenfest.

Worth stating plainly, since it is the terminus of a long thread: the
programme reached a real obstruction rather than running out of ideas. Every
downstream piece works once three dimensions exist -- 1/r^2 by Gauss's law,
clumping by attraction, structure by long range -- and the three itself does
not come from anywhere inside the rules.
"""

import numpy as np
from collections import deque


# ------------------------------------------------------------- structures

def hypertorus(d, L):
    """Periodic d-dimensional hypercubic lattice, degree 2d."""
    N = L ** d
    strides = [L ** i for i in range(d)]
    adj = [set() for _ in range(N)]
    for u in range(N):
        for i in range(d):
            c = (u // strides[i]) % L
            v = u + (((c + 1) % L) - c) * strides[i]
            adj[u].add(v)
            adj[v].add(u)
    return adj, N


def kings_torus(L):
    """2D torus with diagonals: dimension 2, degree 8 -- the degree control."""
    N = L * L
    adj = [set() for _ in range(N)]
    for i in range(L):
        for j in range(L):
            u = i * L + j
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    v = ((i + di) % L) * L + (j + dj) % L
                    adj[u].add(v)
                    adj[v].add(u)
    return adj, N


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


# ------------------------------------------------------------ measurement

def shells_from(adj, s):
    dist = {s: 0}
    q = deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                q.append(w)
    return np.bincount(np.fromiter(dist.values(), int))


def profile(adj, n, srcs=6, seed=1):
    rng = np.random.default_rng(seed)
    acc = {}
    tot, cnt = 0.0, 0
    for _ in range(srcs):
        c = shells_from(adj, int(rng.integers(n)))
        tot += float((np.arange(len(c)) * c).sum())
        cnt += int(c.sum())
        for r in range(len(c)):
            acc[r] = acc.get(r, 0.0) + float(c[r])
    return {k: v / srcs for k, v in acc.items()}, tot / cnt


def d_shells(sh, lo, hi):
    rs = [r for r in sorted(sh) if lo <= r <= hi and sh[r] > 0]
    if len(rs) < 3:
        return float("nan")
    x = np.log(np.array(rs, float))
    y = np.log(np.array([sh[r] for r in rs]))
    return 1.0 + float(np.polyfit(x, y, 1)[0])


# --------------------------------------------------------------- dynamics

class Substrate:
    """Arbitrary seed graph under the close/open conserved dynamics."""

    def __init__(self, adj, n, cap):
        self.n = n
        self.cap = cap
        self.adj = [set(a) for a in adj]
        self.el = []
        self.pos = {}
        for u in range(n):
            for v in self.adj[u]:
                if u < v:
                    self.pos[(u, v)] = len(self.el)
                    self.el.append((u, v))
        self.bank = 0

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

    def close(self, rng):
        for _ in range(40):
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
            self.bank += 1
            return True
        return False

    def open(self, rng):
        if self.bank < 1:
            return False
        for _ in range(40):
            a, b = self.el[int(rng.integers(len(self.el)))]
            if self.adj[a] & self.adj[b]:
                self._del(a, b)
                self.bank -= 1
                return True
        return False


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


def run(adj0, n, deg, hi_window, seed=3):
    U = Substrate(adj0, n, cap=deg + 4)
    rng = np.random.default_rng(seed)
    sh0, md0 = profile(U.adj, n)
    d0 = d_shells(sh0, 2, hi_window)
    rnd = random_regular(n, deg, np.random.default_rng(9))
    _, md_rand = profile(rnd, n)
    marks = sorted({max(1, int(s * n)) for s in MARKS})
    rows = []
    turn = 0
    for t in range(1, marks[-1] + 1):
        for _ in range(8):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
        if t in marks:
            sh, md = profile(U.adj, n)
            rows.append(dict(sweeps=t / n, d=d_shells(sh, 2, hi_window),
                             md=md,
                             phi=(md - md_rand) / (md0 - md_rand)))
    return d0, md0, md_rand, rows


def main():
    print("=" * 78)
    print("1. PRISTINE CALIBRATION at N = 4096  (64^2 = 16^3 = 8^4)")
    print("=" * 78)
    seeds = []
    for d, L in ((2, 64), (3, 16), (4, 8)):
        adj, n = hypertorus(d, L)
        sh, md = profile(adj, n)
        diam = max(sh)
        hi = max(4, diam // 3)
        est = d_shells(sh, 2, hi)
        seeds.append(("%d-torus" % d, adj, n, 2 * d, hi, est))
        print("  %-14s N=%d  degree %2d  diameter %3d  window [2,%2d]  "
              "d_est %.3f" % ("%d-torus" % d, n, 2 * d, diam, hi, est))
    adj, n = kings_torus(64)
    sh, md = profile(adj, n)
    hi = max(4, max(sh) // 3)
    est = d_shells(sh, 2, hi)
    seeds.append(("2-torus+diag", adj, n, 8, hi, est))
    print("  %-14s N=%d  degree %2d  diameter %3d  window [2,%2d]  d_est %.3f"
          % ("2-torus+diag", n, 8, max(sh), hi, est))
    print()
    print("  The last row is the degree control: dimension 2 like the")
    print("  2-torus, degree 8 like the 4-torus.")
    print()

    print("=" * 78)
    print("2. DECAY UNDER THE SAME LOCAL CONSERVED DYNAMICS")
    print("=" * 78)
    summary = []
    for name, adj, n, deg, hi, d0 in seeds:
        d_start, md0, md_rand, rows = run(adj, n, deg, hi)
        print("  %s   pristine d_est %.3f, mean dist %.2f, random %.2f"
              % (name, d_start, md0, md_rand))
        print("     %8s %9s %10s %9s" % ("sweeps", "d_est", "mean dist", "phi"))
        for r in rows:
            print("     %8.2f %9.3f %10.2f %9.3f"
                  % (r["sweeps"], r["d"], r["md"], r["phi"]))
        h = half_life([r["sweeps"] for r in rows], [r["phi"] for r in rows])
        summary.append((name, deg, d_start, h, rows[-1]["d"]))
        print("     half-life %.2f sweeps" % h)
        print()

    print("=" * 78)
    print("3. IS ANY DIMENSION PREFERRED?")
    print("=" * 78)
    print("  A sweep is N moves, but a degree-2d graph has dN edges, so a")
    print("  sweep touches a smaller fraction of a denser graph. The per-edge")
    print("  column is the fair comparison.")
    print()
    print("  %-16s %8s %12s %12s %14s"
          % ("seed", "degree", "pristine d", "half-life", "per edge"))
    for name, deg, d0, h, dl in summary:
        print("  %-16s %8d %12.3f %12.2f %14.3f"
              % (name, deg, d0, h, h * 2.0 / deg))
    print()
    print("  Rows 1 and 4: same dimension (2), degree 4 vs 8.")
    print("  Rows 3 and 4: same degree (8), dimension 4 vs 2.")
    print("  Whichever pair differs more is what the dynamics responds to.")
    print()
    print("  Late-time d_est above (0.107, -1.420) is meaningless: once the")
    print("  metric is gone the diameter collapses and no scaling window")
    print("  exists. Only the early trajectory and the half-life are usable.")
    print()


if __name__ == "__main__":
    main()
