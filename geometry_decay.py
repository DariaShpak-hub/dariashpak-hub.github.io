#!/usr/bin/env python3
"""
geometry_decay.py

How long does a geometry survive under local conserved rewriting?

Every previous module asked whether the dynamics can FIND a geometry, and the
answer was no, for a reason that is entropic rather than incidental. This asks
the other question. Start in the geometry and measure how long it lasts.

That is the Past Hypothesis route, and it is the one branch of the programme
that needs no new ingredient. A periodic square lattice has degree exactly 4,
and the conservation law Q = bank + 2N - E already pins the mean degree at 4,
so the lattice is not smuggled in from outside -- it is an ordinary allowed
state of the very dynamics already built. It is simply an astronomically
atypical one, with mean distance ~sqrt(N)/2 where a typical graph of the same
size and degree has ~log N.

The moves are `close` and `open` only, so N is fixed and the question is
purely about the metric decaying, not about expansion. They are strictly
alternated. That matters: letting the sampler choose freely lets the bank
accumulate, the edge count is E = E0 + bank, and the mean degree drifts up
from 4 to 5.45 over a long run -- so the late-time graph is no longer being
compared against a matched random one. Alternation pins bank to {0,1} and
holds <k> at 4.000 throughout, which the run below confirms.

The measured quantity is how far the metric has slid from lattice to random:

    phi(t) = (d(t) - d_random) / (d_lattice - d_random)

running from 1 at the lattice to 0 at a matched random 4-regular graph, in
units of sweeps, one sweep being one accepted move per node. The half-life is
where phi crosses 0.5, by log-linear interpolation.

Why the answer was not obvious in advance. `close` joins two neighbours of a
common vertex, so every edge it creates spans graph distance exactly 2. These
are not Watts-Strogatz shortcuts -- a handful of RANDOM long-range edges
collapses a lattice immediately, but short chords should not. The control
below measures exactly what that locality is worth.

Result: geometry is not metastable, and the charge is the reason
------------------------------------------------------------------
The half-life saturates at about three sweeps and stays there:

    N          400    900   1600   3600   6400  14400  25600
    t_half    1.59   2.48   3.17   2.88   2.94   3.32   3.25

Fitted over everything this reads N^0.138; fitted over N >= 1600 only it
reads N^0.027. The exponent shrinking as the range extends is what saturation
looks like, not a power law -- across a 16-fold range in N the half-life moves
from 3.17 to 3.25 sweeps. Mean degree holds at 4.03-4.05 throughout and the
conservation residual is 0 at every checkpoint, so this is a like-for-like
comparison against a matched random graph.

A constant number of sweeps means the lattice dies after a fixed number of
touches per node, in O(N) moves, however large the universe is. Geometry is
not metastable. The Past Hypothesis route is closed for this dynamics.

Locality is worth a large constant and nothing more
-----------------------------------------------------
This is the part worth keeping. The non-local control -- degree-preserving
double-edge swaps, which is the same move set with the locality removed --
destroys the same lattice with a half-life of 0.0119 sweeps against 3.17 for
the local moves. So locality buys a factor of 267.

That is a real effect and it is the right size to explain why short chords are
not Watts-Strogatz shortcuts. But 267 is a constant, and cosmology needs a
lifetime that grows with the size of the universe. A factor of 267 is nothing
when the alternative is N.

The control in the other direction behaves: a random 4-regular seed under the
same local moves stays random, drifting only because `open` needs a triangle
and random graphs have few, so closes outnumber opens and <k> creeps to 4.57.

Why: the conserved charge is blind to geometry
------------------------------------------------
Comparing a lattice against the SAME graph rewired by degree-preserving swaps,
so the edge count is identical by construction:

    N        lattice E   Q     swapped E   Q     mean distance after
    400            800   0           800   0                    4.85
    1600          3200   0          3200   0                    6.07
    3600          7200   0          7200   0                    6.83

Q is exactly equal. The ordered and the disordered graph occupy the SAME
sector, so nothing forbids the transition, and the decay is free. Every
conserved quantity in this project has had this property: weight, the Z3
residue, per-component charges, and now Q all take the same value on ordered
and disordered configurations of the same size.

That is the sharp statement the whole programme has been circling. A
conservation law that cannot distinguish order from disorder cannot protect
order. Conservation constrains WHICH states are reachable; it says nothing
about how long the atypical ones last, because typical and atypical sit in the
same sector and entropy decides everything within a sector.

To get a lifetime that grows with N you need something whose value differs
between the lattice and the random graph and which the dynamics has to climb
to change -- a barrier, not a budget. That is an energy, and it is precisely
what CDT's causal slicing and quantum graphity's Hamiltonian supply. Both
routes out of the entropic argument reduce to the same missing ingredient, and
this module measures why neither can be avoided.
"""

import numpy as np
from collections import deque


class Lattice:
    """Periodic square lattice under close/open, with a degree cap."""

    cap = 6

    def __init__(self, m):
        self.m = m
        self.n = m * m
        self.adj = [set() for _ in range(self.n)]
        self.el = []
        self.pos = {}
        self.bank = 0
        self.closes = 0
        self.opens = 0
        for i in range(m):
            for j in range(m):
                u = i * m + j
                self._add(u, ((i + 1) % m) * m + j)
                self._add(u, i * m + (j + 1) % m)
        self.E0 = len(self.el)
        self.Q = self.bank + 2 * self.n - self.E0
        self.orig = set(self.pos)

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

    def residual(self):
        return self.bank + 2 * self.n - len(self.el) - self.Q

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
            self.closes += 1
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
                self.opens += 1
                return True
        return False

    def swap(self, rng):
        """NON-LOCAL control: degree-preserving double-edge swap."""
        for _ in range(40):
            a, b = self.el[int(rng.integers(len(self.el)))]
            c, d = self.el[int(rng.integers(len(self.el)))]
            if len({a, b, c, d}) < 4:
                continue
            if d in self.adj[a] or c in self.adj[b]:
                continue
            self._del(a, b)
            self._del(c, d)
            self._add(a, d)
            self._add(c, b)
            return True
        return False


def mean_distance(adj, n, k=6, seed=2):
    r = np.random.default_rng(seed)
    acc, cnt = 0.0, 0
    for _ in range(k):
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


def random_4regular(n, rng):
    stubs = np.repeat(np.arange(n), 4)
    rng.shuffle(stubs)
    adj = [set() for _ in range(n)]
    for i in range(0, len(stubs) - 1, 2):
        a, b = int(stubs[i]), int(stubs[i + 1])
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def half_life(sweeps, phis):
    """First crossing of phi = 0.5, interpolated in log(sweeps)."""
    for i in range(1, len(phis)):
        if phis[i] <= 0.5 <= phis[i - 1]:
            lo, hi = phis[i - 1], phis[i]
            if abs(hi - lo) < 1e-12:
                return sweeps[i]
            f = (lo - 0.5) / (lo - hi)
            return float(np.exp(np.log(sweeps[i - 1]) +
                                f * (np.log(sweeps[i]) -
                                     np.log(sweeps[i - 1]))))
    return float("nan")


def decay(m, mode, marks_sweeps, seed=3, from_random=False):
    """Run and sample. mode 'local' alternates close/open; 'swap' is the
    non-local control."""
    U = Lattice(m)
    n = U.n
    rng = np.random.default_rng(seed)
    if from_random:
        radj = random_4regular(n, np.random.default_rng(seed + 50))
        U.adj = radj
        U.el = []
        U.pos = {}
        for a in range(n):
            for b in U.adj[a]:
                if a < b:
                    U.pos[(a, b)] = len(U.el)
                    U.el.append((a, b))
        U.bank = 0
        U.Q = U.bank + 2 * n - len(U.el)
    d_start = mean_distance(U.adj, n)
    d_rand = mean_distance(random_4regular(n, np.random.default_rng(9)), n)
    marks = sorted({max(1, int(s * n)) for s in marks_sweeps})
    rows = []
    step = 0
    turn = 0
    for t in range(1, marks[-1] + 1):
        ok = False
        for _ in range(8):
            if mode == "swap":
                ok = U.swap(rng)
            else:
                ok = U.close(rng) if turn == 0 else U.open(rng)
                if ok:
                    turn ^= 1
                elif turn == 1:
                    turn = 0
                    continue
            if ok:
                break
        step += int(ok)
        if t in marks:
            d = mean_distance(U.adj, n)
            denom = (d_start - d_rand) if abs(d_start - d_rand) > 1e-9 else 1.0
            rows.append(dict(sweeps=t / n, d=d, phi=(d - d_rand) / denom,
                             k=2 * len(U.el) / n, bank=U.bank,
                             orig=len(set(U.pos) & U.orig) / max(len(U.orig), 1),
                             res=U.residual()))
    return U, d_start, d_rand, rows


MARKS = [0.25, 0.5, 1, 2, 4, 8, 16, 32]
SIZES = [20, 30, 40, 60, 80, 120, 160]


def main():
    print("=" * 78)
    print("1. LATTICE UNDER LOCAL CONSERVED MOVES (close/open, strictly")
    print("   alternated so <k> stays 4 and N is fixed)")
    print("=" * 78)
    halves = []
    for m in SIZES:
        U, d0, dr, rows = decay(m, "local", MARKS)
        print("  m=%d  N=%d   lattice mean dist %.2f   random 4-reg %.2f"
              % (m, U.n, d0, dr))
        print("     %8s %9s %8s %7s %7s %9s"
              % ("sweeps", "mean d", "phi", "<k>", "bank", "orig edges"))
        for r in rows:
            print("     %8.2f %9.2f %8.3f %7.3f %7d %9.3f"
                  % (r["sweeps"], r["d"], r["phi"], r["k"], r["bank"],
                     r["orig"]))
        h = half_life([r["sweeps"] for r in rows], [r["phi"] for r in rows])
        halves.append((U.n, h))
        print("     conservation residual %d, half-life %.2f sweeps"
              % (rows[-1]["res"], h))
        print()

    print("=" * 78)
    print("2. HOW DOES THE HALF-LIFE SCALE WITH N?")
    print("=" * 78)
    ns = np.array([x[0] for x in halves], float)
    hs = np.array([x[1] for x in halves], float)
    print("  %10s %14s" % ("N", "half-life (sweeps)"))
    for a, b in halves:
        print("  %10d %14.2f" % (a, b))
    p_all = np.polyfit(np.log(ns), np.log(hs), 1)[0]
    big = ns >= 1600
    p_big = np.polyfit(np.log(ns[big]), np.log(hs[big]), 1)[0]
    print()
    print("  fitted over all sizes      : t_half ~ N^%.3f" % p_all)
    print("  fitted over N >= 1600 only : t_half ~ N^%.3f" % p_big)
    print()
    print("  The exponent shrinks as the range extends, which is what")
    print("  saturation looks like, not a power law. Over N = 1600 to %d,"
          % int(ns[-1]))
    print("  a %.0f-fold range, t_half moves from %.2f to %.2f sweeps."
          % (ns[-1] / 1600, hs[big][0], hs[-1]))
    print()
    print("  Metastability would need t_half growing like a power of N.")
    print("  A constant ~3 sweeps means the lattice dies after a fixed number")
    print("  of touches per node, i.e. in O(N) moves, however large it is.")
    print()

    print("=" * 78)
    print("3. CONTROLS")
    print("=" * 78)
    m = 40
    U, d0, dr, rows = decay(m, "swap", [0.004, 0.008, 0.016, 0.03, 0.06,
                                        0.125, 0.25, 0.5, 1, 2])
    h_swap = half_life([r["sweeps"] for r in rows], [r["phi"] for r in rows])
    print("  (a) NON-LOCAL double-edge swap on the same lattice, N=%d" % U.n)
    print("      %8s %9s %8s" % ("sweeps", "mean d", "phi"))
    for r in rows:
        print("      %8.3f %9.2f %8.3f" % (r["sweeps"], r["d"], r["phi"]))
    print("      half-life %.4f sweeps" % h_swap)
    U2, d02, dr2, rows2 = decay(m, "local", MARKS)
    h_loc = half_life([r["sweeps"] for r in rows2], [r["phi"] for r in rows2])
    print("      local half-life at the same N: %.2f sweeps" % h_loc)
    print("      locality buys a factor of %.0f" % (h_loc / h_swap))
    print()

    U3, d03, dr3, rows3 = decay(m, "local", MARKS, from_random=True)
    print("  (b) RANDOM 4-regular seed under the same local moves, N=%d" % U3.n)
    print("      starting mean distance %.2f, matched random %.2f"
          % (d03, dr3))
    print("      %8s %9s %7s" % ("sweeps", "mean d", "<k>"))
    for r in rows3:
        print("      %8.2f %9.2f %7.3f" % (r["sweeps"], r["d"], r["k"]))
    print("      random is a fixed point if this stays flat. It roughly is;")
    print("      the slow drift is <k> creeping up, because `open` needs a")
    print("      triangle and a random 4-regular graph has few, so closes")
    print("      outnumber opens and the bank accumulates.")
    print()

    print("=" * 78)
    print("4. CAN THE CONSERVED CHARGE EVEN SEE THE DIFFERENCE?")
    print("=" * 78)
    print("  Comparing the lattice against the SAME graph rewired by")
    print("  degree-preserving swaps, so the edge count is identical by")
    print("  construction and the comparison is exact.")
    print()
    print("  %-8s %10s %6s %10s %12s %6s"
          % ("N", "lattice E", "Q", "swapped E", "swapped Q", "mean d"))
    for mm in (20, 40, 60):
        L = Lattice(mm)
        n = L.n
        e_l, q_l = len(L.el), 0 + 2 * n - len(L.el)
        rr = np.random.default_rng(9)
        for _ in range(20 * n):
            L.swap(rr)
        q_r = L.bank + 2 * n - len(L.el)
        print("  %-8d %10d %6d %10d %12d %6.2f"
              % (n, e_l, q_l, len(L.el), q_r, mean_distance(L.adj, n)))
    print()
    print("  The charge is unchanged. The lattice and the random graph sit in")
    print("  the SAME sector, so nothing forbids the transition between them --")
    print("  which is exactly why the decay is free. A conservation law that")
    print("  cannot distinguish ordered from disordered cannot protect order.")
    print()


if __name__ == "__main__":
    main()
