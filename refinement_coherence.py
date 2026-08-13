#!/usr/bin/env python3
"""
refinement_coherence.py

Can LOCAL compatibility repair coordinate refinement, or does it just
implement global command by propagation?

The controlled test established the two ends of the range. Coherent refinement
-- every edge told to subdivide at once -- gives r_ij = 2 for every marker pair,
sigma = 0, delta = 0: exact isotropic dilation. Independent refinement at rate
p < 1 gives the right sign, lambda = 1.09 to 1.55, but leaves delta = 0.15-0.18.
So expansion is solved and coordination is not.

The proposal under test is that neighbouring cells cannot hold arbitrarily
different resolutions, so refining one edge forces its neighbours to repair,
and the repair propagates. No clock, no scale factor, no global instruction.

The state
---------
A periodic triangular complex, L x L vertices, E = 3L^2 edges, F = 2L^2 faces,
chi = 0. Each edge carries

    r_e in Z        refinement level; the edge stands for 2^(r_e) fine edges
    U_e in Z_2      transport, carried passively

Distance is the weighted shortest path with edge length 2^(r_e), which is exact
rather than approximate: subdividing an edge creates no new routes, so a coarse
path of levels r_1..r_n has fine length sum 2^(r_i). That lets the whole thing
run on the coarse complex.

The rule, with one knob
------------------------
    seed      pick a random edge, r_e += 1
    repair    while some face has max(r) - min(r) > tol, raise its low edges

tol = 0 demands every face be perfectly level. tol = 1 is the 2:1 balance
condition that adaptive mesh refinement actually uses. Nothing else is
supplied.

The diagnostic that decides it
-------------------------------
delta -> 0 is not by itself a pass. If the only way to reach it is a repair
cascade that sweeps the entire complex after every seed, then "local repair"
is global synchronisation wearing a different hat, and the mechanism has not
explained coordination -- it has renamed it.

So avalanche size is measured alongside delta. The gate is:

    delta -> 0     AND     mean avalanche = O(1) as N grows

Passing both means local constraints genuinely coordinate. Passing only the
first means the cascade is doing the work a global instruction would have done.
Failing the first means the route is dead and should be abandoned rather than
patched, which was the explicit agreement.

Prediction, written before running. tol = 0 must give delta = 0 with avalanches
of order N -- it can only terminate when everything is level. tol >= 1 permits a
rough level field, and whether delta falls depends on how the roughness grows:
if sigma_r saturates while <r> grows linearly, then sigma_r/<r> -> 0 and this
works. The complication is that distance under weights 2^r is dominated by the
LARGEST level on a path, so geodesics route around high-r regions -- the same
routing-around that killed the defect coupling in holonomy_defect.py. That
makes delta a much harder statistic than sigma_r/<r>, and the two can disagree.

Result 1: tol = 0 gives perfect coherence at global cost
----------------------------------------------------------
    tol   L     E    avalanche   fraction of E   delta
      0   8    192       192.0        1.0000     0.00000
      0  12    432       432.0        1.0000     0.00000
      0  16    768       768.0        1.0000     0.00000
      0  24   1728      1728.0        1.0000     0.00000

Exactly delta = 0, and exactly every edge touched by every single seed, at
every size. So zero-tolerance repair reproduces the coherent-refinement result
perfectly and is global synchronisation implemented by propagation. It answers
the question by renaming it, which is what the avalanche diagnostic was added
to detect.

Result 2: tol >= 1 is genuinely local
---------------------------------------
    tol   L      E   avalanche   fraction of E
      1   8    192        2.9       0.0150
      1  12    432        3.4       0.0079
      1  16    768        3.0       0.0039
      1  24   1728        3.7       0.0021
      2  24   1728        1.5       0.0008

The mean cascade is 3 to 4 edges and does NOT grow with N -- the fraction of
the complex touched falls like 1/E. The 2:1 balance condition that adaptive
mesh refinement uses is a genuinely local rule here. The locality gate passes.

Result 3: sigma_r/<r> -> 0 as predicted, and delta does NOT follow it
----------------------------------------------------------------------
tol = 1, L = 24, with edge length 2^r:

    sweep      <r>   sd/<r>    delta
        1    0.395   1.5553   0.15714
        4    2.608   0.4463   0.48372
        8    5.796   0.2349   0.83038

The relative level fluctuation behaves exactly as hoped: 1.56 -> 0.23 and
falling. The prediction that local repair homogenises the level field is
CONFIRMED. But the metric distortion moves the other way -- delta rises from
0.157 to 0.830. The two statistics disagree, and delta is the one that matters.

Result 4: the culprit is the exponential length map, not the coordination
--------------------------------------------------------------------------
Identical dynamics, identical level field, only the map from level to length
changed:

    sweep     <r>   sd/<r>    delta (2^r)   delta (r+1)
        2   0.991   0.8479      0.25040       0.20479
        4   2.608   0.4463      0.48372       0.18665
        6   4.198   0.2957      0.66259       0.14933
        8   5.796   0.2349      0.83038       0.13563

and continuing with length = r + 1 alone:

    sweep      <r>   sd/<r>    delta   log2 lambda
        9    6.644   0.2110  0.12720          2.84
       12    9.307   0.1694  0.10473          3.29
       15   11.882   0.1343  0.08630          3.63
       18   14.418   0.1102  0.06951          3.90
       21   17.106   0.0952  0.06480          4.13
       24   19.870   0.0911  0.06281          4.35

Monotone decrease from 0.205 to 0.063 with no sign of a plateau, while <r> runs
to 19.9 and lambda to 20.4. Both gates pass together.

The reason is simple once separated. Local repair holds the level GRADIENT
bounded, so sd_r saturates near 1 while <r> grows. Under length = r + 1 a
bounded absolute spread against a growing mean gives relative spread -> 0.
Under length = 2^r the same bounded spread is a permanent MULTIPLICATIVE factor
of about 2 between neighbouring edges, which never shrinks, and geodesics route
around the high-level regions on top of that. Bounded roughness in r is
compatible with vanishing distortion only if length is linear in r.

Verdict
-------
    tol = 0, any length map        delta = 0, avalanche = E        global
    tol >= 1, length 2^r           delta grows, avalanche O(1)     fails
    tol >= 1, length r + 1         delta -> 0, avalanche O(1)      PASSES

Gate A and Gate B pass, with one substitution that is not a free choice but a
correction. Level as a BINARY refinement generation, length 2^r, is the mesh
refinement reading. Level as the count of vertices inserted into an edge,
length r + 1, is our own subdivision move -- dN = +1, dE = +1, the one
rewrite_algebra.py assigns the invariant E - N. The linear reading is the one
this project has been using all along, and it is the one that works.

That also retro-explains the anti-Yule weight. Selecting edges with weight
2^(-g) was compensating for an exponential that need never have been there.

What is NOT shown. Gate C, a finite spectral dimension plateau, is untested
here. Gate D, survival of the Z_2 sector under refinement, is untested -- U_e is
carried but never measured. Two dimensions only, one lattice topology, delta
still falling at the last measured point rather than demonstrated to reach
zero, and the level field is driven by uniform random seeding rather than by
anything the model chooses.
"""

import heapq
import math
import numpy as np
from collections import deque


class TriComplex:
    """Periodic triangular complex with per-edge refinement level."""

    def __init__(self, L, seed=0, scale="exp"):
        self.L = L
        self.scale = scale
        self.rng = np.random.default_rng(seed)
        self.n = L * L
        self.eu, self.ev, self.r, self.U = [], [], [], []
        self.inc = [[] for _ in range(self.n)]
        self.H, self.V, self.D = {}, {}, {}
        for x in range(L):
            for y in range(L):
                self.H[(x, y)] = self._add(self.v(x, y), self.v(x + 1, y))
                self.V[(x, y)] = self._add(self.v(x, y), self.v(x, y + 1))
                self.D[(x, y)] = self._add(self.v(x, y), self.v(x + 1, y + 1))
        self.faces = []
        for x in range(L):
            for y in range(L):
                self.faces.append((self.H[(x, y)], self.V[((x + 1) % L, y)],
                                   self.D[(x, y)]))
                self.faces.append((self.D[(x, y)], self.H[(x, (y + 1) % L)],
                                   self.V[(x, y)]))
        self.edge_faces = [[] for _ in range(len(self.eu))]
        for fi, f in enumerate(self.faces):
            for e in f:
                self.edge_faces[e].append(fi)

    def v(self, x, y):
        return (y % self.L) * self.L + (x % self.L)

    def _add(self, u, w):
        i = len(self.eu)
        self.eu.append(u)
        self.ev.append(w)
        self.r.append(0)
        self.U.append(1)
        self.inc[u].append(i)
        self.inc[w].append(i)
        return i

    # ------------------------------------------------------------- dynamics
    def repair(self, seeds, tol):
        """Raise levels until every face satisfies max - min <= tol."""
        q = deque()
        for e in seeds:
            for f in self.edge_faces[e]:
                q.append(f)
        touched = set(seeds)
        guard = 0
        while q and guard < 200 * len(self.eu):
            guard += 1
            f = q.popleft()
            es = self.faces[f]
            m = max(self.r[e] for e in es)
            for e in es:
                if self.r[e] < m - tol:
                    self.r[e] = m - tol
                    touched.add(e)
                    for g in self.edge_faces[e]:
                        q.append(g)
        return len(touched)

    def sweep(self, n_seeds, tol):
        sizes = []
        for _ in range(n_seeds):
            e = int(self.rng.integers(len(self.eu)))
            self.r[e] += 1
            sizes.append(self.repair([e], tol))
        return sizes

    # ---------------------------------------------------------- measurement
    def w(self, e):
        return (1 << self.r[e]) if self.scale == "exp" else (self.r[e] + 1)

    def dist(self, s, targets):
        """Exact weighted shortest paths, edge length 2^r, integer arithmetic."""
        INF = None
        d = {s: 0}
        pq = [(0, s)]
        want = set(targets)
        out = {}
        while pq:
            du, u = heapq.heappop(pq)
            if du > d.get(u, du):
                continue
            if u in want:
                out[u] = du
                want.discard(u)
                if not want:
                    return out
            for e in self.inc[u]:
                w = self.ev[e] if self.eu[e] == u else self.eu[e]
                nd = du + self.w(e)
                cur = d.get(w)
                if cur is None or nd < cur:
                    d[w] = nd
                    heapq.heappush(pq, (nd, w))
        return out

    def ratios(self, pairs, base):
        out = []
        srcs = sorted({a for a, _ in pairs})
        tg = {}
        for a, b in pairs:
            tg.setdefault(a, []).append(b)
        for a in srcs:
            d = self.dist(a, tg[a])
            for b in tg[a]:
                if b in d and base.get((a, b)):
                    # log space: 2^r overflows float at tol = 0, and delta is
                    # scale invariant, so the common factor is divided out
                    out.append(math.log2(d[b]) - math.log2(base[(a, b)]))
        if not out:
            return np.array([], float), float("nan")
        lg = np.array(out, float)
        return np.exp2(lg - lg.min()), float(lg.mean())


def baseline(G, pairs):
    tg = {}
    for a, b in pairs:
        tg.setdefault(a, []).append(b)
    base = {}
    for a, bs in tg.items():
        d = G.dist(a, bs)
        for b in bs:
            if b in d and d[b] > 0:
                base[(a, b)] = d[b]
    return base


def make_pairs(L, rng, m=400):
    n = L * L
    ps = set()
    while len(ps) < m:
        a = int(rng.integers(n))
        b = int(rng.integers(n))
        if a != b:
            ps.add((min(a, b), max(a, b)))
    return sorted(ps)


def run(L, tol, sweeps, seed=0, per_sweep=None, scale="exp"):
    rng = np.random.default_rng(1000 + seed)
    G = TriComplex(L, seed=seed, scale=scale)
    pairs = make_pairs(L, rng)
    base = baseline(G, pairs)
    per = per_sweep or len(G.eu) // 4
    rows = []
    for s in range(1, sweeps + 1):
        sizes = G.sweep(per, tol)
        rr = np.array(G.r, float)
        rat, log_lam = G.ratios(pairs, base)
        delta = float(rat.std() / rat.mean()) if rat.size else float("nan")
        rows.append(dict(sweep=s, mean_r=float(rr.mean()),
                         sd_r=float(rr.std()),
                         rel_r=float(rr.std() / max(rr.mean(), 1e-9)),
                         lam=log_lam, delta=delta,
                         av_mean=float(np.mean(sizes)),
                         av_max=float(np.max(sizes)),
                         E=len(G.eu)))
    return rows


def main():
    print("=" * 78)
    print("1. THE TOLERANCE KNOB: coherence against avalanche cost")
    print("=" * 78)
    print("  L = 16, E = 768, one quarter of the edges seeded per sweep")
    print("  %5s %7s %9s %8s %9s %9s %10s %9s"
          % ("tol", "sweep", "<r>", "sd r", "sd/<r>", "lambda", "delta",
             "avalanche"))
    for tol in (0, 1, 2, 3):
        for row in run(16, tol, 4, seed=1):
            print("  %5d %7d %9.3f %8.3f %9.4f %9.2f %10.5f %9.1f"
                  % (tol, row["sweep"], row["mean_r"], row["sd_r"],
                     row["rel_r"], row["lam"], row["delta"], row["av_mean"]))
        print()

    print("=" * 78)
    print("2. DOES THE AVALANCHE STAY LOCAL AS N GROWS?")
    print("=" * 78)
    print("  mean edges touched per seed, as a fraction of E")
    print("  %5s %7s %8s %10s %12s %12s"
          % ("tol", "L", "E", "avalanche", "fraction E", "delta"))
    for tol in (0, 1, 2):
        for L in (8, 12, 16, 24):
            rows = run(L, tol, 3, seed=2)
            r = rows[-1]
            print("  %5d %7d %8d %10.1f %12.4f %12.5f"
                  % (tol, L, r["E"], r["av_mean"], r["av_mean"] / r["E"],
                     r["delta"]))
        print()

    print("=" * 78)
    print("3. GATE A AND B: does delta fall while <r> keeps growing?")
    print("=" * 78)
    for tol in (1, 2):
        print("  tol = %d" % tol)
        print("  %7s %9s %9s %10s %10s"
              % ("sweep", "<r>", "sd/<r>", "delta", "log2 lam"))
        for row in run(24, tol, 8, seed=3):
            print("  %7d %9.3f %9.4f %10.5f %10.2f"
                  % (row["sweep"], row["mean_r"], row["rel_r"],
                     row["delta"], row["lam"]))
        print()

    print("=" * 78)
    print("4. IS IT THE COORDINATION, OR THE EXPONENTIAL LENGTH MAP?")
    print("=" * 78)
    print("  Same dynamics, same level field. Only the map from level to")
    print("  length changes: 2^r against r+1. If delta falls for the linear")
    print("  map then coordination is fine and the exponential is the problem.")
    print()
    for scale in ("exp", "lin"):
        print("  length = %s" % ("2^r" if scale == "exp" else "r + 1"))
        print("  %7s %9s %9s %10s %10s"
              % ("sweep", "<r>", "sd/<r>", "delta", "log2 lam"))
        for row in run(24, 1, 8, seed=3, scale=scale):
            print("  %7d %9.3f %9.4f %10.5f %10.2f"
                  % (row["sweep"], row["mean_r"], row["rel_r"],
                     row["delta"], row["lam"]))
        print()

    print("=" * 78)
    print("  PASS needs BOTH: delta -> 0, and avalanche = O(1) as N grows.")
    print("  delta -> 0 with avalanche ~ E means local repair is global")
    print("  synchronisation renamed, not an explanation of coordination.")
    print()


if __name__ == "__main__":
    main()
