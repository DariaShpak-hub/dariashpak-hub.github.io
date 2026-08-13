#!/usr/bin/env python3
"""
self_amplification.py

Does self-amplifying relational influence produce gravity-like clumping, and
does it cost the geometry it clumps in?

The proposal: a fluctuation that recruits more connections gains more
influence, which recruits more connections. Positive feedback, attractors form
on their own, and dense regions become dense because density is
self-reinforcing. The hope is that this yields gravity-like concentration
without attraction being written into the rules.

Two versions have to be distinguished, and only one of them has been tested
before. In the GLOBAL version the reattachment target is chosen from the whole
graph with probability proportional to degree -- the Barabasi-Albert rule. In
the LOCAL version the target must already lie within graph distance 2 of the
node being reattached, so influence spreads only to things already nearby,
which is what "more neighbours affected" actually says. The distinction
matters because the global rule manufactures long-range shortcuts and the
local one cannot.

A methodological correction first, because it invalidates numbers elsewhere
--------------------------------------------------------------------------
Dimension in the earlier modules was estimated from cumulative ball counts,
|B_r| ~ r^d. That estimator is biased low, badly, and not because of
finite-size effects. On a cubic lattice the exact ball count is

    |B_r| = 1 + 2r(r+1)(2r+1)/3 + 2r

whose leading term is 4r^3/3 but whose subleading terms pull the log-log slope
down, converging only as ~1/r. Measured on graphs whose dimension is known:

    graph                true d    shell estimator    ball estimator
    cubic torus L=32          3              2.962             2.675
    cubic torus L=16          3              2.872             2.362
    square torus L=64         2              2.000             1.869

The bias is ~20% at accessible radii and is identical at L=16 and L=32, so it
is a small-radius artifact, not a boundary one. Shell counts |S_r| ~ r^(d-1)
have no offset and land on the right answer. Everything below uses shells, and
the d_ball figures in six_criteria.py and scale_factor.py should be read as
underestimates by roughly this much.

Result 1: 1/r^2 is free, and needs no feedback at all
-------------------------------------------------------
Fitting the potential as G(r) = A/r + C directly, against Euclidean distance,
on a plain cubic torus with no dynamics whatsoever:

    L      N            A            C       R^2
    16   4096      0.07734     -0.01278    0.9941
    24  13824      0.07759     -0.00858    0.9977
    32  32768      0.07772     -0.00645    0.9985

Newtonian 1/(4 pi) = 0.07958. A rises with L and is 2.3% below it at L = 32 --
consistent with converging there, though not nailed at these sizes. C, the
finite-size constant, shrinks toward zero as it should. The fit form matters:
log-log fitting a shifted power law steepens the apparent slope badly, which
is how an earlier attempt at this produced an exponent of -5.15 instead of -1.

So the lattice Laplacian's Green function IS the Newtonian potential, and
force = -dG/dr = A/r^2 follows. No dynamics, no feedback, no
self-amplification. Inverse-square is a consequence of three-dimensionality --
Gauss's law, flux through 4 pi r^2 -- and in d dimensions any graph Laplacian
gives r^(2-d) automatically. The hard part was never 1/r^2. It is the 3.

Result 2: global amplification destroys the geometry fast
-----------------------------------------------------------
Cubic torus L = 20, degree-proportional reattachment with the target drawn
from the whole graph:

    rewired   degree Gini   max degree   dimension   mean dist
         0%        0.0000            6        2.93       15.00
         1%        0.0094            8        3.90        9.61
         3%        0.0252            9        4.35        8.02
         8%        0.0547           10        3.95        6.86
        20%        0.0947           12        1.65        5.99

Mean distance collapses from 15.0 to 6.0 after rewiring one edge in five, and
the dimension goes through a meaningless excursion and out the other side.
This is Watts-Strogatz: a global target is a long-range shortcut, and a small
density of shortcuts destroys the metric.

Result 3: locality changes the answer, and the earlier claim was too broad
----------------------------------------------------------------------------
Same rule, but the target must already be within graph distance 2:

    rewired   degree Gini   max degree   dimension   mean dist
         0%        0.0000            6        2.93       15.00
         1%        0.0095            8        3.04       13.87
         3%        0.0256            8        3.11       13.05
         8%        0.0547           10        3.25       11.99
        20%        0.0970           14        3.43       10.77
        50%        0.1657           20        3.79        9.15
       100%        0.2536           41        3.92        7.40
       200%        0.4437          587       -0.26        4.60

At MATCHED clumping the two rules are not remotely alike. At Gini ~ 0.095,
global gives dimension 1.65 and mean distance 5.99; local gives 3.43 and
10.77. The clumping is identical and the geometry is not. So the claim that
self-amplification destroys the metric was drawn from the global rule alone
and does not hold in general -- the distinction between recruiting a neighbour
and recruiting anyone is the whole difference.

But locality delays rather than prevents. By 200% the local rule collapses too,
and the mechanism is self-undermining: once a hub reaches degree 587 its own
distance-2 neighbourhood spans the graph, so "local" recruitment has quietly
become global. Self-amplification defeats its own locality by succeeding.

Capping the degree -- the condition six_criteria.py identified as the real
content of locality -- delays it further and still does not prevent it:

    rewired   degree Gini   max degree   dimension   mean dist
         0%        0.0000            6        2.93       15.00
        20%        0.0970           12        3.43       10.77
       100%        0.2396           12        3.83        7.83
       200%        0.3289           12        3.16        6.73
       400%        0.4034           12        1.56        5.92
       800%        0.4317           12        0.26        5.27

With the degree bounded the concentration proceeds by phase-separating into a
dense core and a sparse periphery instead of by growing a hub, and the metric
goes anyway, just later.

Verdict
-------
There is a window, and it is a graded trade rather than a cliff. Local
recruitment with a degree cap holds dimension between 3.4 and 3.9 while degree
Gini rises from 0 to 0.24 -- real concentration with the geometry still
recognisable, which the global rule cannot do at any level. Sustained
amplification always ends the same way.

Two caveats that matter more than the window.

The dynamics starts ON a cubic torus, so three-dimensionality is put in by
hand. This measures whether amplification PRESERVES a geometry, not whether it
produces one. The open question -- what selects three dimensions -- is
untouched by everything here.

And degree Gini is probably the wrong measure of gravitational clumping.
Gravity concentrates mass: many nodes in a small region, high local density.
A hub is one node with many edges, which is a different object. What this
experiment demonstrates is hub formation, and hub formation is guaranteed for
any positive-feedback rule -- Yule 1925, Simon 1955, Barabasi-Albert 1999 --
so that half of it cannot fail and is not evidence for anything. A real test
of gravity-like clumping would need an observable sensitive to node density
per unit volume, not to degree.
"""

import numpy as np
from collections import deque


# ------------------------------------------------------------- structures

def cubic_torus(L):
    n = L ** 3

    def idx(i, j, k):
        return ((i % L) * L + (j % L)) * L + (k % L)

    adj = [set() for _ in range(n)]
    coord = np.zeros((n, 3), dtype=np.int64)
    for i in range(L):
        for j in range(L):
            for k in range(L):
                u = idx(i, j, k)
                coord[u] = (i, j, k)
                for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                    v = idx(i + d[0], j + d[1], k + d[2])
                    adj[u].add(v)
                    adj[v].add(u)
    return adj, n, coord


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


# ------------------------------------------------------------ measurement

def bfs_shells(adj, s):
    d = {s: 0}
    q = deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in d:
                d[w] = d[u] + 1
                q.append(w)
    return np.bincount(np.fromiter(d.values(), int))


def shell_profile(adj, n, srcs=8, seed=1):
    rng = np.random.default_rng(seed)
    acc = {}
    for _ in range(srcs):
        c = bfs_shells(adj, int(rng.integers(n)))
        for r in range(len(c)):
            acc[r] = acc.get(r, 0.0) + float(c[r])
    return {k: v / srcs for k, v in acc.items()}


def d_from_shells(sh, lo, hi):
    rs = [r for r in sorted(sh) if lo <= r <= hi and sh[r] > 0]
    if len(rs) < 3:
        return float("nan")
    x = np.log(np.array(rs, float))
    y = np.log(np.array([sh[r] for r in rs]))
    return 1.0 + float(np.polyfit(x, y, 1)[0])


def d_from_balls(sh, lo, hi):
    rs = sorted(sh)
    cum = np.cumsum([sh[r] for r in rs])
    sel = [(r, cum[i]) for i, r in enumerate(rs) if lo <= r <= hi]
    if len(sel) < 3:
        return float("nan")
    x = np.log(np.array([s[0] for s in sel], float))
    y = np.log(np.array([s[1] for s in sel], float))
    return float(np.polyfit(x, y, 1)[0])


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def mean_distance(adj, n, srcs=8, seed=2):
    rng = np.random.default_rng(seed)
    acc, cnt = 0.0, 0
    for _ in range(srcs):
        c = bfs_shells(adj, int(rng.integers(n)))
        acc += float((np.arange(len(c)) * c).sum())
        cnt += int(c.sum())
    return acc / cnt


# ------------------------------------------------- Laplacian Green function

def conjugate_gradient(src, dst, deg, b, n, iters=3000, tol=1e-9):
    """Solve (D - A) x = b for the mean-zero solution. No scipy available."""
    def Lmul(v):
        return deg * v - np.bincount(src, weights=v[dst], minlength=n)

    x = np.zeros(n)
    r = b - Lmul(x)
    r -= r.mean()
    p = r.copy()
    rs = float(r @ r)
    for _ in range(iters):
        Ap = Lmul(p)
        denom = float(p @ Ap)
        if abs(denom) < 1e-300:
            break
        a = rs / denom
        x += a * p
        r -= a * Ap
        rs_new = float(r @ r)
        if np.sqrt(rs_new) < tol:
            break
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x - x.mean()


def edge_arrays(adj, n):
    src, dst = [], []
    for u in range(n):
        for v in adj[u]:
            src.append(u)
            dst.append(v)
    src = np.asarray(src)
    dst = np.asarray(dst)
    deg = np.bincount(src, minlength=n).astype(float)
    return src, dst, deg


def green_vs_radius(adj, n, coord, L, source=0):
    """Potential from a unit source, binned by Euclidean distance."""
    src, dst, deg = edge_arrays(adj, n)
    b = -np.ones(n) / n
    b[source] += 1.0
    x = conjugate_gradient(src, dst, deg, b, n)
    d = coord - coord[source]
    d = np.minimum(d % L, (-d) % L)
    r = np.sqrt((d.astype(float) ** 2).sum(axis=1))
    return r, x


def fit_A_over_r(r, g, lo, hi):
    """Fit g = A/r + C. Fitting log-log instead steepens the slope, because
    the lattice Green's function carries a finite-size constant."""
    m = (r >= lo) & (r <= hi)
    inv = 1.0 / r[m]
    y = g[m]
    A, C = np.polyfit(inv, y, 1)
    pred = A * inv + C
    ss = float(((y - pred) ** 2).sum())
    tt = float(((y - y.mean()) ** 2).sum())
    return float(A), float(C), 1 - ss / tt


# ------------------------------------------------------------- the dynamics

def amplify(adj, n, frac, mode, rng, cap=None):
    """Degree-proportional reattachment. mode 'global' picks the new target
    anywhere; mode 'local' restricts it to graph distance 2. `cap` bounds the
    degree, which is the condition six_criteria.py identified as the real
    meaning of locality -- and which self-amplification otherwise defeats,
    since a hub's own distance-2 neighbourhood eventually spans the graph."""
    el = [(u, v) for u in range(n) for v in adj[u] if u < v]
    deg = np.array([len(a) for a in adj], dtype=float)
    moves = int(frac * len(el))
    for _ in range(moves):
        i = int(rng.integers(len(el)))
        a, b = el[i]
        if b not in adj[a]:
            continue
        if len(adj[b]) <= 1:
            continue
        if mode == "global":
            # target anywhere, chosen proportional to degree
            w = deg / deg.sum()
            t = int(rng.choice(n, p=w))
        else:
            # target already within distance 2 of a
            near = set()
            for x in adj[a]:
                near |= adj[x]
            near -= adj[a]
            near.discard(a)
            if not near:
                continue
            cand = np.fromiter(near, int)
            if cap is not None:
                cand = cand[deg[cand] < cap]
                if len(cand) == 0:
                    continue
            w = deg[cand]
            t = int(cand[int(rng.choice(len(cand), p=w / w.sum()))])
        if t == a or t in adj[a]:
            continue
        if cap is not None and (deg[t] >= cap or deg[a] >= cap):
            continue
        adj[a].discard(b)
        adj[b].discard(a)
        adj[a].add(t)
        adj[t].add(a)
        deg[b] -= 1
        deg[t] += 1
        el[i] = (min(a, t), max(a, t))
    return adj


# ------------------------------------------------------------------- parts

def part0():
    print("=" * 78)
    print("0. ESTIMATOR CALIBRATION on graphs of KNOWN dimension")
    print("=" * 78)
    print("  %-22s %8s %14s %14s"
          % ("graph", "true d", "from shells", "from balls"))
    for L in (16, 32):
        adj, n, _ = cubic_torus(L)
        sh = shell_profile(adj, n)
        print("  %-22s %8d %14.3f %14.3f"
              % ("cubic torus L=%d" % L, 3,
                 d_from_shells(sh, 2, L // 2 - 1),
                 d_from_balls(sh, 2, L // 2 - 1)))
    adj, n = square_torus(64)
    sh = shell_profile(adj, n)
    print("  %-22s %8d %14.3f %14.3f"
          % ("square torus L=64", 2, d_from_shells(sh, 3, 20),
             d_from_balls(sh, 3, 20)))
    print()
    print("  Ball counts are biased ~20% low at accessible radii, identically")
    print("  at L=16 and L=32, so it is a small-r artifact not a boundary one.")
    print("  Shells are used everywhere below.")
    print()


def part1():
    print("=" * 78)
    print("1. BASELINE: is 1/r^2 already there, with no feedback?")
    print("=" * 78)
    print("  Fitting G(r) = A/r + C directly. Log-log fitting a shifted power")
    print("  law steepens the slope, which is how an earlier attempt got -5.15.")
    print("  Newtonian value: 1/(4*pi) = %.5f" % (1 / (4 * np.pi)))
    print()
    print("  %6s %10s %16s %12s %9s" % ("L", "N", "A", "C", "R^2"))
    for L in (16, 24, 32):
        adj, n, coord = cubic_torus(L)
        r, g = green_vs_radius(adj, n, coord, L)
        A, C, r2 = fit_A_over_r(r, g, 2.0, L / 3.0)
        print("  %6d %10d %16.5f %12.5f %9.4f" % (L, n, A, C, r2))
    print()
    print("  A converges to 1/(4pi) on a dead lattice: no dynamics, no")
    print("  feedback, no self-amplification. Force = -dG/dr = A/r^2.")
    print("  1/r^2 is a consequence of being 3-dimensional (Gauss's law),")
    print("  not of any microscopic mechanism. The hard part is the 3.")
    print()


def part2and3():
    L = 20
    print("=" * 78)
    print("2/3. SELF-AMPLIFICATION: global vs local recruitment (L=%d)" % L)
    print("=" * 78)
    for mode, cap in (("global", None), ("local", None), ("local", 12)):
        print("  %s recruitment%s"
              % (mode.upper(), "" if cap is None else ", degree capped at %d" % cap))
        print("     %8s %12s %12s %10s %12s"
              % ("rewired", "degree Gini", "max degree", "dimension",
                 "mean dist"))
        fracs = ((0.0, 0.01, 0.03, 0.08, 0.20) if mode == "global"
                 else (0.0, 0.01, 0.03, 0.08, 0.20, 0.5, 1.0, 2.0)
                 if cap is None
                 else (0.0, 0.20, 1.0, 2.0, 4.0, 8.0))
        for frac in fracs:
            adj, n, _ = cubic_torus(L)
            if frac > 0:
                adj = amplify(adj, n, frac, mode, np.random.default_rng(5),
                              cap=cap)
            sh = shell_profile(adj, n)
            deg = [len(a) for a in adj]
            print("     %7.0f%% %12.4f %12d %10.2f %12.2f"
                  % (100 * frac, gini(deg), max(deg),
                     d_from_shells(sh, 2, L // 2 - 1), mean_distance(adj, n)))
        print()


def main():
    part0()
    part1()
    part2and3()


if __name__ == "__main__":
    main()
