#!/usr/bin/env python3
"""
scale_factor.py

Can a scale factor be read off the graph itself, without assuming a dimension
and without a clock?

A scale factor is not "the graph got bigger". It is a very specific claim:
that there exist markers whose mutual separations all grow by the SAME factor,
so that one number a suffices to describe every distance at once. Cosmology
gets this from comoving coordinates -- labels that do not move while proper
distance between them stretches uniformly. A graph has the same ingredients:
subdivide never deletes a node, so the seed nodes persist forever and are
exactly comoving markers, and graph distance between them is proper distance.

    a_ij  =  d(i, j) / d_0(i, j)        for each pair of seed nodes

If a scale factor exists, the a_ij all coincide. If they scatter, expansion is
inhomogeneous and no single a describes it. That is the first test, and it is
a pure graph property -- no schedule, no rates, no clock.

The second test is dimension, measured rather than assumed, two independent
ways that have to agree:

    d_ball   from local volume growth,  |B_r(x)| ~ r^d      (static, intrinsic)
    d_exp    from the expansion itself, a ~ N^(1/d)         (kinematic)

For a d-dimensional expanding space these are the same d, because volume goes
as a^d. If they disagree, the graph is not behaving like an expanding space of
any dimension -- it is just accumulating nodes.

Using N as the clock is not smuggling time: it is the causal-set convention
that volume IS time, and it is the only monotone the growth supplies. But it
does mean a ~ N^(1/d) has no content on its own -- it is close to a tautology.
The content is entirely in whether d_exp agrees with d_ball, and in whether
the a_ij agree with each other. Both can fail, and both do.

Three growth rules are compared, all starting from small seeds:

    subdivide     {x,y} -> {x,z},{z,y}                     pure refinement
    subdivide+    the same, plus a shortcut from z          refinement + links
    1->3 move     insert a node inside a triangle           the 2D DT move

Result 1: pure subdivide gives a real scale factor, of a 1D universe
--------------------------------------------------------------------
The expansion is real and the construction works. Averaged over 8 histories
on the 6-cycle, mean a goes 347, 700, 1399, 2796, 5593 across N = 2500 to
40000: exactly linear, a ~ N^1.000.

But it is not homogeneous, and it never becomes homogeneous:

    N        2500    5000   10000   20000   40000
    spread   75.8%   75.5%   75.6%   75.6%   75.4%
    min/max  0.068   0.068   0.068   0.069   0.070

The closest pair of markers expands 15 times less than the furthest, and that
ratio is frozen -- identical to three digits over a 16-fold increase in N. The
cause is a Yule process: an original edge already cut into k pieces has k
chances of being picked next, so chain lengths grow proportionally to
themselves and their shares follow a Dirichlet distribution, which does not
concentrate. Whatever allotment the first few dozen steps happen to make is
permanent. This is a fixed-shear expansion, closer to Kasner than to FRW, and
no amount of running fixes it because there is no mechanism that could.

The dimension is consistent, and it is one:

    d_exp  = 1.00     from a ~ N^1.000
    d_ball = 0.99     converged, drifting only 0.99 -> 1.00 over the range
    mean distance ~ N^1.000, power-law residual 0.0000 against 0.4995 for log

Subdividing edges cannot produce anything but threads. It is a refinement: it
changes resolution, never topology. Whatever dimension the seed had is the
dimension you keep.

Result 2: the 1->3 move grows volume with no expansion whatsoever
------------------------------------------------------------------
The manifestly two-dimensional move never deletes an edge, so the four K4
markers stay mutually adjacent forever:

    a = 1.000 at every N, at every pair, exactly, with 0.0% spread.

Volume grows 10000-fold; comoving distance does not move at all. This is the
configuration-space-versus-space distinction made measurable rather than
rhetorical, and it is the sharpest single number here.

The graph is also far too compact for a dimension to exist: mean distance goes
only 5.13 -> 6.22 while N goes 2500 -> 40000, diameter 10 -> 13, and no
scaling window survives, so d_ball is unmeasurable. Power-law and logarithmic
fits are statistically indistinguishable (residuals 0.0019 against 0.0021),
and the nominal d_H = 14.84 the power fit returns is meaningless.

Worth stating because it corrects an obvious guess: this is NOT the famous
d_H = 4 of two-dimensional dynamical triangulations. The 1->3 move alone
generates *stacked* triangulations -- Apollonian networks, which are
small-world -- not uniform random ones. Getting the d_H = 4 ensemble requires
flip moves to equilibrate, which this growth rule does not have. Growing a
triangulation and sampling a triangulation are different things.

Result 3: shortcuts give a fake dimension, and the naive test passes it
------------------------------------------------------------------------
This one nearly went in as a success. Adding one shortcut per birth gives

    d_exp  = 3.76      from a ~ N^0.266
    d_ball = 2.52, 2.30, 2.75, 2.80, 3.07   over N = 2500 .. 40000

and a first pass called those consistent at about 3 -- a three-dimensional
emergent space, which is exactly the result one wants and should therefore
distrust. It is an artifact. Two checks kill it:

    mean distance fits a logarithm better than a power law
        (residual 0.0245 against 0.0420)
    d_ball never converges, drifting steadily upward 2.52 -> 3.07

Both say small-world. |B_r| grows exponentially, so there is no power law to
fit, and the log-log slope simply tracks the fitting window, which grows with
N. A dimension that has not stopped moving is not a dimension. Measuring at
one size would have reported d = 3 with a straight face.

Verdict
-------
Subdivide does produce a scale factor. The comoving construction works, the
markers persist, no clock is needed, and a(N) is a genuine graph observable --
that much of the programme is sound. But the universe it produces is
one-dimensional and permanently anisotropic at 75%, and no amount of
subdividing changes either, because refinement preserves topology by
construction.

The two obvious repairs fail in opposite directions. A move that adds no
adjacencies (1->3) gives unbounded volume and zero expansion. A move that adds
random adjacencies (shortcut) gives expansion with no dimension at all. What
is missing in both cases is the same thing: new adjacencies have to be local
in the geometry that already exists, and neither rule has any notion of
locality to respect.

So dimension is not emerging from growth here. It is inherited from the seed,
or destroyed. Any rule that will give d = 3 has to know about 3 somehow, which
is the obstruction CDT hit and solved only by imposing causal structure by
hand. That is a substantive obstruction, not a missing detail.

Consequence for the cosmology: w_eff still has nothing to be a prediction
about. Not because the clock is arbitrary -- the comoving construction removes
that objection -- but because there is no a(N) whose dimension was not simply
put in at the seed.
"""

import numpy as np
from collections import defaultdict, deque


# ----------------------------------------------------------------- growth

def grow_subdivide(seed_edges, n_seed, target, rng, shortcut=0.0):
    """{x,y} -> {x,z},{z,y}. With prob `shortcut`, also link z to a
    neighbour of x, which adds an adjacency refinement never would."""
    adj = defaultdict(set)
    edges = list(seed_edges)
    for x, y in edges:
        adj[x].add(y)
        adj[y].add(x)
    n = n_seed
    while n < target:
        k = int(rng.integers(len(edges)))
        x, y = edges[k]
        z = n
        n += 1
        adj[x].discard(y)
        adj[y].discard(x)
        adj[x].add(z); adj[z].add(x)
        adj[y].add(z); adj[z].add(y)
        edges[k] = (x, z)
        edges.append((z, y))
        if shortcut and rng.random() < shortcut:
            cand = [w for w in adj[x] if w != z]
            if cand:
                w = cand[int(rng.integers(len(cand)))]
                if w not in adj[z]:
                    adj[z].add(w); adj[w].add(z)
                    edges.append((z, w))
    return adj, n


def grow_triangulate(target, rng):
    """The 2D dynamical-triangulation 1->3 move, seeded on the K4 sphere."""
    adj = defaultdict(set)
    for a in range(4):
        for b in range(4):
            if a != b:
                adj[a].add(b)
    tris = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    n = 4
    while n < target:
        k = int(rng.integers(len(tris)))
        a, b, c = tris[k]
        v = n
        n += 1
        for u in (a, b, c):
            adj[u].add(v); adj[v].add(u)
        tris[k] = (a, b, v)
        tris.append((a, c, v))
        tris.append((b, c, v))
    return adj, n


# ------------------------------------------------------------ measurement

def bfs(adj, src, limit=None):
    """Distances from src. Returns dict node -> distance."""
    dist = {src: 0}
    q = deque([src])
    while q:
        u = q.popleft()
        if limit is not None and dist[u] >= limit:
            continue
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def geometry(adj, rng, sources=12):
    """Ball-growth dimension, mean distance and diameter.

    The fitting window is scaled to the graph rather than fixed: a fixed
    window ran past the diameter on the compact graphs and returned a
    NEGATIVE dimension, because the cumulative ball saturates at N and the
    log-log slope goes to zero. The window is now [2, rmax/3], and if that
    leaves fewer than four radii the dimension is reported as unmeasurable,
    which for a small-world graph is the correct answer rather than a
    failure -- there is no power law to fit.
    """
    tot = defaultdict(float)
    rmax = 0
    dsum = 0.0
    dcnt = 0
    keys = list(adj)
    for _ in range(sources):
        s = keys[int(rng.integers(len(keys)))]
        d = bfs(adj, s)
        vals = np.fromiter(d.values(), int)
        rmax = max(rmax, int(vals.max()))
        dsum += float(vals.sum())
        dcnt += len(vals)
        cum = np.cumsum(np.bincount(vals))
        for r in range(1, len(cum)):
            tot[r] += cum[r]
    mean_d = dsum / dcnt
    hi = rmax // 3
    rs = [r for r in sorted(tot) if 2 <= r <= hi]
    if len(rs) < 4:
        return float("nan"), mean_d, rmax
    x = np.log(np.array(rs, float))
    y = np.log(np.array([tot[r] / sources for r in rs]))
    return float(np.polyfit(x, y, 1)[0]), mean_d, rmax


def scaling_exponent(ns, vals):
    """Fit vals ~ N^p and also vals ~ c*ln N; return p and which fits better."""
    ns = np.asarray(ns, float)
    vals = np.asarray(vals, float)
    p, b = np.polyfit(np.log(ns), np.log(vals), 1)
    res_pow = float(np.sum((np.log(vals) - (p * np.log(ns) + b)) ** 2))
    c, b2 = np.polyfit(np.log(ns), vals, 1)
    res_log = float(np.sum(((vals - (c * np.log(ns) + b2)) / vals.mean()) ** 2))
    return float(p), res_pow, res_log


def comoving(adj, markers, base):
    """a_ij = d(i,j)/d_0(i,j) over all marker pairs."""
    out = []
    for i in markers:
        d = bfs(adj, i)
        for j in markers:
            if j > i and (i, j) in base and base[(i, j)] > 0:
                out.append(d[j] / base[(i, j)])
    return np.array(out)


def baseline(adj, markers):
    base = {}
    for i in markers:
        d = bfs(adj, i)
        for j in markers:
            if j > i:
                base[(i, j)] = d[j]
    return base


# ----------------------------------------------------------------- drivers

def cycle_seed(k):
    return [(i, (i + 1) % k) for i in range(k)], k


SE, NS = cycle_seed(6)
MARK6 = list(range(6))
MARK4 = list(range(4))


def seed_adjacency(edges):
    a = defaultdict(set)
    for x, y in edges:
        a[x].add(y)
        a[y].add(x)
    return a


BASE6 = baseline(seed_adjacency(SE), MARK6)
BASE4 = baseline(seed_adjacency([(i, j) for i in range(4)
                                 for j in range(4) if i < j]), MARK4)

RULES = [
    ("pure subdivide", MARK6, BASE6,
     lambda N, s: grow_subdivide(SE, NS, N, np.random.default_rng(s))),
    ("subdivide + shortcut", MARK6, BASE6,
     lambda N, s: grow_subdivide(SE, NS, N, np.random.default_rng(s),
                                 shortcut=1.0)),
    ("1->3 triangle move", MARK4, BASE4,
     lambda N, s: grow_triangulate(N, np.random.default_rng(s))),
]

CK = [2500, 5000, 10000, 20000, 40000]


def part_comoving(seeds=8):
    print("=" * 74)
    print("1. COMOVING TEST: do all marker pairs expand by the same factor?")
    print("=" * 74)
    print("  spread is the coefficient of variation of a_ij across marker")
    print("  pairs, averaged over %d independent growth histories." % seeds)
    print()
    for name, mark, base, build in RULES:
        print("  " + name)
        print("  %10s %12s %12s %10s" % ("N", "mean a", "spread", "min/max"))
        for N in CK:
            means, cvs, ratio = [], [], []
            for s in range(seeds):
                adj, n = build(N, 100 + s)
                a = comoving(adj, mark, base)
                means.append(a.mean())
                if a.mean() > 0:
                    cvs.append(a.std() / a.mean())
                    ratio.append(a.min() / a.max())
            print("  %10d %12.2f %11.1f%% %10.3f"
                  % (N, float(np.mean(means)), 100 * float(np.mean(cvs)),
                     float(np.mean(ratio))))
        print()


def part_dimension():
    print("=" * 74)
    print("2. DIMENSION: measured two ways, which have to agree")
    print("=" * 74)
    print("  d_ball from |B_r| ~ r^d, window [2, rmax/3]; nan means the graph")
    print("  is too compact for any scaling window to exist.")
    print()
    summary = {}
    for name, mark, base, build in RULES:
        print("  " + name)
        print("  %10s %9s %10s %9s %12s"
              % ("N", "d_ball", "mean dist", "diam", "a(N)"))
        ns, meand, avals, dbs = [], [], [], []
        for N in CK:
            adj, n = build(N, 100)
            db, md, rmax = geometry(adj, np.random.default_rng(3))
            a = comoving(adj, mark, base).mean()
            ns.append(n); meand.append(md); avals.append(a); dbs.append(db)
            print("  %10d %9.2f %10.3f %9d %12.3f"
                  % (n, db, md, rmax, a))
        p_d, rp, rl = scaling_exponent(ns, meand)
        # A dimension only means something if it has stopped moving. A d_ball
        # that drifts upward with N is the fitting window growing, not a
        # dimension -- the failure mode that makes small-world graphs look
        # low-dimensional if you only measure at one size.
        fin = [d for d in dbs if not np.isnan(d)]
        drift = (max(fin) - min(fin)) if fin else float("nan")
        stable = bool(fin) and drift < 0.15 * max(1.0, float(np.mean(fin)))
        summary[name] = (ns, avals, p_d, rp, rl,
                         float(np.mean(fin)) if fin else float("nan"), stable)
        print("     mean distance ~ N^%.3f   (power-law resid %.4f, "
              "log-law resid %.4f)" % (p_d, rp, rl))
        print("     -> %s" % ("power law, d_H = %.2f" % (1 / p_d) if rp < rl
                              else "logarithmic: small-world, no finite d_H"))
        if fin:
            print("     d_ball drift over the range: %.2f -> %.2f  (%s)"
                  % (fin[0], fin[-1],
                     "converged" if stable else "STILL DRIFTING, not a dimension"))
        print()

    print("=" * 74)
    print("3. DOES d_exp AGREE WITH d_ball?")
    print("=" * 74)
    print("  %-22s %8s %8s  %s" % ("rule", "d_exp", "d_ball", "verdict"))
    for name, (ns, avals, p_d, rp, rl, db, stable) in summary.items():
        if max(avals) - min(avals) < 1e-9:
            print("  %-22s %8s %8.2f  no expansion: a is constant"
                  % (name, "n/a", db))
            continue
        p_a, _, _ = scaling_exponent(ns, avals)
        d_exp = 1 / p_a if p_a != 0 else float("inf")
        near = (not np.isnan(db)) and abs(d_exp - db) < 0.25 * max(1.0, db)
        if not stable:
            v = "numbers agree but d_ball never converged - no dimension to agree on"
        elif rl < rp:
            v = "small-world: mean distance is logarithmic, not a power law"
        elif near:
            v = "consistent: a genuine %.0f-dimensional expansion" % round(db)
        else:
            v = "DISAGREE - not a d-dimensional space"
        print("  %-22s %8.2f %8.2f  %s" % (name, d_exp, db, v))
    print()


def main():
    part_comoving()
    part_dimension()


if __name__ == "__main__":
    main()
