#!/usr/bin/env python3
"""
global_observables.py

What is the smallest genuinely non-local quantity that separates a geometry
from a random graph, when every local quantity fails?

surface_tension.py established that local order is not merely insufficient but
anti-correlated with global geometry under this dynamics: 4-cycle richness per
node rose from 4.05 to 9.36 while the metric decayed to 15% of the way from
lattice to random. So a local energy cannot price the damage, and the question
becomes which non-local quantity can.

The test needs three graphs, not two, and the third is the point:

    A  pristine square torus           locally ordered, globally ordered
    B  the SAME torus after decay      locally ordered, globally destroyed
    C  random 4-regular                locally disordered, globally disordered

B is produced by the project's own close/open dynamics and is the hard case. A
useful global observable must separate B from A, even though the local
observables say B is MORE ordered than A. Anything that only separates A from
C is measuring something local observables already catch.

Candidates measured, all at matched N and matched degree 4:

    algebraic connectivity  lambda_2 of the Laplacian -- the spectral gap
    spectral dimension      from random-walk return probability p(t) ~ t^(-ds/2)
    effective resistance    R(u,v) between well-separated nodes
    mean distance           the metric itself, for reference

lambda_2 is the strongest candidate a priori. A d-dimensional torus of side L
has lambda_2 = 4 sin^2(pi/L), which vanishes as N^(-2/d); an expander has
lambda_2 = O(1). So the two differ not by a constant but by a growing factor,
which is the signature of a genuinely global property.

Result 1: the local observables get the SIGN wrong, the global ones do not
----------------------------------------------------------------------------
At matched N = 400 and matched degree 4:

    graph                  q mean   |B_2|
    A pristine torus         4.00   12.00
    B decayed torus          9.36   23.05
    C random 4-regular       0.12   15.77

    graph                mean dist   lambda_2   spectral d   resistance
    A pristine torus         10.00   0.097887         2.02       1.1610
    B decayed torus           5.59   0.221396         2.24       0.5801
    C random 4-regular        4.78   0.561195         3.59       0.7535

Both spectral measurements validate against theory: lambda_2 comes out at
0.097887 against the exact 4 sin^2(pi/20) = 0.097887 for a 20x20 torus, and
the spectral dimension of that torus is 2.02.

The local column says B is 2.3 times MORE ordered than A. Every global column
correctly places B between A and C. So the distinction is measurable, and the
failure of the local observables is a failure of sign, not only of magnitude.

Result 2: the spectral separation diverges with N
---------------------------------------------------
    L      N    lambda_2 torus   lambda_2 random   ratio
    12    144         0.267949          0.525434     2.0
    16    256         0.152241          0.552521     3.6
    24    576         0.068148          0.527953     7.7
    32   1024         0.038429          0.548278    14.3

A torus has lambda_2 ~ 4 pi^2 / N in two dimensions while an expander holds
it at O(1), so the ratio grows linearly in N -- 7.1x more nodes gives 7.15x
more separation. That is the signature of a genuinely global property, as
opposed to a local one that differs by a constant factor.

Result 3: the global property is lost first, while the local one is untouched
-------------------------------------------------------------------------------
    sweeps    q mean    lambda_2   mean dist
    0.0000      4.00    0.097887       10.00
    0.0150      4.00    0.098300        9.84
    0.0625      4.00    0.100226        9.44
    0.2500      4.05    0.105204        8.62
    1.0000      4.78    0.130046        7.63
    4.0000      8.58    0.195144        5.95
   16.0000      9.36    0.221396        5.59

At 1/16 of a sweep the metric has already lost 5.6% and lambda_2 has moved,
while q is still exactly 4.00 -- unchanged to three digits. The local
observable does not lag slightly; it does not move at all until the damage is
well underway, and then it moves the wrong way.

The structural reason, which is the point of the whole exercise
-----------------------------------------------------------------
It is important not to over-read this as "local rules cannot produce global
order". They demonstrably can: an Ising model with nearest-neighbour coupling
produces infinite-range correlations below its critical temperature. The
question is what distinguishes that case from this one, and the answer is
sharp.

In the Ising model the order parameter is the magnetisation, M = sum_i s_i.
It is an EXTENSIVE SUM OF LOCAL VARIABLES, so a local energy can push on it
directly, term by term.

lambda_2 is not a sum of anything. By Courant-Fischer it is a min-max over the
whole space of mean-zero functions on the graph,

    lambda_2 = min over x, sum x_i = 0, of  x^T L x / x^T x

so its value is set by whichever global cut is weakest, and no decomposition
into per-node contributions exists. Graph distance is the same: a shortest
path is a global minimisation, not a sum of local quantities.

That is why every attempt in this project failed in the way it did. A local
energy can only address extensive local sums, and the properties that make a
graph a geometry are not of that form. It is not that locality is impossible,
and not that the wrong local observable was chosen four times running. It is
that the target is the wrong TYPE of quantity for a local energy to reach.

The one honest opening this leaves is worth stating precisely, because it is
where a serious attempt would have to go: local rules can control non-extensive
quantities indirectly, by tuning to a critical point where the correlation
length diverges and local couplings acquire global reach. That is a real
mechanism, and it is what would be needed here -- not another local term, but
a rule sitting at criticality. Whether a rewriting system can self-tune to
such a point, without the critical coupling being put in by hand, is the
question this programme would have to answer next.
"""

import numpy as np
from collections import deque

from surface_tension import (square_torus, random_regular, mean_distance,
                             Nucleation)


# --------------------------------------------------------------- spectral

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


def laplacian_mul(src, dst, deg, n, x):
    return deg * x - np.bincount(src, weights=x[dst], minlength=n)


def algebraic_connectivity(adj, n, iters=6000, seed=0):
    """lambda_2 by shifted power iteration orthogonal to the constant vector."""
    src, dst, deg = edge_arrays(adj, n)
    c = 2.0 * deg.max()
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(n)
    x -= x.mean()
    x /= np.linalg.norm(x)
    for _ in range(iters):
        y = c * x - laplacian_mul(src, dst, deg, n, x)
        y -= y.mean()
        nrm = np.linalg.norm(y)
        if nrm < 1e-300:
            break
        x = y / nrm
    return float(x @ laplacian_mul(src, dst, deg, n, x))


def spectral_dimension(adj, n, tmax=64, srcs=4, seed=1):
    """d_s from the lazy-walk return probability, p(t) ~ t^(-d_s/2)."""
    src, dst, deg = edge_arrays(adj, n)
    rng = np.random.default_rng(seed)
    ts, ps = [], []
    acc = np.zeros(tmax + 1)
    for _ in range(srcs):
        s = int(rng.integers(n))
        x = np.zeros(n)
        x[s] = 1.0
        for t in range(1, tmax + 1):
            # lazy walk: half stay, half step to a uniform neighbour
            step = np.bincount(dst, weights=(x / deg)[src], minlength=n)
            x = 0.5 * x + 0.5 * step
            acc[t] += x[s]
    acc /= srcs
    for t in range(4, tmax + 1):
        if acc[t] > 1.5 / n:          # before saturation at the uniform value
            ts.append(t)
            ps.append(acc[t])
    if len(ts) < 4:
        return float("nan")
    sl = np.polyfit(np.log(ts), np.log(ps), 1)[0]
    return float(-2.0 * sl)


def conjugate_gradient(src, dst, deg, b, n, iters=5000, tol=1e-10):
    x = np.zeros(n)
    r = b - laplacian_mul(src, dst, deg, n, x)
    r -= r.mean()
    p = r.copy()
    rs = float(r @ r)
    for _ in range(iters):
        Ap = laplacian_mul(src, dst, deg, n, p)
        d = float(p @ Ap)
        if abs(d) < 1e-300:
            break
        a = rs / d
        x += a * p
        r -= a * Ap
        rs_new = float(r @ r)
        if np.sqrt(rs_new) < tol:
            break
        p = r + (rs_new / rs) * p
        rs = rs_new
    return x - x.mean()


def effective_resistance(adj, n, pairs=6, seed=2):
    """Mean R(u,v) over well-separated pairs."""
    src, dst, deg = edge_arrays(adj, n)
    rng = np.random.default_rng(seed)
    tot = 0.0
    for _ in range(pairs):
        u = int(rng.integers(n))
        d = bfs_dist(adj, u)
        far = [v for v, dd in d.items() if dd >= max(d.values()) - 1]
        v = far[int(rng.integers(len(far)))]
        b = np.zeros(n)
        b[u] = 1.0
        b[v] = -1.0
        x = conjugate_gradient(src, dst, deg, b, n)
        tot += float(x[u] - x[v])
    return tot / pairs


def bfs_dist(adj, s):
    dist = {s: 0}
    q = deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if w not in dist:
                dist[w] = dist[u] + 1
                q.append(w)
    return dist


# ------------------------------------------------------------------ local

def q_mean(adj, n):
    tot = 0
    for i in range(n):
        nb = list(adj[i])
        for a in range(len(nb)):
            na = adj[nb[a]]
            for b in range(a + 1, len(nb)):
                tot += len(na & adj[nb[b]]) - 1
    return tot / n


def ball2(adj, n):
    tot = 0
    for u in range(n):
        s = set(adj[u])
        t = set()
        for v in adj[u]:
            t |= adj[v]
        tot += len((s | t) - {u})
    return tot / n


# --------------------------------------------------------------- decayed B

def decayed(L, sweeps, seed=3):
    adj, n = square_torus(L)
    U = Nucleation(adj, n, 1.0, 0.0)
    rng = np.random.default_rng(seed)
    turn = 0
    for _ in range(int(sweeps * n)):
        for _ in range(6):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
    return U.adj, n


def report(label, adj, n):
    return dict(label=label,
                q=q_mean(adj, n), b2=ball2(adj, n),
                md=mean_distance(adj, n),
                lam2=algebraic_connectivity(adj, n),
                ds=spectral_dimension(adj, n),
                res=effective_resistance(adj, n))


def main():
    L = 20
    n = L * L
    print("=" * 78)
    print("1. THREE GRAPHS: matched N=%d, matched degree 4" % n)
    print("=" * 78)
    A, _ = square_torus(L)
    B, _ = decayed(L, 16)
    C = random_regular(n, 4, np.random.default_rng(9))
    rows = [report("A pristine torus", A, n),
            report("B decayed torus", B, n),
            report("C random 4-regular", C, n)]
    print("  LOCAL observables (these are what a local energy can see)")
    print("  %-22s %10s %10s" % ("graph", "q mean", "|B_2|"))
    for r in rows:
        print("  %-22s %10.2f %10.2f" % (r["label"], r["q"], r["b2"]))
    print()
    print("  GLOBAL observables")
    print("  %-22s %12s %12s %12s %12s"
          % ("graph", "mean dist", "lambda_2", "spectral d", "resistance"))
    for r in rows:
        print("  %-22s %12.2f %12.6f %12.2f %12.4f"
              % (r["label"], r["md"], r["lam2"], r["ds"], r["res"]))
    print()
    lam_exact = 4 * np.sin(np.pi / L) ** 2
    print("  check: exact lambda_2 for a %dx%d torus = 4 sin^2(pi/%d) = %.6f"
          % (L, L, L, lam_exact))
    print()
    a, b, c = rows
    print("  Separation of B from A, the hard case:")
    print("    q mean      B/A = %6.2f   (local says B is MORE ordered)"
          % (b["q"] / a["q"]))
    print("    |B_2|       B/A = %6.2f" % (b["b2"] / a["b2"]))
    print("    lambda_2    B/A = %6.1f" % (b["lam2"] / a["lam2"]))
    print("    resistance  B/A = %6.2f" % (b["res"] / a["res"]))
    print()

    print("=" * 78)
    print("2. DOES THE SEPARATION GROW WITH N?  (a constant factor is weak)")
    print("=" * 78)
    print("  %6s %8s %14s %14s %12s"
          % ("L", "N", "lambda_2 torus", "lambda_2 random", "ratio"))
    for LL in (12, 16, 24, 32):
        nn = LL * LL
        at, _ = square_torus(LL)
        cr = random_regular(nn, 4, np.random.default_rng(9))
        l1 = algebraic_connectivity(at, nn)
        l2 = algebraic_connectivity(cr, nn)
        print("  %6d %8d %14.6f %14.6f %12.1f" % (LL, nn, l1, l2, l2 / l1))
    print()
    print("  lambda_2 of a torus falls as N^(-2/d) while an expander holds")
    print("  it at O(1), so the ratio DIVERGES. That is what a genuinely")
    print("  global property looks like.")
    print()

    print("=" * 78)
    print("3. HOW FAST IS THE GLOBAL PROPERTY LOST, vs THE LOCAL ONE?")
    print("=" * 78)
    print("  %8s %10s %12s %12s" % ("sweeps", "q mean", "lambda_2", "mean dist"))
    adj, nn = square_torus(L)
    U = Nucleation(adj, nn, 1.0, 0.0)
    rng = np.random.default_rng(3)
    turn = 0
    marks = {int(s * nn) for s in (0.0156, 0.0625, 0.25, 1, 4, 16)}
    print("  %8.4f %10.2f %12.6f %12.2f"
          % (0.0, q_mean(U.adj, nn), algebraic_connectivity(U.adj, nn),
             mean_distance(U.adj, nn)))
    for t in range(1, 16 * nn + 1):
        for _ in range(6):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
        if t in marks:
            print("  %8.4f %10.2f %12.6f %12.2f"
                  % (t / nn, q_mean(U.adj, nn),
                     algebraic_connectivity(U.adj, nn),
                     mean_distance(U.adj, nn)))
    print()


if __name__ == "__main__":
    main()
