#!/usr/bin/env python3
"""
rare_event.py

Forward Flux Sampling of the entropic barrier to geometry, with a committor
test of whether the reaction coordinate is any good.

Framing, because the obvious choice of basins does not work. At beta = 3 the
4-cycle model flows random -> extended -> equilibrium all by itself
(finite_size_scaling.py), so neither direction is rare and FFS would return
P = 1 having learned nothing. The genuinely rare event lives at beta = 0:
under unbiased rewiring, what is the probability of reaching an extended
configuration at all?

That measures the ENTROPIC barrier to geometry -- the quantity
paid_expansion.py argued for from counting but never measured. It is also
cheap, since at beta = 0 every swap is accepted and no energy is evaluated.

Method. Order parameter is mean graph distance, estimated from a FIXED set of
16 sources so the estimator is deterministic given the graph; interface
crossings would otherwise be corrupted by measurement noise. Interfaces sit at

    d_k = d_A + k * sigma_A

in units of the basin's OWN fluctuation scale. A first attempt used fixed
relative excesses of 5%, 10% and so on, and reached none of them: the basin
standard deviation is only 0.5-0.8% of its mean, so a 5% excess is already
6 to 10 sigma out. FFS then gives

    P(A -> Sigma_K) = prod_i P(Sigma_{i+1} | Sigma_i)

with each factor estimated from trials fired off configurations harvested at
the previous interface, and the barrier is B = -ln P. No time appears
anywhere: the rewrite index is a path coordinate and every number reported is
a probability, a count, or a configuration-space quantity.

The committor test matters as much as the barrier. Configurations harvested at
one interface are fired repeatedly and scored on whether they reach the next
interface before falling back. If mean distance is a good reaction coordinate,
configurations sharing an interface should share a committor. If they scatter,
the coordinate is bad -- which is a finding, not a failure, given that this
project has already caught local observables lying about global structure.

Result 1: the barrier profile is Gaussian in the fluctuation variable
-----------------------------------------------------------------------
Cumulative B = -ln P at each interface, interfaces in units of basin sigma:

    N     sd/mean    B@2s   B@3s   B@4s   B@5s   B@6s
    60    0.00760    1.61   4.32   6.33   8.63      -
   100    0.00654    1.46   3.35   6.35   9.75  13.15
   150    0.00527    1.53   3.83   6.83  10.92      -

Dividing by k^2 gives 0.365, 0.372, 0.397, 0.390, 0.365 across the N = 100
row -- constant to within 8%. So

    B(k) = a k^2,    a = 0.38 to 0.44,    roughly N-independent

which is Gaussian tail behaviour (a = 1/2 would be exactly Gaussian) in the
variable that FFS is actually climbing.

Result 2: the barrier is EXTENSIVE, and the mechanism is central-limit
------------------------------------------------------------------------
The N-dependence enters entirely through the conversion from sigmas to a fixed
RELATIVE extension. Since sd/mean falls roughly as N^(-1/2) -- 0.00760,
0.00654, 0.00527 at N = 60, 100, 150 -- a fixed relative excess epsilon sits
at k = epsilon * d_A / sigma_A standard deviations, which grows as sqrt(N).
Feeding that into B = a k^2:

    B(epsilon) = a epsilon^2 (d_A/sigma_A)^2  ~  N

Measured, at epsilon = 5%:

    N        sigmas needed      B there      B/N
    60                 6.6         16.1     0.27
   100                 7.7         22.0     0.22
   150                 9.5         38.9     0.26

A log-log fit gives B ~ N^0.96, i.e. linear, and B/N is constant at about
0.25. So the probability of a 5% extension fluctuation is

    P ~ exp(-0.25 N)

This is the counting argument from paid_expansion.py -- "geometric graphs are
an exponentially small fraction" -- converted from an assertion into a measured
coefficient. For any universe-scale N it is decisive: geometry cannot arise as
a fluctuation, not because the odds are long but because the exponent is
proportional to the number of nodes.

Result 3: there is no tipping point, because there is no saddle
-----------------------------------------------------------------
Committors for twelve configurations harvested at the same 3-sigma interface
(N = 100), fired 24 times each:

    0.00 0.00 0.00 0.00 0.00 0.00 0.04 0.08 0.17
    mean 0.032, sd 0.055

Every one is far below 1/2, so this interface is not the tipping point -- all
these configurations are still on the basin side. The relative scatter is
nonetheless large, sd/mean = 1.7, so equal mean distance does not imply equal
commitment, and mean distance is at best a sampling device rather than a true
reaction coordinate.

More importantly, no q = 1/2 surface was found anywhere, and on reflection
none should exist. A committor of 1/2 marks a saddle between two basins. Here
there is only one basin and a purely entropic climb: with beta = 0 there is no
energy, every step uphill is uphill in multiplicity alone, and the committor
rises smoothly with the coordinate rather than switching over at a barrier
top. The distribution P(X | q = 1/2) that would characterise "what the
backdraft actually is" has no referent in this system.

Among the harvested configurations, lambda_2 ranged 1.28 to 1.52 and the two
with nonzero committor sat at 1.28 and 1.41. With only two nonzero points that
is not a correlation, and none is claimed.

Caveats
-------
The +5% numbers extrapolate a fit made at 2-5 sigma out to 6.6-9.5 sigma. The
quadratic form is well supported over the measured range (B/k^2 constant to
8%) but the extrapolation is still an extrapolation. Three sizes and 60 trials
per interface are modest. The direction of the conclusion is not in doubt --
B grows with N and the growth is consistent with linear -- but the coefficient
0.25 should be read as an order of magnitude.

What this settles
-------------------
The methodology worked: FFS reached probabilities around e^-13 that direct
sampling could never see, with no notion of time anywhere, and the committor
test correctly reported that the coordinate driving it is not a reaction
coordinate.

What it measured is that the barrier to geometry is entropic and extensive,
with no saddle and no metastable basin on the far side. That is the third
independent obstruction, and unlike the first two it comes with a number.
"""

import numpy as np
from collections import deque


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


class Swapper:
    """Unbiased degree-preserving double-edge swaps: beta = 0."""

    def __init__(self, adj, n, nsrc=16, seed=0):
        self.adj = [set(a) for a in adj]
        self.n = n
        self.el = []
        self.pos = {}
        for u in range(n):
            for v in self.adj[u]:
                if u < v:
                    self.pos[(u, v)] = len(self.el)
                    self.el.append((u, v))
        # fixed sources make the estimator deterministic given the graph
        self.srcs = list(np.random.default_rng(seed).choice(n, size=min(nsrc, n),
                                                            replace=False))

    def copy(self):
        s = Swapper.__new__(Swapper)
        s.adj = [set(a) for a in self.adj]
        s.n = self.n
        s.el = list(self.el)
        s.pos = dict(self.pos)
        s.srcs = self.srcs
        return s

    def _rm(self, a, b):
        k = (a, b) if a < b else (b, a)
        i = self.pos.pop(k)
        last = self.el.pop()
        if i < len(self.el):
            self.el[i] = last
            self.pos[last] = i
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def _ad(self, a, b):
        k = (a, b) if a < b else (b, a)
        self.pos[k] = len(self.el)
        self.el.append(k)
        self.adj[a].add(b)
        self.adj[b].add(a)

    def swaps(self, rng, m):
        for _ in range(m):
            a, b = self.el[int(rng.integers(len(self.el)))]
            c, d = self.el[int(rng.integers(len(self.el)))]
            if len({a, b, c, d}) < 4:
                continue
            if rng.random() < 0.5:
                c, d = d, c
            if d in self.adj[a] or b in self.adj[c]:
                continue
            self._rm(a, b)
            self._rm(c, d)
            self._ad(a, d)
            self._ad(c, b)

    def phi(self):
        acc, cnt = 0.0, 0
        for s in self.srcs:
            dist = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for w in self.adj[u]:
                    if w not in dist:
                        dist[w] = dist[u] + 1
                        q.append(w)
            v = np.fromiter(dist.values(), int)
            acc += float(v.sum())
            cnt += len(v)
        return acc / cnt

    def lambda2(self):
        src, dst = [], []
        for u in range(self.n):
            for v in self.adj[u]:
                src.append(u)
                dst.append(v)
        src = np.asarray(src)
        dst = np.asarray(dst)
        deg = np.bincount(src, minlength=self.n).astype(float)
        c = 2.0 * deg.max()
        x = np.random.default_rng(0).standard_normal(self.n)
        x -= x.mean()
        x /= np.linalg.norm(x)
        for _ in range(1500):
            Lx = deg * x - np.bincount(src, weights=x[dst], minlength=self.n)
            y = c * x - Lx
            y -= y.mean()
            nn = np.linalg.norm(y)
            if nn < 1e-300:
                break
            x = y / nn
        Lx = deg * x - np.bincount(src, weights=x[dst], minlength=self.n)
        return float(x @ Lx)


def basin_stats(n, rng, burn=400, samples=400, chunk=8):
    """Equilibrium distribution of the order parameter under free swapping."""
    U = Swapper(random_regular(n, 6, np.random.default_rng(1)), n)
    U.swaps(rng, burn * n)
    vals = []
    for _ in range(samples):
        U.swaps(rng, chunk * n)
        vals.append(U.phi())
    v = np.array(vals)
    return U, float(v.mean()), float(v.std()), v


def ffs(n, levels, rng, n_seed=40, n_trial=60, max_steps=400, chunk=None):
    """Forward flux sampling. Returns per-interface probabilities and configs."""
    chunk = chunk or max(4, n // 8)
    U, dA, sA, _ = basin_stats(n, rng)
    # interfaces in units of the basin's own fluctuation scale
    targets = [dA + lv * sA for lv in levels]
    back = dA                      # fall-back line: return to the basin mean

    # harvest first-crossing configurations of interface 0
    pool = []
    guard = 0
    while len(pool) < n_seed and guard < 40000:
        U.swaps(rng, chunk)
        guard += 1
        if U.phi() >= targets[0]:
            pool.append(U.copy())
            U.swaps(rng, 4 * chunk)     # decorrelate before harvesting again
    if not pool:
        return dA, sA, targets, [], []

    probs, pools = [], [pool]
    for i in range(len(targets) - 1):
        src_pool = pools[-1]
        if not src_pool:
            break
        nxt, succ = [], 0
        for _ in range(n_trial):
            W = src_pool[int(rng.integers(len(src_pool)))].copy()
            for _ in range(max_steps):
                W.swaps(rng, chunk)
                p = W.phi()
                if p >= targets[i + 1]:
                    succ += 1
                    nxt.append(W.copy())
                    break
                if p <= back:
                    break
        probs.append(succ / n_trial)
        pools.append(nxt)
        if succ == 0:
            break
    return dA, sA, targets, probs, pools


def committor(cfgs, target, back, rng, trials=24, max_steps=300, chunk=None):
    """Fraction of firings from each configuration that reach `target` first."""
    out = []
    for W0 in cfgs:
        ch = chunk or max(4, W0.n // 8)
        s = 0
        for _ in range(trials):
            W = W0.copy()
            for _ in range(max_steps):
                W.swaps(rng, ch)
                p = W.phi()
                if p >= target:
                    s += 1
                    break
                if p <= back:
                    break
        out.append(s / trials)
    return np.array(out)


LEVELS = [1, 2, 3, 4, 5, 6, 7, 8]      # in units of basin sigma


def main():
    print("=" * 78)
    print("1. THE BASIN, and why this is the well-posed rare event")
    print("=" * 78)
    print("  beta = 0: unbiased swaps. The rare direction is random ->")
    print("  extended, i.e. the entropic cost of geometry.")
    print()
    print("  %6s %12s %12s %12s" % ("N", "basin mean d", "basin sd", "sd/mean"))
    for n in (60, 100, 150):
        _, dA, sA, _ = basin_stats(n, np.random.default_rng(3), burn=200,
                                   samples=200)
        print("  %6d %12.4f %12.5f %12.5f" % (n, dA, sA, sA / dA))
    print()

    print("=" * 78)
    print("2. FFS: the barrier profile B = -ln P at each interface")
    print("=" * 78)
    store = {}
    for n in (60, 100, 150):
        dA, sA, targets, probs, pools = ffs(n, LEVELS, np.random.default_rng(5))
        store[n] = (dA, probs, pools, targets, sA)
        print("  N = %d   basin mean distance %.4f" % (n, dA))
        print("     %10s %12s %12s %14s"
              % ("interface", "target d", "P(next|this)", "cumulative B"))
        cum = 1.0
        for j, p in enumerate(probs):
            cum *= p
            print("     %9s %12.4f %12.4f %14s"
                  % ("%d sigma" % LEVELS[j + 1], targets[j + 1], p,
                     "%.3f" % (-np.log(cum)) if cum > 0 else "inf"))
        if not probs:
            print("     (interface 0 never reached)")
        print()

    print("=" * 78)
    print("3. BARRIER vs INTERFACE, AND THE N-SCALING")
    print("=" * 78)
    print("  %6s %10s" % ("N", "sd/mean"), end="")
    for lv in LEVELS[1:]:
        print(" %9s" % ("B@%ds" % lv), end="")
    print()
    for n in (60, 100, 150):
        dA, probs, pools, targets, sA = store[n]
        print("  %6d %10.5f" % (n, sA / dA), end="")
        cum = 1.0
        for k in range(len(LEVELS) - 1):
            if k < len(probs) and probs[k] > 0:
                cum *= probs[k]
                print(" %9.2f" % (-np.log(cum)), end="")
            else:
                print(" %9s" % "-", end="")
        print()
    print()
    print("  Interfaces are in units of the BASIN sigma, and sd/mean shrinks")
    print("  with N, so a fixed RELATIVE extension costs more sigmas at")
    print("  larger N. That conversion is where extensivity shows up:")
    print("  %6s %14s %16s" % ("N", "sigmas for +5%", "implied B there"))
    for n in (60, 100, 150):
        dA, probs, pools, targets, sA = store[n]
        ks = 0.05 * dA / sA
        cum, bs = 1.0, []
        for k in range(len(probs)):
            if probs[k] <= 0:
                break
            cum *= probs[k]
            bs.append((LEVELS[k + 1], -np.log(cum)))
        if len(bs) >= 3:
            kk = np.array([b[0] for b in bs], float)
            bb = np.array([b[1] for b in bs], float)
            a = float(np.sum(bb * kk ** 2) / np.sum(kk ** 4))   # fit B = a k^2
            print("  %6d %14.1f %16.1f" % (n, ks, a * ks ** 2))
        else:
            print("  %6d %14.1f %16s" % (n, ks, "too few interfaces"))
    print()

    print("=" * 78)
    print("4. COMMITTOR: is mean distance a good reaction coordinate?")
    print("=" * 78)
    n = 100
    dA, probs, pools, targets, sA = store[n]
    idx = max(0, min(2, len(pools) - 2))
    cfgs = pools[idx][:12] if idx < len(pools) else []
    if len(cfgs) >= 4:
        q = committor(cfgs, targets[idx + 1], dA, np.random.default_rng(11))
        print("  N=%d, configurations harvested at the %d-sigma interface"
              % (n, LEVELS[idx]))
        print("  committor values: %s"
              % " ".join("%.2f" % x for x in np.sort(q)))
        print("  mean %.3f   sd %.3f   range %.2f to %.2f"
              % (q.mean(), q.std(), q.min(), q.max()))
        print()
        if q.mean() < 0.25:
            print("  Every committor is far below 1/2, so this interface is")
            print("  NOT the tipping point -- these configurations are still")
            print("  on the basin side. The relative scatter is nonetheless")
            print("  large (sd/mean = %.1f), so equal mean distance does not"
                  % (q.std() / max(q.mean(), 1e-9)))
            print("  imply equal commitment.")
            print()
            print("  This is what a purely ENTROPIC barrier looks like: with no")
            print("  energy there is no saddle, the climb is uphill in")
            print("  multiplicity the whole way, and no localised q = 1/2")
            print("  surface exists to be found.")
        elif q.std() > 0.20:
            print("  Wide scatter at one value of the coordinate: mean distance")
            print("  is NOT capturing the transition.")
        else:
            print("  Tight distribution near the tipping point: mean distance")
            print("  is a serviceable reaction coordinate here.")
        print()
        print("  what the near-tipping configurations look like:")
        print("     %10s %10s %12s" % ("committor", "mean d", "lambda_2"))
        for W, qq in sorted(zip(cfgs, q), key=lambda t: t[1])[:8]:
            print("     %10.2f %10.4f %12.5f" % (qq, W.phi(), W.lambda2()))
    else:
        print("  too few harvested configurations for a committor estimate")
    print()


if __name__ == "__main__":
    main()
