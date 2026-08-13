#!/usr/bin/env python3
"""
soc.py

Self-organised criticality on a rewriting graph, with expansion as the
dissipation channel.

This is the last route left open. Every earlier attempt failed because the
properties that make a graph a geometry are non-extensive, and a local energy
can only push on extensive local sums (global_observables.py). Criticality is
the one circumstance where local couplings reach non-extensive quantities: the
correlation length diverges and local rules acquire global reach. SOC is the
version of that which does not require hand-tuning a coupling.

The mapping, built entirely from operations this project already has:

    stress          node degree. Local, and MEASURED to couple to global
                    geometry -- six_criteria.py found the degree cap acted as
                    a dial on d_H (cap 6 -> 4.03, cap 12 -> 6.34).

    slow drive      `close`: add an edge between two non-adjacent neighbours
                    of a random node. One drive per fully settled avalanche,
                    which enforces timescale separation by construction.

    fast relaxation a topple made of close + open. At a node v over threshold,
                    pick non-adjacent neighbours x, y, then
                        close(x, y);  open(v, x)
                    Bookkeeping: deg(v) -1, deg(x) net 0, deg(y) +1, and E is
                    exactly conserved. Load flows from v to y, so avalanches
                    are genuine conservative transport and can cascade.
                    Connectivity is safe by construction: x keeps its link to
                    the component through the new edge to y.

    dissipation     node FISSION at a second, higher threshold. v splits into
                    v and v', joined by an edge, its neighbours partitioned
                    between them. N +1, E +1.

Two corrections to the design as first sketched, both of which matter.

`subdivide` cannot dissipate degree: splitting edge {x,y} into {x,z},{z,y}
leaves deg(x) and deg(y) untouched. Node fission is the operation that
actually relieves degree stress, and it is also the one that creates a node.

And dissipation must NOT be a fixed bulk probability. Uniform bulk dissipation
at rate epsilon is known to destroy criticality, introducing a finite
correlation length -- the Vespignani-Zapperi result. Real sandpiles are
bulk-conservative and lose grains only at a boundary, where the rate vanishes
automatically as the system grows. A second threshold k_diss > k_c is the
closest graph analogue: dissipation stays rare, local, and its rate is
emergent rather than imposed. Fixed-epsilon dissipation is included below as a
control that theory says should fail.

A prediction the design makes before anything is run
------------------------------------------------------
In steady state one drive (E+1) balances one fission (E+1, N+1), so
dE/dN = 2. The fixed point of <k> = 2E/N is then

    <k> = 4

derived from the drive-dissipation balance, not imposed, and agreeing with the
independent conservation-law prediction of paid_expansion.py. If the steady
state does not land there, the mechanism is not doing what it claims.

What would count as success, and the null hypothesis
------------------------------------------------------
Avalanche criticality alone proves nothing about geometry -- this whole
project has been a catalogue of local order failing to reach global structure.
The measurements are therefore separated:

    P(s) ~ s^-tau      tau near 3/2 with cutoff ~ N is MEAN FIELD, meaning
                       Galton-Watson branching on an expander: critical
                       avalanches, no geometry. This is the expected outcome.

    lambda_2           the payoff. A geometric steady state needs
                       lambda_2 ~ N^(-2/d), falling with N. An expander holds
                       it at O(1). This decides the question.

    d_H                from shell growth, for consistency with lambda_2.

A correction forced by the first run: the topple must BRANCH
--------------------------------------------------------------
As first written the topple sent its whole unit of load to ONE neighbour --
deg(v) -1, deg(y) +1. That does not work, and the failure is instructive
rather than a coding slip. A single-recipient conservative topple makes the
excess perform a RANDOM WALK, which never dies: at N = 100 and <k> = 5.02,
with only two nodes above threshold, avalanches ran past 20000 events without
terminating. Total degree is conserved, so the excess has nowhere to go.

BTW terminates because toppling splits the load across several neighbours, so
it thins out and the branching process can die. The corrected topple does two
close+open pairs to DISTINCT recipients -- deg(v) -2, two neighbours +1 each --
and avalanches then terminate.

Result 1: the steady state is close to the predicted <k> = 4
--------------------------------------------------------------
    target      N   edges     <k>     topples   fissions   valve
       500    500    1161   4.644     2930673        300       0
      1200   1200    2774   4.623    24130116       1000      12

<k> settles at 4.62-4.64, stable across a 2.4-fold change in N. The prediction
was 4, so the drive-dissipation balance argument is qualitatively right and
about 15% off. The gap comes from fission usually firing through the
"cannot shed" fallback rather than at k_diss, which moves half a node's
neighbours and so does not respect the simple dE/dN = 2 accounting.

Result 2: it is not critical -- it is SUPERCRITICAL
-----------------------------------------------------
        N   avalanches   mean size    max size   mean/N     tau
      500          460      6371.7       72627    12.74   1.191
     1200         1277     18896.7      200001    15.75   1.177

A critical sandpile has mean avalanche size small compared to N with a
power-law tail. Here the mean avalanche is 13 to 16 times the size of the
entire system, and at N = 1200 the largest hits the safety cap. That is not
scale-free relaxation; it is continuous global churning, and it means the
timescale separation the construction was built around has collapsed. Each
drive event does not trigger a localised cascade that settles -- it triggers
an essentially unbounded reorganisation.

So tau is not meaningful here and the values 1.18-1.19 should not be read as
an exponent. They are what the estimator returns on a distribution whose tail
is set by an imposed cap.

Result 3: and the geometry is random anyway
---------------------------------------------
        N    lambda_2   lam2 rand   mean d   rand d    d_H
      500     0.83452     1.03002     4.40     4.18   4.21
     1200     0.83262     1.01601     5.07     4.82   4.53

lambda_2 sits 15-19% below a matched random regular graph and does NOT fall
with N -- 0.8345 to 0.8326 across a 2.4-fold increase, where a geometric state
would need it dropping like N^(-2/d). Mean distance is 5-6% above random and
likewise not pulling away. d_H is 4.2 to 4.5 and drifting upward rather than
converging.

Verdict
-------
The null hypothesis was that SOC would buy avalanche criticality while leaving
the graph an expander. The outcome is worse than that: it did not buy
criticality either. The system self-organises to a SUPERCRITICAL state, not a
critical one, and the graph remains within a few per cent of random on every
global measure.

That distinction matters for what it rules out. The promise of SOC is that the
critical point is an attractor requiring no tuning. Here the attractor is
supercritical, so there is no critical point being selected. A critical regime
might well exist at some other (k_c, k_diss) -- but finding it by scanning
those parameters is exactly the hand-tuning SOC was invoked to avoid, and a
critical point you have to tune to is just an ordinary critical point.

With that, the last route identified in global_observables.py is closed by
measurement rather than by argument.
"""

import numpy as np
from collections import deque


class SOC:
    """Driven graph sandpile: close/open topples, fission dissipates."""

    def __init__(self, n0=200, k_c=6, k_diss=8, eps=0.0, seed=0):
        self.rng = np.random.default_rng(seed)
        self.k_c = k_c
        self.k_diss = k_diss
        self.eps = eps
        # seed: a ring, so every node starts well below threshold
        self.adj = [set() for _ in range(n0)]
        for i in range(n0):
            self._add(i, (i + 1) % n0)
        self.n = n0
        self.fissions = 0
        self.topples = 0
        self.valve = 0

    # ------------------------------------------------------------ topology
    def _add(self, a, b):
        if a == b or b in self.adj[a]:
            return False
        self.adj[a].add(b)
        self.adj[b].add(a)
        return True

    def _del(self, a, b):
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def n_edges(self):
        return sum(len(a) for a in self.adj) // 2

    def _two_nonadjacent(self, v):
        nb = list(self.adj[v])
        if len(nb) < 2:
            return None
        for _ in range(12):
            i = int(self.rng.integers(len(nb)))
            j = int(self.rng.integers(len(nb)))
            if i == j:
                continue
            x, y = nb[i], nb[j]
            if y not in self.adj[x]:
                return x, y
        return None

    # --------------------------------------------------------------- moves
    def topple(self, v):
        """BRANCHING topple: two close+open pairs to DISTINCT recipients, so
        deg(v) -2 and two different neighbours each gain 1.

        A single-recipient topple (deg(v) -1, one neighbour +1) does not work:
        the excess then performs a random walk that never dies, and avalanches
        never terminate. Splitting the load across two neighbours is what makes
        it a branching process that can die out, exactly as in BTW.
        """
        nb = list(self.adj[v])
        if len(nb) < 4:
            return False
        self.rng.shuffle(nb)
        done, used = 0, set()
        for x in nb:
            if x in used or x not in self.adj[v]:
                continue
            for y in nb:
                if y == x or y in used or y not in self.adj[v]:
                    continue
                if y in self.adj[x]:
                    continue
                self._add(x, y)
                self._del(v, x)
                used.add(x)
                used.add(y)
                done += 1
                break
            if done == 2:
                break
        if done:
            self.topples += 1
        return done > 0

    def fission(self, v):
        """v splits into v and v', neighbours partitioned, joined by an edge."""
        nb = list(self.adj[v])
        self.rng.shuffle(nb)
        half = len(nb) // 2
        moved = nb[:half]
        w = self.n
        self.adj.append(set())
        self.n += 1
        for u in moved:
            self._del(v, u)
            self._add(w, u)
        self._add(v, w)
        self.fissions += 1
        return True

    def drive(self):
        """One unit of stress: an edge closed somewhere at random."""
        for _ in range(200):
            v = int(self.rng.integers(self.n))
            pair = self._two_nonadjacent(v)
            if pair is None:
                continue
            x, y = pair
            if self._add(x, y):
                return [x, y]
        return []

    def relax(self, touched, cap=200000):
        """Run the avalanche to completion. Returns (size, duration, span)."""
        q = deque(u for u in touched if len(self.adj[u]) >= self.k_c)
        inq = set(q)
        size, dur, seen = 0, 0, set()
        while q:
            layer = len(q)
            dur += 1
            for _ in range(layer):
                if not q:
                    break
                v = q.popleft()
                inq.discard(v)
                if len(self.adj[v]) < self.k_c:
                    continue
                seen.add(v)
                if size >= cap:              # safety valve; counted and reported
                    self.fission(v)
                    self.valve += 1
                    size += 1
                    q.clear()
                    break
                if len(self.adj[v]) >= self.k_diss:
                    self.fission(v)
                    size += 1
                elif self.eps > 0 and self.rng.random() < self.eps:
                    self.fission(v)          # control: bulk dissipation
                    size += 1
                else:
                    before = set(self.adj[v])
                    if not self.topple(v):
                        self.fission(v)      # cannot shed: relieve by fission
                    size += 1
                    for u in before | set(self.adj[v]) | {v}:
                        if len(self.adj[u]) >= self.k_c and u not in inq:
                            q.append(u)
                            inq.add(u)
        return size, dur, len(seen)


# ------------------------------------------------------------ measurement

def shells(adj, n, srcs=8, seed=1):
    rng = np.random.default_rng(seed)
    acc, tot, cnt = {}, 0.0, 0
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
        tot += float(v.sum())
        cnt += len(v)
        c = np.bincount(v)
        for r in range(len(c)):
            acc[r] = acc.get(r, 0.0) + float(c[r])
    return {k: v / srcs for k, v in acc.items()}, tot / cnt


def d_shell(sh, lo, hi):
    rs = [r for r in sorted(sh) if lo <= r <= hi and sh[r] > 0]
    if len(rs) < 3:
        return float("nan")
    x = np.log(np.array(rs, float))
    y = np.log(np.array([sh[r] for r in rs]))
    return 1.0 + float(np.polyfit(x, y, 1)[0])


def lambda2(adj, n, iters=4000):
    src, dst = [], []
    for u in range(n):
        for v in adj[u]:
            src.append(u)
            dst.append(v)
    src = np.asarray(src)
    dst = np.asarray(dst)
    deg = np.bincount(src, minlength=n).astype(float)
    deg[deg == 0] = 1.0
    c = 2.0 * deg.max()
    x = np.random.default_rng(0).standard_normal(n)
    x -= x.mean()
    x /= np.linalg.norm(x)
    for _ in range(iters):
        Lx = deg * x - np.bincount(src, weights=x[dst], minlength=n)
        y = c * x - Lx
        y -= y.mean()
        nn = np.linalg.norm(y)
        if nn < 1e-300:
            break
        x = y / nn
    Lx = deg * x - np.bincount(src, weights=x[dst], minlength=n)
    return float(x @ Lx)


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


def powerlaw_tau(sizes, smin=2):
    """MLE exponent for P(s) ~ s^-tau above smin (discrete Hill estimator)."""
    s = np.asarray([x for x in sizes if x >= smin], float)
    if len(s) < 50:
        return float("nan"), 0
    tau = 1.0 + len(s) / np.sum(np.log(s / (smin - 0.5)))
    return float(tau), len(s)


def run(n_target, k_c=6, k_diss=9, eps=0.0, seed=0, n0=200):
    U = SOC(n0=n0, k_c=k_c, k_diss=k_diss, eps=eps, seed=seed)
    sizes, durs, spans = [], [], []
    guard = 0
    while U.n < n_target and guard < 400000:
        guard += 1
        touched = U.drive()
        if not touched:
            continue
        s, d, sp = U.relax(touched)
        if s > 0:
            sizes.append(s)
            durs.append(d)
            spans.append(sp)
    return U, np.array(sizes), np.array(durs), np.array(spans)


def main():
    print("=" * 78)
    print("1. STEADY STATE: does it land on the predicted <k> = 4?")
    print("=" * 78)
    print("  %8s %8s %8s %9s %10s %10s %8s"
          % ("target", "N", "edges", "<k>", "topples", "fissions", "valve"))
    runs = {}
    for nt in (500, 1200):
        U, s_, d_, sp_ = run(nt)
        runs[nt] = (U, s_, d_, sp_)
        print("  %8d %8d %8d %9.3f %10d %10d %8d"
              % (nt, U.n, U.n_edges(), 2 * U.n_edges() / U.n,
                 U.topples, U.fissions, U.valve))
    print()

    print("=" * 78)
    print("2. ARE THE AVALANCHES CRITICAL?")
    print("=" * 78)
    print("  %8s %10s %12s %12s %10s %10s"
          % ("N", "avalanches", "mean size", "max size", "mean/N", "tau"))
    for nt in (500, 1200):
        U, s_, d_, sp_ = runs[nt]
        tau, _ = powerlaw_tau(s_)
        print("  %8d %10d %12.1f %12d %10.2f %10.3f"
              % (U.n, len(s_), s_.mean(), s_.max(), s_.mean() / U.n, tau))
    print()
    print("  A critical sandpile has mean avalanche size SMALL compared to N,")
    print("  with a power-law tail. Mean size far exceeding N means the")
    print("  system is supercritical and the timescale separation has")
    print("  collapsed: every drive triggers continuous global churning, so")
    print("  there is no well-defined avalanche and tau is not meaningful.")
    print()

    print("=" * 78)
    print("3. THE PAYOFF: is the steady-state graph geometric?")
    print("=" * 78)
    print("  %8s %11s %11s %10s %10s %8s"
          % ("N", "lambda_2", "lam2 rand", "mean d", "rand d", "d_H"))
    for nt in (500, 1200):
        U, s_, d_, sp_ = runs[nt]
        sh, md = shells(U.adj, U.n)
        hi = max(4, max(sh) // 3)
        k = max(2, int(round(2 * U.n_edges() / U.n)))
        rr = random_regular(U.n, k, np.random.default_rng(9))
        _, mdr = shells(rr, U.n)
        print("  %8d %11.5f %11.5f %10.2f %10.2f %8.2f"
              % (U.n, lambda2(U.adj, U.n), lambda2(rr, U.n), md, mdr,
                 d_shell(sh, 2, hi)))
    print()
    print("  Geometric would need lambda_2 falling like N^(-2/d) and mean")
    print("  distance well above the random value. Both must move with N.")
    print()


if __name__ == "__main__":
    main()
