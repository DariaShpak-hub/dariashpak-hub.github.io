#!/usr/bin/env python3
"""
paid_expansion.py

Node births paid for out of a conserved budget, so the growth rate is derived
from structure instead of stipulated -- and then the reason it still fails to
produce a space, measured rather than guessed.

In scale_factor.py subdivide is free. Nodes appear because the rule says so,
which is not an explanation of why space grows, and the growth rate is
whatever the schedule imposes. Here nothing is free. Three local moves on an
ordinary graph trade against a single budget:

    close      two non-adjacent neighbours of v become adjacent    bank +1
    subdivide  {x,y} -> {x,z},{z,y}, a node is born                bank -1
    open       an edge sitting in a triangle is removed            bank -1

with bank >= 0. Clumping -- closing a triangle, the only move that creates an
adjacency -- is the sole source of budget, and node birth is one of the two
sinks. Every node that is ever born was paid for by a clump. All three moves
are local, and all three preserve connectivity: closure and subdivision
obviously, and removing an edge from a triangle leaves its endpoints joined by
the third vertex.

The conserved quantity is exact and, unusually, has teeth:

    Q  =  bank + 2N - E     is invariant under all three moves

Since the bank stays O(1), conservation PINS the mean degree at 4. That is not
imposed anywhere; it is forced by the accounting, and it is the first
quantitative prediction in this project that comes from a conservation law
rather than from a measurement.

Result 1: the budget works, and the degree prediction holds
-------------------------------------------------------------
Growth is unbounded and the conservation residual is exactly 0 at every
checkpoint. Mean degree converges to the predicted 4 from below:

    events     5000    20000    80000
    N          1188     4878    19414
    <k>       3.992    3.998    3.999

The condition matters, though. The prediction is <k> = 4 + 2(bank - Q)/N, so
it holds only while the bank stays O(1). Weight closure heavily enough
(3:1:1 below) and clumping outruns both sinks, the bank accumulates to ~4086,
and <k> drifts to 6.085. The budget only pins the degree while it is actually
binding.

Result 2: the exchange rate is NOT derived -- this is the failure
------------------------------------------------------------------
At equal propensities the rate is clean, half a node per closure:

    births per close = 0.4989,  opens per close = 0.5010

and it is reparameterisation-invariant by construction, being a ratio of event
counts with no clock in it. That much was the point of building it this way.

But it is not invariant under reweighting the move propensities:

    close:sub:open      1:1:1   2:1:1   1:2:1   1:1:2   3:1:1   1:3:3
    births per close   0.4989  0.4927  0.6758  0.3324  0.3245  0.4974

It ranges over 0.32 to 0.68 depending on a choice that is pure schedule. So
the exchange rate is not a property of the structure after all -- it inherits
the propensities, exactly the reparameterisation-survives / reweighting-fails
split found in schedule_invariance.py. The conservation law constrains the
budget to balance; it does not decide how the budget is spent, and the split
between the two sinks is what G actually is.

This is the claim the module was built to establish, and it does not hold.
What conservation buys is the degree law, which is real, and the guarantee
that births are paid for, which is also real. It does not buy a derived G.

Result 3: it is a small-world graph, not a space
--------------------------------------------------
Mean distance fits a logarithm better than a power law, over a 65-fold range
in N, for both natural ways of choosing where to close:

    v degree-weighted   mean dist 5.89 -> 8.93 while N: 1188 -> 77880
    v uniform           mean dist 5.96 -> 8.06 while N: 1208 -> 78455

Diameter grows like log N and the maximum degree runs away -- 160 and 605
respectively -- even though the mean is pinned at 4. Hubs form.

Tightening locality does not fix it, it stops the universe. Restricting
closure to pairs of degree-2 nodes, the freshest nodes subdivision makes,
deadlocks at N = 8 and stays there for 320000 events: no legal close means no
budget, no budget means no births.

So there is a trilemma, and all three corners have now been measured:

    no closure at all         -> 1-dimensional threads   (scale_factor.py)
    closure, loosely local    -> small-world, no dimension
    closure, strictly local   -> deadlock, no expansion

Result 4: the obstruction is entropy, and that is measurable
--------------------------------------------------------------
The natural next move is to hunt for a fourth rule. That would be fitting.
Instead, compare the model's output against the two reference ensembles at
MATCHED N and matched mean degree 4:

    graph on N = 19414 nodes, <k> = 4      mean dist    diameter
    the model's output                          8.04          26
    random 4-regular graph                      8.33          11
    periodic square lattice                    69.50         138

On mean distance the model is indistinguishable from a random graph of the
same size and degree, 8.04 against 8.33, and both are nowhere near the
lattice's 69.50. (The model's diameter is longer than the random graph's, 26
against 11, because subdivision leaves dangling chains; mean distance is the
robust statistic and it lands on random.) Nothing about the budget, the
locality, or the conservation law moved the result away from typical.

That is the actual obstruction, and it is not specific to these rules.
Geometric graphs are an exponentially small fraction of graphs at fixed N and
degree -- a lattice has diameter ~sqrt(N) where a typical graph has ~log N --
so any dynamics that samples configurations by availability, which is what a
random local rewrite does, runs away from geometry for the same reason a gas
fills a box. Entropy is not on geometry's side.

This is why lattice quantum-gravity programmes need more than local moves.
CDT imposes a causal slicing by hand; quantum graphity puts in an energy
functional and has to COOL through a phase transition to reach its geometric
phase. Both are ways of paying the entropic cost of being low-dimensional. A
purely entropic model of the kind built here cannot get one for free, and the
preceding four rule sets are four instances of that single fact.

Consequence, stated against what the module set out to do. Making births paid
for was supposed to deliver two things: a conservation law with real content,
and a growth rate derived from structure. It delivers the first and not the
second. The law is exact, it makes births impossible without clumping, and it
predicts the mean degree -- a genuine derived number. But G still comes from
the propensities, so the budget disciplines expansion without explaining its
rate.

And the dimension question is not merely unsolved, it is blocked, for a reason
this project's own entropy machinery predicts. Four rule sets have now been
tried; all four land on typical graphs, because typical is what entropy
selects.
"""

import numpy as np
from collections import deque


class Universe:
    """Simple graph with the conserved budget Q = bank + 2N - E."""

    def __init__(self, n=6, mode="edge"):
        self.adj = [set() for _ in range(n)]
        self.el = []
        self.pos = {}
        self.n = n
        self.bank = 0
        self.mode = mode
        self.cnt = {"close": 0, "subdivide": 0, "open": 0}
        for i in range(n):
            self._add(i, (i + 1) % n)
        self.Q = self.bank + 2 * self.n - len(self.el)

    # ------------------------------------------------------------ plumbing
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

    # --------------------------------------------------------------- moves
    def close(self, rng):
        """Clumping: the only move that creates an adjacency. Pays 1."""
        for _ in range(24):
            if self.mode == "node":
                v = int(rng.integers(self.n))
            else:
                a, b = self.el[int(rng.integers(len(self.el)))]
                v = a if rng.random() < 0.5 else b
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            x = nb[int(rng.integers(len(nb)))]
            y = nb[int(rng.integers(len(nb)))]
            if x == y or y in self.adj[x]:
                continue
            if self.mode == "strict" and not (len(self.adj[x]) <= 2
                                              and len(self.adj[y]) <= 2):
                continue
            self._add(x, y)
            self.bank += 1
            self.cnt["close"] += 1
            return True
        return False

    def subdivide(self, rng):
        """A node is born in the middle of an edge. Costs 1."""
        if self.bank < 1:
            return False
        a, b = self.el[int(rng.integers(len(self.el)))]
        self.adj.append(set())
        z = self.n
        self.n += 1
        self._del(a, b)
        self._add(a, z)
        self._add(z, b)
        self.bank -= 1
        self.cnt["subdivide"] += 1
        return True

    def open(self, rng):
        """An edge in a triangle is removed; connectivity survives. Costs 1."""
        if self.bank < 1:
            return False
        for _ in range(24):
            a, b = self.el[int(rng.integers(len(self.el)))]
            if self.adj[a] & self.adj[b]:
                self._del(a, b)
                self.bank -= 1
                self.cnt["open"] += 1
                return True
        return False

    # --------------------------------------------------------- measurement
    def bfs(self, s):
        d = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in self.adj[u]:
                if w not in d:
                    d[w] = d[u] + 1
                    q.append(w)
        return d


def distances(adj_of, n, rng, k=8):
    """Mean distance and diameter from k random sources, via BFS."""
    tot, cnt, dm = 0.0, 0, 0
    for _ in range(k):
        s = int(rng.integers(n))
        d = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in adj_of(u):
                if w not in d:
                    d[w] = d[u] + 1
                    q.append(w)
        v = np.fromiter(d.values(), int)
        tot += float(v.sum())
        cnt += len(v)
        dm = max(dm, int(v.max()))
    return tot / cnt, dm


def evolve(mode, events, seed, weights=(1, 1, 1)):
    rng = np.random.default_rng(seed)
    U = Universe(mode=mode)
    ops = ["close", "subdivide", "open"]
    w = np.array(weights, float)
    w = w / w.sum()
    flat = abs(w - w[0]).max() < 1e-12
    for _ in range(events):
        for _ in range(6):
            i = int(rng.integers(3)) if flat else int(rng.choice(3, p=w))
            if getattr(U, ops[i])(rng):
                break
    return U


def fit(ns, vals):
    ns, vals = np.asarray(ns, float), np.asarray(vals, float)
    p, b = np.polyfit(np.log(ns), np.log(vals), 1)
    rp = float(np.sum((np.log(vals) - (p * np.log(ns) + b)) ** 2))
    c, b2 = np.polyfit(np.log(ns), vals, 1)
    rl = float(np.sum(((vals - (c * np.log(ns) + b2)) / vals.mean()) ** 2))
    return float(p), rp, rl


# ------------------------------------------------------- reference graphs

def random_regular(n, k, rng):
    """Configuration model; a few collisions are dropped, so <k> is ~k."""
    stubs = np.repeat(np.arange(n), k)
    rng.shuffle(stubs)
    adj = [set() for _ in range(n)]
    for i in range(0, len(stubs) - 1, 2):
        a, b = int(stubs[i]), int(stubs[i + 1])
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def square_lattice(n):
    """Periodic square lattice on the largest m*m <= n. Degree exactly 4."""
    m = int(np.sqrt(n))
    adj = [set() for _ in range(m * m)]
    for i in range(m):
        for j in range(m):
            u = i * m + j
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                adj[u].add(((i + di) % m) * m + (j + dj) % m)
    return adj, m * m


# ------------------------------------------------------------------ parts

EVENTS = [5000, 20000, 80000]


def part1():
    print("=" * 76)
    print("1. DOES PAID EXPANSION RUN, AND IS THE DEGREE PREDICTION RIGHT?")
    print("=" * 76)
    print("  conservation law Q = bank + 2N - E  =>  <k> -> 4 once bank is O(1)")
    print()
    print("  %8s %8s %8s %7s %8s %10s %9s %9s"
          % ("events", "N", "E", "bank", "<k>", "residual", "closes", "births"))
    for ev in EVENTS:
        U = evolve("edge", ev, 11)
        print("  %8d %8d %8d %7d %8.3f %10d %9d %9d"
              % (ev, U.n, len(U.el), U.bank, 2 * len(U.el) / U.n,
                 U.residual(), U.cnt["close"], U.cnt["subdivide"]))
    print()


def part2():
    print("=" * 76)
    print("2. THE DERIVED EXCHANGE RATE, AND WHETHER IT SURVIVES THE SCHEDULE")
    print("=" * 76)
    print("  %-18s %9s %14s %14s %9s"
          % ("close:sub:open", "N", "births/close", "opens/close", "<k>"))
    for wts in [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2), (3, 1, 1),
                (1, 3, 3)]:
        U = evolve("edge", 20000, 11, wts)
        c = max(U.cnt["close"], 1)
        print("  %-18s %9d %14.4f %14.4f %9.3f"
              % (":".join(map(str, wts)), U.n, U.cnt["subdivide"] / c,
                 U.cnt["open"] / c, 2 * len(U.el) / U.n))
    print()
    print("  a ratio of event counts is reparameterisation-invariant by")
    print("  construction. But it moves from 0.32 to 0.68 under reweighting,")
    print("  so G is NOT a property of the structure -- it inherits the")
    print("  propensities. Same split as schedule_invariance.py.")
    print("  Note 3:1:1 also breaks the degree law: closure outruns both")
    print("  sinks, the bank accumulates, and <k> drifts off 4.")
    print()


def part3():
    print("=" * 76)
    print("3. IS IT A SPACE?")
    print("=" * 76)
    for mode, label in (("edge", "v degree-weighted"),
                        ("node", "v uniform over nodes"),
                        ("strict", "v strictly local (deg-2 pairs only)")):
        print("  %s" % label)
        ns, ms = [], []
        for ev in EVENTS + [320000]:
            U = evolve(mode, ev, 11)
            md, dm = distances(lambda u: U.adj[u], U.n,
                               np.random.default_rng(2))
            ns.append(U.n)
            ms.append(md)
            print("     ev %7d  N %7d  mean dist %7.2f  diam %4d  max deg %5d"
                  % (ev, U.n, md, dm, max(len(a) for a in U.adj)))
        if max(ns) - min(ns) < 2:
            print("     -> FROZEN: no legal close, so no budget, so no births")
        else:
            p, rp, rl = fit(ns, ms)
            print("     -> mean dist ~ N^%.3f  (pow resid %.5f, log resid "
                  "%.5f): %s" % (p, rp, rl,
                                 "power law, d_H = %.2f" % (1 / p) if rp < rl
                                 else "LOGARITHMIC - small-world"))
        print()


def part4():
    print("=" * 76)
    print("4. IS THE RESULT EVEN ATYPICAL?  matched N and matched mean degree")
    print("=" * 76)
    U = evolve("edge", 80000, 11)
    rng = np.random.default_rng(4)
    md, dm = distances(lambda u: U.adj[u], U.n, np.random.default_rng(2))
    radj = random_regular(U.n, 4, rng)
    rmd, rdm = distances(lambda u: radj[u], U.n, np.random.default_rng(2))
    ladj, ln = square_lattice(U.n)
    lmd, ldm = distances(lambda u: ladj[u], ln, np.random.default_rng(2))
    print("  N = %d, mean degree 4" % U.n)
    print("  %-28s %12s %10s" % ("graph", "mean dist", "diameter"))
    print("  %-28s %12.2f %10d" % ("the model's output", md, dm))
    print("  %-28s %12.2f %10d" % ("random 4-regular graph", rmd, rdm))
    print("  %-28s %12.2f %10d" % ("periodic square lattice", lmd, ldm))
    print()
    print("  the model lands on the random graph, not the lattice.")
    print("  sqrt(N) = %.0f against log(N) = %.1f -- geometry is the atypical"
          % (np.sqrt(U.n), np.log(U.n)))
    print("  option, and nothing in the budget or the locality pays for it.")
    print()


def main():
    part1()
    part2()
    part3()
    part4()


if __name__ == "__main__":
    main()
