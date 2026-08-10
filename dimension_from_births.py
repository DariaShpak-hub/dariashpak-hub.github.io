#!/usr/bin/env python3
"""
dimension_from_births.py

Testing d = m + 1: does the number of edges a node birth creates fix the
dimension?

The claim, derived rather than fitted. Let Q = B + cN - E be invariant and let
a birth have dN = +1, dE = +m, dB = -1. Then

    dQ = -1 + c - m = 0    =>    c = m + 1

and since E = B + cN - Q with the bank bounded, the mean degree is

    <k> = 2E/N -> 2c = 2(m + 1)

There is a cleaner way to see the +1. Over a long run the bank returns to
where it started, so closes must fund births one for one. Each birth therefore
costs the graph m edges of its own PLUS one edge from the close that paid for
it: E/N -> m + 1. The +1 is the price of the node.

Separately, a k-colour commuting lattice -- the Cayley graph of Z^k, which
coloured_memory.py verified gives d = k -- has exactly 2k edges per node.
Setting the two degrees equal,

    2(m + 1) = 2k    =>    d = k = m + 1

So the prediction is that the birth rule fixes the dimension. Our earlier
budget used c = 2, hence <k> = 4, hence a TWO-dimensional target, which is why
4 kept appearing in paid_expansion.py and soc.py.

What this can and cannot show
-------------------------------
The arithmetic half is a genuine prediction and is tested directly below:
vary m, and <k> should follow 2(m+1) with the conservation residual staying 0.

The geometric half is where the caution belongs. Getting the right degree is
necessary and nowhere near sufficient -- a random 6-regular graph also has
degree 6 and is an expander with no dimension at all, as global_observables.py
measured. So the test below deliberately measures d as well, and the
expectation is that the uncoloured budget hits the degree and MISSES the
dimension. That would locate exactly what the colour structure is for: the
arithmetic buys the degree, the commutation buys the geometry, and they agree
only when c = k.

Result 1: the arithmetic is exact
-----------------------------------
    m   c=m+1      N   residual    <k>   predicted    E/N       d
    1     2      3311      0      4.000      4       2.000    6.10
    2     3     20000      0      5.999      6       2.999    6.84
    3     4     20000      0      7.998      8       3.999    6.69
    4     5     20000      0      9.996     10       4.998    7.14

Conservation residual is 0 everywhere and the mean degree lands on 2(m+1) to
four significant figures. E/N lands on m+1 exactly, confirming the accounting:
m edges from the birth, one from the close that funded it.

So the budget relation c = m + 1 and the degree relation <k> = 2c are both
right. That part of the derivation stands.

Result 2: the dimension does not follow, and the bridge is WRONG
-----------------------------------------------------------------
The measured dimension in part 1 is 6.10, 6.84, 6.69, 7.14 -- flat, and nowhere
near m+1 = 2, 3, 4, 5. Setting the degree does not set the geometry. That much
was expected: a random 6-regular graph also has degree 6 and no dimension.

The unexpected part is in the colour model:

    k      N      <k>   predicted    E/N      d
    2    30000   3.894      4       1.947   1.84
    3    30001   3.947      6       1.974   2.84
    4    30000   3.114      8       1.557   3.81

The dimension is right -- 1.84, 2.84, 3.81 for k = 2, 3, 4, as
coloured_memory.py found. But the DEGREE is wrong. It sits near 4 regardless
of k, against the predicted 2k, and E/N stays near 2 instead of tracking k.

That kills the bridge. The argument was: the budget forces degree 2c, a
k-colour lattice has degree 2k, therefore c = k and the birth rule fixes the
dimension. The middle step is false for a GROWN lattice. Only a pristine,
complete Z^k has degree 2k; the grown complex is a sparse spanning subgraph of
it -- edges are created lazily, only when a square demands them -- and a sparse
subgraph of Z^k keeps the dimension k while having far fewer edges. This is
familiar behaviour: a supercritical percolation cluster on Z^d also has
Hausdorff dimension d at much lower degree.

Verdict: I proposed d = m + 1 last turn and it does not hold
--------------------------------------------------------------
What survives is the arithmetic: c = m + 1 exactly, <k> = 2(m+1) exactly,
residual 0. What fails is the inference from degree to dimension, in both
directions at once -- the budgeted model has the right degree and the wrong
dimension, and the coloured model has the right dimension and the wrong degree.

Dimension in this framework comes from the commutation structure alone and is
independent of the degree. The budget therefore does not select it, and "why
three dimensions" is NOT answered by "because a birth creates two edges".

The standing position is unchanged from coloured_memory.py: the colour count
sets the dimension, and the colour count is imposed. The one thing this rules
out is a whole class of answers -- anything that tries to fix the dimension
through the edge or degree accounting, since that accounting is measurably
decoupled from the geometry.
"""

import numpy as np
from collections import deque

from coloured_memory import ColouredComplex


class BudgetGrowth:
    """Q = B + (m+1)N - E invariant; births attach to m nearby nodes."""

    def __init__(self, m, n0=12, seed=0):
        self.m = m
        self.c = m + 1
        self.rng = np.random.default_rng(seed)
        self.adj = [set() for _ in range(n0)]
        self.n = n0
        for i in range(n0):
            self._add(i, (i + 1) % n0)
            self._add(i, (i + 2) % n0)
        self.bank = 0
        self.Q = self.bank + self.c * self.n - self.n_edges()
        self.counts = {"close": 0, "birth": 0}

    def _add(self, u, v):
        if u != v:
            self.adj[u].add(v)
            self.adj[v].add(u)

    def n_edges(self):
        return sum(len(a) for a in self.adj) // 2

    def residual(self):
        return self.bank + self.c * self.n - self.n_edges() - self.Q

    def mean_degree(self):
        return 2 * self.n_edges() / self.n

    def close(self):
        """Join two nodes at distance 2. dE=+1, dB=+1."""
        for _ in range(60):
            v = int(self.rng.integers(self.n))
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            a = nb[int(self.rng.integers(len(nb)))]
            b = nb[int(self.rng.integers(len(nb)))]
            if a == b or b in self.adj[a]:
                continue
            self._add(a, b)
            self.bank += 1
            self.counts["close"] += 1
            return True
        return False

    def birth(self):
        """A node is born attached to m mutually-nearby nodes. dN=+1, dE=+m."""
        if self.bank < 1:
            return False
        for _ in range(60):
            v = int(self.rng.integers(self.n))
            pool = set(self.adj[v]) | {v}
            for u in list(self.adj[v]):
                pool |= self.adj[u]
            pool = list(pool)
            if len(pool) < self.m:
                continue
            pick = self.rng.choice(len(pool), size=self.m, replace=False)
            targets = [pool[int(i)] for i in pick]
            w = self.n
            self.adj.append(set())
            self.n += 1
            for t in targets:
                self._add(w, t)
            self.bank -= 1
            self.counts["birth"] += 1
            return True
        return False

    def grow(self, target):
        while self.n < target:
            if self.bank == 0:
                if not self.close():
                    break
            else:
                if not self.birth():
                    if not self.close():
                        break
        return self


def shells_dimension(adj, n, srcs=8, seed=3):
    rng = np.random.default_rng(seed)
    acc, tot, cnt = {}, 0.0, 0
    for _ in range(srcs):
        s = int(rng.integers(n))
        d = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in d:
                    d[w] = d[u] + 1
                    q.append(w)
        v = np.fromiter(d.values(), int)
        tot += float(v.sum())
        cnt += len(v)
        c = np.bincount(v)
        for r in range(len(c)):
            acc[r] = acc.get(r, 0.0) + float(c[r])
    sh = {r: x / srcs for r, x in acc.items()}
    hi = max(4, max(sh) // 3)
    rs = [r for r in sorted(sh) if 2 <= r <= hi and sh[r] > 0]
    if len(rs) < 3:
        return float("nan"), tot / cnt
    x = np.log(np.array(rs, float))
    y = np.log(np.array([sh[r] for r in rs]))
    return 1.0 + float(np.polyfit(x, y, 1)[0]), tot / cnt


def main():
    print("=" * 78)
    print("1. THE ARITHMETIC HALF: does <k> follow 2(m+1)?")
    print("=" * 78)
    print("  %5s %5s %8s %9s %10s %11s %9s %10s"
          % ("m", "c=m+1", "N", "residual", "<k>", "predicted", "E/N", "d"))
    for m in (1, 2, 3, 4):
        G = BudgetGrowth(m, seed=1).grow(20000)
        d, md = shells_dimension(G.adj, G.n)
        print("  %5d %5d %8d %9d %10.3f %11d %9.3f %10.2f"
              % (m, m + 1, G.n, G.residual(), G.mean_degree(),
                 2 * (m + 1), G.n_edges() / G.n, d))
    print()
    print("  E/N should land on m+1: m edges from the birth plus one from the")
    print("  close that funded it.")
    print()

    print("=" * 78)
    print("2. THE GEOMETRIC HALF: does the colour model match those degrees?")
    print("=" * 78)
    print("  %5s %9s %10s %11s %9s %10s"
          % ("k", "N", "<k>", "predicted", "E/N", "d"))
    for k in (2, 3, 4):
        C = ColouredComplex(k, seed=1).grow(30000, 1.0)
        adj = C.adjacency()
        n = len(adj)
        e = sum(len(a) for a in adj) // 2
        d, md = shells_dimension(adj, n)
        print("  %5d %9d %10.3f %11d %9.3f %10.2f"
              % (k, n, 2 * e / n, 2 * k, e / n, d))
    print()

    print("=" * 78)
    print("3. VERDICT")
    print("=" * 78)
    print("  If column <k> tracks 2(m+1) in part 1 and 2k in part 2, the")
    print("  arithmetic is confirmed and c = k = m+1 is consistent.")
    print("  If d tracks the prediction in part 2 but NOT in part 1, then the")
    print("  budget fixes the degree while only the colours fix the geometry.")
    print()


if __name__ == "__main__":
    main()
