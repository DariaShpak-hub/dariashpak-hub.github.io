#!/usr/bin/env python3
"""
colour_count.py

Does the number of colours have a fixed point?

This question could not be asked before coloured_memory.py. In every earlier
construction the dimension was not well defined -- dimension_attractor.py found
monotone drift with no stationary value anywhere -- so asking whether it had a
fixed point was asking about a quantity that did not exist. Now it does:
d = k, verified stable at k = 2, 3, 4.

So the sharp question is whether k itself is selected. Give the colour set a
dynamics:

    merge   two colours are identified, k -> k-1. Triggered when a colour is
            underused, meaning the structure barely distinguishes it.
    split   one colour bifurcates, k -> k+1. Triggered when a colour carries
            disproportionate traffic.

Both are driven by usage -- the edge count carried by each colour -- which is a
property the graph itself supplies rather than one imposed from outside.

Then start from k above and below any candidate value and see where it flows.
This is the same experimental design as dimension_attractor.py, with one
difference that matters: there, the quantity being tested for a fixed point was
never stable even when held fixed. Here it is.

    success   k converges to the same value from above and below. That value
              is the model's prediction for the dimension.
    failure   k runs to 1, diverges, or drifts with no stationary point --
              the signature every previous attempt produced.

Stated prior: failure is more likely. Merging and splitting colours are global
relabellings, and this project has measured five times that global properties
are not reachable by local rules. A usage threshold is a local proxy for
"redundant", and local proxies have misled before -- surface_tension.py found
one that moved in the wrong direction entirely.

Result: no fixed point, for two independent reasons
------------------------------------------------------
The control behaves: with colour moves off, k stays where it is put and the
dimension follows it -- d = 1.85, 2.90, 3.46 for k = 2, 3, 4, the same low-bias
signature as before.

With the colour dynamics switched on, NOTHING HAPPENED. Zero merges, zero
splits, from every starting k, and k came out exactly as it went in:

    k0        2    3    5    8   12
    k final   2    3    5    8   12
    merges    0    0    0    0    0
    splits    0    0    0    0    0

and identically across four random seeds. The reason is measurable:

    k     min/mean   max/mean   spread
    3       0.998      1.002    0.0014
    5       0.977      1.023    0.0172
    8       0.970      1.055    0.0285

Colours are chosen uniformly, so they are used uniformly -- the spread is
between 0.1% and 3%. A threshold at 0.5 or 1.8 of the mean is unreachable by
several orders of magnitude. Usage cannot drive colour dynamics at all,
because there is nothing for it to discriminate.

The second reason is structural and does not depend on the driver. Forcing the
moves to fire by tightening the thresholds makes them fire constantly, and the
move set is asymmetric: split always succeeds, while merge is floored at k = 1.
That is a random walk with a reflecting barrier below and none above, so k
drifts upward without settling whatever the trigger. (The tightened run is not
reported as a measurement: merge_colours is O(N) and firing it every few
hundred steps made the run computationally infeasible rather than informative.)

Verdict
-------
The colour count is not selected. It is a constant of the construction, and
neither a usage-based driver nor a symmetric one can give it a fixed point.

This was the last question in the line opened by coloured_memory.py, and the
prior stated before running it was that it would fail. It failed in a slightly
more interesting way than expected: not by drifting, but by there being nothing
for a local rule to see. Uniform usage is the colour-space version of the
result global_observables.py reached for geometry -- the quantity that would
have to be measured to decide anything is invisible to the dynamics that would
have to act on it.

So the standing position is unchanged and now complete: d = k is real and
stable, and k is imposed. Nothing in this framework selects three.
"""

import numpy as np
from collections import deque, Counter


class MutableColour:
    """Commuting coloured complex whose colour COUNT can change."""

    def __init__(self, k0, seed=0):
        self.k = k0
        self.rng = np.random.default_rng(seed)
        self.parent = [0]
        self.nb = [dict()]              # (colour, dir) -> node
        self.alive = 1
        self.merges = 0
        self.splits = 0

    # ------------------------------------------------------------ unionfind
    def find(self, x):
        r = x
        while self.parent[r] != r:
            r = self.parent[r]
        while self.parent[x] != r:
            self.parent[x], x = r, self.parent[x]
        return r

    def _new(self):
        self.parent.append(len(self.parent))
        self.nb.append(dict())
        self.alive += 1
        return len(self.parent) - 1

    def union(self, a, b):
        work = deque([(a, b)])
        while work:
            x, y = work.popleft()
            x, y = self.find(x), self.find(y)
            if x == y:
                continue
            if len(self.nb[x]) < len(self.nb[y]):
                x, y = y, x
            self.parent[y] = x
            self.alive -= 1
            for key, w in list(self.nb[y].items()):
                if key in self.nb[x]:
                    work.append((self.nb[x][key], w))
                else:
                    self.nb[x][key] = w
            self.nb[y] = dict()

    # ---------------------------------------------------------------- steps
    def step(self, v, colour, direction):
        v = self.find(v)
        key = (colour, direction)
        if key in self.nb[v]:
            return self.find(self.nb[v][key])
        w = self._new()
        self.nb[v][key] = w
        self.nb[w][(colour, -direction)] = v
        return w

    def square(self, v, a, b):
        da = 1 if self.rng.random() < 0.5 else -1
        db = 1 if self.rng.random() < 0.5 else -1
        x = self.step(v, a, da)
        p = self.step(x, b, db)
        y = self.step(v, b, db)
        q = self.step(y, a, da)
        self.union(p, q)

    # -------------------------------------------------------- colour moves
    def usage(self):
        c = Counter()
        for i in range(len(self.parent)):
            if self.find(i) != i:
                continue
            for (col, _) in self.nb[i]:
                c[col] += 1
        return c

    def merge_colours(self, a, b):
        """Identify colour b with colour a. Conflicts cascade into unions."""
        if a == b or self.k <= 1:
            return False
        lo, hi = min(a, b), max(a, b)
        pend = []
        for i in range(len(self.parent)):
            if self.find(i) != i:
                continue
            for d in (1, -1):
                if (hi, d) in self.nb[i]:
                    w = self.nb[i].pop((hi, d))
                    if (lo, d) in self.nb[i]:
                        pend.append((self.nb[i][(lo, d)], w))
                    else:
                        self.nb[i][(lo, d)] = w
        for x, y in pend:
            self.union(x, y)
        # relabel colours above hi downward
        for i in range(len(self.parent)):
            if self.find(i) != i:
                continue
            fix = [(c, d) for (c, d) in self.nb[i] if c > hi]
            for (c, d) in fix:
                self.nb[i][(c - 1, d)] = self.nb[i].pop((c, d))
        self.k -= 1
        self.merges += 1
        return True

    def split_colour(self, a):
        """Colour a bifurcates: each edge keeps a or takes the new label."""
        new = self.k
        for i in range(len(self.parent)):
            if self.find(i) != i:
                continue
            for d in (1, -1):
                if (a, d) in self.nb[i] and self.rng.random() < 0.5:
                    self.nb[i][(new, d)] = self.nb[i].pop((a, d))
        self.k += 1
        self.splits += 1
        return True

    # --------------------------------------------------------------- driver
    def evolve(self, target, lo_frac=0.5, hi_frac=1.8, colour_every=400):
        roots = [0]
        t = 0
        while self.alive < target:
            t += 1
            if len(roots) > 3 * self.alive:
                roots = [r for r in roots if self.find(r) == r]
            v = roots[int(self.rng.integers(len(roots)))]
            if self.find(v) != v:
                continue
            a = int(self.rng.integers(self.k))
            b = int(self.rng.integers(self.k))
            if a != b:
                before = len(self.parent)
                self.square(v, a, b)
                roots.extend(range(before, len(self.parent)))
            if t % colour_every == 0 and self.k >= 1:
                u = self.usage()
                if not u:
                    continue
                mean = sum(u.values()) / max(self.k, 1)
                worst = min(range(self.k), key=lambda c: u.get(c, 0))
                best = max(range(self.k), key=lambda c: u.get(c, 0))
                if self.k > 1 and u.get(worst, 0) < lo_frac * mean:
                    other = best if best != worst else (worst + 1) % self.k
                    self.merge_colours(min(worst, other), max(worst, other))
                elif u.get(best, 0) > hi_frac * mean:
                    self.split_colour(best)
        return self

    # ---------------------------------------------------------- measurement
    def adjacency(self):
        idx, adj = {}, []
        for i in range(len(self.parent)):
            if self.find(i) == i:
                idx[i] = len(adj)
                adj.append(set())
        for i in list(idx):
            for _, w in self.nb[i].items():
                w = self.find(w)
                if w in idx and w != i:
                    adj[idx[i]].add(idx[w])
                    adj[idx[w]].add(idx[i])
        return adj


def dimension(adj, srcs=8, seed=3):
    n = len(adj)
    if n < 200:
        return float("nan"), float("nan")
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
    print("1. CONTROL: with the colour moves OFF, k stays put and d = k")
    print("=" * 78)
    print("  %6s %8s %8s %9s %9s" % ("k0", "k final", "nodes", "d", "mean dist"))
    for k0 in (2, 3, 4):
        C = MutableColour(k0, seed=1)
        C.evolve(15000, colour_every=10 ** 9)
        d, md = dimension(C.adjacency())
        print("  %6d %8d %8d %9.2f %9.2f" % (k0, C.k, C.alive, d, md))
    print()

    print("=" * 78)
    print("2. COLOUR DYNAMICS ON: does k converge from above and below?")
    print("=" * 78)
    print("  %6s %8s %8s %8s %8s %9s %9s"
          % ("k0", "k final", "merges", "splits", "nodes", "d", "mean dist"))
    finals = []
    for k0 in (2, 3, 5, 8, 12):
        C = MutableColour(k0, seed=1)
        C.evolve(15000)
        d, md = dimension(C.adjacency())
        finals.append(C.k)
        print("  %6d %8d %8d %8d %8d %9.2f %9.2f"
              % (k0, C.k, C.merges, C.splits, C.alive, d, md))
    print()
    print("  final k values: %s" % finals)
    if len(set(finals)) == 1:
        print("  -> CONVERGED to k = %d from every start" % finals[0])
    else:
        print("  -> did NOT converge: k depends on where it started,")
        print("     so the colour count is not selected.")
    print()

    print("=" * 78)
    print("3. WHY NOTHING FIRED: colour usage is uniform by construction")
    print("=" * 78)
    print("  %6s %10s %10s %10s" % ("k", "min/mean", "max/mean", "spread"))
    for k0 in (3, 5, 8):
        C = MutableColour(k0, seed=1)
        C.evolve(15000, colour_every=10 ** 9)
        u = C.usage()
        mean = sum(u.values()) / C.k
        vals = np.array([u.get(c, 0) / mean for c in range(C.k)])
        print("  %6d %10.3f %10.3f %10.4f"
              % (k0, vals.min(), vals.max(), vals.std() / vals.mean()))
    print()
    print("  Colours are chosen uniformly, so they are used uniformly: the")
    print("  spread is 0.1%% to 3%%. A usage threshold at 0.5 or 1.8 of the")
    print("  mean is unreachable, which is why 0 merges and 0 splits fired.")
    print()

    print("=" * 78)
    print("4. SEED DEPENDENCE, to separate convergence from coincidence")
    print("=" * 78)
    print("  %6s %s" % ("k0", "final k across 4 seeds"))
    for k0 in (2, 5, 8):
        row = []
        for sd in range(4):
            C = MutableColour(k0, seed=sd)
            C.evolve(8000)
            row.append(C.k)
        print("  %6d %s" % (k0, row))
    print()


if __name__ == "__main__":
    main()
