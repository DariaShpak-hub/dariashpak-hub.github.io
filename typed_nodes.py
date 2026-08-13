#!/usr/bin/env python3
"""
typed_nodes.py

Dimension from node types: three kinds of thing, therefore three kinds of
connection, therefore three dimensions?

coloured_memory.py put colours on EDGES and found d = k -- but k was a free
parameter, so the dimension was a dial. This puts colours on NODES instead,
which is the proposal that particles come in types. The edge colour is then not
stipulated: it is DERIVED as the unordered pair of the two endpoint types.

    q node types  ->  q(q-1)/2 edge colours

    q     2    3    4    5
    cols  1    3    6   10

And that counting has a fixed point. Demand that the number of ways things can
connect equals the number of kinds of thing -- that the type set can name its
own relations without needing more names than it has:

    q(q-1)/2 = q    =>    q^2 = 3q    =>    q = 3

Unique nonzero solution. Three is the only count where connection-types and
thing-types agree. That is a selection argument rather than a dial, which is
what every previous route in this project failed to supply. It is still an
adopted principle, not a theorem -- but it is a much better class of
assumption than "let k = 3".

Why the square relation is the WRONG one here
-----------------------------------------------
Working out what the structure actually is changes the construction. From a
type-a node the available colours are {a,b} and {a,c}. Take a step and the type
changes, so the options change with it. Trace the square:

    a --{a,b}--> b --{b,c}--> c        ends on type c
    a --{a,c}--> c --{c,b}--> b        ends on type b

Different types, so the two paths cannot be identified. The square does not
close, and the commutation relation that coloured_memory.py enforced is simply
inapplicable. What closes instead is

    a --> b --> c --> a

a TRIANGLE in type space. The generators available depend on where you are, so
this is a groupoid rather than a group, and the natural relations are loops
around the type triangle, not squares.

Two candidate relations are therefore tested, both type-consistent:

    triangle   a->b->c->a returns to the node it started from
    hexagon    a->b->c->a and a->c->b->a agree with each other

The prediction to beat is d = 3 at q = 3. The honest alternatives include
d = q-1, since a node only ever has q-1 colours available to it, and a
non-integer value if the structure is fractal rather than manifold-like.

Result: the fixed point is real arithmetic, but it does not give d = 3
------------------------------------------------------------------------
Triangle relation -- total collapse:

    q         2      3      4      5
    nodes     1      3      4      5

Demanding that a->b->c->a return to its starting node kills the structure
outright. The universe collapses to exactly q nodes, one per type. The relation
is too strong: it makes the loop the identity, so there is nowhere to go.

Hexagon relation -- it grows, and q = 3 gives a thread:

    q   colours   nodes      d    mean dist   degree
    2      1          1    nan        nan       nan
    3      3       4022   0.89     438.95      3.98
    4      6      20002   3.26      25.30      2.84
    5     10      20002   3.54      25.06      2.48

At q = 3, the value the fixed point selects, the measured dimension is 0.89 --
a thread. Mean distance is 439 across 4022 nodes, about N/9, which is what a
long chain looks like. Not 3, and not even 2.

The q = 4 and q = 5 rows read 3.26 and 3.54, but those carry the familiar
signature of no scaling range rather than of real geometry: mean distance 25 at
N = 20000 is small-world, and the value drifting toward ~3.5 is exactly where
the shell estimator saturates when there is nothing to measure.

Why it fails, precisely
-------------------------
The fixed point counts TOTAL colours. The dimension is set by the colours
available AT A NODE, and those are different numbers.

    total colours available anywhere   q(q-1)/2
    colours available at one node      q - 1

A node of type a can only use {a,b} and {a,c}. The third colour {b,c} exists in
the system but is invisible from a. So at q = 3 there are three colours globally
and only two locally, and two local directions cannot support three dimensions.
The measured 0.89 is even below that, because the hexagon relation constrains
the two remaining directions against each other.

Which makes the two requirements provably incompatible:

    fixed point       q(q-1)/2 = q     =>  q = 3,  local directions = 2
    d = 3 needs       q - 1 = 3        =>  q = 4,  total colours = 6 != 4

There is no q satisfying both. The self-consistency condition that picks out 3
and the condition that gives three dimensions select different numbers, and
they cannot be reconciled by any choice of q.

Verdict
-------
The idea was worth testing and the arithmetic behind it is real -- q = 3 is
genuinely the unique fixed point of q(q-1)/2 = q, and that remains the only
non-arbitrary "3" anywhere in this project. But the geometry does not follow
from it. Three types of thing give three kinds of connection globally and only
two usable directions locally, and it is the local count that sets the
dimension.

So node-typing does not select the dimension either. What it does do is
sharpen the target: any mechanism that is going to work has to make the
LOCALLY AVAILABLE direction count come out at 3, not the global label count.
"""

import numpy as np
from collections import deque


class TypedComplex:
    """Nodes carry a type; an edge's colour is the pair of endpoint types."""

    def __init__(self, q, seed=0):
        self.q = q
        self.rng = np.random.default_rng(seed)
        self.parent = [0]
        self.typ = [0]
        self.nb = [dict()]          # (other_type, direction) -> node
        self.alive = 1
        self.type_clash = 0

    def find(self, x):
        r = x
        while self.parent[r] != r:
            r = self.parent[r]
        while self.parent[x] != r:
            self.parent[x], x = r, self.parent[x]
        return r

    def _new(self, t):
        self.parent.append(len(self.parent))
        self.typ.append(t)
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
            if self.typ[x] != self.typ[y]:
                self.type_clash += 1        # must never happen; counted
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

    def step(self, v, target_type, direction):
        """From v (type a), follow the {a,target} edge; create if absent."""
        v = self.find(v)
        a = self.typ[v]
        if target_type == a:
            return v
        key = (target_type, direction)
        if key in self.nb[v]:
            return self.find(self.nb[v][key])
        w = self._new(target_type)
        self.nb[v][key] = w
        self.nb[w][(a, -direction)] = v
        return w

    def relation(self, v, b, c, mode):
        """Close a loop around the type triangle a -> b -> c -> a."""
        v = self.find(v)
        a = self.typ[v]
        if len({a, b, c}) < 3:
            return
        d1 = 1 if self.rng.random() < 0.5 else -1
        d2 = 1 if self.rng.random() < 0.5 else -1
        d3 = 1 if self.rng.random() < 0.5 else -1
        x = self.step(v, b, d1)
        y = self.step(x, c, d2)
        z = self.step(y, a, d3)
        if mode == "triangle":
            self.union(z, v)
        else:                                   # hexagon: the other order too
            p = self.step(v, c, d1)
            r = self.step(p, b, d2)
            s = self.step(r, a, d3)
            self.union(z, s)

    def grow(self, target, mode):
        roots = [0]
        guard = 0
        while self.alive < target and guard < 3_000_000:
            guard += 1
            if len(roots) > 3 * self.alive:
                roots = [r for r in roots if self.find(r) == r]
            v = roots[int(self.rng.integers(len(roots)))]
            if self.find(v) != v:
                continue
            b = int(self.rng.integers(self.q))
            c = int(self.rng.integers(self.q))
            before = len(self.parent)
            self.relation(v, b, c, mode)
            roots.extend(range(before, len(self.parent)))
        return self

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


def measure(adj, srcs=8, seed=3):
    n = len(adj)
    if n < 100:
        return float("nan"), float("nan"), n, float("nan")
    rng = np.random.default_rng(seed)
    acc, tot, cnt = {}, 0.0, 0
    for _ in range(srcs):
        s = int(rng.integers(n))
        d = {s: 0}
        qq = deque([s])
        while qq:
            u = qq.popleft()
            for w in adj[u]:
                if w not in d:
                    d[w] = d[u] + 1
                    qq.append(w)
        v = np.fromiter(d.values(), int)
        tot += float(v.sum())
        cnt += len(v)
        c = np.bincount(v)
        for r in range(len(c)):
            acc[r] = acc.get(r, 0.0) + float(c[r])
    sh = {r: x / srcs for r, x in acc.items()}
    hi = max(4, max(sh) // 3)
    rs = [r for r in sorted(sh) if 2 <= r <= hi and sh[r] > 0]
    deg = sum(len(a) for a in adj) / n
    if len(rs) < 3:
        return float("nan"), tot / cnt, n, deg
    x = np.log(np.array(rs, float))
    y = np.log(np.array([sh[r] for r in rs]))
    return 1.0 + float(np.polyfit(x, y, 1)[0]), tot / cnt, n, deg


def main():
    print("=" * 78)
    print("THE FIXED POINT")
    print("=" * 78)
    print("  %5s %14s %10s" % ("q", "q(q-1)/2 colours", "equal?"))
    for q in range(2, 7):
        cols = q * (q - 1) // 2
        print("  %5d %14d %10s" % (q, cols, "YES" if cols == q else ""))
    print()
    print("  q = 3 is the unique nonzero solution of q(q-1)/2 = q.")
    print()

    for mode in ("triangle", "hexagon"):
        print("=" * 78)
        print("RELATION: %s" % mode.upper())
        print("=" * 78)
        print("  %5s %8s %9s %10s %10s %10s %9s"
              % ("q", "colours", "nodes", "d", "mean dist", "degree", "clashes"))
        for q in (2, 3, 4, 5):
            T = TypedComplex(q, seed=1).grow(20000, mode)
            adj = T.adjacency()
            d, md, n, deg = measure(adj)
            print("  %5d %8d %9d %10.2f %10.2f %10.2f %9d"
                  % (q, q * (q - 1) // 2, n, d, md, deg, T.type_clash))
        print()

    print("=" * 78)
    print("  d = 3 at q = 3 would mean the fixed point selects the dimension.")
    print("  d = q-1 would mean the dimension counts the colours available AT")
    print("  a node, not the colours in total, and q = 3 would give d = 2.")
    print()


if __name__ == "__main__":
    main()
