#!/usr/bin/env python3
"""
shear_damping.py

Can the frozen shear be damped, and can a LOCAL rule do it?

scale_factor.py found that pure subdivision expands but never isotropises:
the spread of a_ij across marker pairs sat at 75.4% and was identical to three
digits over a 16-fold range in N. The cause is a Yule process -- an edge
already cut into k pieces has k chances to be cut next -- so chain-length
shares follow a Dirichlet distribution that never concentrates.

Ellis and van Elst's 1+3 covariant formulation names the missing ingredient.
Their shear propagation equation reads

    sigma-dot<ab> - grad<a u-dot b> = -(2/3) Theta sigma^ab + ... - (E^ab - pi^ab/2)

and the first term on the right is the whole story: EXPANSION DAMPS SHEAR, at
a rate set by the expansion scalar Theta. In vacuum Bianchi I that integrates
to sigma ~ a^-3, so shear dies as the cube of the scale factor for free. Pure
subdivision has expansion and no damping term at all.

It also supplies a target. From their eq. (59) and the CBR bound quoted below
it, (sigma^2 - omega^2)/H^2 < 1e-3 today, so the observed universe has

    sigma/H  <~  0.03

against this model's 0.75. The discrete analogue used here is

    spread = std(a_ij) / mean(a_ij)      over marker pairs

which plays the role of sigma/Theta. It is not the same object -- a_ij is a
ratio of distances, not a rate -- so 0.03 is a benchmark rather than a
like-for-like comparison, but the question "does it decay" is exactly right.

Three growth rules
--------------------
    uniform     pick an edge uniformly. The baseline, and the Yule process.

    local       each edge carries a generation counter g, incremented on
                subdivision, and is picked with weight 2^-g. This is the
                anti-Yule rule and it is genuinely local: g is per-edge
                storage updated only where the rewrite happens, exactly like
                a provenance tag. It works because chain weight is CONSERVED
                under subdivision -- one edge of weight 2^-g becomes two of
                weight 2^-(g+1) -- so every original edge keeps total weight
                1 forever and therefore grows at an equal rate.

    global      compute every a_ij, take the least-expanded marker pair, and
                subdivide an edge on a shortest path between them. This is
                shear damping by construction and is included only to show
                what is achievable; it needs all-pairs distances, so it is
                not a local rule and is not proposed as one.

The honest tension, stated before the result: in GR the -(2/3) Theta sigma
term is not put in by hand, it falls out of the Ricci identity. Here any
damping is stipulated. The question is only whether a rule with purely local
bookkeeping can produce it, or whether damping is intrinsically global like
everything else this project has run into.

Result: a local rule damps the shear, at the central-limit rate
------------------------------------------------------------------
Averaged over independent growth histories, since single Dirichlet draws vary
enormously -- the uniform baseline ranges from 37% to 95% across seeds, which
is why one history gave 95% here and the 8-history average in scale_factor.py
gave 75%.

    uniform (5 histories)
        N        2500     5000    10000    20000    40000
        spread   67.3%    67.5%    67.7%    67.7%    67.7%      ~ N^0.00
        min/max  0.088    0.086    0.084    0.081    0.080

    local, weight 2^-g (3 histories)
        N        2500     5000    10000    20000    40000
        spread    2.8%     1.9%     1.5%     1.3%     0.7%      ~ N^-0.47
        min/max  0.902    0.932    0.938    0.953    0.974

    global control (2 histories)
        N         200      500     1000     2000
        spread    1.0%     0.4%     0.2%     0.1%               ~ N^-0.97

The uniform rule is flat to two decimal places in the exponent -- the Yule
freeze, reproduced exactly. The local rule decays as N^-0.47, which is the
central-limit rate 1/sqrt(N) expected when every original edge grows at the
same rate, and it sits BELOW the 3% CBR benchmark at every checkpoint, ending
at 0.7%. The min/max ratio rises from 0.902 to 0.974, so it is not just the
variance shrinking: every marker pair converges on the same expansion factor.

The global control reaches N^-0.97, twice as fast, which is the difference
between correcting the drift statistically and steering it directly.

What actually did the work
----------------------------
Not a force, and not an energy. The fix was to remove a rich-get-richer bias.
Under uniform picking, a chain of k edges has k chances to grow, so shares
follow a Dirichlet distribution that never concentrates. Weighting by 2^-g
makes each original chain's total weight invariant under subdivision, so every
chain grows at an equal rate and the shares concentrate as 1/sqrt(N).

That is worth stating plainly because of what it changes elsewhere. The rewrite
RULE is identical in both rows above. Only the measure on matches differs, and
it moves an observable from 67.7% frozen to 0.7% decaying. This is the sharpest
demonstration in the project of the point that the measure P(G->G') is physics
rather than bookkeeping.

It also stands in an interesting relation to the conservation results. Every
conserved quantity tested earlier -- weight, the Z3 residue, per-component
charges, Q = bank + 2N - E -- was blind to order, taking identical values on a
lattice and on its randomised twin, which is why none of them could protect
geometry. The quantity conserved here, total weight per lineage, is attached to
provenance rather than to the state, so it is not an invariant of the graph at
all. It is a property of the measure. That is precisely why it can do work the
state invariants could not.

What this does NOT fix
------------------------
The dimension. This is still pure subdivision on a 6-cycle, so the result is an
isotropic, homogeneous, ONE-dimensional universe. Against the six criteria:
isotropy and homogenisation now pass for the first time, and scale consistency
still returns d_exp = d_ball = d_H = 1.00.

So the honest position is that two of the three open criteria fell to a change
of measure, and the one that did not is the one every other route also failed:
nothing here selects a dimension.

And the damping is stipulated. In GR the -(2/3) Theta sigma term follows from
the Ricci identity; here the 2^-g weighting was chosen to cancel a Yule bias
that was already understood. The demonstration is that a rule with purely local
bookkeeping CAN isotropise, not that one does so unbidden.
"""

import numpy as np
from collections import deque


def cycle_seed(k):
    return [(i, (i + 1) % k) for i in range(k)], k


class Grower:
    """Subdivision with per-edge generation counters."""

    def __init__(self, seed_edges, n):
        self.adj = [set() for _ in range(n)]
        self.edges = []          # (u, v)
        self.gen = []            # generation of each edge
        self.n = n
        for u, v in seed_edges:
            self._add(u, v, 0)
        self.markers = list(range(n))

    def _add(self, u, v, g):
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.edges.append((u, v))
        self.gen.append(g)

    def _remove_idx(self, i):
        u, v = self.edges[i]
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        last_e, last_g = self.edges.pop(), self.gen.pop()
        if i < len(self.edges):
            self.edges[i] = last_e
            self.gen[i] = last_g

    def subdivide(self, i):
        u, v = self.edges[i]
        g = self.gen[i]
        self._remove_idx(i)
        z = self.n
        self.adj.append(set())
        self.n += 1
        self._add(u, z, g + 1)
        self._add(z, v, g + 1)
        return z

    def bfs(self, s):
        dist = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in self.adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        return dist

    def a_values(self, base):
        out, pairs = [], []
        for i in self.markers:
            d = self.bfs(i)
            for j in self.markers:
                if j > i and base.get((i, j)):
                    out.append(d[j] / base[(i, j)])
                    pairs.append((i, j))
        return np.array(out, float), pairs

    def path_edges(self, a, b):
        """Indices of edges on one shortest path from a to b."""
        d = self.bfs(a)
        if b not in d:
            return []
        path, cur = [b], b
        while cur != a:
            for y in self.adj[cur]:
                if d.get(y, 1 << 30) == d[cur] - 1:
                    path.append(y)
                    cur = y
                    break
            else:
                return []
        want = set()
        for k in range(len(path) - 1):
            x, y = path[k], path[k + 1]
            want.add((x, y) if x < y else (y, x))
        return [i for i, (u, v) in enumerate(self.edges)
                if ((u, v) if u < v else (v, u)) in want]


def baseline_distances(seed_edges, n):
    adj = [set() for _ in range(n)]
    for u, v in seed_edges:
        adj[u].add(v)
        adj[v].add(u)
    base = {}
    for i in range(n):
        dist = {i: 0}
        q = deque([i])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        for j in range(n):
            if j > i:
                base[(i, j)] = dist[j]
    return base


def run(rule, target_n, seed=7, checkpoints=(2500, 5000, 10000, 20000, 40000)):
    se, ns = cycle_seed(6)
    base = baseline_distances(se, ns)
    G = Grower(se, ns)
    rng = np.random.default_rng(seed)
    marks = sorted(c for c in checkpoints if c <= target_n)
    out = []
    while G.n < target_n:
        if rule == "uniform":
            i = int(rng.integers(len(G.edges)))
        elif rule == "local":
            w = np.exp2(-np.asarray(G.gen, dtype=np.float64))
            c = np.cumsum(w)
            i = int(np.searchsorted(c, rng.random() * c[-1]))
            i = min(i, len(G.edges) - 1)
        elif rule == "global":
            a, pairs = G.a_values(base)
            p = pairs[int(np.argmin(a))]
            cand = G.path_edges(p[0], p[1])
            i = (cand[int(rng.integers(len(cand)))] if cand
                 else int(rng.integers(len(G.edges))))
        else:
            raise ValueError(rule)
        G.subdivide(i)
        if marks and G.n >= marks[0]:
            m = marks.pop(0)
            a, _ = G.a_values(base)
            out.append((G.n, float(a.mean()), float(a.std() / a.mean()),
                        float(a.min() / a.max())))
    return out


def run_avg(rule, target_n, seeds, checkpoints):
    """Average over independent growth histories; single Dirichlet draws are
    wildly variable, so a one-seed baseline is not a fair comparison."""
    acc = {}
    for sd in seeds:
        for n, m, sp, r in run(rule, target_n, seed=sd,
                               checkpoints=checkpoints):
            acc.setdefault(n, []).append((m, sp, r))
    out = []
    for n in sorted(acc):
        v = np.array(acc[n])
        out.append((n, v[:, 0].mean(), v[:, 1].mean(), v[:, 1].min(),
                    v[:, 1].max(), v[:, 2].mean()))
    return out


def report(rule, rows):
    print("  %s" % rule)
    print("     %8s %12s %10s %16s %10s"
          % ("N", "mean a", "spread", "spread range", "min/max"))
    for n, m, sp, lo, hi, r in rows:
        print("     %8d %12.2f %9.1f%% %7.1f%%-%6.1f%% %10.3f"
              % (n, m, 100 * sp, 100 * lo, 100 * hi, r))
    if len(rows) >= 3:
        ns = np.array([r[0] for r in rows], float)
        sp = np.array([r[2] for r in rows], float)
        slope = np.polyfit(np.log(ns), np.log(sp), 1)[0]
        print("     -> spread ~ N^%.2f" % slope)
    print()


def main():
    print("=" * 78)
    print("SHEAR DAMPING: does spread/mean decay, and can a local rule do it?")
    print("=" * 78)
    print("  spread = std(a_ij)/mean(a_ij) over marker pairs, the discrete")
    print("  analogue of sigma/Theta. CBR benchmark: sigma/H <~ 0.03.")
    print()

    CK = (2500, 5000, 10000, 20000, 40000)
    report("uniform (5 histories)", run_avg("uniform", 40000, range(5), CK))
    report("local, weight 2^-g (3 histories)",
           run_avg("local", 40000, range(3), CK))

    print("  global shear-damping control (needs all-pairs distances, so it")
    print("  is NOT a local rule -- included only to show achievability):")
    report("global", run_avg("global", 2000, range(2),
                             (200, 500, 1000, 2000)))

    print("=" * 78)
    print("  A decaying spread means expansion is isotropising the graph.")
    print("  A flat spread at ~75% reproduces the Yule freeze.")
    print()


if __name__ == "__main__":
    main()
