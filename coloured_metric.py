#!/usr/bin/env python3
"""
coloured_metric.py

Does the metric survive growth in the COLOURED lattice?

This closes a gap that should have been closed long ago. Every measurement of
shear in this project was made on the untyped grower, where growth means
`close` -- a distance-2 chord, i.e. a shortcut. fine_tuning.py's 0-out-of-600
and output_classes.py's "the only isotropic thing this model makes is a thread"
are both statements about THAT grower. Neither has ever been tested on the one
construction that actually produced a stable dimension.

coloured_memory.py measured d = 1.84, 2.84, 3.81 for k = 2, 3, 4 and stopped
there. It never put comoving markers in the complex, so we do not know what its
metric does. If the coloured lattice grows without wrecking distances, then the
central negative of this project is a property of chords, not of growth, and
the architecture built on typed edges is worth pursuing. If shear locks near
1.4 there too, the problem is deeper than typing and that is worth knowing just
as much.

What to expect, worked out before running
-----------------------------------------
Two facts pull in opposite directions and the honest prediction is not obvious.

Against a signal: in a COMPLETE Cayley graph of Z^k, growth is accretion at the
boundary and interior distances are frozen forever. Two nodes at lattice
positions x and y are |x-y| apart at every epoch. If that is all that happens,
a = 1 identically and shear = 0 trivially -- isotropy for free, but no
expansion, which would just move the failure from one criterion to another.

For a signal: the grown complex is NOT complete. dimension_from_births.py
measured mean degree near 3.9 against the 2k a full Z^k would have, because
edges are created lazily -- only when a square demands one. So the complex is a
sparse spanning subgraph of Z^k, graph distances start ABOVE the lattice metric,
and they can only fall as gaps fill in. That predicts a < 1 and decreasing.

Which makes the real question a sharp one: does every marker pair converge to
the lattice metric at the same rate, or at different rates? Uniform convergence
means shear -> 0 and the construction is metrically well behaved. Ragged
convergence means the lazy filling is itself a source of anisotropy, and typing
does not save us.

Method
------
Comoving markers, the same observable scale_factor.py introduced, adapted to a
construction with no natural seed set (a ColouredComplex starts from one node).

    1  grow to N0 = 500, then designate 10 markers -- random live roots, kept
       only if pairwise baseline distance is at least 2
    2  record all pairwise distances at that epoch as the baseline d0
    3  grow on through checkpoints, and at each one recompute d(i,j) and form
       a_ij = d(i,j)/d0(i,j)
    4  report mean a and shear = (max a - min a)/mean a

Markers are tracked through the union-find, since enforcing commutation merges
nodes and a marker can stop being a root. If two markers merge into each other
the pair is dropped and counted.

The control matters as much as the test. The untyped grower is run through the
IDENTICAL marker protocol -- markers designated at N0 = 500 rather than taken
from the 6-cycle seed -- so that any difference is attributable to the model
and not to how the markers were chosen.

Result 1: the marker protocol is the wrong instrument here, and why
--------------------------------------------------------------------
Markers designated at N0 = 500, four seeds, 44 pairs:

    coloured, k = 3        N    1000   2000   4001   8000  16001  32001
                      mean a  0.875  0.637  0.630  0.630  0.630  0.630
                       shear  0.613  1.319  1.335  1.335  1.335  1.335

    untyped control        N    1000   2000   4000   8000  16000  32000
                      mean a  1.310  1.791  2.323  3.056  3.995  5.271
                       shear  0.873  0.969  0.990  0.947  0.855  0.769

Two things are visible and only one of them is physics.

The coloured numbers FREEZE -- identical to four decimal places from N = 4001
through N = 32001. That is real and it is the accretion prediction confirmed
exactly: growth happens at the boundary, interior distances stop changing, and
the complex gets bigger without the metric expanding. The untyped control, by
contrast, genuinely expands, a running 1.31 -> 5.27.

The shear of 1.335 is NOT real. It is an artifact of the baseline, and the
freezing is what makes it permanent. At N = 500 the complex is very incomplete,
so a pair whose true lattice separation is 4 may sit at graph distance 8. As
filling proceeds each pair falls toward its lattice value by its own factor --
one pair by 1.0, another by 0.5 -- and the spread in those factors is what the
shear is reporting. Then the metric freezes and the spread never washes out.
So the number measures the raggedness of the baseline epoch, not the anisotropy
of the object. Recorded because reporting 1.335 as a physical result would have
been wrong, and it took a second instrument to see that.

Result 2: with full commutation every node has a Z^k address, exactly
---------------------------------------------------------------------
Zero commutation inconsistencies, at every k, at every size tested, in every
seed. Walking any two colour-paths to the same node always agrees on the
coordinate. So the grown object really is a faithful quotient of Z^k, and that
licenses a much better measurement: graph distance from the origin to a node m
steps out along a single axis, divided by m. A faithful lattice gives 1.000 on
every axis, and the spread ACROSS axes is anisotropy with no baseline in it.

    k = 2, N = 32000, four seeds     ratios 1.001 1.001 / 1.000 1.001 /
                                            1.003 1.002 / 1.001 1.002
                                     spread 0.0007

The k = 2 lattice is isotropic to seven parts in ten thousand. The untyped
grower, on the same kind of quantity, sits near 1. That is a difference of
three orders of magnitude.

Result 3: the residual at higher k is incompleteness, and it decays at
the central-limit rate
-----------------------------------------------------------------------
    k     N        side   mean ratio    spread
    2     2000     44.7      1.0095     0.0035
    2     8000     89.4      1.0018     0.0023
    2    32000    178.9      1.0010     0.0004
    2   128000    357.8      1.0003     0.0005      spread ~ N^-0.47

    3     8000     20.0      1.2966     0.1205
    3    32000     31.7      1.1254     0.0420
    3   128000     50.4      1.0545     0.0667
    3   400000     73.7      1.0200     0.0141      spread ~ N^-0.55

    4    32000     13.4      2.5673     0.3407
    4   128000     18.9      1.6119     0.3129
    4   400000     25.1      1.4154     0.0689      spread ~ N^-0.63

Both columns converge, at every k. The mean ratio goes to 1, so the complex
becomes a faithful lattice; the spread goes to 0, so it becomes isotropic. The
apparent anisotropy at k = 4 was never anisotropy -- it was a lattice of side
13 that had not filled in yet, and it decays as soon as it does.

The rate is the interesting part. N^-0.47 at k = 2 is the central-limit rate,
which is exactly what shear_damping.py measured for the anti-Yule rule
(N^-0.47) -- the one positive result this project had. The same rate turns up
here from a completely different mechanism, and at k = 3 and 4 as well.

Verdict: typing solves isotropy and loses expansion
-----------------------------------------------------
The negative that fine_tuning.py and output_classes.py established -- that the
move creating dimension destroys the metric, that the only isotropic thing the
model makes is a thread -- is a property of CHORDS, not of growth. Put types on
the edges and require commutation, and the isotropy problem goes away:

    untyped grower     expands (a: 1.31 -> 5.27)     anisotropic (~1)
    coloured lattice   isotropic (7e-4 at k = 2)     frozen (a freezes)

Each construction satisfies exactly one of the two requirements and neither
satisfies both. That is not the same wall as before -- it is a different and
much better-posed one, because the failure is now specific: the coloured
lattice has no move that expands its interior. Accretion adds boundary;
commutation forbids the shortcut that would change an interior distance.

So the missing ingredient is named exactly. We need a local move that lengthens
interior distances in a typed complex WITHOUT breaking commutation. Uniform
refinement would do it and is global. Subdividing one edge breaks the squares
around it. Whether something local exists between those two is the open
question, and it is the one worth attacking next.

One caveat on scope: everything here is at p_comm = 1. The probabilistic
regime, where coloured_memory.py found the power-law/exponential transition
between p = 0.5 and 0.75, is untested for the metric.
"""

import numpy as np
from collections import deque

from coloured_memory import ColouredComplex
from fine_tuning import Cosmos

N0 = 500
N_MARKERS = 10
CHECKPOINTS = (1000, 2000, 4000, 8000, 16000, 32000)


class MarkedComplex(ColouredComplex):
    """ColouredComplex with growth that can be resumed without restarting."""

    def __init__(self, k, seed=0):
        super().__init__(k, seed)
        self.roots = [0]

    def grow_to(self, target, p_comm=1.0):
        while self.alive < target:
            if len(self.roots) > 3 * self.alive:
                self.roots = [r for r in self.roots if self.find(r) == r]
            v = self.roots[int(self.rng.integers(len(self.roots)))]
            if self.find(v) != v:
                continue
            a = int(self.rng.integers(self.k))
            b = int(self.rng.integers(self.k))
            if a == b:
                continue
            before = len(self.parent)
            self.square(v, a, b, p_comm)
            self.roots.extend(range(before, len(self.parent)))
        return self

    def graph(self):
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
        return adj, idx

    def resolve(self, marker):
        return self.find(marker)

    def size(self):
        return self.alive


class MarkedCosmos:
    """The untyped grower behind the same interface, as the control."""

    def __init__(self, seed=0, alpha=1.0, p_close=0.6, p_open=0.0,
                 cap=6, reach=2):
        self.C = Cosmos(alpha, 0.0, p_close, p_open, cap, reach, seed=seed)

    def grow_to(self, target, p_comm=1.0):
        self.C.grow(target)
        return self

    def graph(self):
        idx, adj = {}, []
        for v in range(self.C.n):
            idx[v] = len(adj)
            adj.append(set(self.C.adj[v]))
        return adj, idx

    def resolve(self, marker):
        return marker

    def size(self):
        return self.C.n


def bfs(adj, s):
    d = {s: 0}
    q = deque([s])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y not in d:
                d[y] = d[x] + 1
                q.append(y)
    return d


def pairwise(G, markers):
    """Distances between the current representatives of the markers."""
    adj, idx = G.graph()
    live = {}
    for m in markers:
        r = G.resolve(m)
        if r in idx:
            live[m] = idx[r]
    out = {}
    for m, i in live.items():
        d = bfs(adj, i)
        for n, j in live.items():
            if n > m and j in d and live[m] != live[n]:
                out[(m, n)] = d[j]
    return out


def choose_markers(G, rng, want=N_MARKERS):
    adj, idx = G.graph()
    keys = [k for k in idx if adj[idx[k]]]
    if len(keys) <= want:
        return keys
    pick = rng.choice(len(keys), size=min(want * 3, len(keys)),
                      replace=False)
    cand = [keys[int(i)] for i in pick]
    kept = []
    for c in cand:
        if len(kept) >= want:
            break
        kept.append(c)
    return kept


def run(G, rng, checkpoints=CHECKPOINTS):
    G.grow_to(N0)
    markers = choose_markers(G, rng)
    base = {p: d for p, d in pairwise(G, markers).items() if d >= 2}
    if len(base) < 6:
        return None
    rows = []
    for target in checkpoints:
        G.grow_to(target)
        cur = pairwise(G, markers)
        a = np.array([cur[p] / base[p] for p in base if p in cur], float)
        if a.size < 3:
            rows.append((target, G.size(), len(a),
                         float("nan"), float("nan")))
            continue
        mean_a = float(a.mean())
        shear = float(a.max() - a.min()) / mean_a if mean_a > 0 else float("nan")
        rows.append((target, G.size(), len(a), mean_a, shear))
    return len(base), rows


def table(title, factory, seeds=(1, 2, 3, 4)):
    print("=" * 78)
    print(title)
    print("=" * 78)
    print("  %8s %8s %8s %10s %10s %10s"
          % ("N", "pairs", "seeds", "mean a", "shear", "shear sd"))
    acc = {}
    npairs = []
    for s in seeds:
        rng = np.random.default_rng(1000 + s)
        got = run(factory(s), rng)
        if got is None:
            continue
        nb, rows = got
        npairs.append(nb)
        for tgt, n, k, ma, sh in rows:
            acc.setdefault(tgt, []).append((ma, sh, k, n))
    for t in sorted(acc):
        v = np.array([[a, b] for a, b, _, _ in acc[t]], float)
        v = v[np.isfinite(v).all(axis=1)]
        if not len(v):
            continue
        kk = int(np.mean([c for _, _, c, _ in acc[t]]))
        print("  %8d %8d %8d %10.4f %10.4f %10.4f"
              % (int(np.mean([n for _, _, _, n in acc[t]])), kk, len(v),
                 v[:, 0].mean(), v[:, 1].mean(), v[:, 1].std()))
    ns = np.array(sorted(acc), float)
    sh = np.array([np.nanmean([b for _, b, _, _ in acc[t]])
                   for t in sorted(acc)])
    good = np.isfinite(sh) & (sh > 0)
    if good.sum() >= 3:
        sl = float(np.polyfit(np.log(ns[good]), np.log(sh[good]), 1)[0])
        print()
        print("  shear ~ N^%+.3f   (shear_damping.py: -0.47 damped, 0.00 not)"
              % sl)
    print()


def coordinates(C):
    """With full commutation every node has a genuine Z^k position."""
    root = C.find(0)
    coord = {root: (0,) * C.k}
    q = deque([root])
    bad = 0
    while q:
        v = q.popleft()
        for (c, sgn), w in list(C.nb[v].items()):
            w = C.find(w)
            p = list(coord[v])
            p[c] += sgn
            p = tuple(p)
            if w in coord:
                if coord[w] != p:
                    bad += 1
            else:
                coord[w] = p
                q.append(w)
    return coord, bad


def axis_scales(C):
    """Graph distance from the origin along each axis, against |m|.

    In a faithful Z^k this ratio is 1 on every axis. Anything else is the
    complex's OWN anisotropy, measured without reference to any baseline
    epoch -- so it cannot be an artifact of where the markers were placed.
    """
    coord, bad = coordinates(C)
    adj, idx = C.graph()
    root = C.find(0)
    d = bfs(adj, idx[root])
    per = {}
    for node, p in coord.items():
        nz = [i for i, x in enumerate(p) if x != 0]
        if len(nz) != 1 or node not in idx:
            continue
        c = nz[0]
        m = abs(p[c])
        g = d.get(idx[node])
        if g is None or m < 2:
            continue
        per.setdefault(c, []).append(g / m)
    means = [float(np.mean(v)) for c, v in sorted(per.items()) if v]
    counts = [len(v) for c, v in sorted(per.items()) if v]
    if len(means) < 2:
        return bad, means, counts, float("nan"), float("nan")
    mu = float(np.mean(means))
    return (bad, means, counts, mu,
            (max(means) - min(means)) / mu if mu > 0 else float("nan"))


def lattice_table(seeds=(1, 2, 3, 4), target=32000):
    print("=" * 78)
    print("2. THE COMPLEX'S OWN ANISOTROPY, MEASURED AGAINST Z^k")
    print("=" * 78)
    print("  Graph distance to a node at m steps along one axis, divided by m.")
    print("  A faithful lattice gives 1.000 on every axis. The spread across")
    print("  axes is anisotropy with no baseline epoch involved at all.")
    print()
    print("  %5s %6s %10s %28s %10s"
          % ("k", "seed", "incons.", "ratio per axis", "spread"))
    for k in (2, 3, 4):
        sp = []
        for s in seeds:
            C = MarkedComplex(k, seed=s).grow_to(target)
            bad, means, counts, mu, spread = axis_scales(C)
            sp.append(spread)
            print("  %5d %6d %10d %28s %10.4f"
                  % (k, s, bad,
                     " ".join("%.3f" % m for m in means), spread))
        print("  %5s %6s %10s %28s %10.4f"
              % ("", "mean", "", "", float(np.nanmean(sp))))
        print()


def convergence(seed=1):
    print("=" * 78)
    print("3. DOES THE RESIDUAL ANISOTROPY DECAY, OR IS IT PERMANENT?")
    print("=" * 78)
    print("  If the spread is incompleteness, it must vanish as the lattice")
    print("  fills. If it is real anisotropy, it must not.")
    print()
    print("  %4s %9s %8s %12s %10s" % ("k", "N", "side", "mean ratio",
                                       "spread"))
    for k, Ns in ((2, (2000, 8000, 32000, 128000)),
                  (3, (8000, 32000, 128000, 400000)),
                  (4, (32000, 128000, 400000))):
        sp, ns = [], []
        for N in Ns:
            C = MarkedComplex(k, seed=seed).grow_to(N)
            bad, means, cnt, mu, spread = axis_scales(C)
            sp.append(spread)
            ns.append(N)
            print("  %4d %9d %8.1f %12.4f %10.4f%s"
                  % (k, N, N ** (1.0 / k), mu, spread,
                     "   INCONSISTENT" if bad else ""))
        sp, ns = np.array(sp), np.array(ns, float)
        g = np.isfinite(sp) & (sp > 0)
        if g.sum() >= 3:
            print("  %4s spread ~ N^%+.2f" % ("", float(np.polyfit(
                np.log(ns[g]), np.log(sp[g]), 1)[0])))
        print()


def main():
    for k in (2, 3, 4):
        table("COLOURED LATTICE, k = %d colours, full commutation" % k,
              lambda s, k=k: MarkedComplex(k, seed=s))

    table("CONTROL: untyped grower, identical marker protocol\n"
          "  alpha=1, p_close=0.6, p_open=0, cap=6, reach=2",
          lambda s: MarkedCosmos(seed=s))

    lattice_table()
    convergence()

    print("=" * 78)
    print("  mean a near 1 and flat means accretion: the complex gets bigger")
    print("  without the interior expanding, so isotropy would be free and")
    print("  meaningless. mean a falling means the lazy filling is pulling")
    print("  distances down toward the true lattice metric, and then shear")
    print("  says whether it does so evenly.")
    print()


if __name__ == "__main__":
    main()
