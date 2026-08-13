#!/usr/bin/env python3
"""
six_criteria.py

The six conditions a growth rule would have to meet, made operational, and
every rule built so far scored against them.

    1  conservation     births cost something, exactly
    2  locality         a node is created near the structure paying for it
    3  isotropy         no direction permanently privileged
    4  scale consistency  d_exp agrees with d_ball
    5  homogenization   initial differences do not freeze into permanent shear
    6  derived rate     G comes from the state, not the schedule

Criterion 2 needs restating before it can be scored, and that turns out to be
the whole story. Read as "new adjacencies join nearby nodes", it is satisfied
by rules that fail anyway: `close` joins two neighbours of a common vertex and
`shortcut` joins z to a neighbour of x, so both act at graph distance exactly
2, and both produce small-world graphs. Distance-locality is not the
discriminator.

The discriminator is BOUNDED LOCAL CAPACITY. A d-dimensional space has bounded
volume in a ball of fixed radius; these graphs did not, running to maximum
degree 605 while the mean stayed pinned at 4. Locality has to mean a bound on
how much structure can pile up at one node, not a bound on how far an edge
reaches.

That restatement is testable, and testing it is the main result here.

Result: capping the degree buys a real dimension, and it is tunable
---------------------------------------------------------------------
Adding a single condition to `close` -- refuse if either endpoint is already
at the cap -- changes the phase of the model. Over a 64-fold range in N, with
mean distance measured from 8 sources at each of 4 checkpoints:

    cap    mean distance ~ N^p    power resid   log resid    verdict
     6            N^0.248            0.00149     0.01283     power, d_H = 4.03
    12            N^0.158            0.00169     0.00711     power, d_H = 6.34
    none          N^0.061            0.00093     0.00135     indistinguishable

At N = 156271 the three sit at mean distance 27.66, 13.82 and 8.35, against
10.23 for a random 4-regular graph and 197.50 for a periodic square lattice.
The uncapped model is a random graph. The capped one is not: it is displaced
from random by a factor of 2.7 and follows a clean power law.

So bounded degree is the missing ingredient, and the cap acts as a dial on the
emergent Hausdorff dimension -- tighter cap, lower dimension. Cap 4 freezes
the universe at N = 42, the same starvation that killed the strictly-local
rule, so there is a floor.

Two honest qualifications. d_H = 4.03 is suspiciously exactly the
branched-polymer value, which is also what 2D dynamical triangulations give,
so the likely reading is that this lands in the branched-polymer universality
class rather than on a manifold -- a finite dimension is not the same as a
smooth space. And the dimension is set by the cap, which is a parameter, so it
is tuned rather than derived. What is genuinely new is that a finite dimension
exists at all, which no previous rule in this project achieved.

Scorecard
---------
                          free      free +      paid       paid
                          subdiv    chord       uncapped   cap 6
    1 conservation        fail      fail        PASS       PASS
    2 locality            PASS      fail        fail       PASS
    3 isotropy            fail      fail        fail       fail
    4 scale consistency   PASS      unsettled   unsettled  unsettled
    5 homogenization      fail      fail        fail       fail
    6 derived rate        n/a       fail        fail       fail

Best is the capped paid model with two clean passes, one unsettled, three
failures. Free subdivide's pass on 4 is real but hollow -- d_exp, d_ball and
d_H all equal 1.00 because it makes threads, and consistency is easy when
there is only one dimension available.

Criterion 4 is scored "unsettled" rather than failed wherever d_ball has not
stopped moving. For the capped model d_exp = 4.53 and d_H = 4.03 agree
reasonably, but d_ball is still climbing, 2.62 -> 3.06 over the range, so it
may be heading for the same value and simply lagging at these sizes. Calling
that a failure would overstate what was measured; calling it a pass would
overstate it in the other direction.

Criteria 3 and 5 fail for EVERY rule, and that is the most informative line in
the table. Anisotropy sits at 40-51% and inhomogeneity at 8-15%, and neither
decays -- across four rules, four checkpoints and a 64-fold range in N, the
trend is flat every time. Nothing here homogenizes.

One caveat on the absolute levels: with only six markers at base distances 1
to 3, some spread is expected simply because those separations differ, so 45%
should not be read as 45% of physical shear. The trend is what carries the
result, and a genuinely homogenizing dynamics would drive it down whatever its
starting level. None does.

Criterion 6 fails everywhere it applies. The split of the budget between the
two sinks tracks the propensities, so G is inherited from the schedule, not
from the state.

What this says about the programme
-----------------------------------
The six conditions are not independent, and the audit shows they trade against
each other rather than accumulating. Locality is buyable -- cap the degree.
Conservation is buyable -- meter the births. Buying both at once works, and
buys a finite dimension as well, which is genuinely new. But isotropy and
homogenization resist every rule tried, and the reason is the one measured in
paid_expansion.py: a homogeneous configuration is a vanishing fraction of the
configurations available, so a dynamics that samples by availability will not
find one and will not stay in one if placed there. Capping the degree
restricts WHICH typical graph gets sampled; it does not make the dynamics stop
sampling typical graphs.

So the honest position is that conditions 1, 2 and a finite 4 are achievable
together, and 3, 5, 6 are not blocked by any detail of these particular rules
but by the entropic argument. Getting them would need what CDT and quantum
graphity both need: something that penalises non-geometric configurations, so
the dynamics is no longer purely entropic. That is a different kind of model
-- one with an energy, not just a budget.
"""

import numpy as np
from collections import deque

from paid_expansion import Universe, random_regular, square_lattice


# ------------------------------------------------------------------ rules

class FreeSubdivide(Universe):
    """scale_factor.py's rule: births are free, nothing is conserved."""

    def subdivide(self, rng):
        a, b = self.el[int(rng.integers(len(self.el)))]
        self.adj.append(set())
        z = self.n
        self.n += 1
        self._del(a, b)
        self._add(a, z)
        self._add(z, b)
        self.cnt["subdivide"] += 1
        return True

    def close(self, rng):
        return False

    def open(self, rng):
        return False


class FreeShortcut(FreeSubdivide):
    """Free subdivision plus a distance-2 chord: local by distance, not by
    capacity."""

    def close(self, rng):
        for _ in range(24):
            a, b = self.el[int(rng.integers(len(self.el)))]
            v = a if rng.random() < 0.5 else b
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            x = nb[int(rng.integers(len(nb)))]
            y = nb[int(rng.integers(len(nb)))]
            if x == y or y in self.adj[x]:
                continue
            self._add(x, y)
            self.cnt["close"] += 1
            return True
        return False


class Capped(Universe):
    """The paid model with bounded local capacity."""

    cap = 10**9

    def close(self, rng):
        for _ in range(30):
            a, b = self.el[int(rng.integers(len(self.el)))]
            v = a if rng.random() < 0.5 else b
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            x = nb[int(rng.integers(len(nb)))]
            y = nb[int(rng.integers(len(nb)))]
            if x == y or y in self.adj[x]:
                continue
            if len(self.adj[x]) >= self.cap or len(self.adj[y]) >= self.cap:
                continue
            self._add(x, y)
            self.bank += 1
            self.cnt["close"] += 1
            return True
        return False


def make(kind, cap=None):
    if kind == "free-subdivide":
        return FreeSubdivide()
    if kind == "free-shortcut":
        return FreeShortcut()
    U = Capped()
    U.cap = cap if cap else 10**9
    return U


# ------------------------------------------------------------ measurement

def probe(U, rng, sources=8):
    """Mean distance, diameter and the ball-growth dimension in one sweep."""
    tot = defaultdict = {}
    acc, cnt, rmax = 0.0, 0, 0
    cum_by_r = {}
    for _ in range(sources):
        s = int(rng.integers(U.n))
        d = U.bfs(s)
        v = np.fromiter(d.values(), int)
        acc += float(v.sum())
        cnt += len(v)
        rmax = max(rmax, int(v.max()))
        cum = np.cumsum(np.bincount(v))
        for r in range(1, len(cum)):
            cum_by_r[r] = cum_by_r.get(r, 0.0) + float(cum[r])
    hi = rmax // 3
    rs = [r for r in sorted(cum_by_r) if 2 <= r <= hi]
    if len(rs) >= 4:
        x = np.log(np.array(rs, float))
        y = np.log(np.array([cum_by_r[r] / sources for r in rs]))
        d_ball = float(np.polyfit(x, y, 1)[0])
    else:
        d_ball = float("nan")
    return acc / cnt, rmax, d_ball


def shear(U, base):
    """Per-marker anisotropy and between-marker inhomogeneity of a_ij."""
    rows = []
    for i in range(6):
        d = U.bfs(i)
        rows.append([d[j] / base[(min(i, j), max(i, j))]
                     for j in range(6) if j != i])
    rows = np.array(rows, float)
    return (float(rows.mean()),
            float(np.mean(rows.std(axis=1) / rows.mean(axis=1))),
            float(rows.mean(axis=1).std() / rows.mean()))


def run(kind, cap, checkpoints, seed=11, weights=None):
    rng = np.random.default_rng(seed)
    U = make(kind, cap)
    ops = ["close", "subdivide", "open"]
    w = None
    if weights:
        w = np.array(weights, float)
        w = w / w.sum()
    base = {}
    d0 = {v: U.bfs(v) for v in range(6)}
    for i in range(6):
        for j in range(i + 1, 6):
            base[(i, j)] = d0[i][j]
    out = []
    for t in range(1, max(checkpoints) + 1):
        for _ in range(6):
            i = int(rng.integers(3)) if w is None else int(rng.choice(3, p=w))
            if getattr(U, ops[i])(rng):
                break
        if t in checkpoints:
            md, dm, db = probe(U, np.random.default_rng(2))
            a, an, ih = shear(U, base)
            out.append(dict(ev=t, N=U.n, meand=md, diam=dm, d_ball=db,
                            a=a, aniso=an, inhom=ih,
                            maxdeg=max(len(s) for s in U.adj),
                            k=2 * len(U.el) / U.n,
                            res=U.residual(), cnt=dict(U.cnt)))
    return U, out


def fit(ns, vals):
    ns, vals = np.asarray(ns, float), np.asarray(vals, float)
    p, b = np.polyfit(np.log(ns), np.log(vals), 1)
    rp = float(np.sum((np.log(vals) - (p * np.log(ns) + b)) ** 2))
    c, b2 = np.polyfit(np.log(ns), vals, 1)
    rl = float(np.sum(((vals - (c * np.log(ns) + b2)) / vals.mean()) ** 2))
    return float(p), rp, rl


CK = [10000, 40000, 160000, 640000]
RULES = [("free-subdivide", None, "free subdivide"),
         ("free-shortcut", None, "free subdivide + chord"),
         ("paid", None, "paid, degree unbounded"),
         ("paid", 6, "paid, degree capped at 6")]


def main():
    print("=" * 78)
    print("A.  THE CAP SCAN: does bounded local capacity buy a dimension?")
    print("=" * 78)
    store = {}
    for kind, cap, label in RULES:
        U, rows = run(kind, cap, CK)
        store[label] = (U, rows)
        print("  %s" % label)
        print("     %8s %8s %8s %6s %8s %8s %7s %6s"
              % ("events", "N", "meand", "diam", "d_ball", "a", "maxdeg", "<k>"))
        for r in rows:
            print("     %8d %8d %8.2f %6d %8.2f %8.2f %7d %6.2f"
                  % (r["ev"], r["N"], r["meand"], r["diam"], r["d_ball"],
                     r["a"], r["maxdeg"], r["k"]))
        ns = [r["N"] for r in rows]
        if max(ns) - min(ns) < 2:
            print("     -> FROZEN at N = %d" % ns[-1])
        else:
            p, rp, rl = fit(ns, [r["meand"] for r in rows])
            print("     -> meand ~ N^%.3f  (pow %.5f, log %.5f): %s"
                  % (p, rp, rl,
                     "POWER LAW, d_H = %.2f" % (1 / p) if rp < rl
                     else "logarithmic, no finite dimension"))
        print()

    U, rows = store["paid, degree capped at 6"]
    n = rows[-1]["N"]
    rng = np.random.default_rng(4)
    radj = random_regular(n, 4, rng)
    ladj, ln = square_lattice(n)

    def md_of(adj_of, size):
        acc, cnt = 0.0, 0
        rr = np.random.default_rng(2)
        for _ in range(8):
            s = int(rr.integers(size))
            d = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for w_ in adj_of(u):
                    if w_ not in d:
                        d[w_] = d[u] + 1
                        q.append(w_)
            v = np.fromiter(d.values(), int)
            acc += float(v.sum())
            cnt += len(v)
        return acc / cnt

    print("  at N = %d, mean degree 4:" % n)
    print("     %-28s %10s" % ("graph", "mean dist"))
    for lbl in ("paid, degree unbounded", "paid, degree capped at 6"):
        print("     %-28s %10.2f" % (lbl, store[lbl][1][-1]["meand"]))
    print("     %-28s %10.2f" % ("random 4-regular graph",
                                 md_of(lambda u: radj[u], n)))
    print("     %-28s %10.2f" % ("periodic square lattice",
                                 md_of(lambda u: ladj[u], ln)))
    print()

    print("=" * 78)
    print("B.  THE SIX CRITERIA, SCORED")
    print("=" * 78)
    for kind, cap, label in RULES:
        U, rows = store[label]
        ns = [r["N"] for r in rows]
        frozen = max(ns) - min(ns) < 2
        p_md, rp, rl = fit(ns, [r["meand"] for r in rows])
        p_a, _, _ = fit(ns, [max(r["a"], 1e-9) for r in rows])
        d_exp = 1 / p_a if p_a > 1e-6 else float("inf")
        dbs = [r["d_ball"] for r in rows if not np.isnan(r["d_ball"])]
        an = [r["aniso"] for r in rows]
        ih = [r["inhom"] for r in rows]

        # 6: does births/close move when the propensities move?
        rates = []
        for wts in ((1, 1, 1), (1, 2, 1), (1, 1, 2)):
            Uw, rw = run(kind, cap, [40000], weights=wts)
            c = max(rw[-1]["cnt"].get("close", 0), 1)
            rates.append(rw[-1]["cnt"].get("subdivide", 0) / c)
        spread6 = (max(rates) - min(rates)) / max(np.mean(rates), 1e-9)

        print("  %s" % label)
        print("     1 conservation     %s   (residual %d)"
              % ("PASS" if r_ok(U) else "fail", rows[-1]["res"]))
        print("     2 locality         %s   (max degree %d -> %d)"
              % ("PASS" if rows[-1]["maxdeg"] <= 2 * rows[0]["maxdeg"]
                 else "fail", rows[0]["maxdeg"], rows[-1]["maxdeg"]))
        print("     3 isotropy         %s   (anisotropy %.0f%% -> %.0f%%)"
              % ("PASS" if an[-1] < 0.5 * an[0] else "fail",
                 100 * an[0], 100 * an[-1]))
        # d_ball only counts once it has stopped moving; a drifting value is
        # the fitting window growing, not a dimension.
        settled = bool(dbs) and (max(dbs) - min(dbs)) < 0.15 * float(np.mean(dbs))
        ok4 = bool(dbs) and settled and rp < rl \
            and abs(d_exp - dbs[-1]) < 0.25 * dbs[-1]
        print("     4 scale consistency %s  (d_exp %.2f, d_ball %.2f%s, d_H %s)"
              % ("PASS" if ok4 else ("fail" if settled or not dbs
                                     else "unsettled"),
                 d_exp, dbs[-1] if dbs else float("nan"),
                 "" if settled else " still climbing from %.2f" % dbs[0],
                 "%.2f" % (1 / p_md) if rp < rl else "none"))
        print("     5 homogenization   %s   (inhomogeneity %.0f%% -> %.0f%%)"
              % ("PASS" if ih[-1] < 0.5 * ih[0] else "fail",
                 100 * ih[0], 100 * ih[-1]))
        if rows[-1]["cnt"].get("close", 0) == 0:
            print("     6 derived rate     n/a    (no paying move: births are "
                  "free, so there is no rate to derive)")
        else:
            print("     6 derived rate     %s   (births/close spans %.0f%% "
                  "under reweighting)"
                  % ("PASS" if spread6 < 0.1 else "fail", 100 * spread6))
        if frozen:
            print("     (frozen: no growth at all)")
        print()


def r_ok(U):
    return U.residual() == 0 and isinstance(U, Capped)


if __name__ == "__main__":
    main()
