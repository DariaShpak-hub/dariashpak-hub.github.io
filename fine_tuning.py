#!/usr/bin/env python3
"""
fine_tuning.py

The fine-tuning question, asked where we own the measure.

The standard argument -- "vary a constant and the universe stops working,
therefore the constants are tuned" -- has a hole that is usually waved past:
you cannot call a value improbable without a probability distribution over the
alternatives, and physics has none. There is no known range for Lambda, no
prior on alpha, nothing. The "122 digits of tuning" number quietly assumes a
uniform measure on a parameter whose allowed range nobody knows.

Here that hole is closed by construction. We wrote the rules, so we know
exactly what the parameters are and what values they can take. That makes this
the one version of the question that is actually answerable, and it lets us
test the two methodological complaints directly:

    1  ONE-AT-A-TIME IS THE WRONG METHOD. Varying one constant while holding
       the others fixed measures the width of a line through the habitable
       region, not its volume. Adams' joint scan over gravity, fine structure
       and nuclear rates found roughly a quarter of the explored space still
       makes stars, against the near-zero the one-at-a-time literature
       reports. Does the same gap appear here?

    2  THE ANSWER DEPENDS ON THE PRIOR EVEN WHEN YOU OWN THE MODEL. Reweight
       the same samples under a log-uniform prior instead of a uniform one and
       the habitable fraction should move. If it moves a lot, then "the
       universe is tuned to one part in X" is not a fact about the universe,
       it is a fact about the measure you chose -- and that holds even here,
       where we invented the parameters ourselves.

The parameter space
-------------------
Five knobs on the grower this project has been building all along.

    alpha     anti-Yule exponent. Subdivision picks an edge with weight
              2^(-alpha*g), g = the edge's generation. alpha = 0 is uniform
              selection, which self_amplification.py showed runs away into
              anisotropy; alpha = 1 is what shear_damping.py used to get shear
              from 67% down to 0.7% at the central-limit rate. Range [0, 3].

    p_close   relative rate of the close move -- a chord that adds an edge and
              funds node creation. model.py's central result is that this is
              also the move that destroys the metric. Range [0, 1].

    p_open    relative rate of the open move -- deleting an edge that sits in
              a triangle, so connectivity survives. Range [0, 1].

    cap       maximum node degree. six_criteria.py showed this sets d_H
              directly (cap 6 -> 4.03, cap 12 -> 6.34), so it is the knob most
              likely to look "tuned". Integers 3..12.

    reach     how far a close can reach: the chord joins the two ends of a
              random walk of this many steps. reach = 2 is the common-neighbour
              chord used everywhere earlier. Integers 2..5.

A sixth structural choice -- what generation a close-created edge is born with
-- is fixed to the local mean of its neighbourhood, so a new chord inherits the
scale of where it appears rather than resetting the clock. It is a choice, and
it is flagged as one.

What counts as habitable
------------------------
Five criteria, fixed before the scan was run and not adjusted afterwards. Each
is something this project already measured on real runs, so none of them was
invented to make a number come out.

    expands     mean scale factor a >= 2. The relational metric actually grew,
                not just the node count. (scale_factor.py)
    extended    shell counts fit a power law better than an exponential, and
                mean distance >= 5. Rules out both tree-like blowup and
                small-world collapse. (coloured_memory.py, global_observables.py)
    dimension   1.8 <= d <= 4.5 on the shell estimator, which reads about 0.15
                low at these sizes. Rules out threads and expanders.
    isotropic   shear (max a - min a)/mean a <= 0.30. The real bound is
                sigma/H < 4.7e-11, unreachable at N = 2000; 0.30 is the
                threshold that separates shear_damping.py's two rules, which
                sat at 67% and 0.7%.
    stable      max degree <= 3 * mean degree. No hub blowup.

All five must hold. N = 2000, single history per point unless stated.

Result 1: at the criteria as written, NOTHING in the box works
---------------------------------------------------------------
600 points drawn uniformly from the six-knob box, two histories each. The
habitable fraction as the isotropy threshold tau is relaxed:

    tau        0.05  0.10  0.20  0.30  0.50  0.80  1.00  1.20  1.50  2.00  inf
    seed A     .000  .000  .000  .000  .000  .002  .010  .035  .117  .213  .232
    seed B     .000  .000  .000  .000  .000  .000  .002  .022  .110  .230  .250

At the pre-registered tau = 0.30 the answer is 0 out of 600, twice over. The
model is not merely tuned, it is dead everywhere -- and it stays dead until the
isotropy requirement has been loosened by a factor of five.

Two pilot runs say this is structural rather than a small-N artifact. Pure
subdivision at alpha = 1 reproduces shear_damping.py exactly -- shear 0.144,
0.132, 0.188, 0.093, 0.065 at N = 500 to 8000 -- but sits at d = 1.00 the whole
way. Switch the closes on and shear locks near 1.4 and stops moving: 1.33,
1.56, 1.39, 1.35, 1.35, 1.53 at N = 500 to 16000. The isotropy that
subdivision buys is destroyed by the move that creates dimension, and the
damage does not decay with size.

That is model.py's tradeoff again, and this scan is its sharpest statement so
far: it is not a matter of finding the right parameters, because there are no
right parameters. Even with the isotropy criterion deleted entirely (tau = inf)
only 23% of the box survives the other four.

Everything below runs at tau = 1.50, the tightest grid value leaving enough of
the box alive to measure. That is a deliberately weakened criterion and every
number after this point inherits the weakening.

Result 2: the criterion chosen decides the answer
---------------------------------------------------
Failure counts (a point can fail several criteria at once):

                  tau = 0.30            tau = 1.50
    isotropic     589   98.2%           245   46.2%
    dimension     252   42.0%           252   47.5%
    stable        214   35.7%           214   40.4%
    extended      177   29.5%           177   33.4%
    expands         0    0.0%             0    0.0%

Two things worth admitting. First, "expands" never fires -- I wrote a criterion
that does no work at all, which is exactly the failure mode I flagged in
advance and then walked into. Second, moving one threshold turns isotropy from
the sole cause of death into one of four roughly equal causes. The headline
number is a property of the criteria, not of the model.

Result 3: which knobs constrain -- and it is not the ones expected
-------------------------------------------------------------------
Habitability rate by knob value, read straight off the joint samples:

    knob        low third   middle   high third
    cap            0.215     0.090     0.014
    p_close        0.040     0.107     0.193
    reach          0.109     0.156     0.082
    p_open         0.135     0.110     0.095
    beta           0.130     0.092     0.117
    alpha          0.090     0.125     0.125

Only two knobs really bite: the degree cap, strongly and monotonically
downward, and the close rate, strongly upward. The two anti-Yule exponents --
alpha, which shear_damping.py showed was the whole isotropy mechanism, and
beta, its analogue for chord placement -- are almost flat. Once closes are in
play, the mechanism that damps shear stops mattering.

Result 4: one-at-a-time understates the volume by 40 to 70 times
------------------------------------------------------------------
Three habitable anchors, each axis swept with the other five held fixed, each
grid point scored as the fraction of four seeds that pass:

    anchor                                   naive product   joint/naive
    a=1.74 b=0.69 pc=0.84 po=0.37 cap=5 r=3      0.00170          68.5
    a=1.55 b=1.20 pc=0.97 po=0.41 cap=3 r=4      0.00171          68.4
    a=1.06 b=2.36 pc=0.73 po=0.78 cap=4 r=3      0.00293          39.9

    measured joint fraction                      0.11667

If the habitable set were a box, the product of the axial fractions would equal
the joint fraction exactly. It is off by nearly two orders of magnitude, and in
the direction that says the set is a curved ridge rather than a box: the
constraints trade against each other, so a step that ruins the universe along
one axis can be paid for by a step along another. One-at-a-time variation walks
straight off the ridge and reports a cliff, while the ridge itself is long.

This is the Adams result reproduced in miniature, and it is the single
strongest technical objection to the fine-tuning argument as usually stated.

Result 5: the prior moves the answer by 3.6x
----------------------------------------------
Same 600 samples, reweighted:

    uniform                                        0.1167
    log-uniform in alpha, beta, p_close, p_open    0.0431
    1/sqrt in alpha only                           0.1098
    linear tilt toward large p_close               0.1549

A factor of 3.6 between the extremes, from nothing but a change of measure on
parameters we invented ourselves. No evidence can settle which of these is
correct, because there is nothing for the prior to be correct about.

Result 6: the box matters, but only through knobs that constrain
------------------------------------------------------------------
    alpha<=3,  cap<=12, reach<=5     0.1167   (600 samples)
    alpha<=9,  cap<=12, reach<=5     0.1100   (200 samples)
    alpha<=3,  cap<=30, reach<=5     0.0450   (200 samples)
    alpha<=3,  cap<=12, reach<=10    0.0900   (200 samples)

I expected widening alpha to shrink the fraction and it did not -- 0.117 to
0.110, no effect, because alpha is not a constraining knob. Widening the degree
cap from 12 to 30 cuts the fraction by 2.6x, because cap is. So the objection
is real but narrower than the usual statement of it: an unbounded parameter
only dilutes the tuning number if habitability actually depends on it. Nothing
in the model fixes cap <= 12, and nothing ever will.

Verdict
-------
Two separate results, and they point opposite ways.

The physics result is negative and hard. Under criteria fixed before the scan,
zero of 600 parameter settings produce a viable universe, and the reason is
always the same one: the move that creates dimension destroys isotropy, at a
level that does not decay from N = 500 to N = 16000. No knob in this model
fixes that, including the one built specifically to fix it.

The philosophical result is that the fine-tuning argument's machinery does not
survive a case where the answer is checkable. Three numbers that ought to be
bookkeeping turn out to carry the conclusion: the criterion (a factor of
infinity here -- 0 becomes 0.117 on relaxing one threshold), the method
(one-at-a-time understates the volume 40-70x), and the prior (3.6x). None of
this shows the real universe is untuned. It shows that "tuned to one part in X"
is not a measurement unless all three are pinned down, and in physics none of
them are.

Caveats. N = 2000, one seed topology (a 6-cycle), the shell estimator reads
about 0.15 low, and at tau = 1.50 individual verdicts flip between seeds 15% of
the time -- the fractions are stable across seeds (0.1167 vs 0.1100) but single
points near the boundary are not.
"""

import numpy as np
from collections import deque

ALPHA_RANGE = (0.0, 3.0)
PCLOSE_RANGE = (0.0, 1.0)
POPEN_RANGE = (0.0, 1.0)
CAP_RANGE = (3, 12)
REACH_RANGE = (2, 5)

REFERENCE = dict(alpha=1.0, p_close=0.15, p_open=0.15, cap=6, reach=2)

CRITERIA = ("expands", "extended", "dimension", "isotropic", "stable")


class Cosmos:
    """The project's grower with every rule exposed as a parameter."""

    def __init__(self, alpha, beta, p_close, p_open, cap, reach,
                 n0=6, seed=0):
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.p_close = float(p_close)
        self.p_open = float(p_open)
        self.cap = int(cap)
        self.reach = int(reach)
        self.rng = np.random.default_rng(seed)

        self.adj = [set() for _ in range(n0)]
        self.einc = [set() for _ in range(n0)]   # node -> live edge indices
        self.eu, self.ev, self.eg = [], [], []
        self.live = []
        self.bucket = {}            # generation -> list of edge indices
        self.n_live = 0
        self.pos = []               # edge index -> position in its bucket
        self.n = n0
        for i in range(n0):
            self._add(i, (i + 1) % n0, 0)
        self.markers = list(range(n0))
        self.base = self._baseline()

    # ------------------------------------------------------------- topology
    def _add(self, u, v, g):
        i = len(self.eu)
        self.eu.append(u)
        self.ev.append(v)
        self.eg.append(g)
        self.live.append(True)
        b = self.bucket.setdefault(g, [])
        self.pos.append(len(b))
        b.append(i)
        self.n_live += 1
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.einc[u].add(i)
        self.einc[v].add(i)
        return i

    def _kill(self, i):
        g = self.eg[i]
        b = self.bucket[g]
        p = self.pos[i]
        last = b.pop()
        if p < len(b):
            b[p] = last
            self.pos[last] = p
        self.live[i] = False
        self.n_live -= 1
        self.adj[self.eu[i]].discard(self.ev[i])
        self.adj[self.ev[i]].discard(self.eu[i])
        self.einc[self.eu[i]].discard(i)
        self.einc[self.ev[i]].discard(i)

    def _random_live_edge(self):
        if self.n_live == 0:
            return -1
        r = int(self.rng.integers(self.n_live))
        for b in self.bucket.values():
            if r < len(b):
                return b[r]
            r -= len(b)
        return -1

    def _weighted_edge(self, expo):
        """Pick an edge with weight 2^(-expo*g), exactly and in O(#gens)."""
        gens = [g for g, b in self.bucket.items() if b]
        if not gens:
            return -1
        ga = np.array(gens, float)
        w = np.exp2(-expo * ga) * np.array(
            [len(self.bucket[g]) for g in gens], float)
        tot = w.sum()
        if not np.isfinite(tot) or tot <= 0:
            return self._random_live_edge()
        g = gens[int(np.searchsorted(np.cumsum(w), self.rng.random() * tot))]
        b = self.bucket[g]
        return b[int(self.rng.integers(len(b)))]

    # ---------------------------------------------------------------- moves
    def subdivide(self):
        i = self._weighted_edge(self.alpha)
        if i < 0:
            return False
        u, v, g = self.eu[i], self.ev[i], self.eg[i]
        self._kill(i)
        z = self.n
        self.adj.append(set())
        self.einc.append(set())
        self.n += 1
        self._add(u, z, g + 1)
        self._add(z, v, g + 1)
        return True

    def _walk(self, x, steps):
        prev, cur = -1, x
        for _ in range(steps):
            nb = [y for y in self.adj[cur] if y != prev]
            if not nb:
                nb = list(self.adj[cur])
            if not nb:
                return -1
            prev, cur = cur, nb[int(self.rng.integers(len(nb)))]
        return cur

    def close(self):
        """Chord joining the ends of a `reach`-step walk. dE = +1."""
        for _ in range(30):
            i = self._weighted_edge(self.beta)
            if i < 0:
                return False
            x = self.eu[i] if self.rng.random() < 0.5 else self.ev[i]
            y = self._walk(x, self.reach)
            if y < 0 or y == x or y in self.adj[x]:
                continue
            if len(self.adj[x]) >= self.cap or len(self.adj[y]) >= self.cap:
                continue
            self._add(x, y, self._chord_gen(x, y))
            return True
        return False

    def _chord_gen(self, a, b):
        """A new chord inherits the generation scale of where it appears."""
        gs = []
        for x in (a, b):
            for i in self.inc(x):
                gs.append(self.eg[i])
        return int(round(float(np.mean(gs)))) if gs else 0

    def inc(self, x):
        """Live edge indices incident to x. O(deg) via the incidence index."""
        return self.einc[x]

    def open(self):
        """Delete an edge that sits in a triangle, so the graph stays whole."""
        for _ in range(30):
            i = self._random_live_edge()
            if i < 0:
                return False
            u, v = self.eu[i], self.ev[i]
            if self.adj[u] & self.adj[v]:
                self._kill(i)
                return True
        return False

    # ----------------------------------------------------------------- loop
    def grow(self, target):
        w = np.array([1.0, self.p_close, self.p_open])
        w = w / w.sum()
        budget = 40 * target
        steps = 0
        while self.n < target and steps < budget:
            steps += 1
            r = self.rng.random()
            if r < w[0]:
                self.subdivide()
            elif r < w[0] + w[1]:
                self.close()
            else:
                self.open()
        self.stalled = self.n < target
        return self

    # ---------------------------------------------------------- observables
    def bfs(self, s):
        d = {s: 0}
        q = deque([s])
        while q:
            x = q.popleft()
            for y in self.adj[x]:
                if y not in d:
                    d[y] = d[x] + 1
                    q.append(y)
        return d

    def _baseline(self):
        base = {}
        for i in self.markers:
            d = self.bfs(i)
            for j in self.markers:
                if j > i:
                    base[(i, j)] = d[j]
        return base

    def measure(self):
        n_e = sum(1 for f in self.live if f)
        deg = [len(a) for a in self.adj]
        mean_k = 2.0 * n_e / self.n
        max_k = max(deg) if deg else 0

        a = []
        for i in self.markers:
            d = self.bfs(i)
            for j in self.markers:
                if j > i and self.base.get((i, j)) and j in d:
                    a.append(d[j] / self.base[(i, j)])
        a = np.array(a, float)
        mean_a = float(a.mean()) if a.size else float("nan")
        shear = (float(a.max() - a.min()) / mean_a
                 if a.size and mean_a > 0 else float("nan"))

        acc, tot, cnt = {}, 0.0, 0
        keys = [v for v in range(self.n) if self.adj[v]]
        for _ in range(6):
            s = keys[int(self.rng.integers(len(keys)))]
            v = np.fromiter(self.bfs(s).values(), int)
            tot += float(v.sum())
            cnt += len(v)
            c = np.bincount(v)
            for r in range(len(c)):
                acc[r] = acc.get(r, 0.0) + float(c[r])
        sh = {r: x / 6.0 for r, x in acc.items()}
        mean_dist = tot / cnt if cnt else float("nan")
        hi = max(4, max(sh) // 3)
        rs = [r for r in sorted(sh) if 2 <= r <= hi and sh[r] > 0]
        if len(rs) < 3:
            d_shell, power = float("nan"), False
        else:
            xr = np.array(rs, float)
            y = np.log(np.array([sh[r] for r in rs]))
            fp = np.polyfit(np.log(xr), y, 1)
            fe = np.polyfit(xr, y, 1)
            res_p = float(np.sum((y - np.polyval(fp, np.log(xr))) ** 2))
            res_e = float(np.sum((y - np.polyval(fe, xr)) ** 2))
            d_shell = 1.0 + float(fp[0])
            power = res_p < res_e

        reach_all = len(self.bfs(0)) == self.n
        return dict(N=self.n, E=n_e, mean_k=mean_k, max_k=max_k,
                    mean_a=mean_a, shear=shear, d=d_shell,
                    mean_dist=mean_dist, power=power, whole=reach_all,
                    stalled=getattr(self, "stalled", False))



# ---------------------------------------------------------------- criteria

TAU_REGISTERED = 0.30
TAU_GRID = (0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 1e9)


def verdict(m, tau=TAU_REGISTERED):
    """The five criteria. Only the isotropy threshold tau is variable."""
    return {
        "expands": bool(m["mean_a"] >= 2.0) and not m["stalled"],
        "extended": bool(m["power"] and m["mean_dist"] >= 5.0 and m["whole"]),
        "dimension": bool(np.isfinite(m["d"]) and 1.8 <= m["d"] <= 4.5),
        "isotropic": bool(np.isfinite(m["shear"]) and m["shear"] <= tau),
        "stable": bool(m["max_k"] <= 3.0 * m["mean_k"]),
    }


def habitable(m, tau=TAU_REGISTERED):
    return all(verdict(m, tau).values())


# ------------------------------------------------------------ the box

ALPHA_RANGE = (0.0, 3.0)
BETA_RANGE = (0.0, 3.0)
PCLOSE_RANGE = (0.0, 1.0)
POPEN_RANGE = (0.0, 1.0)
CAP_RANGE = (3, 12)
REACH_RANGE = (2, 5)

KNOBS = ("alpha", "beta", "p_close", "p_open", "cap", "reach")

AXIS = dict(
    alpha=[0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0],
    beta=[0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0],
    p_close=[0.0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0],
    p_open=[0.0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0],
    cap=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    reach=[2, 3, 4, 5],
)

REFERENCE = dict(alpha=1.0, beta=0.0, p_close=0.15, p_open=0.15,
                 cap=6, reach=2)


def draw(rng, alpha_hi=ALPHA_RANGE[1]):
    return dict(alpha=float(rng.uniform(0.0, alpha_hi)),
                beta=float(rng.uniform(*BETA_RANGE)),
                p_close=float(rng.uniform(*PCLOSE_RANGE)),
                p_open=float(rng.uniform(*POPEN_RANGE)),
                cap=int(rng.integers(CAP_RANGE[0], CAP_RANGE[1] + 1)),
                reach=int(rng.integers(REACH_RANGE[0], REACH_RANGE[1] + 1)))


def run_point(p, seed, target):
    return Cosmos(p["alpha"], p["beta"], p["p_close"], p["p_open"],
                  p["cap"], p["reach"], seed=seed).grow(target).measure()


def unit(p, alpha_hi=ALPHA_RANGE[1]):
    """Box coordinates in [0,1]^6, so distances mean something."""
    return np.array([p["alpha"] / alpha_hi, p["beta"] / 3.0,
                     p["p_close"], p["p_open"],
                     (p["cap"] - 3) / 9.0, (p["reach"] - 2) / 3.0])


# ------------------------------------------------------------------ report

def main(target=2000, n_joint=600, n_wide=200, seeds=(11, 23)):
    print("=" * 78)
    print("0. THE PRE-REGISTERED CRITERIA, AT THE REFERENCE POINT")
    print("=" * 78)
    m = run_point(REFERENCE, seeds[0], target)
    ok = verdict(m)
    print("  " + "  ".join("%s=%s" % (k, REFERENCE[k]) for k in KNOBS))
    print("  N=%d  <k>=%.2f  a=%.1f  shear=%.3f  d=%.2f  dist=%.1f  %s"
          % (m["N"], m["mean_k"], m["mean_a"], m["shear"], m["d"],
             m["mean_dist"], "power" if m["power"] else "EXP"))
    print("  " + "  ".join("%s:%s" % (c, "ok" if v else "FAIL")
                           for c, v in ok.items()))
    print()

    print("=" * 78)
    print("1. JOINT SCAN: %d points drawn uniformly from the 6-knob box"
          % n_joint)
    print("=" * 78)
    rng = np.random.default_rng(5)
    pts = [draw(rng) for _ in range(n_joint)]
    ms = [run_point(p, seeds[0], target) for p in pts]
    ms2 = [run_point(p, seeds[1], target) for p in pts]

    print("  habitable fraction as a function of the isotropy threshold tau")
    print("  (the other four criteria are held exactly as pre-registered)")
    print()
    print("  %10s %12s %12s %10s" % ("tau", "seed A", "seed B", "flip rate"))
    curve = {}
    for tau in TAU_GRID:
        h1 = np.array([habitable(x, tau) for x in ms], bool)
        h2 = np.array([habitable(x, tau) for x in ms2], bool)
        curve[tau] = h1
        label = "no limit" if tau > 100 else "%.2f" % tau
        star = "  <- pre-registered" if abs(tau - TAU_REGISTERED) < 1e-9 else ""
        print("  %10s %11.4f %12.4f %10.3f%s"
              % (label, h1.mean(), h2.mean(), float((h1 != h2).mean()), star))
    print()

    # operating threshold: the tightest one that leaves a set worth measuring
    tau_op = next((t for t in TAU_GRID if curve[t].mean() >= 0.05), TAU_GRID[-1])
    hab = curve[tau_op]
    print("  Operating threshold for the rest of this run: tau = %.2f"
          % tau_op)
    print("  (the tightest grid value leaving at least 5% of the box alive)")
    print()

    print("=" * 78)
    print("2. WHICH CRITERION DOES THE WORK")
    print("=" * 78)
    for tau, name in ((TAU_REGISTERED, "pre-registered tau"),
                      (tau_op, "operating tau")):
        h = curve[tau]
        nf = int((~h).sum())
        fails = {c: 0 for c in verdict(ms[0]).keys()}
        for x in ms:
            v = verdict(x, tau)
            if not all(v.values()):
                for c, good in v.items():
                    if not good:
                        fails[c] += 1
        print("  %s = %.2f -- %d failures" % (name, tau, nf))
        for c in sorted(fails, key=lambda k: -fails[k]):
            print("     %-11s %5d   %5.1f%%"
                  % (c, fails[c], 100.0 * fails[c] / max(nf, 1)))
        print()

    if not hab.any():
        print("  Nothing habitable anywhere in the box at any threshold.")
        return

    print("=" * 78)
    print("3. WHICH KNOBS ACTUALLY CONSTRAIN   (tau = %.2f)" % tau_op)
    print("=" * 78)
    print("  Habitability rate by knob value, straight off the joint samples.")
    print("  A knob that constrains shows a strong gradient down its column.")
    print()
    both = np.concatenate([hab, np.array([habitable(x, tau_op) for x in ms2])])
    allp = pts + pts
    print("  %-9s %22s %22s %22s" % ("knob", "low third", "middle", "high third"))
    for name in KNOBS:
        vals = np.array([p[name] for p in allp], float)
        lo, hi = np.quantile(vals, [1 / 3, 2 / 3])
        cells = []
        for sel in (vals <= lo, (vals > lo) & (vals <= hi), vals > hi):
            if sel.sum() == 0:
                cells.append("%22s" % "-")
            else:
                cells.append("%13s %8.3f"
                             % ("[%.2f,%.2f]" % (vals[sel].min(),
                                                 vals[sel].max()),
                                float(both[sel].mean())))
        print("  %-9s %s" % (name, " ".join(cells)))
    print()

    print("=" * 78)
    print("4. ONE-AT-A-TIME vs JOINT   (at tau = %.2f)" % tau_op)
    print("=" * 78)
    print("  Habitability here is marginal -- the seed flip rate at this tau")
    print("  is %.0f%% -- so an axis point is scored as the FRACTION of seeds"
          % (100 * float((hab != np.array([habitable(x, tau_op)
                                           for x in ms2])).mean())))
    print("  that pass, not as unanimity. That makes the axis estimator the")
    print("  same kind of average as the joint one, so the ratio is honest.")
    print("  Three anchors, because where you stand changes the answer -- and")
    print("  the real argument always stands at our own universe.")
    print()
    U = np.array([unit(p) for p in pts])
    centre = U[hab].mean(axis=0)
    order = [i for i in np.argsort(((U - centre) ** 2).sum(axis=1))
             if hab[i] and habitable(ms2[i], tau_op)]
    anchors = [pts[i] for i in order[:3]]
    aseeds = (11, 23, 41, 59)
    naives = []
    for a_i, anchor in enumerate(anchors):
        print("  anchor %d:  %s" % (a_i + 1, "  ".join(
            "%s=%s" % (k, ("%.2f" % anchor[k])
                       if isinstance(anchor[k], float) else anchor[k])
            for k in KNOBS)))
        naive = 1.0
        row = []
        for name in KNOBS:
            ps = []
            for v in AXIS[name]:
                p = dict(anchor)
                p[name] = v
                ps.append(np.mean([habitable(run_point(p, s, target), tau_op)
                                   for s in aseeds]))
            frac = float(np.mean(ps))
            naive *= frac
            row.append("%s %.2f" % (name, frac))
        naives.append(naive)
        print("     " + "   ".join(row))
        print("     naive product %.5f" % naive)
    print()
    print("  measured joint fraction          %.5f" % hab.mean())
    print("  naive products                   %s"
          % ", ".join("%.5f" % x for x in naives))
    good = [x for x in naives if x > 0]
    if good:
        print("  joint / naive                    %s"
              % ", ".join("%.1f" % (hab.mean() / x) for x in good))
    print()

    print("=" * 78)
    print("5. THE PRIOR: same samples, different measure   (tau = %.2f)"
          % tau_op)
    print("=" * 78)

    def reweight(f):
        w = np.array([f(p) for p in pts], float)
        return float((w * hab).sum() / w.sum())

    print("  %-44s %s" % ("prior over the box", "habitable"))
    print("  %-44s %.4f" % ("uniform", hab.mean()))
    print("  %-44s %.4f" % (
        "log-uniform in alpha, beta, p_close, p_open",
        reweight(lambda p: 1.0 / ((p["alpha"] + 0.01) * (p["beta"] + 0.01)
                                  * (p["p_close"] + 0.01)
                                  * (p["p_open"] + 0.01)))))
    print("  %-44s %.4f" % (
        "1/sqrt in alpha only",
        reweight(lambda p: 1.0 / np.sqrt(p["alpha"] + 0.01))))
    print("  %-44s %.4f" % (
        "linear tilt toward large p_close",
        reweight(lambda p: p["p_close"] + 0.01)))
    print()

    print("=" * 78)
    print("6. THE BOX: the knobs with no principled outer edge")
    print("=" * 78)
    print("  alpha, cap and reach are unbounded above in principle. Nothing")
    print("  in the model says cap <= 12. Widening tests whether the fraction")
    print("  is a property of the model or of where we drew the box.")
    print()
    rng2 = np.random.default_rng(77)

    def wide(alpha_hi, cap_hi, reach_hi, n):
        out = []
        for _ in range(n):
            p = draw(rng2)
            p["alpha"] = float(rng2.uniform(0.0, alpha_hi))
            p["cap"] = int(rng2.integers(3, cap_hi + 1))
            p["reach"] = int(rng2.integers(2, reach_hi + 1))
            out.append(habitable(run_point(p, seeds[0], target), tau_op))
        return float(np.mean(out))

    print("  %-34s %s" % ("box", "habitable"))
    print("  %-34s %.4f   (%d samples)"
          % ("alpha<=3, cap<=12, reach<=5", hab.mean(), n_joint))
    print("  %-34s %.4f   (%d samples)"
          % ("alpha<=9, cap<=12, reach<=5", wide(9.0, 12, 5, n_wide), n_wide))
    print("  %-34s %.4f   (%d samples)"
          % ("alpha<=3, cap<=30, reach<=5", wide(3.0, 30, 5, n_wide), n_wide))
    print("  %-34s %.4f   (%d samples)"
          % ("alpha<=3, cap<=12, reach<=10", wide(3.0, 12, 10, n_wide), n_wide))
    print()


if __name__ == "__main__":
    main()
