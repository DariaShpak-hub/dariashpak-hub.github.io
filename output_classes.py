#!/usr/bin/env python3
"""
output_classes.py

The fine-tuning question inverted: stop specifying a target, and ask what the
model actually makes.

fine_tuning.py measured how many universes hit a target I chose, and the
answer turned out to depend more on the target than on the model -- 0 out of
600 at the isotropy threshold I fixed in advance, 0.117 after relaxing one
number, 0.23 with that criterion deleted. That is the whole "fine-tuned for
WHAT" problem in one line: the target was mine, and it did most of the work.

So this drops the target entirely. Same 600 parameter points, same runs, but
now nothing is scored as pass or fail. Instead the universes are clustered by
what they measurably ARE -- dimension, mean distance, degree, shear, scale
factor, degree spread -- and the question becomes:

    does the output have structure of its own, or is it a smear?

The distinction matters. If distinct families fall out, those are outcomes the
model selects, and they are targets in the only non-circular sense available:
picked by the dynamics rather than by me. If instead it is one continuous blob,
then there is nothing the model "aims at", every target is imposed from
outside, and the fine-tuning framing has no purchase here at all.

Not begging the question
------------------------
k-means always returns k clusters, so finding clusters proves nothing on its
own. Three guards:

    gap statistic   Tibshirani's test. Compare within-cluster scatter to the
                    same quantity on uniform reference data filling the same
                    PCA-aligned box. Real structure beats the reference; a
                    smear does not. The test is run from k = 1, so it is
                    allowed to answer "no clusters".

    silhouette      how well separated the chosen partition actually is.
                    Below about 0.25 means the labels are decorative.

    predictability  if a family is real it should be reachable -- some region
                    of parameter space should produce it. Tested by leave-one-
                    out nearest-neighbour classification from parameters to
                    class, against a label-shuffled baseline.

Degenerate runs -- the ones where the shell estimator has no scaling range and
returns nan -- are held out as their own class rather than dropped, because
"produced nothing measurable" is itself an outcome.

No scipy or sklearn in this environment, so k-means, the gap statistic, the
silhouette and the PCA are all written out directly.

Result 1: there are no families. It is one continuum.
------------------------------------------------------
All 600 runs were measurable; none degenerated. The gap statistic never turns
over -- it climbs from k = 1 to the edge of the grid and is still climbing:

    k                        1      2      3      4      5      6      7      8
    3 well-separated blobs  0.05   0.00   1.54   1.38   1.31   1.29   1.28   1.28
    uniform, no structure   1.14   0.99   0.99   0.99   0.99   0.99   1.00   1.00
    this model's output     1.04   1.06   1.18   1.28   1.36   1.42   1.46   1.48

The controls make the reading unambiguous. Three genuine clusters produce a
peak at exactly k = 3 and then a decline. Structureless uniform data produces a
flat line. Our output does neither: it keeps improving with every extra cluster,
which is what a continuum does, because a continuum can always be cut again.

Two notes against myself. The absolute gap values are inflated for all three
datasets -- the uniform control should read near 0 and reads near 1.0 -- because
a rotated bounding box encloses more volume than the data inside it. Only the
shape of the curve is trustworthy, which is why the controls are here. And
Tibshirani's "first k satisfying the rule" prescription is not used, because on
the positive control it fires at k = 1 on data that visibly has three clusters.
The peak location is the robust statistic; the rule is not.

A second test, independent of any of that: bin the leading principal component
into 24 bins. Zero are empty. The largest gap between neighbouring universes is
0.197 against a mean spacing of 0.014 -- a ratio of 14 in 600 points, which is
ordinary order-statistic scatter, not a hole. The output has no seams in it.

Result 2: what the continuum is
---------------------------------
Cutting it into 8 bins anyway, purely to describe it, and ordering those bins
by the leading component:

    bin    PC1      d   mean dist    <k>   shear   mean a
     7   -4.19   1.02      418       2.01    0.21     275
     2   -2.52   1.09      240       2.09    0.89     158
     4   -2.44   1.27      141       2.15    1.00      98
     3   -1.23   1.59       54       2.31    1.62      36
     1   -0.17   2.14       28       2.68    1.45      18
     5    0.74   2.76       14       2.83    1.82       8
     6    1.72   3.45        9       3.37    1.48       5
     0    2.55   3.93        7       3.42    1.45       4

Dimension, mean distance, mean degree and scale factor are each monotone across
all eight bins, with |correlation| to PC1 of 0.95, 0.98, 0.94 and 0.98. One
axis carries 60% of the variance and every one of those four observables is
essentially a function of position along it.

So the model makes exactly one thing: a dial running from a THREAD -- d = 1,
mean distance 418, enormous scale factor -- to a SMALL-WORLD graph -- d = 3.9,
mean distance 7, barely expanded. Nothing else. There is no separate family of
outputs anywhere in the box.

And the dial is driven by two knobs. Regressing PC1 on all six, standardised:

    p_close   +1.508      alone, R^2 = 0.578
    reach     +0.986      alone, R^2 = 0.265
    cap       +0.374      alone, R^2 = 0.016
    p_open    -0.149
    beta      +0.020
    alpha     -0.048      alone, R^2 = 0.002
    all six    R^2 = 0.841

The close rate and the chord length set where on the dial a universe lands, and
they explain 84% of it. The two anti-Yule exponents contribute nothing
measurable. Position on the dial is also predictable in the other direction:
leave-one-out nearest neighbour from parameters to bin gets 0.502, against 0.147
for shuffled labels and 0.217 for the majority class.

Result 3: why fine_tuning.py found zero, in one line
------------------------------------------------------
Look at the shear column. It is 0.21 in bin 7 and between 0.89 and 1.82 in
every other bin -- and bin 7 is the one at d = 1.02.

The only isotropic thing this model produces is a thread.

That is why nothing in the box passed. Shear is not controlled by the axis the
model actually varies: its correlation with PC1 is only 0.53, and the degree
spread's is -0.09. Isotropy lives mainly on PC2, an independent direction
carrying 19% of the variance. So moving along the dial that creates dimension
does not buy isotropy back, at any setting, because it is not the same
direction. fine_tuning.py's 0-out-of-600 was not bad luck in the sampling; it
is a statement about which directions in observable space the rules can reach.

Result 4: three is not a place
--------------------------------
d passes through 3 on the way from 1 to 3.9. It is an intermediate value on a
smooth dial, with nothing marking it -- no plateau, no density peak, no bin
boundary. If a universe here is three-dimensional it is because someone chose
p_close near 0.7, not because the model prefers it.

That is the same negative as dimension_attractor.py and coloured_memory.py,
but reached without assuming a target, which makes it stronger: it is not that
we failed to find the mechanism selecting 3, it is that the output space has no
feature at 3 to select.

Verdict, and what it says about "fine-tuned for what"
-------------------------------------------------------
Asked without a target, the model has no targets. Its entire output is a
one-parameter family, indexed mostly by the close rate, running smoothly from
thread to small-world with no gaps, no clusters and no preferred point. Every
"kind of universe" in fine_tuning.py was a window I drew on that dial.

Which answers the question that prompted this. "Fine-tuned for what" has no
internal answer here -- not because the answer is hard, but because there is
nothing in the output distribution for the "what" to refer to. The targets were
never in the physics. They were in the observer, and in this case the observer
was me, writing five criteria before pressing run.

The physics that survives is Result 3, and it is worth more than the philosophy:
isotropy and dimension are close to orthogonal directions in this model's
output space, so no amount of moving along the one reachable axis will ever
deliver both.
"""

import os
import numpy as np

from fine_tuning import Cosmos, draw, run_point, KNOBS

CACHE = os.path.join(
    os.environ.get("SCRATCH", "/tmp"), "output_classes_600.npz")

FEATURES = ("d", "log dist", "<k>", "log shear", "log a", "spread", "power")


def collect(n=600, seed=11, target=2000, rng_seed=5):
    """The same 600 points fine_tuning.py drew, measured and kept."""
    rng = np.random.default_rng(rng_seed)
    pts = [draw(rng) for _ in range(n)]
    ms = [run_point(p, seed, target) for p in pts]
    return pts, ms


def featurise(ms):
    rows, ok = [], []
    for m in ms:
        good = (np.isfinite(m["d"]) and np.isfinite(m["shear"])
                and m["shear"] > 0 and m["mean_dist"] > 0)
        ok.append(good)
        if not good:
            rows.append([np.nan] * len(FEATURES))
            continue
        rows.append([m["d"],
                     np.log10(m["mean_dist"]),
                     m["mean_k"],
                     np.log10(m["shear"]),
                     np.log10(max(m["mean_a"], 1e-6)),
                     m["max_k"] / m["mean_k"],
                     1.0 if m["power"] else 0.0])
    return np.array(rows, float), np.array(ok, bool)


def zscore(X):
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def kmeans(X, k, rng, restarts=12, iters=120):
    best, best_w = None, np.inf
    n = len(X)
    for _ in range(restarts):
        # k-means++ seeding
        cent = [X[rng.integers(n)]]
        for _ in range(k - 1):
            d2 = np.min(((X[:, None, :] - np.array(cent)[None]) ** 2
                         ).sum(axis=2), axis=1)
            tot = d2.sum()
            if tot <= 0:
                cent.append(X[rng.integers(n)])
                continue
            cent.append(X[int(np.searchsorted(np.cumsum(d2),
                                              rng.random() * tot))])
        C = np.array(cent, float)
        lab = np.zeros(n, int)
        for _ in range(iters):
            d2 = ((X[:, None, :] - C[None]) ** 2).sum(axis=2)
            new = np.argmin(d2, axis=1)
            if np.array_equal(new, lab):
                break
            lab = new
            for j in range(k):
                if (lab == j).any():
                    C[j] = X[lab == j].mean(axis=0)
        w = float(((X - C[lab]) ** 2).sum())
        if w < best_w:
            best_w, best = w, (lab.copy(), C.copy())
    return best[0], best[1], best_w


def gap_statistic(X, ks, rng, B=14):
    """Tibshirani: log W_k against uniform reference in the PCA-aligned box."""
    Xc = X - X.mean(axis=0)
    _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    Z = Xc @ Vt.T
    lo, hi = Z.min(axis=0), Z.max(axis=0)
    out = {}
    for k in ks:
        _, _, w = kmeans(X, k, rng)
        logw = np.log(max(w, 1e-300))
        refs = []
        for _ in range(B):
            R = rng.uniform(lo, hi, size=Z.shape) @ Vt + X.mean(axis=0)
            _, _, wr = kmeans(R, k, rng, restarts=4)
            refs.append(np.log(max(wr, 1e-300)))
        refs = np.array(refs)
        out[k] = (float(refs.mean() - logw),
                  float(refs.std() * np.sqrt(1.0 + 1.0 / B)))
    return out


def silhouette(X, lab):
    n = len(X)
    ks = np.unique(lab)
    if len(ks) < 2:
        return float("nan")
    D = np.sqrt(np.maximum(((X[:, None, :] - X[None]) ** 2).sum(axis=2), 0))
    s = np.zeros(n)
    for i in range(n):
        same = lab == lab[i]
        same[i] = False
        a = D[i, same].mean() if same.any() else 0.0
        b = min((D[i, lab == c].mean() for c in ks if c != lab[i]),
                default=0.0)
        s[i] = 0.0 if max(a, b) == 0 else (b - a) / max(a, b)
    return float(s.mean())


def loo_nn_accuracy(P, lab):
    """Can the parameters predict the class? Leave-one-out 1-NN."""
    D = ((P[:, None, :] - P[None]) ** 2).sum(axis=2)
    np.fill_diagonal(D, np.inf)
    return float((lab[np.argmin(D, axis=1)] == lab).mean())


def controls(n, dim, rng, ks, B=10):
    """Positive and negative control for the gap test, same n and dim."""
    out = {}
    G = np.vstack([rng.normal(c * 6.0, 1.0, size=(n // 3, dim))
                   for c in range(3)])
    out["3 well-separated blobs"] = gap_statistic(zscore(G), ks, rng, B=B)
    U = rng.uniform(0.0, 1.0, size=(n, dim))
    out["uniform, no structure"] = gap_statistic(zscore(U), ks, rng, B=B)
    return out


def main():
    print("=" * 78)
    print("0. THE SAMPLE")
    print("=" * 78)
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        pts, ms = list(z["pts"]), list(z["ms"])
        print("  loaded %d runs from cache" % len(ms))
    else:
        pts, ms = collect()
        np.savez(CACHE, pts=np.array(pts, dtype=object),
                 ms=np.array(ms, dtype=object))
        print("  ran %d universes at N = 2000" % len(ms))

    X_all, ok = featurise(ms)
    print("  %d measurable, %d degenerate (no scaling range -- held out as"
          " their own class)" % (int(ok.sum()), int((~ok).sum())))
    print()

    X = zscore(X_all[ok])
    P = np.array([[p[k] for k in KNOBS] for p in pts], float)[ok]
    P = zscore(P)
    rng = np.random.default_rng(19)

    print("=" * 78)
    print("1. IS THERE STRUCTURE AT ALL?")
    print("=" * 78)
    print("  Gap statistic against uniform reference data in the same box.")
    print("  Run from k = 1, so it is allowed to answer 'no clusters'.")
    print()
    ks = list(range(1, 9))
    gs = gap_statistic(X, ks, rng)
    print("  %4s %10s %10s %s" % ("k", "gap", "s(k)", ""))
    pick = max(gs, key=lambda k: gs[k][0])
    for k in ks:
        g, s = gs[k]
        print("  %4d %10.4f %10.4f%s"
              % (k, g, s, "  <- peak" if k == pick else ""))
    print()
    print("  Read the SHAPE, not the value. The absolute gap is inflated for")
    print("  every dataset here because a rotated bounding box has more volume")
    print("  than the data it encloses, which is a known weakness of the")
    print("  PCA-box reference in high dimension -- see the uniform control,")
    print("  which sits near 1.0 instead of near 0. What the shape says is")
    print("  clean: a peak means clusters, a flat line means no structure, and")
    print("  a curve still climbing at the edge of the grid means a continuum.")
    print("  (Tibshirani's first-k rule is not used: the control below shows")
    print("  it misfires at k = 1 on data with three genuine clusters.)")
    print()

    print("  Controls, same sample size and dimension:")
    ctl = controls(len(X), X.shape[1], rng, ks)
    print("     %-24s %s" % ("", "  ".join("%5d" % k for k in ks)))
    for name, g in ctl.items():
        pk = max(g, key=lambda k: g[k][0])
        print("     %-24s %s   peak at k = %d"
              % (name, "  ".join("%5.2f" % g[k][0] for k in ks), pk))
    print("     %-24s %s   peak at k = %d%s"
          % ("this model's output", "  ".join("%5.2f" % gs[k][0] for k in ks),
             pick, "  (grid edge: no turnover)" if pick == ks[-1] else ""))
    print()
    print("  A gap curve that keeps climbing and never turns over is what a")
    print("  CONTINUUM looks like: there is always another way to cut it, so")
    print("  more clusters always fit better. Clusters make the curve level")
    print("  off, as the positive control shows.")
    print()

    print("=" * 78)
    print("2. THE PARTITION")
    print("=" * 78)
    lab, C, _ = kmeans(X, max(pick, 2), rng)
    sil = silhouette(X, lab)
    print("  k = %d   silhouette = %.3f   (below ~0.25 means decorative)"
          % (max(pick, 2), sil))
    print()
    mu = X_all[ok].mean(axis=0)
    sd = X_all[ok].std(axis=0)
    sd[sd == 0] = 1.0
    print("  %5s %6s %s" % ("class", "n", "".join("%11s" % f
                                                  for f in FEATURES)))
    for j in range(max(pick, 2)):
        sel = lab == j
        vals = C[j] * sd + mu
        print("  %5d %6d %s" % (j, int(sel.sum()),
                                "".join("%11.2f" % v for v in vals)))
    if (~ok).any():
        print("  %5s %6d %s" % ("degen", int((~ok).sum()),
                                "%11s" % "no scaling range"))
    print()

    print("=" * 78)
    print("3. ARE THE CLASSES REACHABLE FROM PARAMETER SPACE?")
    print("=" * 78)
    acc = loo_nn_accuracy(P, lab)
    base = np.mean([loo_nn_accuracy(P, rng.permutation(lab))
                    for _ in range(20)])
    freq = np.bincount(lab).max() / len(lab)
    print("  leave-one-out 1-NN, parameters -> class   %.3f" % acc)
    print("  same with the labels shuffled             %.3f" % base)
    print("  majority-class baseline                   %.3f" % freq)
    print()
    print("  %5s %s" % ("class", "".join("%10s" % k for k in KNOBS)))
    Praw = np.array([[p[k] for k in KNOBS] for p in pts], float)[ok]
    for j in range(max(pick, 2)):
        sel = lab == j
        print("  %5d %s" % (j, "".join("%10.2f" % v
                                       for v in Praw[sel].mean(axis=0))))
    print()

    print("=" * 78)
    print("4. WHAT SEPARATES THEM")
    print("=" * 78)
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    for c in range(3):
        print("  PC%d (%.0f%% of variance)" % (c + 1, 100 * var[c]))
        order = np.argsort(-np.abs(Vt[c]))
        print("     " + "   ".join("%s %+.2f" % (FEATURES[i], Vt[c][i])
                                   for i in order[:4]))
    print()

    print("=" * 78)
    print("5. IS IT ONE CONTINUUM OR SEVERAL FAMILIES?")
    print("=" * 78)
    t = Xc @ Vt[0]
    print("  If the classes are bins on a single axis, then ordering them by")
    print("  PC1 orders every observable at once, with no crossings.")
    print()
    means = np.array([X_all[ok][lab == j].mean(axis=0)
                      for j in range(max(pick, 2))])
    tm = np.array([t[lab == j].mean() for j in range(max(pick, 2))])
    o = np.argsort(tm)
    print("  %5s %8s %s" % ("class", "PC1", "".join("%11s" % f
                                                    for f in FEATURES)))
    for j in o:
        print("  %5d %8.2f %s" % (j, tm[j],
                                  "".join("%11.2f" % v for v in means[j])))
    print()
    for i, f in enumerate(FEATURES):
        v = means[o, i]
        mono = np.all(np.diff(v) >= 0) or np.all(np.diff(v) <= 0)
        r = float(np.corrcoef(t, X_all[ok][:, i])[0, 1])
        print("     %-10s monotone across classes: %-5s   corr with PC1 %+.3f"
              % (f, "yes" if mono else "no", r))
    print()

    print("  A continuum should also have no holes in it. PC1, 24 bins:")
    print()
    cnt, edges = np.histogram(t, bins=24)
    for c, e in zip(cnt, edges):
        print("     %7.2f  %-40s %d" % (e, "#" * int(40 * c / cnt.max()), c))
    st = np.sort(t)
    sp = np.diff(st)
    print()
    print("     largest gap between neighbouring universes  %.3f" % sp.max())
    print("     mean spacing                                %.3f" % sp.mean())
    print("     ratio                                       %.1f"
          % (sp.max() / sp.mean()))
    print("     empty bins                                  %d of 24"
          % int((cnt == 0).sum()))
    print()

    print("  And which knob drives it. Least squares of PC1 on the six knobs,")
    print("  all standardised, so the coefficients are directly comparable:")
    print()
    A = np.column_stack([P, np.ones(len(P))])
    coef, *_ = np.linalg.lstsq(A, t, rcond=None)
    pred = A @ coef
    r2 = 1.0 - float(((t - pred) ** 2).sum() / ((t - t.mean()) ** 2).sum())
    for i, k in enumerate(KNOBS):
        print("     %-10s %+.3f" % (k, coef[i]))
    print("     %-10s %.3f" % ("R^2", r2))
    print()
    single = {}
    for i, k in enumerate(KNOBS):
        B2 = np.column_stack([P[:, i], np.ones(len(P))])
        c2, *_ = np.linalg.lstsq(B2, t, rcond=None)
        p2 = B2 @ c2
        single[k] = 1.0 - float(((t - p2) ** 2).sum()
                                / ((t - t.mean()) ** 2).sum())
    print("  R^2 from each knob alone:")
    for k in sorted(single, key=lambda x: -single[x]):
        print("     %-10s %.3f" % (k, single[k]))
    print()


if __name__ == "__main__":
    main()
