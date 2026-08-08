#!/usr/bin/env python3
"""
lambda_gradient.py

Does volume entropy correlate with the entropic-gradient observables?

The two surviving dimensionless quantities of this investigation measure
apparently different things:

    lambda   -- volume entropy, lim (1/k) ln|B_k|, the growth rate of
                accessible state space. A property of how the grammar
                expands.
    gradient -- the static entropic ordering: g(M), its correlation with
                ln Omega(M), and the exact entropy gap Delta S_B. A property
                of how multiplicity is distributed.

If expansion and thermodynamic ordering are one phenomenon, as an
entropic-cosmology reading would need, these should covary across grammars.
If they are independent knobs, they should not.

Method
------
Every non-empty subset of the five available rules is a grammar. For each,
the depth-k ball around a fixed seed is enumerated under a common state cap,
and both observables are measured *on that same truncated ball* so the
truncation enters both consistently. Only edges internal to the ball are used
for the gradient, since Omega is undefined for targets outside it.

Truncation is the main caveat: for grammars with unbounded state spaces the
ball is not the whole system, and the exact "global gradient sums to zero"
identity holds only up to boundary terms.

Result: the question is not answerable with this rule set
---------------------------------------------------------
Enumerating to a fixed depth gives corr(lambda, entropy gap) = +0.98, but the
correlation is spurious. At fixed depth a faster-growing grammar mechanically
reaches more states, so corr(lambda, ln N) = +0.989, and the entropy gap is
known independently to scale as ~0.67 ln N. Controlling for ln N drops the
partial correlation to -0.23, and with collinearity that high the partial is
itself unstable.

Re-running at matched state count does not rescue it. The eight grammars that
reach N ~ 1000 fall into two clusters, the correlation changes sign
(r = -0.93 with Spearman only -0.37, so it is not monotone), and within the
cluster that actually varies in lambda (0.507 to 0.593, a 17% spread) the gap
is flat to 0.6%. Several entries are also exact duplicates -- adding 'smooth'
to a grammar containing 'contract' changes nothing -- so the effective sample
is around five, not nineteen.

The flat gap across an appreciable lambda range is weak evidence that the two
observables are independent, consistent with the earlier structural finding
that entropy production and volume growth are separately dialable. But it is
weak. Settling this needs a family of genuinely distinct grammars large
enough to break the lambda-N collinearity, which five rules cannot supply.
"""

import argparse
import math
from collections import defaultdict
from itertools import combinations

from rewrite_lab import Hypergraph, continuation_graph, path_seed
from rewrite_research import COMPOSE, DECOMPOSE, CONTRACT
from expansion import SUBDIVIDE, SMOOTH
from rewrite_thermo import COARSE

ALL_RULES = [COMPOSE, DECOMPOSE, CONTRACT, SUBDIVIDE, SMOOTH]


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def explore(seed, rules, max_depth, max_states):
    """Depth-indexed ball plus the induced subgraph on it."""
    seen = {seed.canonical(): seed}
    frontier = [seed]
    vols = [1]

    for _ in range(max_depth):
        nxt = []
        for H in frontier:
            for r in rules:
                for G in r.apply(H):
                    c = G.canonical()
                    if c not in seen:
                        seen[c] = G
                        nxt.append(G)
        vols.append(len(seen))
        frontier = nxt
        if not frontier or len(seen) >= max_states:
            break

    graph = defaultdict(set)
    for c, H in seen.items():
        for r in rules:
            for G in r.apply(H):
                gc = G.canonical()
                if gc in seen:
                    graph[c].add(gc)
    return seen, graph, vols


def volume_entropy(vols):
    """lambda = mean of ln|B_k| - ln|B_{k-1}| over the last third."""
    if len(vols) < 4 or vols[-1] == vols[-2]:
        return 0.0
    d = [math.log(vols[i]) - math.log(vols[i - 1])
         for i in range(1, len(vols)) if vols[i] > vols[i - 1]]
    if not d:
        return 0.0
    tail = d[len(d) // 2:]
    return sum(tail) / len(tail)


def gradient_stats(states, graph, seed, cg_name):
    fn = COARSE[cg_name]
    label = {k: fn(states[k]) for k in states}
    blocks = defaultdict(list)
    for k in states:
        blocks[label[k]].append(k)
    omega = {M: len(v) for M, v in blocks.items()}
    N = len(states)

    acc = defaultdict(float)
    cnt = defaultdict(int)
    for u in states:
        su = math.log(omega[label[u]])
        for v in graph.get(u, ()):
            acc[label[u]] += math.log(omega[label[v]]) - su
            cnt[label[u]] += 1
    g = {M: (acc[M] / cnt[M] if cnt[M] else 0.0) for M in omega}

    eq = sum((omega[M] / N) * math.log(omega[M]) for M in omega)
    start = math.log(omega[label[seed.canonical()]])

    xs = [math.log(omega[M]) for M in omega]
    ys = [g[M] for M in omega]
    return {
        "gap": eq - start,
        "g_seed": g[label[seed.canonical()]],
        "corr_g": pearson(xs, ys),
        "mean_abs_g": sum(abs(v) for v in g.values()) / len(g),
        "macro": len(omega),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--depth", type=int, default=11)
    ap.add_argument("--cap", type=int, default=2500)
    ap.add_argument("--path", type=int, default=4)
    ap.add_argument("--coarse", default="degseq", choices=sorted(COARSE))
    ap.add_argument("--min-states", type=int, default=12)
    args = ap.parse_args(argv)

    seed = path_seed(args.path)
    rows = []

    print("%-26s %6s %7s %8s %8s %8s %8s" % (
        "grammar", "N", "lambda", "gap", "g_seed", "corr_g", "mean|g|"))

    for size in range(1, len(ALL_RULES) + 1):
        for combo in combinations(ALL_RULES, size):
            states, graph, vols = explore(seed, list(combo),
                                          args.depth, args.cap)
            if len(states) < args.min_states:
                continue
            lam = volume_entropy(vols)
            st = gradient_stats(states, graph, seed, args.coarse)
            name = "+".join(r.name[:4] for r in combo)
            print("%-26s %6d %7.4f %+8.4f %+8.4f %+8.4f %8.4f" % (
                name, len(states), lam, st["gap"], st["g_seed"],
                st["corr_g"], st["mean_abs_g"]))
            rows.append((lam, st))

    print()
    lam = [r[0] for r in rows]
    print("across %d grammars, correlation of lambda with:" % len(rows))
    for key, tag in (("gap", "entropy gap dS_B"),
                     ("g_seed", "gradient at seed"),
                     ("corr_g", "corr(g, lnOmega)"),
                     ("mean_abs_g", "mean |g|")):
        ys = [r[1][key] for r in rows]
        ok = [(a, b) for a, b in zip(lam, ys) if b == b]
        if len(ok) < 3:
            continue
        xs2 = [a for a, _ in ok]
        ys2 = [b for _, b in ok]
        print("  %-20s  pearson=%+.4f   spearman=%+.4f  (n=%d)" % (
            tag, pearson(xs2, ys2), spearman(xs2, ys2), len(ok)))


if __name__ == "__main__":
    main()
