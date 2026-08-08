#!/usr/bin/env python3
"""
resolution_optimum.py

Does the arrow peak at a particular coarse-graining resolution?

Earlier scans showed dS_B = 0 at both extremes -- one block gives Omega = N
everywhere, the discrete partition gives Omega = 1 everywhere -- and nonzero
in between, which looked like a window with an optimum somewhere inside it.
This tests whether the optimum is a property of *resolution* or of something
else.

A theorem first. If the prepared state is a uniformly random microstate, it
lands in block M with probability Omega(M)/N, so

    E[S_prep] = sum_M (Omega(M)/N) ln Omega(M) = S_eq

and therefore E[dS_B] = 0 for *every* coarse-graining, at every resolution.
Resolution alone cannot produce an arrow. A positive dS_B means the
preparation is atypical for that particular partition, not that the partition
has the right grain.

So the well-posed question is: holding the seed fixed, which partition makes
it most atypical, and does the answer sit at intermediate resolution?

Three families are compared at each block count k:

    random      partitions into k blocks, drawn uniformly
    isolating   the seed alone in its own block, the rest split at random
    natural     arity, degree sequence, Weisfeiler-Leman -- structured
                coarse-grainings, plotted at whatever k they happen to have

If the natural coarse-grainings beat random at their own k, the arrow comes
from structure. If the isolating family beats everything and does so at small
k, then there is no intrinsic optimal resolution at all -- only a most
revealing partition, which is a different claim.

Result: no interior optimum. The earlier "window" reading was wrong.
--------------------------------------------------------------------
Random partitions give mean dS_B of essentially zero at every resolution
(-0.027 to +0.040 over k = 2 to 304), confirming the theorem.

The isolating family falls monotonically with k. On path(8), N = 304:

    k          2      3      4      8     30    152    304
    dS_B    5.695  5.005  4.603  3.765  2.388  0.975  0.568

The maximum sits at the coarsest non-trivial resolution, k = 2, at 5.695
against ln N = 5.717 -- the ceiling. Isolating the seed against a single
lumped block captures the whole available arrow, and nothing peaks in
between.

So the earlier scan showing 0 at both extremes and 3.39 in the middle was not
a peak in resolution. The two zeros are separate degeneracies, Omega = N
constant at one end and Omega = 1 at the other, and everything between them
was tracking structure rather than grain.

What does survive is that structure dominates resolution. Against 400 random
partitions at matched k, arity and degseq both sit at the 100th percentile,
wl0 at the 90th. degseq even beats the isolating construction at comparable k
(2.560 against about 2.2) because it isolates the seed *and* leaves the other
blocks unequal, which raises S_eq.

The honest statement replacing "optimal resolution" is therefore: the arrow
is maximised by a partition that isolates the preparation while leaving the
rest lumped, not by any particular grain. The right level of description is
not a scale but whichever description makes the prepared state atypical --
a statement about relevance, close to a sufficient statistic, and weaker than
an intrinsic scale the system picks out for itself.
"""

import numpy as np
from collections import defaultdict

from rewrite_lab import continuation_graph, path_seed
from rewrite_research import RULE_SETS
from rewrite_thermo import COARSE


def gap(blocks_of, states, seed_key):
    """dS_B = S_eq - S_prep for a labelling function on the state list."""
    blocks = defaultdict(int)
    for s in states:
        blocks[blocks_of[s]] += 1
    n = len(states)
    S_eq = sum((c / n) * np.log(c) for c in blocks.values())
    return S_eq - np.log(blocks[blocks_of[seed_key]])


def random_partition(rng, states, k):
    lab = rng.integers(0, k, size=len(states))
    return {s: int(lab[i]) for i, s in enumerate(states)}


def isolating_partition(rng, states, k, seed_key):
    """Seed alone; everything else split at random into k-1 blocks."""
    others = [s for s in states if s != seed_key]
    lab = {seed_key: 0}
    if k > 1:
        r = rng.integers(1, k, size=len(others))
        for i, s in enumerate(others):
            lab[s] = int(r[i])
    else:
        for s in others:
            lab[s] = 0
    return lab


def main():
    for path_len in (6, 8):
        seed = path_seed(path_len)
        states_d, graph, _ = continuation_graph(
            seed, rules=RULE_SETS["G1-reversible"]["rules"])
        states = sorted(states_d)
        seed_key = seed.canonical()
        n = len(states)
        rng = np.random.default_rng(4)

        print("=" * 70)
        print("path(%d):  N = %d states,  ln N = %.4f" % (path_len, n, np.log(n)))
        print("=" * 70)

        print("%6s %12s %12s %12s" % ("k", "random mean", "random max", "isolating"))
        best_rand = (-1, None)
        for k in sorted(set([2, 3, 4, 6, 8, 12, 20, 30, n // 2, n - 1, n])):
            if k < 1 or k > n:
                continue
            R = 400
            vals = [gap(random_partition(rng, states, k), states, seed_key)
                    for _ in range(R)]
            iso = np.mean([gap(isolating_partition(rng, states, k, seed_key),
                               states, seed_key) for _ in range(60)])
            if max(vals) > best_rand[0]:
                best_rand = (max(vals), k)
            print("%6d %12.5f %12.5f %12.5f"
                  % (k, float(np.mean(vals)), float(max(vals)), float(iso)))

        print()
        print("  natural coarse-grainings (structured):")
        print("  %-10s %6s %12s %14s" % ("name", "k", "dS_B", "vs random@k"))
        for cn in ("trivial", "arity", "degseq", "wl0", "wl1"):
            fn = COARSE[cn]
            lab = {s: fn(states_d[s]) for s in states}
            k_nat = len(set(lab.values()))
            g = gap(lab, states, seed_key)
            comp = [gap(random_partition(rng, states, k_nat), states, seed_key)
                    for _ in range(400)]
            pct = float((np.array(comp) < g).mean())
            print("  %-10s %6d %12.5f %13.1f%%" % (cn, k_nat, g, 100 * pct))
        print()
        print("  E[dS_B] over random partitions should be 0 at every k.")
        print("  theoretical max over ALL partitions = ln N = %.4f,"
              % np.log(n))
        print("  attained by isolating the seed against one other block (k=2).")
        print()


if __name__ == "__main__":
    main()
