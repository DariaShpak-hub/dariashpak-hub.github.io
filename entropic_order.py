#!/usr/bin/env python3
"""
entropic_order.py

Thermodynamic ordering on a continuation graph, formulated without time.

Why this module exists
----------------------
The relaxation experiments in rewrite_thermo.py measure <Delta S_B> per step
of a Markov chain. That imports a background time parameter the rewrite
system does not itself provide: the step counter t is added structure, not
something the grammar supplies. If thermodynamic ordering is to be primitive
rather than derived from time, it must be definable on the static
continuation graph.

It can be. Boltzmann's asymmetry lives in the counting of states, not in the
dynamics, so it survives deleting the dynamics entirely.

The construction
----------------
Given a coarse-graining into macrostates with multiplicities Omega(M), define
for each macrostate the mean entropy difference across one rewrite:

    g(M) = mean over microstates u in M, and rewrites u -> v,
           of  [ ln Omega(macro(v)) - ln Omega(M) ]

No trajectory, no chain, no time -- g is a weighted count of edges in a fixed
graph.

Two facts then hold, and together they are the whole content of the arrow:

1. Summed over all directed edges, the entropy difference is exactly zero on
   a symmetric graph, since every edge is counted once in each direction.
   There is no global orientation. The grammar has no preferred direction.

2. Conditioned on a macrostate, g(M) is strongly negative in ln Omega(M):
   from a low-multiplicity macrostate almost every rewrite leads to a
   higher-multiplicity one, and from a high-multiplicity macrostate the
   reverse.

So "the system evolves toward equilibrium" is not a statement about time. It
is the statement that a *prepared* low-Omega macrostate sits on a positive
entropic gradient, in a state space where high-Omega macrostates vastly
outnumber low-Omega ones. Ordering by S_B is available directly; a time
parameter is one way to traverse that ordering, not a prerequisite for it.
"""

import argparse
import math
from collections import defaultdict

from rewrite_lab import continuation_graph, path_seed
from rewrite_research import GRAMMARS
from rewrite_thermo import COARSE


def macrostates(states, graph, cg_name):
    fn = COARSE[cg_name]
    label = {k: fn(states[k]) for k in states}
    blocks = defaultdict(list)
    for k in states:
        blocks[label[k]].append(k)
    omega = {M: len(v) for M, v in blocks.items()}
    return label, blocks, omega


def gradient(states, graph, label, omega):
    """
    Static entropic gradient per macrostate, and the global sum.

    Returns (g, total, n_edges) where g[M] is the mean change in ln Omega
    across one rewrite leaving a microstate of M, and total is the sum of
    that change over every directed edge in the graph.
    """
    acc = defaultdict(float)
    cnt = defaultdict(int)
    total = 0.0
    n_edges = 0

    for u in states:
        Mu = label[u]
        su = math.log(omega[Mu])
        for v in graph.get(u, ()):
            d = math.log(omega[label[v]]) - su
            acc[Mu] += d
            cnt[Mu] += 1
            total += d
            n_edges += 1

    g = {M: (acc[M] / cnt[M] if cnt[M] else 0.0) for M in omega}
    return g, total, n_edges


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def report(n, rules, cg_name):
    seed = path_seed(n)
    states, graph, _ = continuation_graph(seed, rules=rules)
    label, blocks, omega = macrostates(states, graph, cg_name)
    g, total, n_edges = gradient(states, graph, label, omega)

    # No global orientation: on a symmetric graph every edge contributes its
    # entropy difference once in each direction, so the sum must vanish.
    assert abs(total) < 1e-9 * max(n_edges, 1), \
        f"global entropy sum {total} is not zero -- graph is not symmetric"

    xs = [math.log(omega[M]) for M in omega]
    ys = [g[M] for M in omega]
    r = pearson(xs, ys)

    seed_M = label[seed.canonical()]
    N = len(states)
    up = sum(1 for M in omega if g[M] > 1e-12)
    dn = sum(1 for M in omega if g[M] < -1e-12)

    print(f"N={N:>5}  macro={len(omega):>4}  coarse={cg_name:<8} "
          f"sum_edges(dS)={total:+.2e}  corr(g, lnOmega)={r:+.4f}")
    print(f"        prepared seed macrostate: Omega={omega[seed_M]:>4}  "
          f"lnOmega={math.log(omega[seed_M]):.3f}  g={g[seed_M]:+.4f}")
    print(f"        macrostates with g>0: {up}   g<0: {dn}")

    ordered = sorted(omega, key=lambda M: omega[M])
    lo, hi = ordered[0], ordered[-1]
    print(f"        lowest  Omega={omega[lo]:>4}: g={g[lo]:+.4f}")
    print(f"        highest Omega={omega[hi]:>4}: g={g[hi]:+.4f}")
    print()
    return {"N": N, "corr": r, "g_seed": g[seed_M],
            "omega_seed": omega[seed_M], "total": total}


def entropy_gap(states, graph, cg_name, seed):
    """
    The total arrow, exactly, with no dynamics.

        Delta S_B = sum_M (Omega(M)/N) ln Omega(M)  -  ln Omega(M_0)

    The first term is the equilibrium mean of S_B: under a kernel with uniform
    stationary distribution the equilibrium macrostate weight is exactly
    Pi(M) = Omega(M)/N. The second is the prepared macrostate's entropy. Their
    difference reproduces the measured end-to-end relaxation to machine
    precision, so the entire quantity is fixed by counting -- no trajectory,
    no clock, no simulation.

    Its immediate consequence is that the *sign* of the arrow is not a fact
    about system size. It is negative whenever the prepared macrostate happens
    to sit above the equilibrium mean, and positive when below. What genuinely
    varies with size is not this endpoint difference but the smoothness of the
    path taken to it.
    """
    label, blocks, omega = macrostates(states, graph, cg_name)
    N = len(states)
    eq = sum((omega[M] / N) * math.log(omega[M]) for M in omega)
    start = math.log(omega[label[seed.canonical()]])
    return eq - start, eq, start


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--paths", type=int, nargs="+", default=[4, 5, 6, 7, 8, 9])
    ap.add_argument("--grammar", choices=sorted(GRAMMARS),
                    default="G1-reversible")
    ap.add_argument("--coarse", nargs="+", default=["degseq"],
                    choices=sorted(COARSE))
    args = ap.parse_args(argv)

    rules = GRAMMARS[args.grammar]["rules"]
    print("=" * 72)
    print(f"grammar = {args.grammar}   (static -- no time parameter)")
    print("=" * 72)
    print()

    for cg in args.coarse:
        for n in args.paths:
            report(n, rules, cg)


if __name__ == "__main__":
    main()

