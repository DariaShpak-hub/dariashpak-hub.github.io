#!/usr/bin/env python3
"""
correlation_length.py

Graph correlation length xi, and whether 1/xi^2 supplies a small parameter.

The proposal under test replaces an inserted fundamental scale with three
quantities computed from the rewrite system itself: the entropy gap, an
emergent dimension, and a graph correlation length, giving

    Lambda_eff = (1 / xi^2) F(Delta S_B, d_eff, ...).

Unlike an effective area, this could in principle produce smallness without
tuning: if xi grows with the state space then 1/xi^2 falls as a power law.

xi is measured from the relaxation rate of the reversible chain --
D(p_t || uniform) decays like lambda_2^(2t), and xi = -1/ln lambda_2 -- so it
is a length in rewrites, a pure count.

Measured over N = 7 to 880, xi rises from 4.4 to 26.5 and 1/xi^2 falls from
5.3e-2 to 1.4e-3, so the mechanism does generate suppression. But the ratio
xi / ln N is flat at about 3.9 from N = 41 onward, indicating xi ~ ln N rather
than a power of N. Under a logarithmic law, reaching 1e-122 needs ln N ~ 1e61,
that is N ~ exp(1e61). Even reading the (poorer) power-law fit xi ~ N^0.33 at
face value requires N ~ 1e183. Two decades of N cannot decisively separate
these, but neither derives the smallness: both convert "why is Lambda tiny"
into "why is the state space that large".
"""
import sys, math
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from collections import defaultdict
from rewrite_lab import continuation_graph, path_seed
from rewrite_research import RULE_SETS
from rewrite_thermo import metropolis, evolve, shannon

rules = RULE_SETS["G1-reversible"]["rules"]


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


print("%6s %8s %10s %10s %10s" % ("N", "diameter", "xi_spec", "1/xi^2", "xi/lnN"))
Ns, xis = [], []
for n in range(4, 10):
    seed = path_seed(n)
    states, graph, _ = continuation_graph(seed, rules=rules)
    keys = list(states)
    N = len(keys)
    P = metropolis(graph, keys)

    # diameter via BFS from the seed (graph distance in rewrites)
    dist = {seed.canonical(): 0}
    frontier = [seed.canonical()]
    while frontier:
        nxt = []
        for u in frontier:
            for v in graph.get(u, ()):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    nxt.append(v)
        frontier = nxt
    diam = max(dist.values())

    # relaxation: H(t) = D(p_t || uniform) decays like lambda_2^(2t)
    H = []
    for t, p in evolve(P, {seed.canonical(): 1.0}, 400):
        h = math.log(N) - shannon(p)
        if h <= 1e-14:
            break
        H.append((t, math.log(h)))
    tail = H[len(H) // 2:]
    rate = slope([t for t, _ in tail], [h for _, h in tail])   # = 2 ln lambda2
    lam2 = math.exp(rate / 2.0)
    xi = -1.0 / math.log(lam2) if 0 < lam2 < 1 else float("inf")

    print("%6d %8d %10.3f %10.3e %10.3f" % (N, diam, xi, 1.0 / xi**2, xi / math.log(N)))
    Ns.append(N); xis.append(xi)

beta = slope([math.log(N) for N in Ns], [math.log(x) for x in xis])
print()
print("fit  xi ~ N^beta  with beta = %.4f" % beta)
print("so   1/xi^2 ~ N^%.4f" % (-2 * beta))
need = 122 * math.log(10) / (2 * beta) if beta > 0 else float("inf")
print("to reach 1/xi^2 ~ 1e-122 requires N ~ e^%.1f = 1e%.0f"
      % (need, need / math.log(10)))
