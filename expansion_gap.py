#!/usr/bin/env python3
"""
expansion_gap.py

Does expansion open an entropy gap, or does the state track equilibrium?

The hypothesis is that a specially prepared initial state plus expansion into
an enlarging space of configurations yields a thermodynamic arrow. It has a
rate condition attached that is easy to omit: if the state equilibrates as
fast as the space grows, it stays typical and no gap ever opens. So the
quantity to watch is

    gap_k = ln|B_k| - S(p_k)

the shortfall between the entropy the state has and the entropy available in
the region it can currently reach.

Result: sustained expansion is what matters. A saturating rule set opens a
gap while it is still growing and then closes it again -- G1 peaks at 0.695
nats around k=4 and falls back to 0.457 once the ball stops growing at 41
states. Expanding rule sets never close it: 0.259 -> 2.083 under polynomial
growth, 0.202 -> 3.681 under exponential. The ordering follows the growth
rate exactly.

Two caveats worth keeping attached to this.

The effect is partly kinematic. gap_k >= 0 holds by construction, since p_k
is supported on B_k, and the gap widens whenever the ball grows faster than
a diffusive walk can spread into it. The result is a statement about
relative rates, not a discovery that expansion creates entropy.

And this is a coarse-grained gap throughout. Fine-grained volume is
conserved -- Liouville in the continuous case, and here the distribution
never gains information it did not start with. Nothing releases previously
inaccessible degrees of freedom; what grows is the region the state has
spread over, which is why the arrow remains a property of the coarse-graining.
"""
import sys, math
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from collections import defaultdict
from rewrite_lab import continuation_graph, path_seed, Hypergraph
from rewrite_research import COMPOSE, DECOMPOSE, RULE_SETS
from expansion import SUBDIVIDE, SMOOTH, star

def run(name, seed, rules, depth):
    seen = {seed.canonical(): seed}
    frontier = [seed]
    ball = [1]
    succ = defaultdict(set)
    for _ in range(depth):
        nxt = []
        for H in frontier:
            for r in rules:
                for G in r.apply(H):
                    c = G.canonical()
                    if c not in seen:
                        seen[c] = G; nxt.append(G)
        ball.append(len(seen))
        frontier = nxt
        if not frontier or len(seen) > 60000:
            break
    for c, H in seen.items():
        for r in rules:
            for G in r.apply(H):
                gc = G.canonical()
                if gc in seen:
                    succ[c].add(gc)

    p = {seed.canonical(): 1.0}
    print("%-34s %5s %8s %9s %9s %9s" % (name, "k", "|B_k|", "ln|B_k|", "S(p_k)", "gap"))
    for k in range(1, len(ball)):
        nxt = defaultdict(float)
        for u, m in p.items():
            s = succ.get(u, set())
            if not s:
                nxt[u] += m; continue
            nxt[u] += 0.5 * m
            for v in s:
                nxt[v] += 0.5 * m / len(s)
        p = dict(nxt)
        S = -sum(x*math.log(x) for x in p.values() if x > 0)
        lb = math.log(ball[k])
        if k % 2 == 0 or k == len(ball)-1:
            print("%-34s %5d %8d %9.4f %9.4f %9.4f" % ("", k, ball[k], lb, S, lb - S))
    print()

run("G1 finite (saturates)", path_seed(6), RULE_SETS["G1-reversible"]["rules"], 14)
run("subdivide+smooth (polynomial)", star(3), [SUBDIVIDE, SMOOTH], 14)
run("comp+decomp+subdivide (exp)", path_seed(2), [COMPOSE, DECOMPOSE, SUBDIVIDE], 12)
