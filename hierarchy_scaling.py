#!/usr/bin/env python3
"""
hierarchy_scaling.py

Does the entropic residual shrink with state-space size?

The proposal under test is that a large microscopic scale could yield a tiny
effective cosmological constant without tuning, via

    epsilon(N) = Delta S_B / S_eq  ~  N^-alpha,     alpha > 0.

Because Delta S_B = S_eq - S_0 exactly (see entropic_order.entropy_gap), the
residual is simply

    epsilon = 1 - S_0 / S_eq,

so the entire question is how the prepared macrostate's entropy tracks the
equilibrium entropy as the state space grows. That is directly measurable.

Result: it does not shrink. For coarse-grainings placing the seed in a
singleton macrostate (S_0 = 0) the residual is identically 1 at every size;
for the arity profile it *grows*, fitting alpha = -0.37. Delta S_B itself
diverges logarithmically rather than remaining finite, and S_eq diverges at
the same rate, so the ratio is scale-free.

The structural reason is that epsilon and the arrow are the same parameter.
A small residual requires S_0 -> S_eq, i.e. a prepared state of typical
entropy -- which is precisely the case with no arrow. Driving epsilon to
1e-120 would require the preparation entropy to match the equilibrium entropy
to 120 decimal places, relocating the fine-tuning rather than removing it.
"""
import sys, math
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from collections import defaultdict
from rewrite_lab import continuation_graph, path_seed
from rewrite_research import RULE_SETS
from entropic_order import macrostates, entropy_gap

rules = RULE_SETS["G1-reversible"]["rules"]

def slope(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs)
    return sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sxx if sxx > 0 else float("nan")

for cg in ("arity", "degseq", "wl0"):
    print("=== coarse-graining: %s ===" % cg)
    print("%6s %9s %9s %9s %9s %9s" % ("N", "S_eq", "S_0", "dS_B", "eps", "dS/lnN"))
    Ns, eps, gaps = [], [], []
    for n in range(4, 10):
        seed = path_seed(n)
        states, graph, _ = continuation_graph(seed, rules=rules)
        gap, eq, start = entropy_gap(states, graph, cg, seed)
        N = len(states)
        e = gap/eq if abs(eq) > 1e-15 else float("nan")
        print("%6d %9.4f %9.4f %+9.4f %+9.4f %9.4f" % (
            N, eq, start, gap, e, gap/math.log(N)))
        Ns.append(N); eps.append(e); gaps.append(gap)

    good = [(N, e) for N, e in zip(Ns, eps) if e == e and e > 0]
    if len(good) >= 3:
        a = -slope([math.log(N) for N, _ in good], [math.log(e) for _, e in good])
        print("  fitted alpha in eps ~ N^-alpha : %+.4f   (need alpha>0 to help)" % a)
    b = slope([math.log(N) for N in Ns], gaps)
    print("  dS_B ~ %.4f * ln N" % b)
    print()

# What preparation would make epsilon small?
print("=== achievable epsilon over ALL possible preparations, N=880, degseq ===")
seed = path_seed(9)
states, graph, _ = continuation_graph(seed, rules=rules)
label, blocks, omega = macrostates(states, graph, "degseq")
N = len(states)
eq = sum((omega[M]/N)*math.log(omega[M]) for M in omega)
cand = sorted(((eq - math.log(omega[M]))/eq, omega[M]) for M in omega)
print("  S_eq = %.4f   over %d macrostates" % (eq, len(omega)))
print("  smallest |eps|:")
for e, om in sorted(cand, key=lambda t: abs(t[0]))[:4]:
    print("     eps=%+.4f  at Omega=%5d  (lnOmega=%.4f)" % (e, om, math.log(om)))
print("  largest eps: %+.4f at Omega=%d" % (cand[-1][0], cand[-1][1]))
print("  eps of the actual seed: %+.4f" % ((eq - math.log(omega[label[seed.canonical()]]))/eq))
