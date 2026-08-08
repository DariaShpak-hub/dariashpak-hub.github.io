#!/usr/bin/env python3
"""
area_law_feasibility.py

Can S ~ area be distinguished from S ~ volume at enumerable sizes?

The question is well posed and dimensionless: for a region R of the
hypergraph, does the multiplicity of interior configurations scale with the
boundary |dR| or with the bulk |R|? That is the Bekenstein-Hawking
dichotomy stated combinatorially.

It is not measurable here, and the reason is structural rather than a matter
of more compute. On the sectors that can be enumerated, corr(r, mean|dR|) is
exactly 0.0000: the boundary size is symmetric under exchanging a region with
its complement, so it rises to a peak at r = n/2 and falls back, while the
interior grows monotonically (corr +0.96). Every region at these sizes is
essentially all boundary -- there is no bulk to separate from it.

Separating r^2 from r^3 needs r spanning at least a decade, so n of order 200
nodes, and the number of states in a fixed (nodes, weight) sector grows
super-exponentially in n. Exact enumeration is out of reach by an enormous
margin. Monte Carlo sampling within a sector, using the Metropolis kernel
already built in rewrite_thermo, is the only viable route.

Two caveats worth recording before anyone attempts it. Generic combinatorial
counting with free interior degrees of freedom gives a volume law, so that is
the expected answer. And an area law, if found, would place the system with
gapped local Hamiltonians in condensed matter, which have entanglement area
laws without any gravity -- it would not be Bekenstein-Hawking, which is a
claim about the coefficient A/4l_P^2 and needs a length this framework
cannot supply.
"""
import sys, math, random
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
sys.path.insert(0, "/tmp/claude-0/-home-user-dariashpak-hub-github-io/31ff9562-0d3c-5826-99eb-bd71006a2e6c/scratchpad")
from collections import defaultdict
from itertools import combinations
from closed import enumerate_fixed
from entropic_order import pearson

random.seed(7)

print("region size r vs boundary size |dR|, over sectors we can enumerate")
print("%4s %4s %7s %6s %10s %10s %12s" % (
    "n","w","states","r","mean|R_int|","mean|dR|","corr(r,|dR|)"))

for (n, w) in [(6,6),(7,6)]:
    st = enumerate_fixed(n, w)
    if not st: continue
    rs, bs, ints = [], [], []
    for r in range(1, n):
        tot_b = tot_i = cnt = 0
        for H in list(st.values())[:400]:
            for R in combinations(range(n), r):
                Rs = set(R)
                interior = sum(1 for e in H.edges if e <= Rs)
                cross    = sum(1 for e in H.edges if (e & Rs) and not (e <= Rs))
                tot_i += interior; tot_b += cross; cnt += 1
        if not cnt: continue
        mi, mb = tot_i/cnt, tot_b/cnt
        print("%4d %4d %7d %6d %10.3f %10.3f" % (n, w, len(st), r, mi, mb))
        rs.append(r); bs.append(mb); ints.append(mi)
    print("     corr(r, mean|dR|) = %+.4f      corr(r, mean interior) = %+.4f"
          % (pearson(rs, bs), pearson(rs, ints)))
    print()

print("scale needed to separate r^2 from r^3:")
for dec in (1, 1.5, 2):
    r_max = 10**dec
    print("  r spanning %5.0f  =>  n ~ %5.0f nodes; sector size grows"
          " super-exponentially in n" % (r_max, r_max*2))
