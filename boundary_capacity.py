#!/usr/bin/env python3
"""
boundary_capacity.py

How much of the total entropy can a boundary hold, in principle?

The entropy chain rule is exact:

    S_total = S(boundary) + S(bulk | boundary)

so the two properties usually wanted together are mutually exclusive. If the
boundary determines the bulk, S(bulk|boundary) -> 0 and therefore
S_total <= |boundary| ln 2: the entire state space is boundary-sized. If
instead the state count scales with the bulk, then
S(bulk|boundary) >= S_total - |boundary| ln 2, which is itself proportional
to bulk. "Simple exterior plus a hidden state space large compared to the
boundary plus strong correlations" is not a hard target, it is inconsistent.

Measured below, with the chain rule verified to ~1e-11. The decisive feature
is that boundary capacity and bulk information move in opposite directions:
crossing-edge count peaks at r = n/2 and falls away, so S(boundary) shrinks
exactly where S(bulk|boundary) grows, and the ratio crosses 1. There is no
region size at which the boundary could carry the interior.

This is why the result is not a defect of one rule set. Classical
configuration entropy is extensive by construction -- a product of local
degrees of freedom -- whereas area-law entropy in physics is a property of
entanglement in quantum states. Black hole entropy is striking precisely
because it violates the extensivity a classical bulk would give.

A construction with Omega ~ exp(|boundary|) is possible in principle: a
topological or gauge system whose bulk degrees of freedom are pure gauge.
But then there is no hidden interior to encode, which is the honest
resolution rather than a workaround.
"""
import sys, math
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from collections import defaultdict
from itertools import combinations
from boundary_entropy import enumerate_labeled, connectivity, entropy

def run(n, w):
    print("n=%d w=%d" % (n, w))
    print("%3s %8s %9s %9s %9s %9s %9s" % (
        "r","bdry slots","S(Y)","S(X|Y)","S(X,Y)","chain err","S(X|Y)/S(Y)"))
    for r in range(2, n):
        R = frozenset(range(r))
        joint = defaultdict(lambda: defaultdict(int))
        for edges in enumerate_labeled(n, w):
            X = frozenset(e for e in edges if e <= R)
            Y = frozenset(e for e in edges if (e & R) and not (e <= R))
            joint[Y][X] += 1
        tot = sum(sum(d.values()) for d in joint.values())
        SY = entropy([sum(d.values()) for d in joint.values()])
        SXY = entropy([k for d in joint.values() for k in d.values()])
        SX_Y = sum((sum(d.values())/tot)*entropy(list(d.values()))
                   for d in joint.values())
        # boundary "slots": possible crossing edges
        bslots = sum(1 for s in (2,3) for c in combinations(range(n), s)
                     if (set(c) & R) and not set(c) <= R)
        print("%3d %8d %9.4f %9.4f %9.4f %9.2e %9.3f" % (
            r, bslots, SY, SX_Y, SXY, abs(SXY-(SY+SX_Y)),
            SX_Y/SY if SY > 0 else float('nan')))
    print()

run(6, 5)
run(7, 5)
run(7, 6)
