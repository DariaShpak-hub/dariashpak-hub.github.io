"""Charges on a window CLOSED under compose+decompose: fixed nodes and weight."""
import sys
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from collections import defaultdict
from itertools import combinations
from rewrite_lab import Hypergraph
from rewrite_research import COMPOSE, DECOMPOSE
from charges import components

RULES = [COMPOSE, DECOMPOSE]

def enumerate_fixed(n, w):
    """All canonical hypergraphs covering exactly n nodes with weight w."""
    binary = [frozenset(c) for c in combinations(range(n), 2)]
    ternary = [frozenset(c) for c in combinations(range(n), 3)]
    out = {}
    for t in range(w // 2 + 1):
        b = w - 2 * t
        if b < 0:
            continue
        for tt in combinations(ternary, t):
            for bb in combinations(binary, b):
                edges = set(tt) | set(bb)
                if len(edges) != t + b:
                    continue
                H = Hypergraph(edges)
                if len(H.nodes()) != n:
                    continue
                c = H.canonical()
                if c not in out:
                    out[c] = H
    return out

print("window closed under compose+decompose (fixed nodes n, weight w)")
print("%4s %4s %8s %12s %s" % ("n", "w", "states", "components", "verdict"))
for n in range(3, 7):
    for w in range(2, 7):
        st = enumerate_fixed(n, w)
        if len(st) < 2:
            continue
        comps = components(st, RULES)
        verdict = ("single sector" if len(comps) == 1
                   else "%d sectors -> hidden charge" % len(comps))
        print("%4d %4d %8d %12d %s" % (n, w, len(st), len(comps), verdict))
