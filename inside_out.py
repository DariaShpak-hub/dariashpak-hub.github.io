#!/usr/bin/env python3
"""
inside_out.py

Does sigma map the concentration regime onto the expansion regime?

sigma is edge complementation, sigma(e) = V \\ e, which on five nodes sends
each binary edge to a ternary one and back. A state with b binary and t
ternary edges therefore maps to one with (t, b): sigma reflects the (b, t)
plane across the diagonal b = t.

The two regimes sit on the two sides of that diagonal. Concentration is
motion toward low b and high t -- fewer, larger bindings, the direction
compose drives. Expansion is motion toward high b and low t, the decompose
direction. So the question becomes whether the reflection is a symmetry of

    the multiplicity  Omega(b, t)     -- the entropy landscape
    the rule set                      -- the dynamics

and these can differ. Three things are checked:

    1  is Omega(b,t) = Omega(t,b) exactly?
    2  does the rule set commute with the reflection?
    3  where does the cascade's terminal region map to?

If (1) holds and (2) fails, the state space is symmetric between concentrated
and dispersed while the rules prefer one direction -- which would say the
asymmetry we observe is imposed by the dynamics, not by the space of states.

Result: yes for the state space, no for the dynamics
----------------------------------------------------
The multiplicity is exactly symmetric: max |Omega(b,t) - Omega(t,b)| = 0 over
all 121 cells, since Omega(b,t) = C(10,b) C(10,t). The fully concentrated
state (0,10) and the fully dispersed state (10,0) have identical entropy,
zero, both being unique. The diagonal b = t is simultaneously the fixed set of
the reflection and the entropy maximum, at (5,5) with S_B = 11.059.

The rule set breaks the reflection completely: 0 of 30 compose moves have
their reflection in the set. The obstruction is an arity mismatch rather than
an accident. compose is 2 binary -> 1 ternary, so its reflection is
2 ternary -> 1 binary, while decompose is 1 ternary -> 2 binary. The rule set
is self-inverse but not self-mirror, and no relabelling can identify a 2->1
rule with a 1->2 one.

So the asymmetry between concentration and dispersion is imposed by the
dynamics, not by the space of states: the entropy landscape is perfectly
balanced and the rules pick a side.

Adding "2 ternary -> 1 binary" and "1 binary -> 2 ternary" restores
equivariance -- that is the symmetrisation used in sector_breaking.py, where
it gave 0 equivariance violations and A_M = 0 exactly. Turning the system
inside out is therefore a genuine involution of the model; the bundled rules
simply fail to respect it, and that one arity mismatch is where the whole
directional preference lives.
"""

import numpy as np
from math import comb, log
from itertools import combinations

N_NODES = 5
EDGES = [frozenset(c) for c in combinations(range(N_NODES), 2)] + \
        [frozenset(c) for c in combinations(range(N_NODES), 3)]
IDX = {e: i for i, e in enumerate(EDGES)}
NB = sum(1 for e in EDGES if len(e) == 2)
NT = len(EDGES) - NB
V = frozenset(range(N_NODES))
SIG = [IDX[frozenset(V - e)] for e in EDGES]

COMPOSE = [(i, j, IDX[EDGES[i] | EDGES[j]])
           for i in range(NB) for j in range(i + 1, NB)
           if len(EDGES[i] & EDGES[j]) == 1]
DECOMPOSE = [(k, i, j) for (i, j, k) in COMPOSE]


def omega(b, t):
    return comb(NB, b) * comb(NT, t)


def check_measure():
    print("1. is the multiplicity symmetric under the reflection?")
    worst = 0
    for b in range(NB + 1):
        for t in range(NT + 1):
            worst = max(worst, abs(omega(b, t) - omega(t, b)))
    print("   max |Omega(b,t) - Omega(t,b)| over all %d cells : %d"
          % ((NB + 1) * (NT + 1), worst))
    print("   -> entropy landscape is %s"
          % ("EXACTLY symmetric" if worst == 0 else "not symmetric"))
    print()
    print("   S_B along the reflection axis and off it:")
    print("   %6s %6s %12s %12s %10s" % ("b", "t", "S_B(b,t)", "S_B(t,b)", "diff"))
    for b, t in ((0, 10), (2, 8), (4, 6), (5, 5), (6, 4), (8, 2), (10, 0)):
        print("   %6d %6d %12.5f %12.5f %10.2e"
              % (b, t, log(omega(b, t)), log(omega(t, b)),
                 abs(log(omega(b, t)) - log(omega(t, b)))))
    print()


def check_dynamics():
    print("2. does the rule set commute with the reflection?")
    # a move is a triple (i,j,k): two binaries consumed, one ternary made.
    # its sigma image consumes two ternaries and makes one binary, which is
    # a differently shaped rule. Count how many images land back in the set.
    comp = {(min(i, j), max(i, j), k) for i, j, k in COMPOSE}
    dec = {(k, min(i, j), max(i, j)) for k, i, j in DECOMPOSE}
    img_comp = {(min(SIG[i], SIG[j]), max(SIG[i], SIG[j]), SIG[k])
                for i, j, k in COMPOSE}
    both = comp | dec
    inside = sum(1 for m in img_comp if m in both)
    print("   compose moves                       : %d" % len(comp))
    print("   their reflections landing in the set: %d" % inside)
    print("   -> rule set is %s under the reflection"
          % ("equivariant" if inside == len(img_comp) else "NOT equivariant"))
    print()
    print("   reason: compose is 2 binary -> 1 ternary, so its reflection is")
    print("   2 ternary -> 1 binary, while decompose is 1 ternary -> 2 binary.")
    print("   The arities do not match, so no relabelling can identify them.")
    print()


def check_regimes():
    print("3. where do the regimes map?")
    print("   %-28s %-14s %-14s" % ("state", "(b,t)", "reflection"))
    cases = [("cascade terminal (all ternary)", (0, 10)),
             ("cascade seed (path, 4 binary)", (4, 0)),
             ("entropy ridge", (5, 5)),
             ("all binary", (10, 0))]
    for name, (b, t) in cases:
        print("   %-28s %-14s %-14s  S_B %.4f -> %.4f"
              % (name, "(%d,%d)" % (b, t), "(%d,%d)" % (t, b),
                 log(omega(b, t)), log(omega(t, b))))
    print()
    # where is the entropy maximum, and is the diagonal it?
    best = max(((b, t) for b in range(NB + 1) for t in range(NT + 1)),
               key=lambda x: omega(*x))
    print("   entropy maximum at (b,t) = %s, on the diagonal: %s"
          % (str(best), best[0] == best[1]))
    print("   the diagonal b = t is the fixed set of the reflection,")
    print("   so the two regimes are its two sides, exchanged exactly.")
    print()


if __name__ == "__main__":
    check_measure()
    check_dynamics()
    check_regimes()
