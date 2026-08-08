#!/usr/bin/env python3
"""
effective_pressure.py

Does the released state have negative effective pressure?

With no energy function, the free energy of a purely entropic system is
F = -TS, so

    p_eff = -(dF/dV) = T (dS/dV)

and the sign of the pressure is the sign of dS/dV. Negative pressure requires
entropy to *decrease* as volume grows -- the rubber-band condition, not the
ideal-gas one.

This needs no trajectory and no clock. Take the volume to be the weight
w = b + 2t, which is read off the edge structure rather than declared, and
compute the equilibrium equation of state directly by counting:

    Omega_n(w) = sum over b + 2t = w of  C(C(n,2), b) * C(C(n,3), t)
    S_eq(w)    = ln Omega_n(w)

with C(n,2) binary slots and C(n,3) ternary slots on n nodes. Then read off
where dS/dw changes sign.

The point of scanning n is that a turnover in Omega(w) is what produces
negative pressure, and a turnover exists only because the slots run out. If
the crossover weight grows with system size, the negative-pressure region is
a saturation effect that recedes as the system grows, rather than a property
of the dynamics.

Result: negative pressure is real, and it is half-filling
---------------------------------------------------------
dS/dw is positive below w = 15 and negative above it on five nodes, and the
self-triggered cascade reached <w> = 16.24, so the released state does sit in
the negative-pressure regime with no external agent and no inserted Lambda.

The mechanism is entirely combinatorial. Omega(w) is a product of binomials
and therefore symmetric about half-filling, so the sign change sits exactly at
w_max/2 at every size:

    n            4    5    6    7    8    9   10
    w_max       14   30   55   91  140  204  285
    p<0 from     7   15   28   46   70  102  143

Negative pressure is just the far side of that symmetric peak: entropy falls
past the midpoint because the slots are running out, the same reason a nearly
full lattice gas has fewer configurations than a half full one.

That disqualifies it as a dark-energy analogue. Since w_max grows as roughly
n^3/3, the crossover recedes as the system grows, so at fixed density in the
large-n limit the system is always below half filling and always at positive
pressure. Negative pressure here requires the universe to be running out of
room, which is the opposite of the regime where dark energy dominates.

It does unify two results previously reported as separate. The entropy peak
and the pressure sign change are the same event -- both are Omega(w) turning
over at w_max/2 -- so the entropy pulse and the negative pressure are one
phenomenon, not two.

And p < 0 was never sufficient in any case: acceleration needs
p < -rho c^2 / 3, and rho requires an energy density this framework does not
have.
"""

import math
from math import comb, log


def omega(n, w):
    nb, nt = comb(n, 2), comb(n, 3)
    tot = 0
    for t in range(0, min(nt, w // 2) + 1):
        b = w - 2 * t
        if 0 <= b <= nb:
            tot += comb(nb, b) * comb(nt, t)
    return tot


def equation_of_state(n):
    nb, nt = comb(n, 2), comb(n, 3)
    wmax = nb + 2 * nt
    S = {}
    for w in range(wmax + 1):
        o = omega(n, w)
        if o > 0:
            S[w] = log(o)
    ws = sorted(S)
    # central difference for dS/dw
    dS = {}
    for i in range(1, len(ws) - 1):
        dS[ws[i]] = (S[ws[i + 1]] - S[ws[i - 1]]) / (ws[i + 1] - ws[i - 1])
    cross = None
    for i in range(1, len(ws) - 2):
        a, b = ws[i], ws[i + 1]
        if a in dS and b in dS and dS[a] > 0 >= dS[b]:
            cross = (a, b)
            break
    return S, dS, cross, wmax


def main():
    print("equation of state  S_eq(w) = ln Omega(w),  p_eff ~ dS/dw")
    print()
    print("n=5 in detail (the system used in the cascade runs)")
    S, dS, cross, wmax = equation_of_state(5)
    print("%5s %14s %11s %10s" % ("w", "Omega(w)", "S_eq", "dS/dw"))
    for w in sorted(S):
        if w % 2 == 0 or w in (13, 15, 17):
            mark = ""
            if w in dS:
                mark = "  <-- p>0" if dS[w] > 0 else "  <-- p<0"
            print("%5d %14d %11.5f %10.5f%s" % (
                w, omega(5, w), S[w], dS.get(w, float('nan')), mark))
    print("  sign change between w =", cross, " (max w =", wmax, ")")
    print()

    print("does the negative-pressure region recede as the system grows?")
    print("%4s %9s %9s %14s %14s" % (
        "n", "bin slots", "ter slots", "w_max", "p<0 begins at w"))
    for n in range(4, 11):
        S, dS, cross, wmax = equation_of_state(n)
        print("%4d %9d %9d %14d %14s" % (
            n, comb(n, 2), comb(n, 3), wmax,
            str(cross[1]) if cross else "none"))
    print()
    print("cascade reached <w> = 16.24 at the last measured step,")
    print("against the n=5 sign change above.")


if __name__ == "__main__":
    main()
