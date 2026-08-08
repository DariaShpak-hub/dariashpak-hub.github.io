#!/usr/bin/env python3
"""
effective_cosmology.py

An effective expansion history read off the self-triggered cascade.

Identification: a^3 ~ exp(S_cg), so ln a = S/3, H = d ln a / dt,
q = -1 - Hdot/H^2, and w_eff = -1 - (2/3) Hdot/H^2. Acceleration is q < 0,
equivalently w_eff < -1/3; w_eff < -1 would be phantom, a stricter condition.

Result: the wrong cosmology, in a specific and falsifiable way.

    early    brief acceleration -- q = -0.283 at t=5, min w_eff = -0.569
    then     monotonic deceleration, Hdot < 0 everywhere while H > 0
    t=3697   H crosses zero
    after    contraction

So this predicts acceleration at *early* times and deceleration then
recollapse at late times. The observed universe is the reverse: decelerating
early, accelerating now. Under this identification the model is falsified by
supernovae and BAO, which is the outcome the programme was built to produce.

Three caveats bound how much that falsification means.

The identification a^3 ~ exp(S) is a choice, not a derivation. So is using
rewrite-step as the clock, and w_eff is not reparameterization invariant --
under t -> tau(t) it shifts by (2/3) tau''/(tau' H). The whole expansion
history therefore inherits the clock arbitrariness established earlier, so
what is falsified is this reading, not the framework.

And w_eff diverges near t=3697 (reaching 2.9e6) purely because H -> 0 in the
denominator. That is a coordinate artifact of the H=0 crossing, not physics.
"""
import sys, numpy as np
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from cascade import build, IDXA, BIN_MASK, W, N, shannon, IDX

src, dst, od, live = build((IDXA & BIN_MASK) == 0, range(10, 20))
seed = 0
for e in ((0,1),(1,2),(2,3),(3,4)): seed |= 1 << IDX[frozenset(e)]
p = np.zeros(N); p[seed] = 1.0

S = [shannon(p)]
for t in range(6000):
    w = np.zeros(N); np.divide(p, od, out=w, where=live)
    nxt = np.bincount(dst, weights=0.5*w[src], minlength=N)
    nxt += 0.5*p; nxt[~live] += 0.5*p[~live]; p = nxt
    S.append(shannon(p))
S = np.array(S)

lna = S / 3.0
H = np.gradient(lna)                       # d ln a / dt
Hdot = np.gradient(H)
with np.errstate(divide="ignore", invalid="ignore"):
    q = -1.0 - Hdot / H**2                 # deceleration parameter
    weff = -1.0 - (2.0/3.0) * Hdot / H**2

print("%7s %11s %11s %11s %11s %9s" % ("t","S","H","Hdot","w_eff","q"))
for t in [1,5,20,100,500,1000,2000,3000,3500,3690,3700,4000,5000,5999]:
    print("%7d %11.5f %11.3e %11.3e %11.4f %9.4f" % (
        t, S[t], H[t], Hdot[t], weff[t], q[t]))

pos = H > 1e-12
print()
print("H > 0 for t in [%d, %d]; H changes sign at t = %d"
      % (1, int(np.flatnonzero(pos)[-1]), int(np.flatnonzero(pos)[-1]) + 1))
seg = slice(5, int(np.flatnonzero(pos)[-1]) - 5)
print("while H > 0:  Hdot < 0 everywhere: %s" % bool((Hdot[seg] < 0).all()))
print("              min w_eff = %.4f   max w_eff = %.4f"
      % (weff[seg].min(), weff[seg].max()))
print("              any w_eff < -1 (acceleration)? %s"
      % bool((weff[seg] < -1).any()))
print("              any q < 0 (accelerating)?      %s"
      % bool((q[seg] < 0).any()))
