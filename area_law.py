#!/usr/bin/env python3
"""
area_law.py

Does boundary entropy scale like a surface or like the bulk?

The proposal under test is to drop spatial dimension entirely and define an
effective area from an emergent boundary entropy, A_eff ~ S_boundary, giving
a curvature scale Lambda_eff ~ F(Delta S_B) / A_eff.

Two findings, one dimensional and one structural.

Dimensionally it does not help. S_boundary is a pure number, so A_eff =
l^2 * S_boundary and the proportionality constant carries the [L]^2. This is
exactly Bekenstein-Hawking, where S = A / 4 l_P^2 and the Planck length is
what converts area into entropy. Without an independent l there is no area,
only a dimensionless count.

Structurally, an area law does not avoid dimension -- it encodes it. Writing
the surface exponent as

    s = dln|dB_k| / dln|B_k|,

a volume law is s = 1 (boundary proportional to bulk) while an area law in
dimension d is s = (d-1)/d, which inverts to d = 1/(1-s). The area relation
and the dimension are the same data, so one cannot replace the other.

Measured: the polynomially growing grammars do obey an area law, s = 0.69 to
0.84 rising with seed valence, while the exponentially growing grammar gives
s = 1.04, a volume law with d = infinity. The grammar whose accessible volume
expands exponentially is therefore precisely the one with no area law. Note
the two dimension estimates (from volume growth, and from 1/(1-s)) disagree
at these depths, so the scaling is not yet asymptotic.
"""
import sys, math
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
from expansion import ball_volumes, SUBDIVIDE, SMOOTH, star
from rewrite_research import COMPOSE, DECOMPOSE, GRAMMARS
from rewrite_lab import path_seed


def slope(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    return (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
            if sxx > 0 else float("nan"))


def analyse(name, seed, rules, depth, pred_d=None):
    v = ball_volumes(seed, rules, depth, max_states=300000)
    surf = [v[i] - v[i - 1] for i in range(1, len(v))]
    # drop saturated tail (surface 0) and the first couple of points
    pts = [(v[i], surf[i - 1]) for i in range(1, len(v)) if surf[i - 1] > 0]
    pts = pts[max(1, len(pts) // 3):]
    if len(pts) < 3:
        print("%-34s insufficient data" % name)
        return
    s = slope([math.log(b) for b, _ in pts], [math.log(a) for _, a in pts])
    ks = list(range(1, len(v)))
    d_fit = slope([math.log(k) for k in ks[len(ks) // 3:]],
                  [math.log(v[k]) for k in ks[len(ks) // 3:]])
    # area law with dimension d predicts s = (d-1)/d, i.e. d = 1/(1-s)
    d_from_s = 1.0 / (1.0 - s) if s < 0.999 else float("inf")
    print("%-34s s=%6.3f  d_fit=%6.3f  d_from_s=%7.3f  %s"
          % (name, s, d_fit, d_from_s,
             ("pred d=%d" % pred_d) if pred_d else ""))


print("s = dln|dB| / dln|B|.   s=1 volume law (boundary ~ bulk).")
print("s=(d-1)/d area law with dimension d, so d = 1/(1-s).")
print()
for m in (3, 4, 5, 6):
    analyse("star(%d) subdivide+smooth" % m, star(m),
            [SUBDIVIDE, SMOOTH], 20, pred_d=m - 1)
print()
analyse("path(2) compose+decompose+subdivide", path_seed(2),
        [COMPOSE, DECOMPOSE, SUBDIVIDE], 12)
analyse("path(6) G1 (finite, saturates)", path_seed(6),
        GRAMMARS["G1-reversible"]["rules"], 10)
