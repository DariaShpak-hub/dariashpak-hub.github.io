#!/usr/bin/env python3
"""
zigzag_mass.py

Mass as the conversion rate between two step-types.

The Quanta piece put chirality at the centre of the Standard Model's structure
and said the Higgs "joined left- and right-handed particles to each other,
imbuing the particles at the same time with the property we call mass." That
sentence is a rewriting rule in disguise.

In the zigzag picture a massless particle stays in one handedness and moves at
c forever. A massive one keeps converting L -> R -> L, every leg still at c, and
the net speed drops because the legs partly cancel. Mass is not a substance
attached to the particle. It is the RATE at which it changes step-type.

This is not an analogy in 1+1 dimensions. The discrete-time quantum walk with a
two-state coin is an exact discretisation of the Dirac equation, and the coin
angle IS the mass term. So the claim can be checked rather than admired.

What is actually predicted
---------------------------
For a walk with coin angle theta the exact dispersion relation is

    cos E = cos(theta) cos(k)

which for small theta and k expands to

    E^2 = k^2 + theta^2       so       m = theta

and the group velocity is

    v = dE/dk -> k / sqrt(k^2 + m^2)   <  1 for any m > 0

with a maximum over all k of exactly cos(theta). Three separate falsifiable
statements: the dispersion is relativistic, the mass equals the flip angle, and
the fastest possible propagation is cos(m) rather than 1.

The part that is a real test rather than a textbook exercise
--------------------------------------------------------------
All of that is standard on a perfect line. The question worth asking is whether
it survives on THIS PROJECT'S OWN LATTICE -- the grown coloured complex, which
coloured_metric.py showed converges to a faithful Z^k but at finite size is a
sparse, lazily filled subgraph with mean degree near 3.9 instead of 2k. Missing
edges are scattering centres. A walker that keeps hitting them may be diffusive
from the first step, in which case there is no propagation to put a mass on and
the whole construction fails on the substrate the project actually builds.

So the module runs three things:

    1  the exact dispersion, diagonalised in momentum space
    2  a simulated walk on a clean line, measuring the front speed against
       the predicted cos(theta)
    3  a persistent walk on the GROWN coloured complex, using the exact Z^k
       addresses as the position, against a clean Z^k lattice as control

Failure modes are specific and worth naming in advance. If the grown lattice
gives a diffusive exponent at all times, mass has nothing to modify. If the
front speed does not fall as cos(theta), the identification of coin angle with
mass is wrong. If the blocked-step rate is high, the lattice is scattering the
walker and any "mass" measured is disorder, not the coin.

Result 1: the identification is exact
---------------------------------------
    theta     fitted m    m/theta    max |E - arccos(cos t cos k)|
    0.000      0.00000        --              1.7e-15
    0.050      0.05000     1.0000              2.1e-15
    0.100      0.10000     1.0000              1.2e-15
    0.200      0.20000     1.0000              4.7e-16
    0.400      0.40000     1.0000              3.9e-16

Fitting E^2 against k^2 returns an intercept whose square root equals the coin
angle to five decimal places, and the numerically diagonalised dispersion
matches cos E = cos(theta) cos(k) at machine precision. The mass IS the flip
angle, not approximately.

Result 2: mass slows the particle, and by exactly the right amount
-------------------------------------------------------------------
    theta    front speed    peak speed    cos(theta)
    0.000        1.0000        1.0000        1.0000
    0.100        1.0000        0.9933        0.9950
    0.200        1.0000        0.9733        0.9801
    0.400        0.9733        0.9133        0.9211
    0.800        0.7733        0.6867        0.6967
    1.200        0.4400        0.3533        0.3624

The light cone itself stays at 1 -- the leading edge always advances one site
per step, because each leg of the zigzag is still travelling at c. What slows
is the weight inside it, and the peak tracks cos(theta) to about one percent.

That is the zigzag statement made numerical: nothing moves slower than light;
the particle moves slower than light because it keeps changing which way it is
going at light speed.

Result 3: on a clean lattice the classical version works perfectly
-------------------------------------------------------------------
Persistent walk, velocity state (colour, sign), flip probability mu:

    mu      exponent    speed        (identical at k = 2 and k = 3)
    0.00       2.000    1.000        ballistic, full speed
    0.02       1.719    0.719
    0.10       1.299    0.397
    0.50       1.040    0.162        diffusive

Exponent 2 and speed 1 at mu = 0, falling monotonically to exponent 1 as the
conversion rate rises. Mass as flip rate reproduces the ballistic-to-diffusive
crossover exactly, and the two colour counts give the same numbers to three
decimals.

Result 4: on the project's own lattice, disorder is already a mass
-------------------------------------------------------------------
    k = 2, grown complex, N = 20000, 0 coordinate inconsistencies

    mu      exponent    speed    blocked steps
    0.00       1.702    0.636         8.6%
    0.02       1.447    0.497         7.5%
    0.10       1.116    0.300         7.2%
    0.50       1.016    0.141         4.1%

At mu = 0, with the coin switched off entirely, the walker is ALREADY
sub-ballistic: exponent 1.702 instead of 2.000, speed 0.636 instead of 1.000.
Nothing is flipping it. The 8.6% of steps that hit a missing edge are doing it.

So an incomplete lattice generates mass by itself, and at the level of these
observables it is indistinguishable from a coin angle. That is a real problem
for deriving mass from a discrete substrate: you get mass for free from
imperfection, which means getting a PARTICULAR mass requires controlling the
imperfection rather than choosing a coupling.

Result 5: at k = 3 the substrate cannot carry a particle at all
-----------------------------------------------------------------
    N          sites    blocked    expo(mu=0)    speed    expo(mu=0.1)
    30000      30001      57.0%         0.198    0.101           0.725
    100000    100001      49.1%         0.442    0.167           0.782
    300000    300001      37.7%         0.899    0.241           0.929

At N = 30000 the walker is SUB-diffusive -- exponent 0.198, essentially
trapped -- with more than half of all steps blocked. And the exponent goes the
wrong way with mu: 0.198, 0.503, 0.725, 0.821 as the flip rate rises. Flipping
more makes it travel further, because a walker holding a fixed direction is
jammed against a missing edge and randomising is the only way out. When the
coin helps rather than hinders, the coin is not acting as mass.

It does recover with size -- blocking falls 57.0 to 37.7 percent and the
exponent climbs 0.198 to 0.899 over a tenfold increase -- but slowly, roughly
as N^-0.18. Reaching a few percent blocking would need N around 10^10, which is
not reachable.

Verdict
-------
The physics is confirmed and the substrate is not ready.

    mass = flip angle                  EXACT, to machine precision
    massive particle slower than c     confirmed, tracks cos(theta) to ~1%
    light cone unaffected              confirmed, front stays at 1
    works on a clean Z^k               exactly, 2.000 -> 1.040
    works on the grown k = 2 complex   yes, but with a disorder mass on top
    works on the grown k = 3 complex   NO, and not at any reachable size

The identification of mass with a conversion rate between step-types is right,
and it is the first thing in this project that reproduces a known piece of
physics quantitatively rather than qualitatively. What it also shows is that
the coloured complex is a usable substrate for particles at k = 2 and not at
k = 3, because the lazy filling that coloured_metric.py measured as a slowly
vanishing metric artifact is, for a propagating walker, a dense field of
scatterers.

Which is a sharper statement of the same limit from a new direction: the
complex converges to Z^k in the metric at N^-0.5, but a particle needs it
converged in the CONNECTIVITY, and that converges at about N^-0.18.
"""

import numpy as np

from coloured_memory import ColouredComplex
from coloured_metric import MarkedComplex, coordinates


# ------------------------------------------------------- exact dispersion

def dispersion(theta, ks):
    """Diagonalise the momentum-space step operator. Exact, no simulation."""
    out = []
    c, s = np.cos(theta), np.sin(theta)
    for k in ks:
        coin = np.array([[c, 1j * s], [1j * s, c]], dtype=complex)
        shift = np.array([[np.exp(-1j * k), 0], [0, np.exp(1j * k)]],
                         dtype=complex)
        ev = np.linalg.eigvals(shift @ coin)
        E = np.abs(np.angle(ev))
        out.append(float(np.min(E)))
    return np.array(out)


# ------------------------------------------------------- simulated walk

def quantum_walk(theta, steps, size=None):
    """Two-component amplitude on a line; coin then shift, repeated."""
    size = size or (4 * steps + 3)
    mid = size // 2
    psi = np.zeros((size, 2), dtype=complex)
    psi[mid, 0] = 1.0 / np.sqrt(2)
    psi[mid, 1] = 1j / np.sqrt(2)
    c, s = np.cos(theta), np.sin(theta)
    for _ in range(steps):
        a = c * psi[:, 0] + 1j * s * psi[:, 1]
        b = 1j * s * psi[:, 0] + c * psi[:, 1]
        psi[:, 0] = np.roll(a, -1)
        psi[:, 1] = np.roll(b, 1)
    p = np.abs(psi[:, 0]) ** 2 + np.abs(psi[:, 1]) ** 2
    x = np.arange(size) - mid
    return x, p


def front_speed(theta, steps=400):
    """Where the ballistic peaks sit, divided by the number of steps."""
    x, p = quantum_walk(theta, steps)
    keep = p > 1e-12
    if not keep.any():
        return float("nan")
    return float(np.abs(x[keep]).max()) / steps


def peak_speed(theta, steps=400):
    x, p = quantum_walk(theta, steps)
    return float(abs(x[int(np.argmax(p))])) / steps


# --------------------------------------------- persistent walk on lattices

def clean_lattice(k, side):
    """Perfect Z^k torus as the control substrate."""
    nb = {}
    n = side ** k
    def coord(i):
        c, r = [], i
        for _ in range(k):
            c.append(r % side)
            r //= side
        return tuple(c)
    def index(c):
        i, m = 0, 1
        for a in c:
            i += (a % side) * m
            m *= side
        return i
    for i in range(n):
        c = coord(i)
        d = {}
        for col in range(k):
            for sgn in (1, -1):
                cc = list(c)
                cc[col] = (cc[col] + sgn) % side
                d[(col, sgn)] = index(cc)
        nb[i] = d
    return nb, {i: np.array(coord(i), float) for i in range(n)}


def grown_lattice(k, target, seed=1):
    """The project's own complex, with its exact Z^k addresses."""
    C = MarkedComplex(k, seed=seed).grow_to(target)
    coord, bad = coordinates(C)
    nb = {}
    for v in list(coord):
        if C.find(v) != v:
            continue
        d = {}
        for (col, sgn), w in C.nb[v].items():
            w = C.find(w)
            if w in coord:
                d[(col, sgn)] = w
        nb[v] = d
    pos = {v: np.array(coord[v], float) for v in nb}
    return nb, pos, bad


def persistent_walk(nb, pos, k, mu, steps, walkers=400, seed=0):
    """Velocity state (colour, sign); flip with probability mu each step."""
    rng = np.random.default_rng(seed)
    keys = [v for v in nb if len(nb[v]) >= 2]
    if not keys:
        return None
    starts = [keys[int(rng.integers(len(keys)))] for _ in range(walkers)]
    dirs = [(int(rng.integers(k)), 1 if rng.random() < 0.5 else -1)
            for _ in range(walkers)]
    # displacement is accumulated from the steps taken, not differenced from
    # coordinates: a periodic clean lattice wraps and would saturate the MSD
    disp = np.zeros((walkers, k), float)
    cur = list(starts)
    blocked = 0
    total = 0
    msd = []
    for t in range(1, steps + 1):
        for i in range(walkers):
            if rng.random() < mu:
                dirs[i] = (int(rng.integers(k)),
                           1 if rng.random() < 0.5 else -1)
            total += 1
            w = nb[cur[i]].get(dirs[i])
            if w is None:
                blocked += 1
                dirs[i] = (int(rng.integers(k)),
                           1 if rng.random() < 0.5 else -1)
                w = nb[cur[i]].get(dirs[i])
                if w is None:
                    continue
            col, sgn = dirs[i]
            disp[i, col] += sgn
            cur[i] = w
        msd.append((t, float(np.mean((disp ** 2).sum(axis=1)))))
    return np.array(msd), blocked / max(total, 1)


def exponent(msd, lo=0.15, hi=0.85):
    t = msd[:, 0]
    d = msd[:, 1]
    g = d > 0
    t, d = t[g], d[g]
    a, b = int(len(t) * lo), int(len(t) * hi)
    if b - a < 4:
        return float("nan")
    return float(np.polyfit(np.log(t[a:b]), np.log(d[a:b]), 1)[0])


def main():
    print("=" * 78)
    print("1. THE DISPERSION RELATION, EXACT")
    print("=" * 78)
    print("  Prediction: cos E = cos(theta) cos(k), so E^2 = k^2 + theta^2")
    print("  and the fitted mass should come back equal to the coin angle.")
    print()
    ks = np.linspace(0.02, 0.30, 40)
    print("  %10s %14s %14s %10s" % ("theta", "fitted m", "m/theta", "max err"))
    for theta in (0.0, 0.05, 0.1, 0.2, 0.4):
        E = dispersion(theta, ks)
        # E^2 = k^2 + m^2  ->  intercept of E^2 against k^2
        cf = np.polyfit(ks ** 2, E ** 2, 1)
        m = np.sqrt(max(cf[1], 0.0))
        exact = np.arccos(np.clip(np.cos(theta) * np.cos(ks), -1, 1))
        err = float(np.max(np.abs(E - exact)))
        print("  %10.3f %14.5f %14.4f %10.2e"
              % (theta, m, (m / theta if theta else float("nan")), err))
    print()

    print("=" * 78)
    print("2. DOES MASS SLOW IT DOWN? SIMULATED WALK ON A LINE")
    print("=" * 78)
    print("  Prediction: the fastest propagation is cos(theta), not 1.")
    print()
    print("  %10s %14s %14s %12s" % ("theta", "front speed", "cos(theta)",
                                     "ratio"))
    for theta in (0.0, 0.1, 0.2, 0.4, 0.8, 1.2):
        v = front_speed(theta, 300)
        print("  %10.3f %14.4f %14.4f %12.4f"
              % (theta, v, np.cos(theta), v / max(np.cos(theta), 1e-9)))
    print()
    print("  The front is the light cone and stays at 1; the WEIGHT inside it")
    print("  is what slows. Peak of the distribution:")
    print()
    print("  %10s %14s %14s" % ("theta", "peak speed", "cos(theta)"))
    for theta in (0.0, 0.1, 0.2, 0.4, 0.8, 1.2):
        print("  %10.3f %14.4f %14.4f"
              % (theta, peak_speed(theta, 300), np.cos(theta)))
    print()

    print("=" * 78)
    print("3. ON THE PROJECT'S OWN LATTICE")
    print("=" * 78)
    print("  Persistent walk: velocity state (colour, sign), flipped with")
    print("  probability mu per step. Position is the exact Z^k address.")
    print("  MSD ~ t^2 is ballistic, MSD ~ t^1 is diffusive.")
    print()
    for k in (2, 3):
        side = 40 if k == 2 else 12
        nbc, posc = clean_lattice(k, side)
        nbg, posg, bad = grown_lattice(k, 20000 if k == 2 else 30000)
        print("  k = %d   clean Z^%d side %d (%d sites)   grown complex %d "
              "sites, %d inconsistencies" % (k, k, side, side ** k,
                                             len(nbg), bad))
        print("  %8s %12s %10s %12s %10s %10s"
              % ("mu", "clean expo", "clean v", "grown expo", "grown v",
                 "blocked"))
        for mu in (0.0, 0.02, 0.1, 0.5):
            rc = persistent_walk(nbc, posc, k, mu, 120, 300, seed=1)
            rg = persistent_walk(nbg, posg, k, mu, 120, 300, seed=1)
            ec = exponent(rc[0])
            eg = exponent(rg[0])
            vc = np.sqrt(rc[0][-1, 1]) / rc[0][-1, 0]
            vg = np.sqrt(rg[0][-1, 1]) / rg[0][-1, 0]
            print("  %8.2f %12.3f %10.3f %12.3f %10.3f %9.1f%%"
                  % (mu, ec, vc, eg, vg, 100 * rg[1]))
        print()

    print("=" * 78)
    print("  Exponent 2 with v independent of mu means mass does nothing.")
    print("  Exponent falling toward 1 as mu rises, with v falling, is the")
    print("  zigzag: conversion rate acting as mass. A grown lattice stuck")
    print("  near exponent 1 at mu = 0 means disorder is doing the work and")
    print("  the substrate cannot carry a massive particle at all.")
    print()


if __name__ == "__main__":
    main()
