#!/usr/bin/env python3
"""
finite_size_scaling.py

Is there a real phase transition in a global geometric observable, or only a
crossover?

This is the experiment that lets the whole programme say no. A backdraft needs
a metastable basin, a barrier, and a loss of metastability at a critical point.
If instead every quantity varies smoothly, the characteristic scale stays
finite, and the susceptibility peak does not grow with N, then there is no
transition and no backdraft mechanism -- whatever else the model does.

Setup. The energy stays LOCAL (H = -(4-cycle count), as in geometric_phase.py)
while the order parameter is GLOBAL, which is the split global_observables.py
argued for: local energies can only push on extensive local sums, so if a
transition in a non-extensive quantity exists at all it has to show up here.
Control parameter is beta. Moves are degree-preserving double-edge swaps on
6-regular graphs, so N and the degree are fixed throughout and the comparison
across beta is like-for-like.

Order parameter, normalised so it is intensive:

    phi = mean graph distance / mean distance of a random 6-regular graph

phi = 1 is indistinguishable from random, phi > 1 is extended.

Two diagnostics, in the order that matters:

    hysteresis     cooling from random against heating from a cubic lattice.
                   If the branches do not meet, the system is either at a
                   first-order transition or simply not equilibrated, and no
                   scaling analysis is valid until that is settled.

    susceptibility chi = N * Var(phi) at each beta. At a genuine continuous
                   transition the peak DIVERGES with N. At a crossover it
                   saturates. This is the decisive measurement.

Sizes are perfect cubes -- 216, 512 -- so a cubic lattice is available as the
ordered starting configuration at every size.

Result: no equilibrium extended phase exists, so there is nothing to
transition into
------------------------------------------------------------------------
The susceptibility scan the module was built around is NOT reported, because
running it revealed it would have been meaningless. What killed it is the
finding that replaces it.

At short times the two branches disagree badly -- gaps up to 0.244 at N = 216
-- and the order parameter falls with N at fixed beta (phi = 1.405, 1.355,
1.174 at beta = 3 for N = 216, 512, 1000), which is the signature of larger
systems equilibrating more slowly rather than of a transition sharpening. A
susceptibility peak measured under those conditions would be an artifact of
where each size happened to sit along its relaxation, not a critical point.

Running one trajectory properly long, at N = 216 and beta = 3, shows why:

    sweeps      phi   4cyc/node   tri/node   mean dist   diam
       100   1.1712       4.653      0.060       3.756      7
       500   1.4599       9.384      0.037       4.682      8
      1500   2.2454      14.000      0.060       7.201     15
      3500   2.5082      15.972      0.139       8.044     20
      7500   1.0502      16.574      0.144       3.368     11
     15500   1.0014      16.690      0.144       3.212      7

The energy improves MONOTONICALLY throughout -- four-cycles per node climb
from 4.65 to 16.69 and never turn round -- while the extension rises to
phi = 2.51 and then collapses back to 1.0014, which is the random-graph value
to three digits.

So the extended structure is a long-lived TRANSIENT, peaking near 3500 sweeps,
on the road to an equilibrium that is locally dense and globally
indistinguishable from random. It is the "locally ordered, globally destroyed"
configuration that global_observables.py had to introduce by hand as case B --
here the dynamics produces it as its own equilibrium.

A correction to geometric_phase.py
------------------------------------
That module reported d_H = 4.03 and |B_2| = 24.05 against the cubic lattice's
24.00, and flagged the regime as "a long-lived non-equilibrium regime, not a
thermodynamic phase". The flag was right and is now confirmed quantitatively:
those measurements were taken around 100 to 3200 sweeps, inside the rising
transient, and the state does not survive.

But its other claim was wrong. It said the clique ground state was
"dynamically inaccessible" because triangles stayed near zero. Longer running
reaches 16.69 four-cycles per node, which EXCEEDS the 15 of disjoint K_7,
while triangles stay at 0.144. So the optimum is being approached, and the
structure it approaches is bipartite-blob like rather than clique like --
consistent with disjoint K_{6,6}, which would give 18.75. Inaccessibility was
the wrong diagnosis; slowness was the right one.

Verdict
-------
The answer to the question this module asked is no. There is no equilibrium
extended phase in this model, so there is no phase transition into one, and
therefore no backdraft mechanism of the kind the cosmological picture needs.

That is a clean negative and it is the one the design was meant to allow. It
also disposes of the most encouraging measurement the project had produced,
which is the right outcome when the encouraging measurement was taken before
equilibrium.
"""

import numpy as np
from collections import deque

from geometric_phase import (random_regular, cubic_lattice, count4_total,
                             c4_through, both_pairs, mean_distance)


def build_state(adj, n):
    el = [(u, v) for u in range(n) for v in adj[u] if u < v]
    pos = {e: i for i, e in enumerate(el)}
    return el, pos


def sweep(adj, n, el, pos, c4, beta, rng, nsteps):
    """nsteps Metropolis double-edge swaps against H = -(4-cycle count)."""
    def rm(a, b):
        k = (a, b) if a < b else (b, a)
        i = pos.pop(k)
        last = el.pop()
        if i < len(el):
            el[i] = last
            pos[last] = i
        adj[a].discard(b)
        adj[b].discard(a)

    def ad(a, b):
        k = (a, b) if a < b else (b, a)
        pos[k] = len(el)
        el.append(k)
        adj[a].add(b)
        adj[b].add(a)

    for _ in range(nsteps):
        a, b = el[int(rng.integers(len(el)))]
        c, d = el[int(rng.integers(len(el)))]
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            c, d = d, c
        if d in adj[a] or b in adj[c]:
            continue
        before = (c4_through(adj, a, b) + c4_through(adj, c, d)
                  - both_pairs(adj, (a, b), (c, d)))
        rm(a, b)
        rm(c, d)
        ad(a, d)
        ad(c, b)
        after = (c4_through(adj, a, d) + c4_through(adj, c, b)
                 - both_pairs(adj, (a, d), (c, b)))
        dH = -(after - before)
        if dH <= 0 or rng.random() < np.exp(-beta * dH):
            c4 += after - before
        else:
            rm(a, d)
            rm(c, b)
            ad(a, b)
            ad(c, d)
    return c4


def phi_of(adj, n, ref):
    return mean_distance(adj, n) / ref


def branch(adj, n, betas, ref, rng, eq_sweeps, samp, gap):
    """Walk beta in the given order, carrying the configuration forward."""
    el, pos = build_state(adj, n)
    c4 = count4_total(adj, n)
    out = []
    for beta in betas:
        c4 = sweep(adj, n, el, pos, c4, beta, rng, int(eq_sweeps * n))
        vals = []
        for _ in range(samp):
            c4 = sweep(adj, n, el, pos, c4, beta, rng, int(gap * n))
            vals.append(phi_of(adj, n, ref))
        v = np.array(vals)
        out.append((beta, float(v.mean()), float(v.var()),
                    c4 / n))
    return out


BETAS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
SIZES = [(6, 216), (8, 512)]


def trajectory(n, beta, chunks, seed=7):
    """Long single run, tracking energy and extension together."""
    from geometric_phase import triangles_per_node
    from self_amplification import shell_profile
    ref = mean_distance(random_regular(n, 6, np.random.default_rng(9)), n)
    A = [set(a) for a in random_regular(n, 6, np.random.default_rng(1))]
    el, pos = build_state(A, n)
    c4 = count4_total(A, n)
    rng = np.random.default_rng(seed)
    rows, tot = [], 0
    for ch in chunks:
        c4 = sweep(A, n, el, pos, c4, beta, rng, int(ch * n))
        tot += ch
        sh = shell_profile(A, n, srcs=6)
        rows.append((tot, phi_of(A, n, ref), c4 / n,
                     triangles_per_node(A, n), mean_distance(A, n), max(sh)))
    return ref, rows


def main():
    print("=" * 78)
    print("1. HYSTERESIS: do cooling and heating branches meet?")
    print("=" * 78)
    for L, n in SIZES:
        ref = mean_distance(random_regular(n, 6, np.random.default_rng(9)), n)
        cold = random_regular(n, 6, np.random.default_rng(1))
        up = branch([set(a) for a in cold], n, BETAS, ref,
                    np.random.default_rng(7), 40, 8, 4)
        lat, _ = cubic_lattice(L)
        down = branch([set(a) for a in lat], n, BETAS[::-1], ref,
                      np.random.default_rng(7), 40, 8, 4)[::-1]
        print("  N = %d   (random-graph reference mean distance %.3f)"
              % (n, ref))
        print("     %6s %14s %14s %10s"
              % ("beta", "phi from rnd", "phi from lat", "gap"))
        for (b1, m1, _, _), (_, m2, _, _) in zip(up, down):
            print("     %6.2f %14.4f %14.4f %10.4f" % (b1, m1, m2, m2 - m1))
        print()
    print("  The branches do approach each other, so the short-time gap is")
    print("  slow equilibration rather than a first-order transition. But 40")
    print("  sweeps is nowhere near enough, as part 2 shows, so any scaling")
    print("  analysis at these times would be measuring a transient.")
    print()

    print("=" * 78)
    print("2. LONG RUN: the extended state is a TRANSIENT")
    print("=" * 78)
    print("  N=216, beta=3.0, single run from a random start.")
    print()
    ref, rows = trajectory(216, 3.0, (100, 400, 1000, 2000, 4000, 8000))
    print("  random-graph reference mean distance %.3f" % ref)
    print("  %8s %8s %11s %10s %10s %8s"
          % ("sweeps", "phi", "4cyc/node", "tri/node", "mean dist", "diam"))
    for r in rows:
        print("  %8d %8.4f %11.3f %10.3f %10.3f %8d" % r)
    print()
    print("  The energy improves MONOTONICALLY (4-cycles per node 4.65 to")
    print("  16.69) while the extension rises to phi = 2.51 and then falls")
    print("  back to 1.00, exactly the random-graph value. The extended")
    print("  structure is a long-lived transient on the way to an")
    print("  equilibrium that is locally dense and globally random.")
    print()
    print("  Note 16.69 four-cycles per node exceeds the 15 of disjoint K_7,")
    print("  and triangles stay at 0.14, so the endpoint is bipartite-blob")
    print("  like rather than clique like.")
    print()

    print("=" * 78)
    print("3. VERDICT")
    print("=" * 78)
    print("  There is no equilibrium extended phase, so there is no phase")
    print("  transition into one, and no backdraft. The susceptibility scan")
    print("  this module was built to run is not reported: at the sweep")
    print("  counts it can afford, every size sits at a different point on a")
    print("  NON-MONOTONIC transient, so its peaks would be artifacts.")
    print()


if __name__ == "__main__":
    main()
