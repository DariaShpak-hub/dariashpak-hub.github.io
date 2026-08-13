#!/usr/bin/env python3
"""
matter_clumping.py

Gravity-like clumping measured as MASS DENSITY, not as hub degree -- and
whether producing it costs the geometry.

self_amplification.py measured concentration by the Gini coefficient of node
degree, and that is the wrong observable. Gravity concentrates mass: many
particles in a small region. A hub is one node with many edges, a different
object entirely, and hub formation is guaranteed for any positive-feedback
rule (Yule 1925, Simon 1955, Barabasi-Albert 1999), so it cannot fail and is
not evidence for anything.

Separating the two needs matter, because in a bare graph "density" and
"connectivity" are the same quantity. So a conserved population of particles
lives on a cubic torus, mean occupancy 20 per node, and two independent
couplings are scanned:

    beta   particles hop preferentially toward neighbours holding more
           particles -- matter attracts matter, an inserted attraction
    alpha  edges are reattached toward particle-dense nodes within graph
           distance 2 -- the geometry responds to matter, which is the
           self-amplifying relational influence being proposed

beta = alpha = 0 is pure diffusion and fixes the null: at occupancy 20 the
density contrast must come out at the Poisson value 1/sqrt(20) = 0.224.

Result: the two effects are orthogonal, and the proposed one is the weaker
--------------------------------------------------------------------------
    beta  alpha   sigma_m   top 1%   max n   dimension   deg Gini
     0.0    0.0     0.233    0.016      35        2.89     0.0000
     0.0    2.0     0.642    0.026      60        4.35     0.3549
     1.0    0.0     1.067    0.055     152        2.89     0.0000
     2.0    0.0     2.981    0.155     346        2.89     0.0000
     4.0    0.0     2.715    0.147     375        2.89     0.0000
     2.0    2.0     4.818    0.344    1113        4.40     0.3375
     4.0    2.0     4.203    0.275    1057        4.39     0.3269
     4.0    4.0     4.116    0.273     690        4.20     0.3256

The control lands on Poisson exactly, 0.233 against 0.224 predicted, so the
measurement is calibrated.

Matter attraction alone (beta = 2, alpha = 0) gives strong concentration --
contrast 13x Poisson, the densest 1% of nodes holding 15.5% of all mass, a
peak node at 346 particles against a mean of 20 -- with the dimension at
2.89, IDENTICAL to the untouched lattice, because the geometry is never
modified.

Geometry feedback alone (beta = 0, alpha = 2) does the opposite. Contrast
reaches only 0.642, under 3x Poisson, while the dimension goes to 4.35 and
the degree Gini to 0.355. It buys the least clumping of any setting tested
and pays for it with the metric.

So mass concentration and metric destruction are not two sides of one
mechanism. They are independent, and the self-amplifying route delivers
mostly the destruction. The earlier reading -- that clumping and geometry
pull against each other -- came from measuring clumping by degree, where the
two are the same thing by construction.

Note also that beta = 4 clumps slightly LESS than beta = 2 (2.715 against
2.981). Very strong hopping bias traps particles in whichever well they first
enter, so the dynamics freezes before the largest wells win.

But the overdensities are not structures
------------------------------------------
Mean pairwise distance among the densest 1% of nodes, against the same number
of randomly chosen nodes on the same torus:

    beta  alpha    densest 1%   random   ratio
     0.0    0.0         5.823    5.941   0.980
     2.0    0.0         5.378    5.941   0.905
     0.0    2.0         6.036    5.941   1.016
     2.0    2.0         5.915    5.941   0.996

At best a 10% deficit. The dense sites are scattered across the torus, not
gathered anywhere. So the strong contrast at beta = 2 is condensation onto
isolated sites, which is the standard zero-range-process behaviour, and not
the extended contiguous overdensity that structure formation means. There are
heavy nodes; there are no halos.

The reason is that a nearest-neighbour hopping bias is a CONTACT interaction.
Each site hoards from its immediate neighbours and the wells never learn about
each other, so they neither merge nor organise. Gravity forms extended
structure because it is long-ranged -- distant mass pulls -- and that long
range is precisely what the 1/r^2 of self_amplification.py comes from, which
in turn comes from Gauss's law in three dimensions.

Which closes the loop rather tidily. Even matter clumping, the one thing here
that looked like it worked, does not produce gravitational structure on its
own. It needs the long-range force, and the long-range force needs the stable
3D metric.

What this does and does not show
----------------------------------
It does not derive gravity. beta is an inserted attraction -- a local hopping
bias put in by hand -- so "insert attraction, get clumping" is not a finding.

What is a finding is the geometry column. Newtonian gravity clumps matter on a
fixed background, and this reproduces exactly that: strong clustering with the
metric untouched, because nothing requires the substrate to change for mass to
concentrate. The proposal was that clumping should emerge from geometry
responding to density. Measured against mass rather than degree, that
mechanism is both weaker at clumping and destructive of the geometry, and
there is no setting where it beats plain attraction on a static lattice.

The open question is unchanged and untouched: what selects three dimensions.
Everything else -- 1/r^2 by Gauss's law, clumping by attraction -- follows
once a stable 3D substrate exists, and nothing in this thread produces one.
"""

import numpy as np

from self_amplification import (cubic_torus, shell_profile, d_from_shells,
                                gini)

OCC = 20


def diffuse(n_arr, nbr, rng, moves, beta):
    """Move particles; with beta > 0 they prefer fuller neighbours."""
    N = len(n_arr)
    p = n_arr / n_arr.sum()
    srcs = rng.choice(N, size=moves, p=p)
    u = rng.random(moves)
    for idx in range(moves):
        s = int(srcs[idx])
        if n_arr[s] <= 0:
            continue
        nb = nbr[s]
        if len(nb) == 0:
            continue
        if beta == 0:
            t = int(nb[int(u[idx] * len(nb))])
        else:
            w = (n_arr[nb] + 1.0) ** beta
            c = np.cumsum(w)
            t = int(nb[int(np.searchsorted(c, u[idx] * c[-1]))])
        n_arr[s] -= 1
        n_arr[t] += 1


def rewire(adj, nbr, n_arr, rng, moves, alpha, cap):
    """Reattach edges toward dense nodes within graph distance 2."""
    N = len(adj)
    for _ in range(moves):
        a = int(rng.integers(N))
        if not adj[a] or len(adj[a]) >= cap:
            continue
        b = list(adj[a])[int(rng.integers(len(adj[a])))]
        if len(adj[b]) <= 1:
            continue
        near = set()
        for x in adj[a]:
            near |= adj[x]
        near -= adj[a]
        near.discard(a)
        near = {t for t in near if len(adj[t]) < cap}
        if not near:
            continue
        cand = np.fromiter(near, int)
        w = (n_arr[cand] + 1.0) ** alpha
        t = int(cand[int(rng.choice(len(cand), p=w / w.sum()))])
        adj[a].discard(b)
        adj[b].discard(a)
        adj[a].add(t)
        adj[t].add(a)
    for i in range(N):
        nbr[i] = np.fromiter(adj[i], int)


def spatial_spread(coord, L, sel, rng, trials=400):
    """Mean pairwise periodic distance among `sel`, and the random baseline.

    Distinguishes a localised overdensity from scattered heavy sites.
    """
    def mpd(idx):
        a = coord[idx]
        tot, cnt = 0.0, 0
        for _ in range(trials):
            i, j = rng.integers(len(idx), size=2)
            if i == j:
                continue
            d = a[i] - a[j]
            d = np.minimum(d % L, (-d) % L)
            tot += float(np.sqrt((d.astype(float) ** 2).sum()))
            cnt += 1
        return tot / max(cnt, 1)

    rand = rng.choice(len(coord), size=len(sel), replace=False)
    return mpd(sel), mpd(rand)


def run(L, beta, alpha, sweeps=80, rw=150, cap=12, seed=3):
    adj, N, coord = cubic_torus(L)
    rng = np.random.default_rng(seed)
    n_arr = np.full(N, OCC, dtype=np.int64)
    nbr = [np.fromiter(a, int) for a in adj]
    for _ in range(sweeps):
        diffuse(n_arr, nbr, rng, 2 * N, beta)
        if alpha > 0:
            rewire(adj, nbr, n_arr, rng, rw, alpha, cap)
    sh = shell_profile(adj, N)
    order = np.argsort(n_arr)[::-1]
    k = max(2, N // 100)
    top = order[:k]
    d_sel, d_rand = spatial_spread(coord, L, top, np.random.default_rng(11))
    return dict(sigma=float(n_arr.std() / n_arr.mean()),
                top1=float(n_arr[top].sum() / n_arr.sum()),
                maxn=int(n_arr.max()),
                dim=d_from_shells(sh, 2, L // 2 - 1),
                dgini=gini([len(a) for a in adj]),
                d_sel=d_sel, d_rand=d_rand)


def main():
    L = 12
    print("=" * 78)
    print("MASS CLUMPING vs METRIC   (cubic torus L=%d, N=%d, occupancy %d)"
          % (L, L ** 3, OCC))
    print("=" * 78)
    print("  beta  = matter attracts matter (inserted hopping bias)")
    print("  alpha = geometry reattaches toward dense nodes (the proposal)")
    print("  Poisson null for the density contrast: 1/sqrt(%d) = %.3f"
          % (OCC, 1 / np.sqrt(OCC)))
    print()
    print("  %6s %6s %9s %8s %8s %10s %9s"
          % ("beta", "alpha", "sigma_m", "top 1%", "max n", "dimension",
             "deg Gini"))
    rows = {}
    for beta, alpha in ((0, 0), (0, 2), (1, 0), (2, 0), (4, 0),
                        (2, 2), (4, 2), (4, 4)):
        r = run(L, beta, alpha)
        rows[(beta, alpha)] = r
        print("  %6.1f %6.1f %9.3f %8.3f %8d %10.2f %9.4f"
              % (beta, alpha, r["sigma"], r["top1"], r["maxn"], r["dim"],
                 r["dgini"]))
    print()
    print("  Control (0,0) lands on Poisson, so the contrast is calibrated.")
    print("  (2,0): 13x Poisson clumping with the dimension UNCHANGED at 2.89.")
    print("  (0,2): under 3x Poisson, and the dimension is gone.")
    print("  Mass concentration and metric damage are independent effects,")
    print("  and the self-amplifying route supplies mostly the damage.")
    print()

    print("=" * 78)
    print("ARE THE OVERDENSITIES LOCALISED, OR JUST SCATTERED HEAVY SITES?")
    print("=" * 78)
    print("  Mean pairwise distance among the densest 1%% of nodes, against")
    print("  the same count of randomly chosen nodes. Lower = clustered.")
    print()
    print("  %6s %6s %16s %16s %10s"
          % ("beta", "alpha", "densest 1%", "random nodes", "ratio"))
    for key in ((0, 0), (2, 0), (0, 2), (2, 2)):
        r = rows[key]
        print("  %6.1f %6.1f %16.3f %16.3f %10.3f"
              % (key[0], key[1], r["d_sel"], r["d_rand"],
                 r["d_sel"] / r["d_rand"]))
    print()


if __name__ == "__main__":
    main()
