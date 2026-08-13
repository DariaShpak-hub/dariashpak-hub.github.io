#!/usr/bin/env python3
"""
coherence_field.py

Can a local coherence interaction pay for geometry?

geometry_decay.py showed a lattice dies in O(N) moves because every conserved
quantity is blind to order. The proposal under test: give each node a
two-state coherence variable, let neighbours prefer to agree,

    H = -J sum_<ij> s_i s_j

and ask whether the lattice half-life stops being O(N) and starts growing.
The mechanism is borrowed from magnetism, not the physics -- what is wanted
is the competition between a local interaction favouring order and entropy
favouring disorder.

There is a necessary condition to check first, and checking it is cheap.
Protection requires the coherence energy to be LOWER on the lattice than on a
random graph of the same size and the same edge count. If the two are equal
the field is exactly as blind as the charge and can protect nothing.

Result 1: in the ordered phase the field is exactly blind
-----------------------------------------------------------
Comparing a square lattice against the same graph rewired by degree-preserving
swaps, so N and E match by construction, with spins equilibrated at each J:

      J      lattice <H>/N   random <H>/N     difference
    0.100         -0.02222       -0.02356       +0.00133
    0.200         -0.09067       -0.08444       -0.00622
    0.300         -0.18267       -0.16933       -0.01333
    0.441         -0.60327       -0.70512       +0.10185
    0.600         -1.16000       -1.16800       +0.00800
    0.900         -1.79200       -1.79200        0.00000
    1.500         -3.00000       -3.00000        0.00000

At large J the difference is exactly zero, to every digit. That is not a
numerical accident, it is forced: once the spins are fully aligned every term
s_i s_j equals 1, so H = -J E, which depends only on the edge count. A
ferromagnet deep in its ordered phase cannot see the shape of the graph it
lives on, only how many bonds it has. The conserved charge Q = bank + 2N - E
was blind for exactly the same reason.

Result 2: near criticality the field has the WRONG sign
---------------------------------------------------------
The one place the two differ appreciably is around J = 0.44, and there the
random graph has the LOWER energy, -0.705 against -0.603. The coherence field
prefers disorder.

This is not bad luck either. The ordering temperature of an Ising model rises
with connectivity: a random 4-regular graph is locally tree-like with
J_c = atanh(1/3) = 0.347, while the square lattice has J_c = 0.441. At any J
between those values the random graph is already ordered while the lattice is
not, so the random graph wins on energy. Expanders order more easily than
lattices, because order spreads faster when neighbourhoods branch.

A caveat on Results 1 and 2, which matters. What actually weights a graph in
the joint Boltzmann distribution is the spin FREE energy F(G) = -ln Z(G), not
the mean energy <H> tabulated above. At large J the distinction vanishes --
the spins are frozen aligned, Z = 2 exp(J E) for either graph, so the free
energies are equal too and the blindness is exact. Near criticality it does
not vanish, and the energy comparison there is suggestive rather than
decisive. That is why the half-life measurement below is the authoritative
test: the joint dynamics samples F automatically, whatever it is.

Result 3: the field changes nothing
-------------------------------------
Half-life of the lattice, in sweeps, against the J = 0 control:

      J          N=400   N=900   N=1600   N=3600
    0.000         1.59    2.48     3.17     2.88     (control)
    0.300         2.04    2.64     2.48     3.00
    0.441         1.76    2.34     2.68     2.93
    0.900         2.21    2.34     2.42     2.83

Every row is flat in N and sits at two to three sweeps. Deviations from the
control run to +-30% with no systematic sign and no trend in either J or N.
The lattice still dies after a fixed number of touches per node, in O(N)
moves, exactly as it did with no coherence field at all.

So the answer to the question this module asked is no. Adding the simplest
local coherence interaction does not change O(N) into anything that grows.

(One confound, reported rather than hidden: at J = 0.9 the mean degree drifts
to 4.36 instead of 4.04, because aligned spins make adding a bond always
favourable and removing one always costly, so closes outrun opens. It does not
rescue anything -- that row is as flat as the others.)

Why the analogy does not transfer
-----------------------------------
Magnetism stabilises SPIN order, and spin order is easiest on precisely the
graphs that are not geometric. Ordering spreads faster when neighbourhoods
branch, so an expander orders at a weaker coupling than a lattice does. A
ferromagnetic term therefore rewards connectivity, which is the opposite of
what a geometry needs.

And in the regime where the interaction is strongest -- deep order -- it
degenerates into a function of the edge count alone and stops seeing the graph
entirely. That is the same failure mode as the conserved charge, arrived at
from a different direction, and the exact zeros in Result 1 are the cleanest
statement of it anywhere in this project.

The contrast with geometric_phase.py is the useful part. The 4-cycle energy
worked, partially, because it is a functional of the graph's own loop
structure -- geometric information, read directly off the graph. An Ising
energy is a functional of a field living ON the graph, and couples to the
graph only through the edge count and, weakly and with the wrong sign, through
the loop corrections near criticality.

The lesson is not "coherence fields cannot work" but something sharper: the
thing that pays for geometry has to be a function of the geometry, not a
function of matter arranged on it. A field added on top can only see the
graph through whatever the field's own physics happens to couple to, and
ferromagnetic coupling sees bonds, not shape.
"""

import numpy as np
from collections import deque

from geometry_decay import Lattice, mean_distance, random_4regular, half_life


class Coherent(Lattice):
    """Lattice carrying an Ising coherence variable on each node."""

    def __init__(self, m, J):
        super().__init__(m)
        self.J = J
        self.s = np.ones(self.n, dtype=np.int64)

    def local_field(self, u):
        return sum(int(self.s[w]) for w in self.adj[u])

    def spin_flip(self, rng):
        u = int(rng.integers(self.n))
        dH = 2.0 * self.J * int(self.s[u]) * self.local_field(u)
        if dH <= 0 or rng.random() < np.exp(-dH):
            self.s[u] = -self.s[u]

    def spin_burn(self, rng, sweeps=40):
        for _ in range(sweeps * self.n):
            self.spin_flip(rng)

    def coherence_energy(self):
        e = 0
        for u in range(self.n):
            for v in self.adj[u]:
                if v > u:
                    e += int(self.s[u]) * int(self.s[v])
        return -self.J * e / self.n

    # graph moves, now Metropolis-accepted against the coherence energy
    def close(self, rng):
        for _ in range(40):
            a, b = self.el[int(rng.integers(len(self.el)))]
            v = a if rng.random() < 0.5 else b
            nb = list(self.adj[v])
            if len(nb) < 2:
                continue
            x = nb[int(rng.integers(len(nb)))]
            y = nb[int(rng.integers(len(nb)))]
            if x == y or y in self.adj[x]:
                continue
            if len(self.adj[x]) >= self.cap or len(self.adj[y]) >= self.cap:
                continue
            dH = -self.J * int(self.s[x]) * int(self.s[y])
            if dH > 0 and rng.random() >= np.exp(-dH):
                continue
            self._add(x, y)
            self.bank += 1
            self.closes += 1
            return True
        return False

    def open(self, rng):
        if self.bank < 1:
            return False
        for _ in range(40):
            a, b = self.el[int(rng.integers(len(self.el)))]
            if not (self.adj[a] & self.adj[b]):
                continue
            dH = self.J * int(self.s[a]) * int(self.s[b])
            if dH > 0 and rng.random() >= np.exp(-dH):
                continue
            self._del(a, b)
            self.bank -= 1
            self.opens += 1
            return True
        return False


def randomized_copy(m, seed=9, per_node=20):
    L = Lattice(m)
    r = np.random.default_rng(seed)
    for _ in range(per_node * L.n):
        L.swap(r)
    return L


def spin_equilibrate(adj, n, J, rng, sweeps=200, start_up=False):
    s = (np.ones(n, dtype=np.int64) if start_up
         else rng.integers(0, 2, n).astype(np.int64) * 2 - 1)
    nb = [np.fromiter(a, int) for a in adj]
    for _ in range(sweeps):
        for u in rng.permutation(n):
            dH = 2.0 * J * int(s[u]) * int(s[nb[u]].sum())
            if dH <= 0 or rng.random() < np.exp(-dH):
                s[u] = -s[u]
    return s


def energy_of(adj, n, s, J):
    e = 0
    for u in range(n):
        for v in adj[u]:
            if v > u:
                e += int(s[u]) * int(s[v])
    return -J * e / n


def decay(m, J, marks_sweeps, seed=3, spins_per_move=6):
    U = Coherent(m, J)
    n = U.n
    rng = np.random.default_rng(seed)
    if J > 0:
        U.spin_burn(rng, sweeps=30)
    d_start = mean_distance(U.adj, n)
    d_rand = mean_distance(random_4regular(n, np.random.default_rng(9)), n)
    marks = sorted({max(1, int(s * n)) for s in marks_sweeps})
    rows = []
    turn = 0
    for t in range(1, marks[-1] + 1):
        for _ in range(8):
            ok = U.close(rng) if turn == 0 else U.open(rng)
            if ok:
                turn ^= 1
                break
            if turn == 1:
                turn = 0
        if J > 0:
            for _ in range(spins_per_move):
                U.spin_flip(rng)
        if t in marks:
            d = mean_distance(U.adj, n)
            rows.append(dict(sweeps=t / n, d=d,
                             phi=(d - d_rand) / (d_start - d_rand),
                             k=2 * len(U.el) / n,
                             mag=float(abs(U.s.mean()))))
    return U, d_start, d_rand, rows


MARKS = [0.25, 0.5, 1, 2, 4, 8, 16, 32]


def main():
    print("=" * 78)
    print("1. NECESSARY CONDITION: can the field tell a lattice from random?")
    print("=" * 78)
    m = 30
    L, R = Lattice(m), randomized_copy(m)
    n = L.n
    print("  N=%d   lattice E=%d   randomized E=%d   (matched by construction)"
          % (n, len(L.el), len(R.el)))
    print()
    print("  %7s %15s %15s %13s %9s"
          % ("J", "lattice <H>/N", "random <H>/N", "difference", "lat |m|"))
    for J in (0.1, 0.2, 0.3, 0.4407, 0.6, 0.9, 1.5):
        up = J > 0.5
        sl = spin_equilibrate(L.adj, n, J, np.random.default_rng(4), start_up=up)
        sr = spin_equilibrate(R.adj, n, J, np.random.default_rng(4), start_up=up)
        el = energy_of(L.adj, n, sl, J)
        er = energy_of(R.adj, n, sr, J)
        print("  %7.3f %15.5f %15.5f %13.5f %9.3f"
              % (J, el, er, el - er, abs(sl.mean())))
    print()
    print("  Zero difference at large J is forced, not numerical: aligned")
    print("  spins give H = -J*E, which sees only the edge count.")
    print("  Positive difference near J=0.44 means the field PREFERS random.")
    print()

    print("=" * 78)
    print("2. DOES THE HALF-LIFE GROW WITH N?")
    print("=" * 78)
    print("  %7s %8s %10s %10s %9s %9s"
          % ("J", "N", "half-life", "final <k>", "final |m|", "vs J=0"))
    base = {}
    for J in (0.0, 0.3, 0.4407, 0.9):
        for mm in (20, 30, 40, 60):
            U, d0, dr, rows = decay(mm, J, MARKS)
            h = half_life([r["sweeps"] for r in rows],
                          [r["phi"] for r in rows])
            if J == 0.0:
                base[U.n] = h
            print("  %7.3f %8d %10.2f %10.3f %9.3f %9s"
                  % (J, U.n, h, rows[-1]["k"], rows[-1]["mag"],
                     "-" if J == 0.0 else "%+.0f%%"
                     % (100 * (h / base[U.n] - 1))))
        print()

    print("=" * 78)
    print("3. VERDICT")
    print("=" * 78)
    print("  Protection would show as the half-life growing with N at some J.")
    print("  It does not: every J gives 2-3 sweeps, flat in N, matching the")
    print("  J=0 control from geometry_decay.py. The coherence field changes")
    print("  nothing.")
    print()
    print("  Why: a ferromagnet rewards connectivity, and expanders order at")
    print("  weaker coupling than lattices do (J_c = 0.347 vs 0.441). Deep in")
    print("  the ordered phase it degenerates to H = -J*E and sees only the")
    print("  bond count. Compare geometric_phase.py, where a 4-cycle energy DID")
    print("  move the system, because it reads the graph's own loop structure")
    print("  instead of a field living on top of it.")
    print()


if __name__ == "__main__":
    main()
