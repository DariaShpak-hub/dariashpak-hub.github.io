#!/usr/bin/env python3
"""
necessity.py

Constructive necessity test: what is forced once relation, alternative and
consistency exist, and nothing else?

The project has been testing machinery it built later rather than the seed it
started from. This goes back to the seed. Three primitives only --

    RELATION      a state says which distinctions stand in relation
    ALTERNATIVE   the same relational exterior admits more than one interior
    CONSISTENCY   a local substitution must remain compatible with what it
                  touches

-- and nothing else. No t, x, d, energy, mass, probability, entropy, particle,
bank, colour, transport.

The claim to be tested is that conservation and multiplicity are not further
assumptions but algebraic shadows of RESTRICTION, and that an effective
interaction follows from counting alone. Every ingredient of the later model
then earns one of three labels:

    LOGICALLY REQUIRED      cannot be removed and still have a system
    MATHEMATICALLY EMERGENT follows with no extra assumption
    ADDED ASSUMPTION        put in by hand, must justify itself

Why this is exact rather than sampled
---------------------------------------
Six distinctions, all 2^15 = 32768 labelled relational states, enumerated in
full. Nothing is estimated: components are found by union-find over the entire
state space, multiplicities are exact integer counts, and the invariant claim
is checked against every state rather than a sample.

The two layers
---------------
    R1  ALTERNATIVE alone. Any single relation may be toggled. This is local
        substitution with no compatibility requirement at all.

    R2  ALTERNATIVE + CONSISTENCY. A substitution must preserve what each
        distinction commits to its exterior. Concretely the double swap:
        relations ab and cd are replaced by ad and cb. Every vertex keeps its
        degree, so no distinction's external commitment changes -- this is
        interface preservation with nothing geometric assumed.

The measurements are the same at both layers, so the difference between the
columns IS the contribution of consistency:

    sectors        connected components of the state space under the moves
    invariant      is there a quantity constant on components, unposited
    multiplicity   the exact count of states per sector
    drift          the stationary distribution of the bare walk
    information    what a later state still says about its sector
    interaction    Omega(m) for m relations between two halves; whether
                   -ln Omega(m) has a minimum, which is an effective force
                   derived from counting and nothing else

Result 1: consistency is what creates everything
--------------------------------------------------
    layer                          sectors   moves/state   deg vector const
    R0  relation only                32768           0.0   (vacuous)
    R1  alternative only                 1          15.0   no
    R2  alternative + consistency     6944           5.6   YES

With alternatives and no compatibility requirement the entire state space is
ONE component. Nothing is conserved, nothing distinguishes where you started,
and the walk entropy is exactly maximal (10.3972 of 10.3972) because every
state has the same 15 moves. Local substitution on its own produces no physics
whatsoever.

Add the single requirement that a substitution preserve what each distinction
commits to its exterior, and the space shatters into 6944 sectors.

Result 2: the invariant is discovered, complete, and exact
-----------------------------------------------------------
    connected components under the consistency move   6944
    distinct degree vectors over all 32768 states     6944

Equal. Not approximately -- the component count and the degree-vector count are
the same integer, and the degree vector is verified constant on every component
across all 32768 states. So the degree vector is a COMPLETE invariant: two
states are mutually reachable if and only if they agree on it.

Nothing posited it. It is the algebraic shadow of restriction, exactly as
predicted. Conservation is not a primitive.

Result 3: multiplicity, drift and persistent information follow
-----------------------------------------------------------------
    sector sizes         min 1, median 2, max 70
    walk entropy R1      10.3972 of max 10.3972   -- no drift at all
    walk entropy R2      10.1639 of max 10.3972   -- genuine drift
    information R1       I = ln(1)    = 0.0000 nats
    information R2       I = ln(6944) = 8.8456 nats, forever

Multiplicity is just counting, and it is wildly unequal once sectors exist: 1
state in the smallest, 70 in the largest. At R1 every state has the same number
of moves so the stationary distribution is uniform and there is no drift; at R2
move counts differ and the walk is measurably biased toward high-multiplicity
states.

The information line is the sharpest. Under R1 a state says nothing about where
it came from -- everything mixes. Under R2, 8.85 nats about the initial sector
survive for all time, with no memory object, no record, no transport variable
and no stored bit. Persistent information is a consequence of restriction.

Result 4: an effective interaction, from counting and nothing else
-------------------------------------------------------------------
Split the six distinctions into two halves and count states by m, the number of
relations across the cut.

    R1, whole space          m   0    1     2     3     4     5  ...
                        Omega   64  576  2304  5376  8064  8064
                    -ln Omega  -4.16 -6.36 -7.74 -8.59 -8.99 -8.99

    R2, largest sector       m   0    2     4     6
                        Omega    1   18    45     6
                    -ln Omega  0.00 -2.89 -3.81 -1.79

-ln Omega has a minimum at m = 4 in both. That is an effective potential with a
preferred separation and a restoring tendency either side of it, and the only
thing behind it is how many arrangements exist. No force was introduced.

The R2 row carries something extra: only EVEN m occurs. Cross-cut parity is
conserved, because a double swap changes the cross count by (a xor b)+(c xor d)
against (a xor d)+(c xor b), which agree mod 2. A Z_2 sector appears without
being added.

But it is not an independent invariant, and saying otherwise would be wrong:
m = sum of degrees over one half, minus twice the internal edges, so m mod 2 is
already a function of the degree vector -- which the component count proved
complete. The parity is a shadow of the shadow.

VERDICT
-------
    ingredient               label                    evidence
    relation                 LOGICALLY REQUIRED       no state space without it
    alternative              LOGICALLY REQUIRED       nothing changes without it
    consistency              ADDED ASSUMPTION         but it generates the rest
    conservation             EMERGENT                 6944 = 6944, exact
    multiplicity             EMERGENT                 counting, 1 to 70
    statistical drift        EMERGENT (R2 only)       10.1639 vs 10.3972
    persistent information   EMERGENT                 8.8456 nats, no record
    effective interaction    EMERGENT                 -ln Omega minimum at m = 4
    parity / Z_2 sector      EMERGENT, NOT NEW        a function of the degrees

    bank B                   ADDED ASSUMPTION         absent here
    colour count k           ADDED ASSUMPTION         absent here
    transport U_e            ADDED ASSUMPTION         absent here
    refinement level r       ADDED ASSUMPTION         absent here
    METRIC                   ABSENT                   nothing here is a distance
    DIMENSION                ABSENT                   nothing here has one

Five things the project spent modules building turn out to be free: conserved
charges, multiplicity, drift, persistent information and an effective
interaction. All five follow from restricting alternatives, and none needs a
mechanism.

The two that do not appear are the two the project has never managed to obtain
either: there is no distance and no dimension anywhere in this construction.
That is consistent rather than disappointing -- it says metric is genuinely a
further structure and not a consequence of relational consistency, which is the
same conclusion the geometry experiments reached from the other direction.

Scope. Six distinctions, one notion of consistency (degree preservation), one
bipartition for the interaction test. The double swap is a natural reading of
interface preservation but not the only one, and a different consistency
condition would give a different invariant -- the claim tested here is that SOME
invariant is forced, not that the degree vector is privileged.
"""

import numpy as np
from itertools import combinations

N = 6
PAIRS = list(combinations(range(N), 2))
NP = len(PAIRS)
BIT = {p: i for i, p in enumerate(PAIRS)}
NSTATES = 1 << NP

A_SIDE, B_SIDE = {0, 1, 2}, {3, 4, 5}
CROSS = [i for i, (u, v) in enumerate(PAIRS)
         if (u in A_SIDE) != (v in A_SIDE)]


def degrees(s):
    d = [0] * N
    for i, (u, v) in enumerate(PAIRS):
        if s >> i & 1:
            d[u] += 1
            d[v] += 1
    return tuple(d)


def toggle_moves(s):
    return [s ^ (1 << i) for i in range(NP)]


def swap_moves(s):
    """ab, cd present and ad, cb absent -> remove ab, cd; add ad, cb."""
    out = []
    present = [PAIRS[i] for i in range(NP) if s >> i & 1]
    for (a, b), (c, d) in combinations(present, 2):
        if len({a, b, c, d}) != 4:
            continue
        for x, y in (((a, d), (c, b)), ((a, c), (b, d))):
            i = BIT.get((min(x), max(x)))
            j = BIT.get((min(y), max(y)))
            if i is None or j is None:
                continue
            if (s >> i & 1) or (s >> j & 1):
                continue
            t = s
            t &= ~(1 << BIT[(a, b)])
            t &= ~(1 << BIT[(c, d)])
            t |= (1 << i) | (1 << j)
            out.append(t)
    return out


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def analyse(movefn, name):
    uf = UF(NSTATES)
    deg = np.zeros(NSTATES, dtype=np.int64)
    for s in range(NSTATES):
        ms = movefn(s)
        deg[s] = len(ms)
        for t in ms:
            uf.union(s, t)
    root = np.array([uf.find(s) for s in range(NSTATES)])
    roots, sizes = np.unique(root, return_counts=True)
    comp = {r: i for i, r in enumerate(roots)}
    lab = np.array([comp[r] for r in root])

    dv = [degrees(s) for s in range(NSTATES)]
    n_degseq = len(set(dv))
    const = True
    seen = {}
    for s in range(NSTATES):
        c = lab[s]
        if c in seen and seen[c] != dv[s]:
            const = False
            break
        seen.setdefault(c, dv[s])

    tot = deg.sum()
    pi = deg / tot if tot else np.zeros(NSTATES)
    live = pi[pi > 0]
    H = float(-(live * np.log(live)).sum())
    return dict(name=name, sectors=len(roots), sizes=sizes, labels=lab,
                degvec_constant=const, n_degseq=n_degseq,
                mean_moves=float(deg.mean()), H=H,
                Hmax=float(np.log(NSTATES)), deg=deg)


def interaction(labels, want_sector):
    """Omega(m) for m relations across the A|B cut, inside one sector."""
    counts = {}
    for s in range(NSTATES):
        if labels[s] != want_sector:
            continue
        m = sum(1 for i in CROSS if s >> i & 1)
        counts[m] = counts.get(m, 0) + 1
    return counts


def main():
    print("=" * 78)
    print("LAYER R0: RELATION ONLY")
    print("=" * 78)
    print("  %d distinctions, %d possible relations, %d states."
          % (N, NP, NSTATES))
    print("  With no alternatives there are no moves, so every state is its")
    print("  own sector. Nothing is conserved because nothing changes, and")
    print("  the question 'what is invariant' has no content yet.")
    print()

    res = {}
    for fn, nm in ((toggle_moves, "R1  alternative only"),
                   (swap_moves, "R2  alternative + consistency")):
        res[nm] = analyse(fn, nm)

    print("=" * 78)
    print("WHAT EACH LAYER PRODUCES")
    print("=" * 78)
    print("  %-32s %10s %10s %14s"
          % ("layer", "sectors", "moves/state", "deg vector const"))
    for nm, r in res.items():
        print("  %-32s %10d %10.1f %14s"
              % (nm, r["sectors"], r["mean_moves"],
                 "YES" if r["degvec_constant"] else "no"))
    print()
    print("  distinct degree vectors over all %d states: %d"
          % (NSTATES, res["R2  alternative + consistency"]["n_degseq"]))
    print()

    for nm, r in res.items():
        print("  %s" % nm)
        print("     sectors               %d" % r["sectors"])
        print("     largest sector        %d states (%.1f%% of the space)"
              % (r["sizes"].max(), 100.0 * r["sizes"].max() / NSTATES))
        print("     multiplicity spread   min %d, median %d, max %d"
              % (r["sizes"].min(), int(np.median(r["sizes"])),
                 r["sizes"].max()))
        print("     walk entropy          %.4f of max %.4f"
              % (r["H"], r["Hmax"]))
        print("     information about the sector that survives forever:")
        print("        I = ln(sectors) = %.4f nat%s"
              % (np.log(r["sectors"]), "" if r["sectors"] == 1 else "s"))
        print()

    print("=" * 78)
    print("EFFECTIVE INTERACTION FROM COUNTING ALONE")
    print("=" * 78)
    print("  Split the six distinctions into two halves. Inside one sector,")
    print("  count states by m, the number of relations across the cut. If")
    print("  -ln Omega(m) has a minimum then there is a preferred m, and that")
    print("  is an effective force with nothing but counting behind it.")
    print()
    for nm, r in res.items():
        big = int(np.argmax(np.bincount(r["labels"])))
        c = interaction(r["labels"], big)
        if not c:
            continue
        print("  %s   (largest sector, %d states)" % (nm, sum(c.values())))
        print("     %4s %10s %12s %10s" % ("m", "Omega(m)", "-ln Omega",
                                           "force"))
        ms = sorted(c)
        V = {m: -np.log(c[m]) for m in ms}
        for m in ms:
            f = (V[m + 1] - V[m]) if (m + 1) in V else float("nan")
            print("     %4d %10d %12.4f %10.4f" % (m, c[m], V[m], f))
        mstar = min(ms, key=lambda m: V[m])
        print("     minimum of -ln Omega at m = %d" % mstar)
        print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    rows = [("relation", "LOGICALLY REQUIRED", "no state space without it"),
            ("alternative", "LOGICALLY REQUIRED", "nothing changes without it"),
            ("consistency", "ADDED ASSUMPTION", "but it generates the rest"),
            ("conservation", "EMERGENT", "%d sectors = %d degree vectors"
             % (res["R2  alternative + consistency"]["sectors"],
                res["R2  alternative + consistency"]["n_degseq"])),
            ("multiplicity", "EMERGENT", "counting alone"),
            ("statistical drift", "EMERGENT (R2 only)", "walk entropy below max"),
            ("persistent information", "EMERGENT", "%.4f nats, no record kept"
             % np.log(res["R2  alternative + consistency"]["sectors"])),
            ("effective interaction", "EMERGENT", "-ln Omega has a minimum"),
            ("bank / colour / transport", "ADDED ASSUMPTION", "absent here"),
            ("metric", "ABSENT", "nothing here is a distance"),
            ("dimension", "ABSENT", "nothing here has one")]
    for a, b, c in rows:
        print("  %-26s %-20s %s" % (a, b, c))
    print()


if __name__ == "__main__":
    main()
