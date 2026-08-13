#!/usr/bin/env python3
"""
rewrite_algebra.py

Stop positing conservation laws. Compute them.

Every invariant in this project so far was guessed and then checked: Q = B + 2N
- E was written down by hand, verified to have residual 0, and only afterwards
understood. That works once. It does not scale, and it quietly biases the
search toward the invariants we already had in mind.

This is the engine version. A rewrite is described only by what it does to the
cell counts,

    delta n = (dN_0, dN_1, dN_2, ...)

those vectors are stacked into a matrix R, and every exact invariant is an
integer solution of

    R q = 0,      giving   Q = q . n

Doing this over the integers rather than the reals matters, and the difference
is not cosmetic. A real null space finds free invariants and nothing else.
Smith normal form finds free invariants AND torsion -- quantities conserved only
modulo some integer. A rule set with an invariant that is exactly preserved mod
2 has a parity sector, which is precisely the structure a defect number would
live in, and a floating-point null space cannot see it at all.

The complex, and why the boundary operator earns its place
------------------------------------------------------------
The state is a chain complex rather than a graph:

    C_0 <- C_1 <- C_2 <- ...     with    d_{k-1} d_k = 0

That identity is checked here, not assumed, because it is what separates a
topological face from an accidental graph cycle -- the distinction this project
has been conflating since the beginning. With it, Betti numbers come from
linear algebra instead of from cycle counting:

    beta_k = dim ker d_k - rank d_{k+1}

Validation before use
----------------------
An engine that discovers laws is worthless if it cannot rediscover the ones we
already know, so every part is run against a known answer first:

    d^2 = 0             on a circle, a sphere and a torus
    Betti numbers       (1,1), (1,0,1), (1,2,1) -- known exactly
    Euler characteristic 0, 2, 0
    invariant discovery must find Q = B + 2N - E from the close/open/subdivide
                        change vectors, having never been told it
    invariant discovery must find chi from the Pachner moves
    invariant discovery must find NOTHING when the rules admit nothing

The last one matters as much as the others. A method that always returns an
invariant is not measuring anything.

The distortion classifier
--------------------------
A candidate rewrite is also scored geometrically before any simulation, which
is what would have caught `close` on day one. Take the patch boundary vertices,
form their distance matrix M before and M' after, and fit the best pure
dilation:

    delta_R = min_lambda ||M' - lambda M||_F / ||M||_F

    lambda > 1, delta small     expansion
    lambda < 1                  a shortcut
    delta large                 an anisotropic tear

Result 1: the topology layer validates exactly
------------------------------------------------
    complex                        counts     d^2=0   Betti (F2)   chi
    circle  (triangle boundary)    [3,3]        yes    [1, 1]        0
    sphere  (tetrahedron bdry)     [4,6,4]      yes    [1, 0, 1]     2
    torus   (3x3 grid)             [9,27,18]    yes    [1, 2, 1]     0

Every Betti number and Euler characteristic is the known answer, and d^2 = 0
holds on all three. So the chain-complex machinery is trustworthy before it is
used on anything unknown -- which is the whole point of running it on objects
whose answers were fixed in advance.

Result 2: it rediscovers our own conservation law, unprompted
--------------------------------------------------------------
    rules                             invariant found       expected
    close/open/subdivide, no bank     none                  none
    close/open/subdivide, with bank   B + 2N - E            Q = B + 2N - E
    subdivision only                  N - E                 the cycle rank
    Pachner 1-3 and 2-2               3V - E, 2V - F        chi in the span: YES
    a rule doubling faces             V, E; torsion [2]     F mod 2

The second line is the one that matters. Q = B + 2N - E took a long time to
arrive at by hand, and the engine returns it from three change vectors with no
knowledge of what a bank is or why it was introduced.

The first line matters just as much: given the same rules WITHOUT a bank the
answer is "none", correctly. A method that always returns an invariant measures
nothing, and this one declines when there is nothing to find. It also explains
the bank retrospectively -- the bank is not a physical object, it is the extra
coordinate that makes the kernel non-empty.

The Pachner row shows a real subtlety. The kernel is two-dimensional and the
returned basis is 3V - E and 2V - F, neither of which looks like chi. But
chi = V - E + F is an integer combination of them, which the span check
confirms. So the engine finds the invariant LATTICE; naming a preferred
generator inside it is a human choice, not a mathematical one.

The last row is the reason for doing this over Z. A real null space sees only V
and E. Smith normal form additionally returns torsion coefficient 2, i.e. F is
conserved modulo 2 -- a parity sector, exactly the structure a defect number
would live in, and completely invisible to floating point.

Result 3: the classifier reproduces months of simulation in milliseconds
-------------------------------------------------------------------------
Best pure dilation of the patch boundary, on a 9-node ring:

    move                       lambda   delta   reads as
    subdivide one edge          1.111   0.157   tear
    subdivide every edge        2.000   0.000   expansion
    close (distance-2 chord)    0.889   0.157   shortcut

This is a retrospective check with some force behind it. holonomy_defect.py
needed a 64-fold budget scan and a control to establish that localized
refinement cannot expand while uniform refinement can; scale_factor.py,
paid_expansion.py and model.py spent far longer establishing that close is a
shortcut. The classifier returns both from three tiny patches, instantly, with
lambda = 0.889 < 1 for close and delta = 0.000 for uniform subdivision.

Had this existed at the start, `close` would have been rejected before any
cosmology was built on it.

Verdict
-------
The change of representation pays for itself immediately. Three things that
were previously hand-work are now mechanical and exact: topology via d^2 = 0
and Betti numbers, conservation via the integer kernel including torsion, and
the shortcut/tear/expansion distinction via the dilation fit.

What this does NOT yet do is enumerate rewrites. It scores rules it is handed.
Automatic enumeration of valid local complex rewrites, and the Pareto search
over the resulting score vectors, is the next piece -- and it is only worth
building now that each individual score has been validated against a known
answer.

One honest note on the second table: the bank rules were entered with dN and dE
transposed on the first run, and the engine returned B - N + 2E, which is the
correct invariant of the rules it was actually given. The engine caught my
typo, not the other way round. That is the argument for the whole approach in
one line.
"""

import numpy as np
from itertools import combinations
from collections import deque


# --------------------------------------------------------------- complexes

class Complex:
    """Simplicial complex: cells[k] is an ordered list of frozensets."""

    def __init__(self, top_simplices):
        self.cells = {}
        for s in top_simplices:
            for k in range(len(s)):
                for f in combinations(sorted(s), k + 1):
                    self.cells.setdefault(k, set()).add(frozenset(f))
        self.cells = {k: sorted(v, key=lambda x: sorted(x))
                      for k, v in self.cells.items()}
        self.index = {k: {c: i for i, c in enumerate(v)}
                      for k, v in self.cells.items()}

    def dim(self):
        return max(self.cells)

    def counts(self):
        return [len(self.cells.get(k, [])) for k in range(self.dim() + 1)]

    def euler(self):
        return sum((-1) ** k * n for k, n in enumerate(self.counts()))

    def boundary(self, k):
        """d_k : C_k -> C_{k-1} as an integer matrix."""
        if k <= 0 or k not in self.cells:
            return np.zeros((0, len(self.cells.get(0, []))), dtype=object)
        rows = len(self.cells[k - 1])
        cols = len(self.cells[k])
        M = [[0] * cols for _ in range(rows)]
        for j, c in enumerate(self.cells[k]):
            vs = sorted(c)
            for i_, v in enumerate(vs):
                face = frozenset(vs[:i_] + vs[i_ + 1:])
                M[self.index[k - 1][face]][j] = (-1) ** i_
        return M


def matmul(A, B):
    if not A or not B:
        return []
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]


def is_zero(M):
    return all(all(x == 0 for x in row) for row in M)


def rank_f2(M):
    """Rank over F_2 by Gaussian elimination."""
    if not M or not M[0]:
        return 0
    A = [[x % 2 for x in row] for row in M]
    rows, cols, r = len(A), len(A[0]), 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c]:
                A[i] = [(a ^ b) for a, b in zip(A[i], A[r])]
        r += 1
    return r


def betti_f2(K):
    out = []
    for k in range(K.dim() + 1):
        nk = len(K.cells.get(k, []))
        rk = rank_f2(K.boundary(k)) if k >= 1 else 0
        rk1 = rank_f2(K.boundary(k + 1)) if (k + 1) in K.cells else 0
        out.append(nk - rk - rk1)
    return out


# ------------------------------------------------------- Smith normal form

def smith(Ain):
    """U A V = D diagonal, over Z. Returns D, V, rank."""
    A = [list(map(int, row)) for row in Ain]
    if not A or not A[0]:
        return [], [], 0
    m, n = len(A), len(A[0])
    V = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def col_op(j1, j2, q):                      # col j1 -= q * col j2
        for i in range(m):
            A[i][j1] -= q * A[i][j2]
        for i in range(n):
            V[i][j1] -= q * V[i][j2]

    def col_swap(j1, j2):
        for i in range(m):
            A[i][j1], A[i][j2] = A[i][j2], A[i][j1]
        for i in range(n):
            V[i][j1], V[i][j2] = V[i][j2], V[i][j1]

    def row_op(i1, i2, q):
        for j in range(n):
            A[i1][j] -= q * A[i2][j]

    def row_swap(i1, i2):
        A[i1], A[i2] = A[i2], A[i1]

    t = 0
    while t < min(m, n):
        piv = None
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < best):
                    best, piv = abs(A[i][j]), (i, j)
        if piv is None:
            break
        row_swap(t, piv[0])
        col_swap(t, piv[1])
        while True:
            done = True
            for i in range(t + 1, m):
                if A[i][t]:
                    row_op(i, t, A[i][t] // A[t][t])
                    if A[i][t]:
                        row_swap(t, i)
                        done = False
            for j in range(t + 1, n):
                if A[t][j]:
                    col_op(j, t, A[t][j] // A[t][t])
                    if A[t][j]:
                        col_swap(t, j)
                        done = False
            if done and all(A[i][t] == 0 for i in range(t + 1, m)) \
                    and all(A[t][j] == 0 for j in range(t + 1, n)):
                break
        t += 1
    D = [A[i][i] if i < min(m, n) else 0 for i in range(min(m, n))]
    rank = sum(1 for d in D if d != 0)
    return D, V, rank


def integer_kernel(R):
    """Basis of {q : R q = 0} over Z, plus the torsion coefficients of R."""
    if not R or not R[0]:
        return [], []
    D, V, rank = smith(R)
    n = len(R[0])
    basis = [[V[i][j] for i in range(n)] for j in range(rank, n)]
    torsion = [d for d in D if d not in (0, 1, -1)]
    return basis, torsion


def _rank_q(M):
    """Rank over Q by fraction-free elimination."""
    A = [list(map(int, r)) for r in M]
    if not A or not A[0]:
        return 0
    rows, cols, r = len(A), len(A[0]), 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c]:
                f, g = A[r][c], A[i][c]
                A[i] = [f * a - g * b for a, b in zip(A[i], A[r])]
        r += 1
    return r


def normalise(q):
    from math import gcd
    g = 0
    for x in q:
        g = gcd(g, abs(x))
    if g > 1:
        q = [x // g for x in q]
    for x in q:
        if x != 0:
            return [-y for y in q] if x < 0 else q
    return q


def show(q, names):
    parts = []
    for c, nm in zip(q, names):
        if c == 0:
            continue
        parts.append("%s%d%s" % ("+" if c > 0 and parts else
                                 ("-" if c < 0 else ""), abs(c), nm)
                     if abs(c) != 1 else
                     "%s%s" % ("+" if c > 0 and parts else
                               ("-" if c < 0 else ""), nm))
    return " ".join(parts) if parts else "0"


# ------------------------------------------------------------- classifier

def distortion(before_adj, after_adj, boundary):
    """Best pure dilation of the patch boundary, and the residual."""
    def dmat(adj, nodes):
        M = np.zeros((len(nodes), len(nodes)))
        for a, i in enumerate(nodes):
            d = {i: 0}
            q = deque([i])
            while q:
                x = q.popleft()
                for y in adj[x]:
                    if y not in d:
                        d[y] = d[x] + 1
                        q.append(y)
            for b, j in enumerate(nodes):
                M[a, b] = d.get(j, np.inf)
        return M
    M = dmat(before_adj, boundary)
    Mp = dmat(after_adj, boundary)
    ok = np.isfinite(M) & np.isfinite(Mp) & (M > 0)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    lam = float((M[ok] * Mp[ok]).sum() / (M[ok] ** 2).sum())
    res = float(np.linalg.norm(Mp[ok] - lam * M[ok])
                / np.linalg.norm(M[ok]))
    return lam, res


# ------------------------------------------------------------------- main

def torus_grid(L=3):
    tris = []
    for x in range(L):
        for y in range(L):
            a = (x % L) * L + (y % L)
            b = ((x + 1) % L) * L + (y % L)
            c = ((x + 1) % L) * L + ((y + 1) % L)
            e = (x % L) * L + ((y + 1) % L)
            tris.append(frozenset((a, b, c)))
            tris.append(frozenset((a, c, e)))
    return tris


def main():
    print("=" * 78)
    print("1. VALIDATION: d^2 = 0, Betti numbers, Euler characteristic")
    print("=" * 78)
    tests = [
        ("circle  (triangle boundary)",
         [frozenset(c) for c in combinations(range(3), 2)], [1, 1], 0),
        ("sphere  (tetrahedron boundary)",
         [frozenset(c) for c in combinations(range(4), 3)], [1, 0, 1], 2),
        ("torus   (3x3 grid)", torus_grid(3), [1, 2, 1], 0),
    ]
    print("  %-32s %10s %8s %16s %8s"
          % ("complex", "counts", "d^2=0", "Betti (F2)", "chi"))
    for name, tops, want_b, want_chi in tests:
        K = Complex(tops)
        ok = all(is_zero(matmul(K.boundary(k - 1), K.boundary(k)))
                 for k in range(2, K.dim() + 1))
        b = betti_f2(K)
        print("  %-32s %10s %8s %16s %8d   %s"
              % (name, K.counts(), "yes" if ok else "NO", b, K.euler(),
                 "ok" if (b == want_b and K.euler() == want_chi) else "MISMATCH"))
    print()

    print("=" * 78)
    print("2. INVARIANT DISCOVERY, against rules whose answer we already know")
    print("=" * 78)

    cases = [
        ("close / open / subdivide, WITHOUT a bank",
         ["N", "E"],
         [[0, 1], [0, -1], [1, 1]],
         "nothing -- this is why the bank was introduced"),
        ("close / open / subdivide, WITH a bank",
         ["B", "N", "E"],
         [[1, 0, 1], [-1, 0, -1], [-1, 1, 1]],
         "Q = B + 2N - E"),
        ("subdivision only",
         ["N", "E"],
         [[1, 1]],
         "E - N (the cycle rank)"),
        ("Pachner 1-3 and 2-2 in 2D",
         ["V", "E", "F"],
         [[1, 3, 2], [0, 0, 0]],
         "chi = V - E + F"),
        ("a rule that doubles faces",
         ["V", "E", "F"],
         [[0, 0, 2]],
         "F mod 2 -- torsion, invisible to a real null space"),
    ]
    named = {"close / open / subdivide, WITH a bank": ("B + 2N - E",
                                                       [1, 2, -1]),
             "Pachner 1-3 and 2-2 in 2D": ("chi = V - E + F", [1, -1, 1])}
    def in_span(basis, target):
        """Is `target` an integer combination of the kernel basis?"""
        if not basis:
            return False
        A = [[b[i] for b in basis] for i in range(len(target))]
        aug = [row + [t] for row, t in zip(A, target)]
        return rank_f2(A) == rank_f2(aug) and _rank_q(A) == _rank_q(aug)

    for name, names, R, expect, *chk in cases:
        basis, tors = integer_kernel(R)
        print("  %s" % name)
        print("     change vectors      %s" % R)
        if basis:
            for q in basis:
                print("     invariant found     %s" % show(normalise(q), names))
        else:
            print("     invariant found     none")
        if tors:
            print("     torsion             %s" % tors)
        if name in named:
            lbl, tgt = named[name]
            print("     %-19s %s"
                  % ("contains " + lbl + "?",
                     "YES" if in_span(basis, tgt) else "NO"))
        print("     expected            %s" % expect)
        print()

    print("=" * 78)
    print("3. THE GEOMETRIC CLASSIFIER ON OUR OWN OLD MOVES")
    print("=" * 78)
    L = 9
    base = {i: {(i - 1) % L, (i + 1) % L} for i in range(L)}
    bnd = [0, 3, 6]

    sub = {i: set(v) for i, v in base.items()}
    sub[0].discard(1); sub[1].discard(0)
    sub[L] = {0, 1}; sub[0].add(L); sub[1].add(L)

    clo = {i: set(v) for i, v in base.items()}
    clo[0].add(2); clo[2].add(0)

    allsub = {}
    nxt = L
    for i in range(L):
        allsub[i] = set()
    for i in range(L):
        j = (i + 1) % L
        allsub[nxt] = {i, j}
        allsub[i].add(nxt)
        allsub[j].add(nxt)
        nxt += 1

    print("  %-34s %10s %10s   %s" % ("move", "lambda", "delta", "reads as"))
    for nm, g in (("subdivide one edge", sub),
                  ("subdivide every edge", allsub),
                  ("close (distance-2 chord)", clo)):
        lam, res = distortion(base, g, bnd)
        tag = ("shortcut" if lam < 0.99 else
               "expansion" if lam > 1.01 and res < 0.1 else
               "tear" if res >= 0.1 else "refinement, metric unchanged")
        print("  %-34s %10.3f %10.3f   %s" % (nm, lam, res, tag))
    print()


if __name__ == "__main__":
    main()
