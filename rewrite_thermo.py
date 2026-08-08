#!/usr/bin/env python3
"""
rewrite_thermo.py

Boltzmann entropy and entropy-production rate for motif coarse-grainings of
a continuation graph.

Two quantities, with very different status on a reversible grammar:

Boltzmann entropy
    S_B(M) = ln Omega(M), where Omega(M) is the number of microstates
    (continuation-graph states) carrying motif label M. Non-trivial, and
    genuinely dependent on the choice of coarse-graining.

Entropy production rate
    sigma = 1/2 sum_{u,v} [pi(u)P(u,v) - pi(v)P(v,u)]
                          * ln[ pi(u)P(u,v) / (pi(v)P(v,u)) ]

    On a grammar whose rules are mutually inverse the continuation graph is
    symmetric, the chain satisfies detailed balance, and sigma = 0 exactly.
    Coarse-graining cannot rescue this. For any partition into blocks A, B,
    pi-weighted lumping gives

        Pi(A)Phat(A,B) = sum_{u in A, v in B} pi(u)P(u,v)
                       = sum_{u in A, v in B} pi(v)P(v,u)      (fine balance)
                       = Pi(B)Phat(B,A)

    so coarse detailed balance is inherited from fine detailed balance, for
    *every* motif definition. A reversible fine chain cannot manufacture
    coarse entropy production. (The converse is the familiar one: coarse-
    graining an irreversible chain *under*-estimates its true production.)

Meaningful entropy production therefore requires a grammar that breaks
detailed balance. The --grammar flag runs the mixed grammar for contrast,
where the chain is absorbing and the relevant object is the quasi-stationary
distribution rather than the stationary one.
"""

import argparse
import hashlib
import math
from collections import defaultdict

from rewrite_lab import Hypergraph, continuation_graph, path_seed
from rewrite_research import GRAMMARS

# ---------- Motif coarse-grainings -------------------------------------


def cg_arity(H):
    """Multiset of edge sizes: the crudest motif census."""
    prof = defaultdict(int)
    for e in H.edges:
        prof[len(e)] += 1
    return tuple(sorted(prof.items()))


def cg_degseq(H):
    """Sorted degree sequence."""
    deg = defaultdict(int)
    for e in H.edges:
        for v in e:
            deg[v] += 1
    return tuple(sorted(deg.values()))


def _digest(obj):
    """
    Stable content hash, used to keep refined colours from growing without
    bound across rounds.

    It must be a *global* function of the colour's content. Renumbering
    colours to small ints per graph would break WL's monotonicity: colour 0
    in one state would be unrelated to colour 0 in another, making signatures
    incomparable across states and letting an extra round merge blocks it
    should only ever split. hashlib is used rather than hash() because the
    latter is salted per process for str/bytes.
    """
    return hashlib.blake2b(repr(obj).encode(), digest_size=8).hexdigest()


def _wl(H, rounds):
    """
    Weisfeiler-Leman colour refinement on the incidence structure.

    Initial colour of a node is its (degree, multiset of incident edge sizes).
    Each round refines by the multiset, over incident edges, of the edge size
    together with the colours of the edge's other members. The state signature
    is the sorted multiset of final node colours -- an isomorphism-invariant
    local motif census, refining monotonically with the number of rounds.
    """
    inc = defaultdict(list)
    for e in H.edges:
        for v in e:
            inc[v].append(e)

    colour = {v: _digest((len(es), tuple(sorted(len(e) for e in es))))
              for v, es in inc.items()}

    for _ in range(rounds):
        nxt = {}
        for v, es in inc.items():
            sig = []
            for e in es:
                others = tuple(sorted(colour[u] for u in e if u != v))
                sig.append((len(e), others))
            nxt[v] = _digest((colour[v], tuple(sorted(sig))))
        colour = nxt

    return tuple(sorted(colour.values()))


def cg_wl0(H):
    return _wl(H, 0)


def cg_wl1(H):
    return _wl(H, 1)


def cg_wl2(H):
    return _wl(H, 2)


def cg_trivial(H):
    """Everything in one block -- the maximal coarse-graining."""
    return ()


COARSE = {
    "arity": cg_arity,
    "degseq": cg_degseq,
    "wl0": cg_wl0,
    "wl1": cg_wl1,
    "wl2": cg_wl2,
    "trivial": cg_trivial,
}

# ---------- Chains -----------------------------------------------------


def transition(graph, keys, lazy=True):
    """Uniform-successor chain, optionally lazy. Sinks absorb."""
    P = defaultdict(dict)
    for u in keys:
        succ = sorted(graph.get(u, ()))
        if not succ:
            P[u][u] = 1.0
            continue
        if lazy:
            P[u][u] = 0.5
            for v in succ:
                P[u][v] = P[u].get(v, 0.0) + 0.5 / len(succ)
        else:
            for v in succ:
                P[u][v] = 1.0 / len(succ)
    return P


def metropolis(graph, keys):
    """
    Metropolis kernel whose stationary distribution is *uniform* over
    microstates -- the microcanonical ensemble.

    This matters for testing S = k ln Omega. Boltzmann's relation presupposes
    equal a priori probability of microstates, but a uniform-successor random
    walk on an undirected graph equilibrates to pi ~ degree, not to uniform.
    Under that walk ln Omega is simply not the equilibrium entropy, and the
    "defect" column measures precisely that mismatch.

    Proposing a uniform neighbour and accepting with probability
    min(1, deg(u)/deg(v)) makes the kernel symmetric, hence doubly stochastic,
    hence uniform-stationary, while remaining reversible. Rejections stay put,
    which also makes the chain aperiodic -- so the bipartite-periodicity
    problem that required a lazy walk earlier disappears on its own.
    """
    deg = {u: len(graph.get(u, ())) for u in keys}
    P = defaultdict(dict)
    for u in keys:
        succ = sorted(graph.get(u, ()))
        if not succ:
            P[u][u] = 1.0
            continue
        stay = 1.0
        for v in succ:
            acc = min(1.0, deg[u] / deg[v]) if deg[v] else 1.0
            p = acc / deg[u]
            P[u][v] = P[u].get(v, 0.0) + p
            stay -= p
        if stay > 1e-15:
            P[u][u] = P[u].get(u, 0.0) + stay
    return P


def evolve(P, p0, steps):
    """Exact forward evolution of the *fine* distribution."""
    p = dict(p0)
    for t in range(steps + 1):
        yield t, p
        nxt = defaultdict(float)
        for u, mass in p.items():
            if mass <= 0:
                continue
            for v, pr in P[u].items():
                nxt[v] += mass * pr
        p = dict(nxt)


def stationary(P, keys, iters=50000, tol=1e-15):
    n = len(keys)
    pi = {k: 1.0 / n for k in keys}
    for _ in range(iters):
        nxt = defaultdict(float)
        for u in keys:
            for v, p in P[u].items():
                nxt[v] += pi[u] * p
        delta = sum(abs(nxt.get(k, 0.0) - pi[k]) for k in keys)
        pi = {k: nxt.get(k, 0.0) for k in keys}
        if delta < tol:
            break
    return pi


def quasi_stationary(P, keys, absorbing, iters=50000, tol=1e-15):
    """
    Quasi-stationary distribution: power iteration on the transient
    sub-stochastic kernel with renormalisation at each step.

    On an absorbing chain the stationary distribution sits entirely on the
    absorbing set and nothing moves, so sigma is trivially zero. The QSD is
    the distribution conditioned on not yet having been absorbed, and is the
    object for which entropy production has content.
    """
    trans = [k for k in keys if k not in absorbing]
    if not trans:
        return {}
    n = len(trans)
    q = {k: 1.0 / n for k in trans}
    tset = set(trans)
    for _ in range(iters):
        nxt = defaultdict(float)
        for u in trans:
            for v, p in P[u].items():
                if v in tset:
                    nxt[v] += q[u] * p
        mass = sum(nxt.values())
        if mass <= 0:
            break
        nxt = {k: nxt.get(k, 0.0) / mass for k in trans}
        delta = sum(abs(nxt[k] - q[k]) for k in trans)
        q = nxt
        if delta < tol:
            break
    return q


# ---------- Thermodynamic quantities -----------------------------------


def shannon(dist):
    return -sum(p * math.log(p) for p in dist.values() if p > 0)


def entropy_production(P, pi, keys):
    """
    Returns (sigma_finite, n_absolutely_irreversible).

    A transition with pi(u)P(u,v) > 0 but P(v,u) = 0 is absolutely
    irreversible and contributes +infinity. Those are counted separately and
    excluded from sigma_finite, which is then the production rate restricted
    to the bidirectional support.
    """
    sigma = 0.0
    hard = 0
    for u in keys:
        for v, puv in P[u].items():
            if u == v:
                continue
            fwd = pi.get(u, 0.0) * puv
            rev = pi.get(v, 0.0) * P.get(v, {}).get(u, 0.0)
            if fwd <= 0 and rev <= 0:
                continue
            if fwd > 0 and rev <= 0:
                hard += 1
                continue
            if rev > 0 and fwd <= 0:
                continue
            sigma += 0.5 * (fwd - rev) * math.log(fwd / rev)
    return sigma, hard


def lump(P, pi, keys, label):
    """pi-weighted lumping onto motif classes. Returns (Pi, Phat, blocks)."""
    blocks = defaultdict(list)
    for k in keys:
        blocks[label[k]].append(k)

    Pi = {M: sum(pi.get(k, 0.0) for k in ks) for M, ks in blocks.items()}
    Phat = defaultdict(dict)
    for M, ks in blocks.items():
        if Pi[M] <= 0:
            continue
        acc = defaultdict(float)
        for u in ks:
            for v, p in P[u].items():
                acc[label[v]] += pi.get(u, 0.0) * p
        for N, f in acc.items():
            Phat[M][N] = f / Pi[M]
    return Pi, Phat, blocks


def report(name, states, graph, cg_names, use_qsd=False):
    keys = list(states)
    P = transition(graph, keys, lazy=True)
    absorbing = {k for k in keys if not graph.get(k)}

    if use_qsd and absorbing:
        pi = quasi_stationary(P, keys, absorbing)
        dist_kind = "quasi-stationary"
        active = list(pi)
    else:
        pi = stationary(P, keys)
        dist_kind = "stationary"
        active = keys

    sigma, hard = entropy_production(P, pi, active)

    print(f"{name}: {len(states)} states  [{dist_kind}]")
    print(f"  fine  S_Gibbs = {shannon(pi):.4f} nats     "
          f"S_B(total) = ln {len(states)} = {math.log(len(states)):.4f}")
    print(f"  fine  sigma   = {sigma:.3e} nats/step"
          + (f"   (+{hard} absolutely irreversible transitions => +inf)"
             if hard else ""))
    print(f"  {'motif':<9} {'macro':>5} {'<S_B>':>8} {'maxS_B':>8} "
          f"{'S(Pi)':>8} {'defect':>8} {'sigma_hat':>11}")

    partitions = {}
    for cn in cg_names:
        fn = COARSE[cn]
        label = {k: fn(states[k]) for k in keys}
        Pi, Phat, blocks = lump(P, pi, active, label)
        partitions[cn] = {k: label[k] for k in active}

        omega = {M: len(ks) for M, ks in blocks.items()}
        mean_SB = sum(Pi.get(M, 0.0) * math.log(omega[M])
                      for M in Pi if omega.get(M))
        max_SB = math.log(max(omega.values())) if omega else 0.0

        # Within-block entropy of pi, and the exact decomposition
        #     S(pi) = S(Pi) + sum_M Pi(M) H(pi | M).
        # <S_B> >= H_within always, with equality iff pi is uniform inside
        # every block. The gap ("defect") therefore measures how badly the
        # equal-a-priori-probability assumption behind Boltzmann entropy
        # fails for this motif class -- i.e. whether it is a legitimate
        # macrostate at all.
        H_within = 0.0
        for M, ks in blocks.items():
            if Pi.get(M, 0.0) <= 0:
                continue
            sub = {k: pi.get(k, 0.0) / Pi[M] for k in ks if pi.get(k, 0.0) > 0}
            H_within += Pi[M] * shannon(sub)
        assert abs((shannon(Pi) + H_within) - shannon(pi)) < 1e-9, \
            "entropy decomposition S(pi) = S(Pi) + <H_within> failed"

        s_hat, hard_hat = entropy_production(Phat, Pi, list(Pi))

        print(f"  {cn:<9} {len(blocks):>5} {mean_SB:>8.4f} {max_SB:>8.4f} "
              f"{shannon(Pi):>8.4f} {mean_SB - H_within:>8.4f} {s_hat:>11.3e}"
              + (f"  (+{hard_hat} inf)" if hard_hat else ""))

    # WL refinement must be monotone in the number of rounds.
    for a, b in (("wl0", "wl1"), ("wl1", "wl2")):
        if a in partitions and b in partitions:
            for u in active:
                for v in active:
                    if partitions[b][u] == partitions[b][v]:
                        assert partitions[a][u] == partitions[a][v], \
                            f"{b} is not a refinement of {a}"
    print()


def relaxation(name, states, graph, cg_name, steps, seed_key=None):
    """
    The gas-in-a-box experiment (section 5): prepare a low-entropy initial
    macrostate and ask whether <Delta S_B> > 0 as it relaxes.

    The fine chain is evolved *exactly* and the macrostate distribution is
    obtained by projection, so nothing here depends on the coarse process
    being Markov -- which in general it is not.

    Note what the arrow does and does not come from. The kernel is reversible
    and has zero steady-state entropy production; the dynamics contains no
    arrow. The arrow comes from preparing a low-Omega initial condition in a
    state space where high-Omega macrostates vastly outnumber low-Omega ones.
    That is Boltzmann's answer, reproduced rather than programmed in.
    """
    keys = list(states)
    P = metropolis(graph, keys)
    N = len(keys)

    # Sanity: the kernel must be symmetric (=> doubly stochastic => uniform
    # stationary => S = ln Omega is the right equilibrium entropy).
    for u in keys:
        for v, p in P[u].items():
            if u != v:
                assert abs(p - P[v].get(u, 0.0)) < 1e-12, "kernel not symmetric"

    fn = COARSE[cg_name]
    label = {k: fn(states[k]) for k in keys}
    blocks = defaultdict(list)
    for k in keys:
        blocks[label[k]].append(k)
    omega = {M: len(ks) for M, ks in blocks.items()}

    if seed_key is None:
        # Prepare the lowest-Omega macrostate: "all the gas on the left".
        M0 = min(omega, key=lambda M: (omega[M], sorted(blocks[M])[0]))
        seed_key = sorted(blocks[M0])[0]

    p0 = {seed_key: 1.0}

    print(f"{name}  motif={cg_name}  N={N}  macrostates={len(blocks)}")
    print(f"  start: Omega={omega[label[seed_key]]}  "
          f"S_B={math.log(omega[label[seed_key]]):.4f}   "
          f"equilibrium S_B(max)=ln {N} = {math.log(N):.4f}")
    print(f"  {'t':>4} {'<S_B>':>9} {'S_fine':>9} {'H=D(p||u)':>10} "
          f"{'dS_B':>10}")

    prev_SB = None
    traj = []
    for t, p in evolve(P, p0, steps):
        Pi = defaultdict(float)
        for k, m in p.items():
            Pi[label[k]] += m
        SB = sum(m * math.log(omega[M]) for M, m in Pi.items() if m > 0)
        S_fine = shannon(p)
        H = math.log(N) - S_fine          # D(p || uniform)
        dSB = None if prev_SB is None else SB - prev_SB
        traj.append((t, SB, S_fine, H, dSB))
        prev_SB = SB

    show = [r for r in traj if r[0] in
            (0, 1, 2, 3, 5, 8, 12, 20, 35, 60, 100, steps)]
    for t, SB, S_fine, H, dSB in show:
        print(f"  {t:>4} {SB:>9.4f} {S_fine:>9.4f} {H:>10.4f} "
              + ("     --   " if dSB is None else f"{dSB:>+10.5f}"))

    # H-theorem: for a doubly stochastic kernel the fine Shannon entropy is
    # non-decreasing, so H = D(p || uniform) must never increase.
    for (t0, _, _, h0, _), (t1, _, _, h1, _) in zip(traj, traj[1:]):
        assert h1 <= h0 + 1e-12, f"H-theorem violated at t={t0}->{t1}"

    total = traj[-1][1] - traj[0][1]
    steps_up = sum(1 for r in traj[1:] if r[4] > 1e-15)
    steps_dn = sum(1 for r in traj[1:] if r[4] < -1e-15)
    print(f"  total <Delta S_B> = {total:+.4f} nats over {steps} steps "
          f"({steps_up} up, {steps_dn} down)")
    print(f"  monotone increase: {steps_dn == 0}")
    print()
    return total, steps_dn


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--paths", type=int, nargs="+", default=[4, 5, 6],
                    metavar="N", help="path seeds (default: 4 5 6 -> 7/16/41)")
    ap.add_argument("--grammar", choices=sorted(GRAMMARS),
                    default="G1-reversible")
    ap.add_argument("--coarse", nargs="+", default=["arity", "degseq",
                                                    "wl1", "wl2", "trivial"],
                    choices=sorted(COARSE))
    ap.add_argument("--qsd", action="store_true",
                    help="use the quasi-stationary distribution (absorbing chains)")
    ap.add_argument("--relax", action="store_true",
                    help="run the relaxation experiment (section 5) instead")
    ap.add_argument("--steps", type=int, default=200,
                    help="relaxation steps (default: 200)")
    args = ap.parse_args(argv)

    rules = GRAMMARS[args.grammar]["rules"]
    print("=" * 72)
    print(f"grammar = {args.grammar}"
          + ("   [relaxation from prepared low-entropy state]"
             if args.relax else ""))
    print("=" * 72)
    print()

    for n in args.paths:
        states, graph, _ = continuation_graph(path_seed(n), rules=rules)
        if args.relax:
            for cn in args.coarse:
                relaxation(f"path({n})", states, graph, cn, args.steps)
        else:
            report(f"path({n})", states, graph, args.coarse, use_qsd=args.qsd)


if __name__ == "__main__":
    main()
