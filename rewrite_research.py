#!/usr/bin/env python3
"""
rewrite_research.py

Experimental instrument for the question: which properties of a rewrite
rule set determine whether non-trivial effective structure emerges at all?

Motivation
----------
A rule set whose rules are mutually inverse produces a *symmetric* continuation
graph. On such a graph several natural observables are degenerate by
construction, not by accident:

    * SCC structure   -- one component, always
    * recoverability  -- identically 1
    * probability flux-- identically 0 (detailed balance)
    * "vacuum"        -- every state qualifies, vacuously

The degeneracy has a single cause: a conserved quantity plus reversible rules
gives no ordering on states. What breaks it is a *Lyapunov function* --
a quantity monotone under every rewrite but not constant. This module supplies
the machinery to look for one and to measure what appears once it exists.

Design
------
Three rule sets are compared on identical seeds, with predictions registered
in advance (see RULE_SETS below):

    G1  compose + decompose            reversible; predicted fully degenerate
    G2  compose only                   strictly monotone; predicted acyclic
    G3  compose + decompose + contract mixed; predicted layered hierarchy

G3 is the interesting one. ``compose``/``decompose`` preserve node count and
are mutually inverse, so they generate reversible dynamics *within* a layer.
``contract`` merges two nodes joined by an edge, strictly decreasing node
count, so it drives irreversible descent *between* layers. That is the
standard physical picture -- fast reversible dynamics within a shell, slow
irreversible relaxation downward, a unique ground state at the bottom -- built
from nothing but rewrite rules.
"""

import argparse
from collections import defaultdict

from rewrite_lab import (
    Hypergraph,
    Rule,
    compose_binary_to_ternary,
    decompose_ternary,
    continuation_graph,
    path_seed,
    weight,
)

# ---------- An irreversible rule ---------------------------------------


def contract_edge(H):
    """
    Merge two nodes joined by a common edge (edge contraction).

    This is coarse-graining expressed as a rewrite. It strictly decreases the
    node count, so node count is a Lyapunov function for any rule set
    containing this rule. Singleton edges produced by the merge are dropped:
    an edge on one node carries no relational content.
    """
    seen = set()

    for e in H.edges:
        verts = sorted(e)
        for i in range(len(verts)):
            for j in range(i + 1, len(verts)):
                x, y = verts[i], verts[j]

                new_edges = set()
                for f in H.edges:
                    g = frozenset(x if v == y else v for v in f)
                    if len(g) >= 2:
                        new_edges.add(g)

                G = Hypergraph(new_edges)
                c = G.canonical()
                if c not in seen:
                    seen.add(c)
                    yield G


COMPOSE = Rule("compose", compose_binary_to_ternary)
DECOMPOSE = Rule("decompose", decompose_ternary)
CONTRACT = Rule("contract", contract_edge)

RULE_SETS = {
    "G1-reversible": {
        "rules": [COMPOSE, DECOMPOSE],
        "prediction": "one SCC, R=1 everywhere, zero flux, no sinks",
    },
    "G2-monotone": {
        "rules": [COMPOSE],
        "prediction": "acyclic, all SCCs singleton, sinks exist, R<1",
    },
    "G3-mixed": {
        "rules": [COMPOSE, DECOMPOSE, CONTRACT],
        "prediction": "layered: non-trivial transient SCCs, one terminal SCC",
    },
}

# ---------- Potentials -------------------------------------------------


def node_count(H):
    return len(H.nodes())


def edge_count(H):
    return len(H.edges)


POTENTIALS = {
    "weight": weight,
    "nodes": node_count,
    "edges": edge_count,
}


def classify_potential(states, graph, fn):
    """
    Classify a candidate potential over the transition relation.

    Returns one of "constant", "lyapunov" (non-increasing along every edge and
    non-constant), or "none".
    """
    values = {k: fn(H) for k, H in states.items()}
    if len(set(values.values())) == 1:
        return "constant"
    for u in graph:
        for v in graph[u]:
            if v in values and values[v] > values[u]:
                return "none"
    return "lyapunov"


# ---------- Graph structure --------------------------------------------


def sccs(graph, nodes):
    """Iterative Tarjan. Returns a list of components (lists of keys)."""
    index, low, onstack = {}, {}, {}
    stack, result, counter = [], [], [0]

    for root in nodes:
        if root in index:
            continue
        index[root] = low[root] = counter[0]
        counter[0] += 1
        stack.append(root)
        onstack[root] = True
        work = [(root, iter(sorted(graph.get(root, ()))))]

        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter[0]
                    counter[0] += 1
                    stack.append(w)
                    onstack[w] = True
                    work.append((w, iter(sorted(graph.get(w, ())))))
                    advanced = True
                    break
                elif onstack.get(w):
                    low[v] = min(low[v], index[w])
            if advanced:
                continue
            work.pop()
            if low[v] == index[v]:
                comp = []
                while True:
                    w = stack.pop()
                    onstack[w] = False
                    comp.append(w)
                    if w == v:
                        break
                result.append(comp)
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
    return result


def terminal_components(graph, comps):
    """Components with no outgoing edge to another component -- the vacua."""
    owner = {}
    for i, c in enumerate(comps):
        for k in c:
            owner[k] = i
    terminal = []
    for i, c in enumerate(comps):
        if all(owner.get(w, i) == i for k in c for w in graph.get(k, ())):
            terminal.append(c)
    return terminal


def is_symmetric(graph):
    return not any(u not in graph.get(v, ()) for u in graph for v in graph[u])


def recoverability(graph, keys):
    """
    R(v) = fraction of out-neighbours w from which v is reachable again.

    This is the local "does the structure reform through the surrounding
    dynamics" measure. Sinks are reported as R = 0.
    """
    # Reverse reachability from each node is expensive; instead compute SCC
    # membership plus forward reachability per component.
    comps = sccs(graph, keys)
    owner = {}
    for i, c in enumerate(comps):
        for k in c:
            owner[k] = i

    # Reachability between components on the condensation DAG.
    dag = defaultdict(set)
    for u in graph:
        for v in graph[u]:
            if owner[u] != owner[v]:
                dag[owner[u]].add(owner[v])

    reach = {}

    def comp_reach(i):
        if i in reach:
            return reach[i]
        reach[i] = set()
        acc = set()
        for j in dag.get(i, ()):
            acc.add(j)
            acc |= comp_reach(j)
        reach[i] = acc
        return acc

    for i in range(len(comps)):
        comp_reach(i)

    out = {}
    for v in keys:
        succ = graph.get(v, set())
        if not succ:
            out[v] = 0.0
            continue
        back = 0
        for w in succ:
            # v is reachable from w iff same component, or v's component is
            # reachable from w's component.
            if owner[w] == owner[v] or owner[v] in reach[owner[w]]:
                back += 1
        out[v] = back / len(succ)
    return out


def stationary_and_flux(graph, keys, iters=20000, tol=1e-14):
    """
    Lazy uniform-successor Markov chain: with probability 1/2 stay, otherwise
    step to a uniformly chosen successor. Sinks absorb.

    The laziness matters. The compose/decompose graph is bipartite -- every
    rewrite flips the parity of the ternary-edge count -- so the plain walk is
    periodic and power iteration oscillates forever without converging. The
    lazy chain is aperiodic and has the same stationary distribution.

    Returns (pi, max_abs_net_flux, total_abs_net_flux) where the net flux on
    an edge is J(u,v) = pi(u)P(u,v) - pi(v)P(v,u). Detailed balance means
    every J is zero. Note that flux is only informative for an *irreducible*
    chain: on an absorbing chain all mass reaches the terminal components and
    the stationary flux is trivially zero.
    """
    n = len(keys)
    if n == 0:
        return {}, 0.0, 0.0
    pi = {k: 1.0 / n for k in keys}

    for _ in range(iters):
        nxt = defaultdict(float)
        for u in keys:
            succ = graph.get(u, set())
            if not succ:
                nxt[u] += pi[u]          # absorbing
            else:
                nxt[u] += 0.5 * pi[u]    # lazy self-loop
                share = 0.5 * pi[u] / len(succ)
                for v in succ:
                    nxt[v] += share
        delta = sum(abs(nxt[k] - pi[k]) for k in keys)
        pi = {k: nxt.get(k, 0.0) for k in keys}
        if delta < tol:
            break

    def P(u, v):
        succ = graph.get(u, set())
        return 0.5 / len(succ) if v in succ else 0.0

    mx = tot = 0.0
    for u in graph:
        for v in graph[u]:
            J = pi.get(u, 0.0) * P(u, v) - pi.get(v, 0.0) * P(v, u)
            mx = max(mx, abs(J))
            tot += abs(J)
    return pi, mx, tot / 2.0


# ---------- Experiment -------------------------------------------------


def run(seed, name, spec, max_states=None):
    states, graph, edge_rules = continuation_graph(
        seed, rules=spec["rules"], max_states=max_states)
    keys = list(states)

    comps = sccs(graph, keys)
    term = terminal_components(graph, comps)
    sizes = sorted((len(c) for c in comps), reverse=True)
    sinks = [k for k in keys if not graph.get(k)]
    R = recoverability(graph, keys)
    pi, mx_flux, tot_flux = stationary_and_flux(graph, keys)

    pots = {p: classify_potential(states, graph, fn)
            for p, fn in POTENTIALS.items()}
    lyap = [p for p, v in pots.items() if v == "lyapunov"]

    mean_R = sum(R.values()) / len(R) if R else 0.0

    print(f"  {name}")
    print(f"    states={len(states)}  transitions={sum(len(v) for v in graph.values())}"
          f"  symmetric={is_symmetric(graph)}")
    print(f"    SCCs={len(comps)}  largest={sizes[0] if sizes else 0}"
          f"  terminal={len(term)}  sinks={len(sinks)}")
    absorbing = len(term) < len(comps)
    flux_note = "" if absorbing else "  (irreducible: flux is meaningful)"
    print(f"    mean R={mean_R:.3f}   max|flux|={mx_flux:.3e}  total|flux|={tot_flux:.3e}"
          f"{flux_note}")
    print(f"    potentials={pots}")
    print(f"    Lyapunov: {lyap if lyap else 'NONE -- no ordering on states'}")
    if len(term) == 1 and len(term[0]) == len(keys):
        print("    vacuum: DEGENERATE -- terminal component is the entire "
              "state space")
    elif term:
        vac = min(term, key=lambda c: (len(c), sorted(c)))
        rep = sorted(vac)[0]
        print(f"    vacuum (smallest of {len(term)} terminal components, "
              f"size {len(vac)}): {rep}")
    print(f"    predicted: {spec['prediction']}")
    print()
    return states, graph


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--path", type=int, default=4, metavar="N",
                    help="seed with a path of N binary edges (default: 4)")
    ap.add_argument("--max-states", type=int, default=None, metavar="M",
                    help="cap distinct states per rule set")
    ap.add_argument("--rules", choices=sorted(RULE_SETS), default=None,
                    help="run only one rule set (default: all)")
    args = ap.parse_args(argv)

    seed = path_seed(args.path)
    print("=" * 68)
    print(f"seed = path({args.path})   nodes={len(seed.nodes())}  weight={weight(seed)}")
    print("=" * 68)
    print()

    names = [args.rules] if args.rules else sorted(RULE_SETS)
    for name in names:
        run(seed, name, RULE_SETS[name], max_states=args.max_states)


if __name__ == "__main__":
    main()
