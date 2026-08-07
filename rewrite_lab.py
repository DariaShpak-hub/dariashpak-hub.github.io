#!/usr/bin/env python3
"""
rewrite_lab.py

Framework for exploring local rewrite grammars on finite unlabeled hypergraphs.

This is a research scaffold rather than a finished physics model.
You provide one or more rewrite rules and the program enumerates the
reachable continuation graph up to a chosen size limit.

Notes on the model
------------------
Both bundled rules preserve the "weight" of a hypergraph
(``sum(len(e) - 1 for e in edges)``): ``compose`` fuses two binary edges
sharing one node into a single ternary edge, and ``decompose`` does the
reverse. Weight is therefore an invariant of the whole continuation graph,
which the self-check under ``__main__`` verifies.
"""

import argparse
from itertools import permutations
from collections import deque, defaultdict

# ---------- Hypergraph -------------------------------------------------


class Hypergraph:
    def __init__(self, edges=()):
        self.edges = frozenset(frozenset(e) for e in edges)
        # Canonical form is expensive (O(n!)); compute it lazily once and
        # cache it, since it is queried from __hash__, __eq__, __repr__ and
        # the search loop.
        self._canonical = None

    def nodes(self):
        out = set()
        for e in self.edges:
            out |= set(e)
        return out

    def canonical(self):
        """
        Canonical representative under node relabeling.

        Returns a hashable, order-independent value (a sorted tuple of
        sorted edge tuples) that is identical for isomorphic hypergraphs.
        The empty hypergraph canonicalizes to the empty tuple.

        Complexity is O(n!) in the number of nodes, so this is only suitable
        for the small graphs a rewrite scaffold explores.
        """
        if self._canonical is not None:
            return self._canonical

        nodes = sorted(self.nodes())
        if not nodes:
            # Always return the same hashable value type. Returning ``self``
            # here (as an earlier version did) caused infinite recursion in
            # __hash__/__repr__ for the empty graph.
            self._canonical = ()
            return self._canonical

        # Relabel onto the fixed alphabet 0..n-1 and take the lexicographically
        # smallest image over all such bijections. Mapping to a fixed alphabet
        # (rather than to a permutation of the graph's own labels) makes the
        # result a true isomorphism invariant: two isomorphic graphs with
        # different surviving labels -- which rules routinely produce -- get the
        # same canonical form, matching the "unlabeled hypergraph" model.
        best = None
        for p in permutations(range(len(nodes))):
            m = dict(zip(nodes, p))
            img = tuple(sorted(tuple(sorted(m[v] for v in e))
                               for e in self.edges))
            if best is None or img < best:
                best = img

        self._canonical = best
        return best

    def __hash__(self):
        return hash(self.canonical())

    def __eq__(self, other):
        if not isinstance(other, Hypergraph):
            return NotImplemented
        return self.canonical() == other.canonical()

    def __repr__(self):
        return str(self.canonical())


# ---------- Rewrite rules ----------------------------------------------


class Rule:
    """
    Rule callback returns iterable of rewritten Hypergraph objects.
    """

    def __init__(self, name, fn):
        self.name = name
        self.fn = fn

    def apply(self, H):
        yield from self.fn(H)


def compose_binary_to_ternary(H):
    """{x,y},{y,z} -> {x,y,z}"""
    edges = list(H.edges)
    seen = set()

    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            a = edges[i]
            b = edges[j]
            if len(a) == len(b) == 2:
                inter = a & b
                if len(inter) != 1:
                    continue
                new = frozenset(a | b)
                if new in H.edges:
                    continue
                new_edges = set(H.edges)
                new_edges.remove(a)
                new_edges.remove(b)
                new_edges.add(new)
                G = Hypergraph(new_edges)
                c = G.canonical()
                if c not in seen:
                    seen.add(c)
                    yield G


def decompose_ternary(H):
    """{x,y,z} -> any two binary edges sharing one node"""
    seen = set()

    for e in list(H.edges):
        if len(e) != 3:
            continue

        verts = list(e)

        for middle in verts:
            others = [v for v in verts if v != middle]
            e1 = frozenset((middle, others[0]))
            e2 = frozenset((middle, others[1]))

            new_edges = set(H.edges)
            new_edges.remove(e)

            if e1 in new_edges or e2 in new_edges:
                continue

            new_edges.add(e1)
            new_edges.add(e2)

            G = Hypergraph(new_edges)
            c = G.canonical()
            if c not in seen:
                seen.add(c)
                yield G


RULES = [
    Rule("compose", compose_binary_to_ternary),
    Rule("decompose", decompose_ternary),
]

# ---------- Search -----------------------------------------------------


def continuation_graph(seed, rules=RULES, max_states=None):
    """
    Enumerate the reachable states from ``seed`` under ``rules``.

    Parameters
    ----------
    seed : Hypergraph
        Starting state.
    rules : list[Rule]
        Rewrite rules to apply at every state.
    max_states : int or None
        Stop once this many distinct states have been visited. ``None``
        (the default) enumerates the full reachable set. When the limit is
        hit, transitions leaving already-visited states are still recorded,
        but no new states are expanded.

    Returns
    -------
    visited : dict
        Maps canonical key -> representative Hypergraph.
    graph : dict
        Adjacency: canonical key -> set of successor canonical keys.
    edge_rules : dict
        Maps (src_key, dst_key) -> set of rule names that produce that
        transition, preserving rule provenance.
    """
    queue = deque([seed])
    visited = {}
    graph = defaultdict(set)
    edge_rules = defaultdict(set)

    while queue:
        H = queue.popleft()
        key = H.canonical()

        if key in visited:
            continue

        if max_states is not None and len(visited) >= max_states:
            # Reached the cap: do not register or expand new states.
            continue

        visited[key] = H

        for rule in rules:
            for G in rule.apply(H):
                gkey = G.canonical()
                graph[key].add(gkey)
                edge_rules[(key, gkey)].add(rule.name)
                if gkey not in visited:
                    queue.append(G)

    return visited, graph, edge_rules


def weight(H):
    return sum(len(e) - 1 for e in H.edges)


def outdegree(graph):
    return {k: len(v) for k, v in graph.items()}


def path_seed(length):
    """Path hypergraph with ``length`` binary edges: {0,1},{1,2},..."""
    return Hypergraph([{i, i + 1} for i in range(length)])


# ---------- CLI --------------------------------------------------------


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument(
        "--path", type=int, default=4, metavar="N",
        help="seed with a path of N binary edges (default: 4)")
    parser.add_argument(
        "--max-states", type=int, default=None, metavar="M",
        help="stop after visiting M distinct states (default: no limit)")
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    seed = path_seed(args.path)

    states, graph, edge_rules = continuation_graph(
        seed, max_states=args.max_states)

    # Weight is invariant under both bundled rules; verify it.
    seed_w = weight(seed)
    assert all(weight(H) == seed_w for H in states.values()), \
        "weight is not invariant under the active rule set"

    print("=" * 60)
    print("Reachable states:", len(states))
    print("Transitions:", sum(len(v) for v in graph.values()))
    print("Invariant weight:", seed_w)
    print("=" * 60)

    for k in sorted(states):
        H = states[k]
        print("State:", k)
        print(" Weight:", weight(H))
        print(" Outdegree:", len(graph[k]))
        print()

    print("Adjacency")
    for k in sorted(graph):
        for dst in sorted(graph[k]):
            rules = ",".join(sorted(edge_rules[(k, dst)]))
            print(k, "--[", rules, "]-->", dst)


if __name__ == "__main__":
    main()
