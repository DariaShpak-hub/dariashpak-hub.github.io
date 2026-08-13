"""What IS the hidden charge? Test candidate invariants against the sectors."""
import sys
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
sys.path.insert(0, "/tmp/claude-0/-home-user-dariashpak-hub-github-io/31ff9562-0d3c-5826-99eb-bd71006a2e6c/scratchpad")
from collections import defaultdict
from closed import enumerate_fixed, RULES
from charges import components

def n_components(H):
    parent = {v: v for v in H.nodes()}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e in H.edges:
        vs = list(e)
        for v in vs[1:]:
            a, b = find(vs[0]), find(v)
            if a != b: parent[a] = b
    return len({find(v) for v in H.nodes()})

def cycle_rank(H):
    # weight - nodes + components  (generalised first Betti number)
    return sum(len(e)-1 for e in H.edges) - len(H.nodes()) + n_components(H)

def max_deg(H):
    d = defaultdict(int)
    for e in H.edges:
        for v in e: d[v] += 1
    return max(d.values()) if d else 0

CANDS = {"n_components": n_components, "cycle_rank": cycle_rank,
         "n_edges": lambda H: len(H.edges), "max_degree": max_deg}

for (n, w) in [(5,4),(5,5),(5,6),(6,4),(6,5),(6,6)]:
    st = enumerate_fixed(n, w)
    comps = components(st, RULES)
    print("n=%d w=%d : %d states, %d sectors" % (n, w, len(st), len(comps)))
    for name, fn in CANDS.items():
        vals = []
        constant = True
        for comp in comps:
            s = {fn(st[c]) for c in comp}
            if len(s) > 1: constant = False
            vals.append(sorted(s))
        if constant:
            flat = [v[0] for v in vals]
            sep = "SEPARATES" if len(set(flat)) == len(flat) else "constant but degenerate"
            print("   %-14s conserved, values=%s  %s" % (name, sorted(flat), sep))
        else:
            print("   %-14s NOT conserved" % name)
    print()
