"""Is the per-component (nodes, weight) multiset the complete charge?"""
import sys
sys.path.insert(0, "/home/user/dariashpak-hub.github.io")
sys.path.insert(0, "/tmp/claude-0/-home-user-dariashpak-hub-github-io/31ff9562-0d3c-5826-99eb-bd71006a2e6c/scratchpad")
from closed import enumerate_fixed, RULES
from charges import components

def comp_profile(H):
    """Multiset of (node count, weight) over connected components."""
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
    nodes, wt = {}, {}
    for v in H.nodes():
        nodes[find(v)] = nodes.get(find(v), 0) + 1
    for e in H.edges:
        r = find(next(iter(e)))
        wt[r] = wt.get(r, 0) + len(e) - 1
    return tuple(sorted((nodes[r], wt.get(r, 0)) for r in nodes))

def n_comp(H):
    return len(comp_profile(H))

print("%3s %3s %7s %8s %9s %9s %s" % (
    "n","w","states","sectors","by n_comp","by profile","complete?"))
for n in range(4, 8):
    for w in range(3, 8):
        st = enumerate_fixed(n, w)
        if len(st) < 2: continue
        comps = components(st, RULES)
        # is each invariant constant on sectors, and does it separate them?
        prof_vals, ncomp_vals = [], []
        ok_prof = ok_ncomp = True
        for c in comps:
            ps = {comp_profile(st[k]) for k in c}
            ns = {n_comp(st[k]) for k in c}
            if len(ps) > 1: ok_prof = False
            if len(ns) > 1: ok_ncomp = False
            prof_vals.append(next(iter(ps))); ncomp_vals.append(next(iter(ns)))
        sep_prof = ok_prof and len(set(prof_vals)) == len(comps)
        sep_ncomp = ok_ncomp and len(set(ncomp_vals)) == len(comps)
        print("%3d %3d %7d %8d %9s %9s %s" % (
            n, w, len(st), len(comps),
            "sep" if sep_ncomp else ("cons" if ok_ncomp else "NOT"),
            "sep" if sep_prof else ("cons" if ok_prof else "NOT"),
            "COMPLETE" if sep_prof else "incomplete"))
