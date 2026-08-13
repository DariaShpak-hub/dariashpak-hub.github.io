# The quantities nobody went back for

Companion to `ORIGINS.md`. If Fork A and Fork B are one operation — removing the
external stage — applied to different variables, then the reassignments were
made **one variable at a time, by hand, over a century**. Simultaneity in 1905.
The time coordinate and the energy zero in 1915. Particle number in 1976, sixty
years late, and only because Unruh happened to ask.

That raises a question the fork analysis does not answer: **which stage-relative
quantities were never reassigned at all?**

Those are places where a contradiction may be sitting undetected — not because
it is subtle, but because nobody looked.

---

## The test

A quantity is **stage-relative** if two people who agree on everything inside
the system can disagree about its value. Something outside must be specified: a
frame, a time origin, an observer, a basis, a partition, a cutoff, an ensemble.

For each such quantity: was it reassigned when the stage was removed, when, and
what happened?

---

## The audit

### Reassigned cleanly — these are the templates

| quantity | outside thing needed | resolution | year |
|---|---|---|---|
| Simultaneity | a frame | became frame-relative | 1905 |
| Acceleration | absolute space (Newton's bucket) | local inertial frame | 1915 |
| Gauge potential A_μ | a gauge choice | A_μ unphysical, **holonomy** physical | 1959 (Aharonov–Bohm) |
| Renormalization scale μ | an arbitrary cutoff | μ-independence, **RG invariants** physical | ~1971 (Wilson) |
| Quantum phase | a phase convention | global phase unphysical, **Berry phase** physical | 1984 |

Note the shared mechanism in the last three. In each case the arbitrary choice
was not eliminated — it was kept, and an **invariant of the choice** was
identified and promoted to the observable. Aharonov–Bohm did it with a loop
integral, Wilson with a flow, Berry with a geometric phase.

### Reassigned catastrophically

| quantity | what happened |
|---|---|
| Energy zero | convention → absolute source. Fork B. 10¹²⁰. |

### Reassigned, sixty years late

| quantity | resolution | year |
|---|---|---|
| Particle number / "what is a particle" | observer-relative; no preferred vacuum in curved spacetime | 1976 (Unruh) |
| Temperature | observer-relative | 1976 |
| The vacuum | observer-relative | 1976 |

Sixty years is the interesting number. The stage was removed in 1915 and the
question "relative to whom is this state empty?" was not asked until 1976. There
is no reason to think 1976 was the last of them.

### Never reassigned

Six, and this is the list that matters.

**1. Rotation.** Newton's bucket makes rotation absolute. Mach argued it should
be relative to the total matter distribution. General relativity **does not
implement Mach's principle** — Gödel's rotating universe (1949) is a solution,
and asymptotically flat spacetimes permit rotation relative to infinity in an
otherwise empty universe. Frame dragging is partial Mach and nobody claims it is
the whole thing. The reassignment was attempted, is incomplete, and stalled.

**2. The coarse-graining that defines entropy.** Entropy in statistical
mechanics is relative to a macrostate partition chosen from outside (Jaynes was
explicit that it measures missing information relative to a description). Then
Bekenstein–Hawking assigns a horizon an **absolute** entropy, A/4. That is Fork
B applied to entropy — a quantity defined only up to a reference choice became
an absolute physical number. The information paradox is what happens where the
two meet.

**3. The preferred basis.** Probability in quantum mechanics is defined relative
to a basis chosen from outside. Decoherence and einselection explain why some
bases are stable — but only given a system/environment split, which is item 4.
The problem was converted, not solved.

**4. The subsystem factorisation — where to cut.** Entanglement entropy depends
on how the Hilbert space is factored, which is an outside choice. In gauge
theories the space **provably does not factorise** across a spatial region
because of the Gauss constraint, and different centre choices give different
entropies (Casini–Huerta–Rosabal 2014). In holography this is the factorisation
problem and it is open.

**5. Local energy density of the gravitational field.** No local tensor exists —
only pseudotensors that depend on coordinates. This is Noether's second theorem
stated physically. Quasi-local proposals (Hawking, Brown–York, Wang–Yau) are
partial and none is accepted.

**6. What counts as "the system" in a gravitating universe.** Asymptotically
flat boundary conditions require an outside. Cosmology has no outside. This is
why there is no S-matrix in de Sitter space, and it is live.

---

## The pattern in the six

Sort the successes against the failures and the division is clean, and it is not
the division anyone would have guessed.

**Everything successfully reassigned was a *labelling* freedom.**
Frame, gauge, phase, scale, coordinate. A group acts on the choice. You can take
a quotient, or find an invariant of the group action and promote it.

**Everything never reassigned is a *partitioning* freedom.**
Where to cut, what to ignore, what to distinguish, what is inside. Rotation
relative to *which* matter. Entropy relative to *which* coarse-graining. Energy
in *which* region.

> There is no group acting on "where to draw the boundary."

Which means the entire mathematical technology physics built for handling
stage-relativity — take the quotient, find the invariant — **does not apply to
half the list.** Gauge theory, the renormalization group, and the theory of
geometric phases are all one technique, and that technique has no purchase on
items 1 through 6.

### The objection, and the refinement

The obvious counterexample: **the renormalization group is a coarse-graining**,
and coarse-graining is a partition. It works beautifully. So why doesn't the
technique transfer?

Because the RG's partition is indexed by **scale**, and scale carries a group —
dilatations. The partition is *generated by* a labelling freedom, which is
exactly why an invariant exists. Cut the modes by momentum shell and there is a
one-parameter family to flow along. Cut a spatial region out of a lattice and
there is no group taking one region to another that acts on the factorisation.

So the refined statement:

> A stage-relative quantity is tractable exactly when the choice is generated by
> a group action that is **not transitive** on the physical states.

Both conditions are load-bearing. If there is no group, there is no invariant —
that is items 1 through 6. If the group is transitive, the quotient is trivial
and carries no information — that is the preferred basis, where U(n) acts
transitively on bases, so quotienting by it destroys everything rather than
resolving anything.

A symmetry that is too large is as useless as none. That is why item 3 is hard
in a way gauge freedom never was, and it is not a matter of effort.

---

## Why the observer branch stalled

The six unreassigned quantities are almost all on the **observer** side of
Fork A. Item 5 is on the time side; the rest are about where an observer draws
lines.

The time side of Fork A got a century of concentrated work — the problem of
time, Wheeler–DeWitt, and all fifteen programmes in `TIMELESS_MODELS.md`. The
observer side got Everett, got decoherence, and largely stopped.

That asymmetry is now explained rather than merely observed. The time side is a
**labelling** problem — reparameterisation is a group action, which is why
people could get traction, and why shape dynamics could prove an exact trade.
The observer side is a **partitioning** problem, and the tools do not exist.

It stalled because it was the half that the available mathematics could not
touch, not because it was less important.

---

## Which of the six this project is already standing in

Three of the six are live in code we have already run, and two of those we have
already measured as decisive.

**Item 2, coarse-graining — measured, and it changes the answer.** Every
dimension estimator we built required choosing a coarse-graining, and shell
versus ball estimators differ by about **20%** on the same structure. The
partition choice is not a detail; it is a fifth of the result.

**Item 6, the boundary — measured, and it flips the sign.** The only persistent
arrow we ever produced came from a receding ceiling on an **open** substrate,
giving a constant 16.6% lag, while every smooth start on a **closed** substrate
gave a *negative* arrow. Same dynamics, opposite conclusion, decided entirely by
what counts as the system.

**Item 4, the factorisation — the one place we got it for free.** In
`necessity.py` the 6944 sectors were not imposed. The partition of the state
space fell out of the move set and then turned out to be exactly the degree
vectors. Nobody chose where to cut.

That last one needs its caveat stated plainly: a partition of a state space is
**not** a tensor factorisation of a Hilbert space. They are different objects and
the resemblance may be superficial. But it is the only instance in this
repository of a partition determined from inside the structure, and it is
testable: does the partition stay intrinsic under a different move set, or was
it an artifact of the one rule I happened to pick?

That is a well-posed experiment on an unreassigned quantity, and it is cheap.

---

## What this list is for

Not a claim that six new contradictions exist. It is a claim about where to
look. The Unruh case is the precedent: a stage-relative quantity sat unexamined
for sixty years, and when someone finally asked, the answer was that the
quantity had no observer-independent value at all.

Five of the six above have never had that question asked with the same
seriousness, and the reason is now identifiable — the technique that worked
everywhere else structurally cannot work on them.
