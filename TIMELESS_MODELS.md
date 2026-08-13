# Timeless models: where each one put the time

A survey of programmes that claim to eliminate fundamental time, sorted by what
they actually did with it. The question is not whether they are good physics —
several are excellent — but whether the time was **removed**, **renamed**, or
**relocated**.

This is a companion to the self-audit of our own modules. The same diagnostic is
applied to everyone, including us.

---

## The diagnostic

Five questions. A programme that answers "no" to all five has genuinely removed
time. Nothing in this survey answers "no" to all five.

1. **Monotonic variable.** Is there a quantity that increases along every
   solution and is used to order states? Then it is a clock, whatever it is
   called.
2. **One-parameter flow.** Is there a one-parameter group of automorphisms —
   unitary, modular, renormalisation — whose parameter labels states? Then the
   structure of time is present under another name.
3. **Asymmetric boundary.** Does the arrow come from a condition imposed at one
   end? Then the direction is stipulated, not derived.
4. **Asymptotics.** Does the framework require "in" and "out" states? Then local
   time was removed and global time was kept.
5. **Derivation order.** Does the construction proceed in steps that are
   *afterwards* declared unphysical? A quotient taken at the end is not the same
   as never having had the structure.

Question 5 is the one that catches the most, and it is the one that caught us.

---

## Class I — Time excluded, then re-derived from within

The equations have no t. Then t reappears, because without it there are no
predictions.

### Wheeler–DeWitt

ĤΨ[h_ij, φ] = 0. Nothing evolves. The universe's wavefunction is a single
timeless object on superspace.

**Where the time went.** Semiclassical WKB time (Banks 1985). Expand in inverse
Planck mass; the heavy gravitational sector carries a phase e^{iS/ħ}, and t is
*defined* by the gradient of S along which the light matter sector obeys a
Schrödinger equation. So time is recovered — but only after picking a single WKB
branch rather than a superposition of them, and only in the regime where the
expansion is valid. Fails 1 and 5.

There is a second, less discussed cost: with no time there is no inner product
and hence no probability interpretation, which is why "the problem of time"
(Isham 1992, Kuchař 1992) is a cluster of problems rather than one.

**Verdict:** excluded, then imported back through an approximation scheme. The
timelessness holds exactly where the theory cannot be used.

### Hartle–Hawking no-boundary

Euclidean path integral over compact geometries with no boundary. There is no
initial moment because there is no edge.

**Where the time went.** Three places. The signature: imaginary time is still a
coordinate direction, and the Euclidean circle's circumference is β = 1/kT —
which means the "no time" formulation is a *thermal* one, and thermal
equilibrium is defined by analyticity in time (the KMS condition). The contour:
the answer depends on which integration contour is chosen, and
Feldbrugge–Lehners–Turok showed the Lorentzian version with the natural contour
is unstable. And the vocabulary: "initial condition" is unavoidable because the
scale factor is being used as the ordering variable.

**Verdict:** excluded from the manifold, retained in the signature and the
boundary choice. Fails 1 and 3.

### Barbour's Platonia and best matching

The most serious attempt on the list. Configuration space (shape space) with no
external time; dynamics from Jacobi's principle, which is genuinely
reparameterisation invariant — the parameter along the curve is arbitrary and
carries no physics. "Time is the measure of change," not a container.

**Where the time went.** Two residues, one small and one large.

Small: a geodesic in shape space is still a **curve** — a one-dimensional
ordered set of configurations. The parameter is arbitrary; the ordering is not.
Reparameterisation invariance removes the units, not the sequence.

Large: the "time capsules" proposal — that we experience time because the
wavefunction concentrates on configurations that contain apparent records — is
a conjecture about the solution, not a result. It requires the state to be
special. And the later Janus-point work with Koslowski and Mercati derives an
arrow from monotonically growing **complexity** away from a central point, which
is a monotonic variable in the sense of question 1.

**Verdict:** the closest anyone gets. Passes 2, 3, 4. Fails 1 in its arrow-
generating extension.

---

## Class II — Time replaced by another time

The most common outcome by a wide margin. A parameter is removed and a
different parameter, with the same formal role, is installed.

### Page–Wootters conditional probabilities

Split the universe into system and clock, impose (Ĥ_S + Ĥ_C)|Ψ⟩ = 0, and define
|ψ(t)⟩_S = ⟨t|Ψ⟩ where |t⟩ is a clock pointer state. Schrödinger evolution falls
out of a static global state. Genuinely elegant.

**Where the time went.** Into the clock's Hamiltonian. For the pointer states to
be distinguishable and the recovered evolution to be correct, the clock's
spectrum must be unbounded — an ideal clock, which no physical system is.
Kuchař (1992) showed the original construction gives wrong two-time correlation
functions; Höhn, Smith and Lock (2019–2021) repaired most of this and showed
equivalence with relational Dirac observables. The repair is real. It does not
change what the clock is.

**Verdict:** time replaced by the pointer of a subsystem stipulated to be a
perfect clock. Fails 1 and 2.

### Thermal time hypothesis (Connes–Rovelli 1994)

Given a state ω on a von Neumann algebra, Tomita–Takesaki construction gives a
canonical one-parameter automorphism group σ_t^ω — the modular flow. The
proposal: physical time **is** this flow, determined by the state rather than by
the geometry. This is the deepest idea in the survey, because it explains *why
this time and not another* and ties time to thermodynamics through the KMS
condition.

**Where the time went.** It is a one-parameter automorphism group with parameter
t. That is what a time is. The hypothesis answers "whose time" — a genuine
advance — and does not touch "whether time". And the state must be chosen; a
different state gives a different flow, so state selection now does the work
that foliation choice did before.

**Verdict:** the clearest possible fail on question 2, and the most valuable
programme on the list anyway. Those are compatible.

### Shape dynamics

Trades refoliation invariance for three-dimensional conformal invariance,
yielding a theory with a genuine global Hamiltonian and a preferred
simultaneity: York time, τ = (2/3)K.

**Where the time went.** Nowhere hidden — this is the honest one. It is a
theorem-like statement of the pattern: you can trade **many-fingered** time for
**one global** time, and the trade is exact. Shape dynamics is the best evidence
in the literature that time is conserved under reformulation.

**Verdict:** explicit. Fails 1 by design and says so.

### Group field theory condensate cosmology

Relational dynamics with a massless scalar field φ as the clock, chosen
precisely because it is monotonic along solutions.

**Verdict:** a clock by construction. Fails 1.

### Rovelli's partial observables

Dynamics as correlations between partial observables, none preferred. No
evolution *in* anything; just a constraint surface and its correlations.
Deflationary and defensible.

**Where the time went.** To extract anything resembling evolution you designate
one variable as the clock, and the correlation only has evolutionary structure
if that variable is monotonic along the trajectory. Generic constrained systems
have no globally monotonic variable (Hájíček's global problem of time), so the
construction works locally and fragments globally.

**Verdict:** honest framework, local clocks. Fails 1 locally, and where it does
not fail 1 it produces no dynamics at all.

### Holographic emergence (AdS/CFT, Ryu–Takayanagi, tensor networks)

Bulk geometry emerges from boundary entanglement; the radial direction emerges
from RG scale.

**Where the time went.** Bulk time does not emerge. It **is** the boundary time.
Space is derived, time is imported unchanged. This gets blurred constantly in
summaries of the programme, and the blur matters: "spacetime is emergent" is
true for three of the four dimensions.

**Verdict:** fails 2 (RG scale is a one-parameter flow) and does not attempt the
time question at all.

---

## Class III — Time replaced by something else, and smuggled elsewhere

### Causal set theory

Kinematics: a locally finite partial order. Malament's theorem gives the
conformal metric from the causal structure; volume from counting; "order +
number = geometry." As a **static** structure this is genuinely time-free — a
partial order is not a time function, and there is no parameter.

**Where the time went.** Into the dynamics. Rideout–Sorkin classical sequential
growth builds the set one element at a time, labelled 1…n, and then imposes
"discrete general covariance": the measure must not depend on the labelling. The
labels vanish from the answers. But the object constructed is a Markov chain,
and a Markov chain is indexed by a step count. Sorkin's own formulation — "the
growth is not in time, time is in the growth" — is the claim, stated clearly,
not a derivation.

**Verdict:** kinematics passes all five. Dynamics fails 5, exactly. This is the
purest example of question 5 in the literature.

Our own `causal_order.py` is a direct test of the same move: the growth poset
our rewriting produces has ordering fraction falling as 1/n (0.0354 → 0.0061 as
n runs 600 → 4801) while Minkowski holds at 0.50. The birth sequence produced a
**tree**, not a causal structure. So in our case the smuggled order did not even
buy the geometry it was smuggled for.

### Causal dynamical triangulations

The most informative empirical result in this entire area, and it points the
wrong way for everyone here.

Euclidean dynamical triangulations — no causal structure imposed — produces
degenerate phases: crumpled with effectively infinite dimension, or branched
polymer with d ≈ 2. Adding a preferred foliation, i.e. **putting causality in by
hand**, produces a phase with an extended four-dimensional de Sitter-like
geometry.

**Verdict:** this is a measurement, not a philosophical position. Somebody asked
"can you get geometry without putting time in?" and the answer came back no.

Our repository reproduced the shape of that answer independently and without
meaning to: every module that produced geometry had smuggled time or a
background lattice, in ten attempts, and the clean modules produced mechanics
only.

### Spin foams

A timeless sum over 2-complexes with no external evolution parameter.

**Where the time went.** A spin foam is a cobordism between an incoming and an
outgoing spin network — a structure with two distinguished ends. Fails 4.

### Amplituhedron and positive geometry

Scattering amplitudes as volumes of a geometric object, with locality and
unitarity emergent rather than assumed. No Lagrangian, no time evolution
anywhere in the formalism.

**Where the time went.** The object computed is an S-matrix, which presupposes
asymptotically free states at t → ∓∞. Local time is removed completely; global
time is retained maximally. Fails 4 as hard as anything can.

### Decoherent histories

A framework for closed-system quantum mechanics with no external observer.

**Where the time went.** A history is a time-ordered sequence of projectors, and
the decoherence functional is defined on those sequences. Time is the
scaffolding the framework is built on. Fails 1 and 5.

### The Past Hypothesis (Boltzmann, Albert, Price)

The arrow is not in the laws; it comes from a low-entropy boundary condition at
one end of the universe.

**Where the time went.** Into the word "past". The hypothesis explains why
entropy increases *toward one end*. It does not explain why that end is called
earlier. Price's objection is exact and has never been answered: the
construction is time-symmetric and the asymmetry is in the labelling.

**Verdict:** fails 3 in the definitional sense. Which is not a criticism of the
physics — it is the physics correctly identifying that the arrow is a boundary
condition.

This one matches our own strongest result. The only thing that ever produced a
persistent arrow in our system was a **receding ceiling** — Ω growing faster
than the system fills it — giving a constant 16.6% lag on an open substrate,
while every smooth start on a closed substrate gave a *negative* arrow. A
boundary condition, not a law. We arrived at the Past Hypothesis from the other
side and by measurement.

### Constructor theory

Statements about which transformations are possible and impossible, rather than
about trajectories. A real shift away from dynamics as the basic mode of
description.

**Where the time went.** A task is an ordered pair — an input attribute and an
output attribute. That arrow is minimal but real. And the programme has not yet
reproduced known dynamics, so the timelessness has not been paid for.

**Verdict:** fails 5 minimally. Unfinished rather than wrong.

### Wolfram model

Hypergraph rewriting; "causal invariance" means the causal graph is independent
of the order in which updates are applied.

**Where the time went.** Time is the update count, explicitly. Causal invariance
is confluence — a statement that the *result* does not depend on the order, not
that there is no order. And branchial space is a space of simultaneous states,
which requires a simultaneity, i.e. a slicing.

**Verdict:** fails 1 and 5 openly. Worth naming because it is frequently
described as the most timeless of the programmes and is the least. Our own
rewriting model inherits this hazard directly and it is why the audit mattered.

### 't Hooft cellular automaton interpretation

Explicit discrete time step. Not a timeless model; included so the list is not
accused of grading on a curve.

---

## Class IV — Genuinely time-free, and what it cost

| structure | passes all five? | what it produces | what it cannot |
|---|---|---|---|
| Causal set **kinematics** (order + number, no growth) | yes | conformal geometry, volume | no dynamics |
| Barbour's shape space, before the arrow question | yes | relational configurations | no arrow |
| Constructor theory tasks | nearly | possibility claims | no known dynamics recovered |
| Our `necessity.py` (relation, alternative, consistency) | yes | conservation, multiplicity, drift, 8.85 nats of persistent information, an effective interaction | no metric, no dimension, no arrow |

The pattern across the table is the survey's actual result:

> **Every structure that is genuinely time-free produces no dynamics and no
> arrow. Every structure that produces dynamics has a parameter.**

---

## The five mechanisms, and how often each is used

| mechanism | used by |
|---|---|
| **Monotonic variable** renamed | GFT scalar clock, York time, complexity (Janus), Wolfram update count, WKB phase |
| **One-parameter flow** — same structure, different word | modular flow (thermal time), unitary group (Page–Wootters), RG scale (holography) |
| **Boundary condition at one end** | Past Hypothesis, no-boundary proposal, our receding ceiling |
| **Asymptotic in/out** | S-matrix, amplituhedron, spin foam cobordism |
| **Derivation order, quotiented afterwards** | causal set sequential growth, decoherent histories, and our own equal-N comparisons |

---

## The deflation this survey has to survive

There is an obvious objection to the headline, and it should be stated before
someone else states it: **"dynamics" may just be *defined* as a one-parameter
family of states.** If so, "no dynamics without a parameter" is analytic, not
empirical, and the survey has discovered a tautology dressed as a pattern.

That objection is partly right, and it sharpens the question rather than
dissolving it. What survives it:

1. The **CDT result is not analytic.** It is a computation, and it says that
   without imposed causality you get a crumpled phase or a branched polymer.
   Nothing about the definition of dynamics predicts *which* degenerate phase.
2. The **shape dynamics trade is not analytic.** That refoliation invariance can
   be exchanged for exactly one global time is a specific structural result, not
   a restatement of what evolution means.
3. **Our own asymmetry is not analytic.** Restriction alone gave conservation,
   multiplicity, statistical drift, persistent information and an effective
   interaction — five pieces of mechanics — with no parameter anywhere. That is
   already a counterexample to "mechanics requires time." What it did not give
   is a **metric, a dimension, or an arrow**. So the true division is not
   dynamics-versus-no-dynamics. It is:

> Mechanics is available without time. Geometry and direction are not — in
> fifteen independent attempts by unrelated communities, including ours.

That is the claim worth pursuing, and it is not a tautology, because nothing in
the definition of "metric" mentions a parameter.

## What would settle it

A single construction that yields a **distance function** and a **preferred
direction** from a structure passing all five diagnostic questions. No entry in
this survey does it. Neither do we.

The negative version — a proof that no such structure exists — would be a real
theorem and would end the programme honestly. This repository has accumulated
evidence pointing that way from several unrelated directions without coming
close to proving it.
