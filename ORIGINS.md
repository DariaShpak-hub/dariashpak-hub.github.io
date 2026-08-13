# Where each contradiction was made

For every disagreement in `MEASUREMENTS.md` and every smuggle in
`TIMELESS_MODELS.md`, the question here is not *what* disagrees but *when the
disagreement became inevitable*. Most of them were built in, decades before the
measurement that exposed them, by a decision that was correct locally and
incompatible globally.

The template is the Minkowski/Hamilton case:

> Minkowski spatialised time, turning it into a frozen coordinate. Hamilton
> operationalised time, turning it into the generator of change. General
> Relativity was built on the first and Quantum Mechanics on the second, which
> is why Wheeler–DeWitt broke.

That is right, and it can be made sharper — sharp enough to be a diagnostic
rather than an analogy. Below, that case first in full, then the same treatment
applied to the rest.

A **fork** is a point where two theories made opposite decisions about the
ontological status of the same thing. A **symptom** is what shows up when the
branches are forced to meet. Most of what gets called an open problem is a
symptom.

---

## Fork A — Something external was required, then abolished

### The worked case: why Wheeler–DeWitt had to break

The chain is tighter than "two different notions of time." Every link is a
theorem or a definition, not an interpretation.

1. **Hamilton (1834).** Time is the parameter of evolution, external to the
   system. H generates translation along it.
2. **Noether (1918), first theorem.** Energy *is defined as* the conserved
   quantity associated with time-translation symmetry. Energy and time are not
   two facts; they are one fact stated twice.
3. **Minkowski (1908) → Einstein (1915).** Time becomes a coordinate, and then
   — via general covariance — a coordinate with **no physical content at all**,
   since any diffeomorphism is a symmetry.
4. **Noether (1918), second theorem.** If the action is invariant under an
   infinite-dimensional group containing translations, the conservation laws
   are *improper*: the divergences vanish identically. Energy conservation in
   GR is not a law, it is an identity — true because of how the theory is
   built, saying nothing about the world.
5. **Pauli (1933), footnote.** If H is bounded below, there is no self-adjoint
   operator conjugate to it. Quantum mechanics **cannot** make time an
   observable. Not "chose not to" — cannot, given that the spectrum is bounded
   below, which is forced by stability.
6. **ADM (1959–62) → DeWitt (1967).** In canonical GR the Hamiltonian is a
   constraint, H ≈ 0, because time translation is a gauge transformation.

Put 1, 5 and 6 together and the collision is a syllogism:

> QM: H generates evolution, and time cannot be an operator.
> GR: H generates gauge.
> Therefore evolution is gauge. Therefore ĤΨ = 0.

Wheeler–DeWitt did not break. It reported, correctly, the consequence of two
prior commitments. Anyone deriving it in 1967 was reading out a result that had
been fixed in 1918.

**The deeper form of the fork.** "Time as parameter" and "time as coordinate" is
the visible layer. Underneath is something more general:

> Hamiltonian mechanics requires an external **stage**. General relativity is
> the theory that abolishes the stage.

That is why the problem has more than one symptom, and why the symptoms look
unrelated until you line them up:

| symptom | the same fork, seen from |
|---|---|
| Problem of time / ĤΨ = 0 | the evolution parameter was external |
| No local energy conservation in GR | energy *was* the time-translation charge (Noether 1 vs Noether 2) |
| No time operator in QM | Pauli — the parameter cannot be internalised |
| The measurement problem | the **observer** was external, for the same structural reason |
| No preferred vacuum (Unruh 1976) | "empty" was defined relative to an external frame |

The measurement problem belongs on this list and is almost never put there. Von
Neumann's two processes — unitary evolution and collapse — need a rule for when
each applies, and the rule mentions "measurement," which the theory does not
define. That is the same shape as needing an external clock: quantum mechanics
structurally requires *something outside the system*, and general relativity is
the theory of the system with no outside.

**Why this fork cannot be resolved by giving up either side.** Both branches are
load-bearing.

- GR's abolition of the stage is not a stylistic preference. Background
  independence is *what made the theory work* — it is the content of the
  equivalence principle and the reason the theory predicts anything Newton
  doesn't.
- QM's external parameter is not a stylistic preference either. Pauli's theorem
  makes it forced, and the premise (spectrum bounded below) is forced by the
  requirement that matter not be able to radiate indefinitely.

So the fork is between two things neither of which can be abandoned cheaply.
This is a much stronger statement than "the theories are incompatible," and it
explains the failure pattern of every attempt at reconciliation: each one gives
up one branch and inherits that branch's original problem.

---

## Fork B — A convention was promoted to an observable

### Λ: the 120 orders of magnitude are not a failed calculation

In every physics before 1915, the **zero of energy was conventional.** Lagrangian
mechanics is invariant under L → L + const. Quantum field theory has this
freedom explicitly: normal ordering discards the vacuum energy, and no
experiment in QFT can see it, because every QFT observable is an energy
*difference*.

Then Einstein put T_μν on the right-hand side. Absolute energy density
gravitates. There is no shift freedom in GR — the zero point is physical.

> A gauge freedom of one theory is a physical observable of the other.

That is the whole vacuum catastrophe. The calculation is not wrong. The
calculation is *meaningless in its own theory* — QFT has no way to define the
quantity GR is asking for, because in QFT that quantity was never observable.
Zel'dovich (1967) was the first to state this quantitatively, half a century
after the fork was made.

This reclassifies the problem. It is usually filed as the worst fine-tuning in
physics. It is more accurately a **well-posedness problem**: one theory is being
asked for a number it does not possess, and 10¹²⁰ is what you get when you
answer anyway using a cutoff you invented.

Note the structural similarity to Fork A: in both, GR takes something that had
been part of the *description* — the coordinate labels, the zero of energy — and
makes it either meaningless or physical. GR is the theory that keeps
renegotiating what is description and what is world, and every renegotiation
left QM on the other side of a line.

---

## Fork C — Two definitions identified while they happened to agree

The most common origin, and the one most likely to be repeated. Pattern: a
quantity gets an **operational** definition (what the apparatus does) and a
**structural** definition (what role it plays in a theory). While precision is
low they coincide. The symbol is shared. Precision improves. They separate. The
shared symbol makes the separation look like a contradiction in the world rather
than in the naming.

### H₀ — 1929 versus 1922

- **Hubble (1929):** H₀ is the slope of a scatter plot of recession velocity
  against distance. An empirical local relation. Fully operational.
- **Friedmann (1922) / Lemaître (1927):** H(t) = ȧ/a is a coefficient in the
  metric of a homogeneous isotropic solution, evaluated at our epoch. Fully
  structural.

These are the same number **only if** the cosmological principle holds and the
model connecting z = 1100 to z = 0 is correct.

What the CMB measures is θ*, an angle, to 0.03% — the best-determined quantity
in cosmology. Converting an angle into an expansion rate requires the sound
horizon r_s and the entire expansion history: H₀ = f(θ*, r_s, ΛCDM). The early
"measurement of H₀" is not a measurement of H₀. It is the value ΛCDM requires.

So the tension has three possible readings and the shared symbol hides which:
a systematic in the ladder, a systematic in the CMB analysis, or the model
between them being wrong. Only the third would be new physics, and only the
third is what the word "tension" implies. The other two are called the same
thing.

This is exactly Fork A's structure repeated: an operationally defined quantity
and a structurally defined quantity, identified while they agreed.

### The neutron lifetime — disappearance versus production

- **Bottle:** counts neutrons that fail to survive. Measures a *disappearance*
  rate.
- **Beam:** counts protons that appear. Measures a *production* rate.

These are the same number only under a **completeness assumption**: that the
neutron decays into exactly the channel we think, with no other exit. Fornal and
Grinberg (2018) pointed out that a ~1% dark decay channel makes both
measurements correct and measuring different things.

So this may not be a contradiction at all. It may be two well-executed
measurements of two different quantities that were given one name because for
sixty years nothing distinguished them. The 5σ would then be a **measurement of
the branching ratio into the unknown channel**, which is the opposite of an
error.

### Inertial mass = gravitational mass → dark matter

The equality of inertial and gravitational mass was, for Newton, an unexplained
empirical coincidence. Einstein promoted it to a principle. Once it is a
principle, an anomalous acceleration has exactly one interpretation available:
unseen source.

The alternative — that the relation between acceleration and source is not what
we wrote — is not ruled out by anything, it is ruled out by the methodological
commitment that **the law is rigid and the content is free.**

That commitment has a track record, and it is one for two on its own most famous
test cases. Le Verrier applied it to Uranus and found Neptune. Le Verrier
applied it to Mercury and found Vulcan, which does not exist; the resolution was
to change the law. The field remembers Neptune.

The empirical case for dark matter is genuinely strong and multi-scale — CMB
peaks, BAO, lensing, the Bullet Cluster — and I am not arguing against it. The
point is narrower: the **residue** is where the fork shows. The radial
acceleration relation is tighter than a separate collisionless substance with
its own dynamics has any right to produce. Under the "add content" branch that
tightness must be a conspiracy of halo and disc. A conspiracy is what a
misidentified definition looks like from inside.

### Dark energy — the same move, one level up

"Dark energy density" is not measured. What is measured is a
distance–redshift relation. Calling the fitted residual an energy density
imports the Fork B assumption that the vacuum has an absolute energy — the very
quantity QFT cannot define. If DESI's evolving w survives, the constant we
inferred was not constant, which is the mildest possible version of "the
identification was wrong."

---

## Fork D — A heuristic was promoted to a criterion

### Naturalness

Chain: Dirac's large numbers hypothesis (1937) → Wilson's RG picture, in which
low-energy parameters are outputs of high-energy physics → 't Hooft's technical
naturalness (1979) → thirty years of BSM model building → the LHC.

Each step was reasonable. The compound is that a heuristic about where to look
became a criterion for what is wrong. And the criterion requires a **probability
measure on parameter space that nobody has**. "This value is unlikely" is not a
statement anyone can currently make.

This is why the hierarchy problem is in `MEASUREMENTS.md` only as a note that it
does not belong there: no measurement disagrees with anything. Strong CP is the
same — θ̄ < 10⁻¹⁰ is a fact, and its being *surprising* is a claim about a prior.

The LHC is often described as having falsified supersymmetry. What it did was
fail to confirm a criterion that was never a measurement. Those are different
enough to matter, and the difference is at the fork, not at the collider.

---

## Contradictions with a structural cause but no conceptual fork

Honesty requires this section, or the pattern above is being fitted to
everything.

**G — the constant with no redundancy.** The scatter in G is not a fork. It has
a structural cause worth naming anyway: **G is the only fundamental constant
that enters no independent phenomenon.** α shows up in atomic spectra, g−2,
the Lamb shift, quantum Hall, and a dozen more, so an error in one route
contradicts the others. G is defined by the very equation being tested, and
every measurement is a force between laboratory masses over centimetres to
metres. There is no cross-check available, ever. A systematic in G is invisible
to consistency by construction, which is a sufficient explanation for why the
scatter is ten times the quoted errors and why it has not converged in two
hundred years.

**α — Cs versus Rb.** Probably plain systematics in one interferometer. Same
method, same quantity, different atom. No fork. I would not bet on new physics
here.

**Lithium-7.** Either stellar depletion or BBN. Hard astrophysics, no
conceptual fork visible. Included so the analysis has a null result in it.

**CMB large-angle anomalies.** One sky, look-elsewhere effects, and no
possibility of an independent measurement. Structurally unresolvable rather than
conceptually forked.

---

## The window

Sort the forks by date and something uncomfortable appears.

| year | what changed status |
|---|---|
| 1908 | time: evolution parameter → coordinate (Minkowski) |
| 1915 | zero of energy: convention → observable (Einstein, T_μν) |
| 1915 | coordinates: labels → pure gauge (general covariance) |
| 1918 | energy conservation: law → identity (Hilbert's conjecture, Noether's second theorem) |
| 1922–29 | H₀: model coefficient *and* empirical slope, named once |
| 1925–27 | time: retained as external parameter (Schrödinger, Heisenberg) |
| 1932 | observer: retained as external (von Neumann) |
| 1933 | time: *cannot* be internalised (Pauli) |

Nearly every foundational contradiction in physics is a fossil of a
twenty-five-year period in which the ontological status of a handful of
quantities was changed in one theory and deliberately preserved in the other.
Nothing since has moved a single one of these. The subsequent century has
produced measurements that expose them at ever-higher precision, and no
reclassification.

And the most important detail: **Hilbert and Noether identified the energy
problem in 1917–18, immediately, at the moment of formulation.** Hilbert
believed no proper energy conservation law exists in GR and wanted a theorem;
Noether's second theorem gave it to him. This was not discovered later by
someone probing for weaknesses. It was known at the start, correctly diagnosed,
and then absorbed as a technicality.

---

## What this means for this project

Three consequences, and one of them is uncomfortable.

**1. The model is on the far side of Fork A by construction, and that is its
strongest property.** A relational structure with no external stage does not
have to *solve* the problem of time; it never creates it. `necessity.py` has no
clock and no observer, and it still produced conservation, multiplicity, drift,
8.85 nats of persistent information and an effective interaction.

**2. The price is exactly what we measured, and now it has a reason.** No
metric, no dimension, no arrow. That is not bad luck. The metric, the dimension
and the arrow were **living on the stage that was abolished.** Fork A says: give
up the external stage and you give up the things defined relative to it. Our ten
failed attempts to get geometry are the local instance of a hundred-year-old
structural fact, and CDT's crumpled-versus-branched-polymer result is the same
fact measured by someone else.

That is a better outcome than a puzzle. The negative results in `LEDGER.md` stop
being a list of things that did not work and become a prediction: *any*
stage-free structure will produce mechanics and fail to produce geometry, and
the ones that appear to succeed will be found to have re-imported a stage. That
is falsifiable, and `TIMELESS_MODELS.md` is fifteen tests of it that all came
out the same way.

**3. The uncomfortable one: Fork C is the trap this project is standing in.**
Every entry in Fork C is an operational definition identified with a structural
one while they agreed. The self-audit already found three of ours — "equal N" as
an unnamed simultaneity, "sweeps" standing in for wall clock, and m = θ where θ
was typed in and then described as derived. Those are the same error as calling
Hubble's slope and Friedmann's coefficient by one symbol. The difference is only
that nobody has built a century on ours yet.

The rule that follows: **never give one name to a quantity the model computes
and a quantity the model was given.** That is the whole of Fork C, and it is the
only one of the four that is fully in our control.

---

## Sources

- [Noether's Theorems and Energy in General Relativity, arXiv:2103.17160](https://arxiv.org/pdf/2103.17160)
- [Emmy Noether on Energy Conservation in General Relativity, arXiv:1912.03269](https://arxiv.org/pdf/1912.03269)
- [Emmy Noether and Her Theorems, Rowe, Annalen der Physik 2024](https://onlinelibrary.wiley.com/doi/full/10.1002/andp.202300479)
- [Noether, 'Invariant Variation Problems', 1918 — Brading](https://www.lms.ac.uk/sites/default/files/files/Events/2018_09%20Brading%20Noether.pdf)
