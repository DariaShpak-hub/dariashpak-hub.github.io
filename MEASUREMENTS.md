# What is actually measured, how well, and where it disagrees

A reference sheet. Nothing in this file comes from our model — it is the
external landscape, assembled so that claims made elsewhere in this repository
can be checked against it.

Three things get confused constantly and are kept apart here:

- **Precision** — how many digits the apparatus returns.
- **Directness** — how many theoretical steps sit between the dial and the
  number that gets quoted. A quantity can be known to 0.1% and still be a fit
  parameter of a model rather than a reading.
- **Reproducibility** — whether two unrelated techniques land in the same place.

A number that is precise, direct and reproduced is knowledge. A number that is
precise and reproduced but *indirect* is knowledge **conditional on the model
that connects them**, and the condition is usually dropped when the number is
quoted.

Numbers below are as of early 2026. Where a value is moving I say so.

---

## Part 1 — The measured things, ranked by directness

### Tier 1: read off an apparatus, no cosmology in between

| quantity | precision | reproduced? |
|---|---|---|
| Second (Cs hyperfine) | exact by definition | — |
| Optical clock frequency ratios | ~10⁻¹⁸ fractional | yes, Al⁺/Yb/Sr agree |
| Electron g−2 | 1.3×10⁻¹³ on g/2 | single group, but stable over 15 yr |
| Rydberg constant | ~1×10⁻¹² | yes |
| h, e, k_B, N_A | exact by SI definition (2019) | — |
| Proton/electron mass ratio | ~1×10⁻¹¹ | yes |
| Muon lifetime → G_F | ~0.5 ppm | yes |
| Neutron lifetime | 0.03% (bottle) | **no — see below** |
| Fine structure constant α | 8×10⁻¹¹ | **no — see below** |
| Newton's G | 2.2×10⁻⁵ (CODATA) | **no — see below** |

Note what that column shows: the four best-measured constants in the table
include three of the project's open contradictions. Precision and agreement are
independent axes.

A specific irony worth holding onto: the electron g−2 is the most precisely
measured quantity in physics and the most precise test of QED — but converting
it into a test requires α, and α is one of the disputed ones. The world's
sharpest confrontation between theory and experiment is currently limited by a
disagreement between two atom-recoil benches.

### Tier 2: measured through a short, well-tested chain

| quantity | precision | notes |
|---|---|---|
| Higgs mass | ~0.1% | direct reconstruction |
| W boson mass | ~10 MeV / 0.012% | needed heavy modelling of QCD; see the CDF story |
| Top quark mass | ~0.3% | definition-dependent at the ~1 GeV level |
| Δm²₂₁, Δm²₃₁ (neutrino) | few % | robust, several experiments |
| Gravitational wave strain | ~10⁻²¹ | GW170817 pinned c_gw to \|c_gw/c − 1\| < ~10⁻¹⁵ |
| PPN γ (Cassini) | γ−1 = (2.1 ± 2.3)×10⁻⁵ | GR passes |
| D/H primordial abundance | ~1% | agrees with BBN |
| CMB angular scale θ* | ~0.03% | the single best-measured cosmological quantity |

### Tier 3: quoted as measurements, actually model-conditional fits

These are precise. They are also outputs of a six-parameter cosmological model
fitted to data, not readings.

| quantity | precision | what it is conditional on |
|---|---|---|
| Ω_c h² (dark matter density) | ~1% | ΛCDM + GR + standard recombination |
| H₀ from CMB | ~0.5% | ΛCDM, via the sound-horizon ruler |
| Age of the universe | ~0.3% | ΛCDM |
| Σm_ν upper limit | — | ΛCDM; the bound moves when w(z) is freed |
| Dark energy density | ~2% | never observed, only inferred from distances |
| Black hole masses (EHT) | ~10% | fitted to a family of simulated images |

The distinction matters most for H₀, because the two sides of the famous tension
are not the same kind of object. One is a distance ladder you can resolve
rung by rung; the other is a parameter in a fit.

---

## Part 2 — Contradictions

### A. Measurement versus measurement

Same quantity, two apparatuses, no theory in dispute. These are the cleanest,
and the hardest to explain away.

**α — fine structure constant.** Berkeley caesium recoil (2018) and LKB rubidium
recoil (2020) both reach ~10⁻¹⁰ and differ by **more than 5σ**. Same method,
different atom. Nobody thinks the atom is the cause; the working assumption is
an unaccounted systematic in one of the two interferometers, likely laser
beam-profile related. Unresolved.

**Neutron lifetime.** Beam average **888.1 ± 2.0 s**; magnetic-bottle average
**877.8 ± 0.3 s**. About **5σ**. Beam counts protons produced; bottle counts
neutrons surviving. Two ways of asking the same question about the same
particle. Forty representatives of every operating experiment met at PSI in
September 2025 specifically about this; BL3 and upgraded UCN traps are being
built to settle it. Unresolved.

**Newton's G.** CODATA quotes 2.2×10⁻⁵ relative uncertainty, but that number is
manufactured: the individual experiments scatter by roughly ten times their own
stated errors, spanning several hundred ppm. This is not a tension between two
camps, it is a field-wide failure of error estimation on the oldest constant in
physics. Unresolved, and quietly so.

**H₀ within the local ladder.** SH0ES (Cepheids) gives ~73; CCHP (TRGB, JAGB)
gives 68.8 ± 1.8 ± 1.3, or 70.4 ± 1.2 ± 1.3 from TRGB alone. Both are
late-universe, both use JWST. So part of the Hubble tension is a *local*
measurement-versus-measurement dispute about distance-ladder calibration, not
only early-versus-late.

**Hadronic vacuum polarisation.** e⁺e⁻ → hadrons cross-sections from CMD-3
disagree with BaBar and KLOE. This is the live disagreement now sitting under
the muon g−2 story — see below.

### B. Measurement versus theory

**The cosmological constant.** Observed vacuum energy density versus the
quantum field theory estimate: between 60 and 120 orders of magnitude,
depending on the cutoff. This is the largest failed prediction in the history
of physics and it has never moved.

**Strong CP.** θ̄ < ~10⁻¹⁰ from the neutron electric dipole moment. The Standard
Model permits any value up to π. Not a contradiction with data — a contradiction
with expectation.

**Lithium-7.** BBN predicts Li/H ≈ 4.7×10⁻¹⁰; the Spite plateau in metal-poor
halo stars shows ~1.6×10⁻¹⁰. A factor of ~3. Every other BBN abundance works to
percent level, which is what makes this one uncomfortable — you cannot blame
the framework without losing deuterium.

**Neutrino mass.** The Standard Model predicts exactly zero. Oscillations show
nonzero. This is a genuine, confirmed, already-patched SM failure, and it is
worth remembering that the patch (right-handed neutrinos or a dimension-5
operator) is a choice, not a measurement.

**Baryon asymmetry.** SM CP violation is too small by roughly 10⁸ to produce the
observed matter excess.

**Muon g−2 — and what happened to it.** Fermilab's final result (June 2025)
reached 127 ppb. The Theory Initiative's 2025 white paper then dropped the
data-driven hadronic vacuum polarisation in favour of lattice QCD, and the
discrepancy **vanished**. But the 5σ did not disappear — it *moved*, from
experiment-versus-theory to lattice-versus-data-driven. This is a template worth
noting: an anomaly can be resolved by changing which measurement the prediction
is built on.

**Naturalness / the hierarchy problem.** Included here because it is usually
listed with the above, and it should not be. There is no measurement that
disagrees with anything. It is an argument that a parameter value is
surprising. Aesthetic dissatisfaction is not data.

### C. Measurement versus observation

Inference through a model chain, set against something you can point an
instrument at. This is where the two categories are most often conflated.

**H₀.** Early-universe ~67–68 versus late-universe ~73, exceeding 5σ, and not
resolved as of early 2026 — JWST confirmed the Cepheid photometry and removed
crowding and dust as explanations, which *deepened* rather than eased it. The
asymmetry to keep in view: the early number is not an observation of the
expansion rate. It is the value ΛCDM requires given the sound horizon. If the
model between the CMB and today is wrong, the "measurement" is a
misattribution, not an error.

**Dark matter.** Inferred gravitationally at every scale — rotation curves,
lensing, CMB acoustic peaks, BAO, cluster collisions — and the inferences agree
with each other beautifully. Direct detection has found nothing; LZ and XENONnT
have pushed WIMP cross-sections down to the neutrino fog. So we have an
extremely well-measured *gravitational effect* and zero non-gravitational
observation of the thing. Both halves of that sentence are true and usually only
one is said at a time.

**Dark energy.** Same structure, and now with an added wrinkle: DESI DR2 (2025)
prefers an evolving equation of state over a constant at roughly 3σ (values
range 2.8–4.2σ depending on which supernova compilation is included). The
significance is contested — several groups argue the CMB/BAO/SN combination is
internally tense and the preference is not robust. If it holds, the constant we
inferred is not constant.

**Neutrino mass, cosmology versus the lab.** DESI+CMB combinations give Σm_ν
bounds tight enough to press against the oscillation floor (0.059 eV for normal
ordering), with the tightest 2σ limits below 0.05 eV. Worse: when the fit is
allowed to run into unphysical territory, the data *prefer negative* effective
neutrino mass, at up to 3σ from oscillation results. A cosmological inference is
disagreeing with a laboratory measurement, and the cosmological side is the
model-dependent one.

**Small-scale structure.** Core-cusp, missing satellites, too-big-to-fail,
planes of satellites. Simulation output versus observed dwarf galaxies.
Substantially eased by baryonic feedback — but feedback is a fitted sub-grid
prescription, so the resolution moved the problem into a tunable layer rather
than removing it.

**Early massive objects.** JWST finds more bright massive galaxies at z > 10
than standard star formation in ΛCDM expected, plus supermassive black holes
too early and the "little red dots". Partly eased by AGN contamination, bursty
star formation and revised dust — the pattern again is a squeeze rather than a
falsification.

**The radial acceleration relation.** Rotation curves obey a very tight
one-parameter regularity linking observed acceleration to baryonic acceleration.
MOND predicts it directly; ΛCDM must produce it as an emergent conspiracy of
halo and disc. The scatter is smaller than the halo-to-halo scatter has any
right to allow. This one gets less attention than it deserves precisely because
the alternative it favours fails elsewhere.

**The cosmic dipole.** Number counts of radio sources and quasars show a dipole
significantly larger than the CMB dipole implies, at 4–5σ. If real it strains
the cosmological principle — the assumption underneath everything in Tier 3.

**CMB large-angle anomalies.** Low quadrupole, hemispherical asymmetry, cold
spot, alignment of low multipoles. All individually weak, all subject to
look-elsewhere, and there is only one sky so they cannot be re-measured.

### D. The base rate — what dissolved

This is the most useful column for calibration, and it is almost never printed
next to the live anomalies.

| anomaly | lifetime | outcome |
|---|---|---|
| OPERA superluminal neutrinos | 2011–2012 | loose fibre connector |
| BICEP2 primordial B-modes | 2014–2015 | galactic dust |
| 750 GeV diphoton | 2015–2016 | fluctuation |
| Proton radius puzzle | 2010–2019 | electronic measurements were wrong; muonic value correct |
| W mass (CDF II, 7σ) | 2022–2024 | CMS and ATLAS agree with the SM; CDF II stands alone |
| Muon g−2 | 2001–2025 | theory input replaced; discrepancy moved to data-vs-lattice |
| S₈ (KiDS) | 2021–2025 | redshift calibration; KiDS-Legacy now agrees with Planck at 0.73σ |

Most 3–5σ anomalies die. The two big ones that survived — neutrino oscillation
and accelerated expansion — survived by being confirmed with an *independent
method* within a few years, not by getting more significant with more data from
the same technique.

Applying that rule to the current list: the neutron lifetime and α are being
attacked by genuinely independent methods and will resolve. H₀ has already
survived a serious independent test (JWST photometry) and is now a dispute about
what kind of object each side is. DESI's w(z) has not yet had its independent
confirmation.

---

## Part 3 — What this implies about "settled physics"

The confident band of tested physics runs from roughly **10⁻¹⁹ m** (LHC) to
**10²⁶ m** (observable horizon), and from about **one second** after the hot
dense phase onward. That is a wide band and it should not be minimised.

But:

- The Planck scale is **sixteen orders of magnitude** past the small end. Every
  quantum gravity programme, including combinatorial ones, is extrapolating
  across a gap larger than the range that has been tested.
- Four of the surviving discrepancies — G, α, the neutron lifetime, H₀ — sit in
  the **middle** of the well-understood regime, not at any frontier. They are
  table-top and near-field, and they are disagreements between apparatus and
  apparatus.
- "Measurement versus measurement" is often filed as a systematics problem, and
  filing it that way stops inquiry. It may well be systematics. Nobody has shown
  which systematic, in any of the four, after decades.

The honest summary is not "physics is in crisis" and not "the foundations are
secure". It is: the framework is extraordinarily well confirmed over a specific
domain, the confirmations mostly test the framework against itself, and the
places where two independent instruments are pointed at the same quantity are
disproportionately where the disagreements are.

---

## Sources

- [Muon g−2 Theory Initiative White Paper (2025)](https://muon-gm2-theory.illinois.edu/white-paper-25/)
- [Muon g−2 Run 2/3 measurement, arXiv:2506.21219](https://arxiv.org/pdf/2506.21219)
- [Muon Experiment Calls It a Wrap, APS Physics 18, 116](https://physics.aps.org/articles/v18/116)
- [CCHP status report, arXiv:2408.06153](https://arxiv.org/abs/2408.06153)
- [The Hubble Tension: A Decade Review](https://iopscience.iop.org/article/10.1088/1674-4527/ae842f)
- [The Hubble tension, CERN Courier](https://cerncourier.com/a/the-hubble-tension/)
- [DESI DR2 results guide](https://www.desi.lbl.gov/2025/03/19/desi-dr2-results-march-19-guide/)
- [Did DESI DR2 truly reveal dynamical dark energy? arXiv:2504.15222](https://arxiv.org/abs/2504.15222)
- [The beam–bottle debate at PSI, CERN Courier](https://cerncourier.com/a/the-beam-bottle-debate-at-psi/)
- [CODATA 2022, arXiv:2409.03787](https://arxiv.org/pdf/2409.03787)
- [α to 81 parts per trillion (LKB)](https://hal.science/hal-03107990/file/main.pdf)
- [KiDS-Legacy cosmic shear, A&A 2025](https://www.aanda.org/articles/aa/full_html/2025/11/aa54908-25/aa54908-25.html)
- [W mass snaps back, CERN Courier](https://cerncourier.com/a/w-mass-snaps-back/)
- [PDG 2025 review: mass and width of the W boson](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-w-mass.pdf)
- [Neutrino cosmology after DESI, arXiv:2407.18047](https://arxiv.org/abs/2407.18047)
