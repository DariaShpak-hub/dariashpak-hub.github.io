# Given the assumption that the world is combinatorial

Everything in this repository is conditional on one assumption nobody has
justified: that the world is a discrete relational structure. This file is the
ledger of what follows from that assumption and what does not, with the
measurement behind each line.

The assumption itself is a **methodological preference**, not a result. No
theorem says geometry must emerge from anything. That it worked for chemistry
and thermodynamics is a track record, not a proof. Every entry below inherits
that caveat.

---

## What follows, with nothing else added

These come from three primitives only — **relation**, **alternative**,
**consistency** — with no time, space, dimension, energy, mass, probability,
entropy, particles or fields. All exact, all enumerated over the complete state
space rather than sampled (`necessity.py`, 6 elements, 32768 states).

| | measurement |
|---|---|
| **Conservation laws** | 6944 sectors = 6944 distinct degree vectors, *equal integers*. The invariant is complete: two states are mutually reachable iff they agree on it. Nothing posited it. |
| **Multiplicity** | Bare counting. Sector sizes run 1 to 70. |
| **Statistical drift** | Walk entropy 10.1639 against a maximum of 10.3972 — present under consistency, exactly absent without it. |
| **Persistent information** | 8.8456 nats surviving forever, with no record, no memory object, no stored bit. |
| **Effective interaction** | −ln Ω over the cross-cut count has a minimum at m = 4. A preferred separation with a restoring tendency, from counting alone. No force introduced. |
| **Parity sectors** | Z₂ structure appears unbidden — though *not* independent: it is a function of the degree vector, which was already complete. |

And the load-bearing negative that makes all of the above meaningful:

**Without the consistency requirement, none of it exists.** One component,
no invariant, walk entropy exactly maximal, **I = 0.0000 nats**. Local
substitution on its own produces no physics whatsoever.

> Conservation is not a primitive. It is the algebraic shadow of restriction.

---

## What follows once colour is imposed

These require typed edges with a commutation relation — a structure that is
*added*, not derived.

| | measurement |
|---|---|
| **Exact coordinates** | Every node has a Zᵏ address, **0 inconsistencies** at every k, size and seed (`coloured_metric.py`) |
| **Isotropy** | Spread 0.0007 at k = 2, N = 32000; decaying as N^−0.47 — the central-limit rate (`coloured_metric.py`) |
| **A dimension** | d = 1.84, 2.84, 3.81 for k = 2, 3, 4 (`coloured_memory.py`) |
| **Mass as a rate** | m = θ to five decimals; peak speed tracks cos θ to ~1%; the light cone stays at 1 (`zigzag_mass.py`) |

The dimension is the colour count, and the colour count is chosen. So this
column buys structure at the price of the thing it was supposed to explain.

---

## What does not follow

| | measurement |
|---|---|
| **A metric** | No distance anywhere in relation + alternative + consistency. Reached independently from the geometry side and the order side. |
| **A preferred dimension** | d is a dial. `output_classes.py` found the output is a single continuum with **no feature at 3** — not a failure to find the mechanism, an absence of anything for a mechanism to select. |
| **The arrow of time** | Absent in four independent places. See below. |
| **Expansion with isotropy** | **0 of 600** parameter settings, twice over, at the pre-registered threshold (`fine_tuning.py`) |
| **Causal structure** | The growth poset is a tree: ordering fraction 0.0354 → 0.0061 as n runs 600 → 4801, falling as 1/n while Minkowski holds at 0.50 (`causal_order.py`) |
| **Particles** | A Z₂ defect exists — gauge-invariant, local, loop intact — but geodesics route around it, so it is metrically invisible (`holonomy_defect.py`) |
| **Mass values** | m = θ, and θ is typed in. Nine Yukawa couplings became nine angles. The hierarchy — a factor of 340,000 electron to top — is untouched. |
| **Gravity** | Nothing. |

---

## The arrow is absent in four independent places

This is the sharpest negative in the project, because the four are unrelated.

1. **The metric.** ds² = −dt² + dx² is exactly invariant under t → −t.
2. **The order.** Reverse every relation in a causal set and you get another
   valid causal set. Orientation is a labelling convention.
3. **The dynamics.** Detailed balance holds; net current across every edge is
   exactly zero. Non-uniform is not directed (`causal_order.py`).
4. **The entropy identity.** E[ΔS_B] = 0 at every resolution. A random rewrite
   has zero expected entropy change.

The only thing that ever produced a persistent arrow was a **receding ceiling** —
Ω growing faster than the system fills it, giving a constant 16.6% lag on an
open substrate while every smooth start on a closed one gave a *negative* arrow
(`gravity_arrow.py`).

That is a boundary condition, not a law.

---

## Three structural results that constrain other people's models

**1. Every layer identification broke something.** The pattern held without
exception:

```
cycle          ≠  face
colour         ≠  memory
missing loop   ≠  curvature
defect         ≠  metric source
spatial parent ≠  causal predecessor
order          ≠  geometry
signature      ≠  orientation
refinement     ≠  expansion
```

Progress came only from *separating* things that had been fused.

**2. Localized refinement cannot produce expansion.** Shortest paths route
around any bounded region made expensive, so dilation saturates no matter how
much refinement is added — measured over a 64-fold budget increase, mean ε
0.012 → 0.017 while a uniform control ran 0.040 → 6.90 (`holonomy_defect.py`).

**3. Disorder generates mass indistinguishably from a coupling.** With the coin
switched off entirely, 8.6% blocked steps produced exponent 1.702 and speed
0.636 — a mass with no coupling anywhere (`zigzag_mass.py`). Any discrete model
claiming to *derive* a mass must first prove it is not measuring its own
lattice defects. And the rates differ: metric convergence goes as N^−0.5 while
connectivity convergence goes as N^−0.18, so a substrate can look geometrically
converged and still be a dense scattering field.

---

## The asymmetry nobody would have predicted

**Order gives mechanics almost for free. It gives geometry not at all.**

Conservation, multiplicity, drift, persistent information and an effective
interaction all fall out of restriction with nothing added. Distance, dimension
and the arrow do not appear under any amount of pressure.

If the combinatorial assumption is right, that asymmetry is the most
informative thing here: it says mechanics is cheap and geometry is expensive,
which is the opposite of how the two are usually ranked.

---

## Instruments left behind

Independent of any interpretation, these work and are calibrated:

- **Invariant discovery** — integer kernel via Smith normal form, including
  torsion. Rediscovered Q = B + 2N − E from three change vectors with no
  knowledge of what a bank is, and correctly returns *nothing* when the rules
  admit nothing (`rewrite_algebra.py`).
- **Rewrite classifier** — best pure dilation of a patch boundary. Reads
  `close` as a shortcut (λ = 0.889) and uniform subdivision as clean expansion
  (λ = 2.000, δ = 0.000), instantly. Would have rejected `close` before any
  cosmology was built on it.
- **Dimension gate** — two independent estimators on order alone, returning
  1.94, 2.96, 3.81 for true dimensions 2, 3, 4 (`causal_growth.py`).
- **Isotropy against Zᵏ** — axis ratios with no baseline epoch, so no marker
  artifact (`coloured_metric.py`).

---

## What would falsify the assumption

Nothing here does. That is worth saying plainly: this is a reconstruction
programme, and reconstruction cannot fail in the way a prediction can.

The nearest thing to a falsification available would be a proof that no
combinatorial structure of the relevant kind admits a metric with a directed
arrow — and this project has assembled evidence pointing that way without
coming close to proving it.

Until then, every line above should be read with its first four words in place.
