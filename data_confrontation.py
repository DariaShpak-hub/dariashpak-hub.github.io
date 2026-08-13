#!/usr/bin/env python3
"""
data_confrontation.py

Every part of the model put against the actual observational number.

The discipline this needs is deciding what is COMPARABLE. The model has no
clock, no length scale and no matter content, so most cosmological observables
have no counterpart in it, and pretending otherwise would be the same error as
fitting w(z) to a schedule-dependent quantity. Only two classes of observable
survive translation:

    dimensionless ratios      shear-to-expansion, density contrast, clustering
    scaling exponents         effective dimension, force law

Everything with a unit attached -- H0, the age, a length -- is excluded by
construction, and is listed as such rather than fudged.

Observational values and their sources
----------------------------------------
    spatial dimension        3; inverse-square gravity verified to ~50 microns
                             in torsion-balance experiments (Eot-Wash), and no
                             deviation at LHC energies

    shear                    (sigma/H)_0 < 4.7e-11 at 95% CL, Saadeh et al.
                             2016, "How isotropic is the Universe?", from
                             Planck. The cruder bound from Ellis & van Elst
                             eq. (59) plus CBR isotropy is ~1e-3; the direct
                             Bianchi-model limit is far tighter and is the one
                             used here.

    CMB density contrast     delta rho/rho ~ 1e-5 at recombination

    homogeneity scale        transition to statistical homogeneity at roughly
                             70-100 h^-1 Mpc

    galaxy clustering        xi(r) = (r/r_0)^-gamma with r_0 ~ 5.4 h^-1 Mpc
                             and gamma ~ 1.8, so xi ~ 18 at 1 h^-1 Mpc:
                             pairs are ~19x more likely than random

    baryon asymmetry         eta = n_B/n_gamma = 6.1e-10

    dark energy              DESI DR2 + CMB + SNe prefer w0 > -1, wa < 0 at
                             2.8-4.2 sigma depending on the SN compilation

    entropy budget           S_now ~ 1e104 k_B, dominated by supermassive
                             black holes (Egan & Lineweaver 2010); the de
                             Sitter horizon bound is ~1e122 k_B, so the
                             universe sits ~18 orders of magnitude BELOW its
                             ceiling

Result
------
Six quantities are comparable. Two pass, four fail, and eight more observables
have no counterpart at all.

    quantity              observed          model                 verdict
    effective dimension   3 (to ~50 um)     1.00 / no selection   FAIL
    shear sigma/H         < 4.7e-11         1.50 / sqrt(N)        PASS
    force law exponent    1/r^2             1/r^2 exactly         PASS (not independent)
    clustering vs random  ~19x at 1 Mpc     1.10x                 FAIL
    density contrast      1e-5 at recomb.   8-15%, non-decaying   FAIL
    baryon asymmetry      6.1e-10           exactly 0             FAIL

The shear pass is the real one and it is quantitative. From shear_damping.py
the local rule gives spread * sqrt(N) = 1.496, constant to 11.9% across the
runs, so spread = 1.50 / sqrt(N). Setting that equal to the Saadeh et al. 2016
bound of 4.7e-11 gives

    N = 1.0e21 nodes

which the observable universe exceeds by a wide margin -- it holds ~1e80
baryons and ~1e185 Planck volumes. Against the cruder Ellis-van Elst bound of
1e-3 the threshold is only N = 2.2e6. So a universe of any realistic size
running this growth rule would look isotropic to well beyond current
sensitivity. That is the first quantitative agreement with data anywhere in
this project.

The force-law pass is real but not independent: 1/r^2 follows from Gauss's law
the moment d = 3 is granted, and d = 3 is exactly what the model fails to
produce. It should be read as "consistent given an assumption the model cannot
supply".

The four failures share one root, already measured: nothing selects a
dimension, so there is no 3D substrate; without a 3D substrate there is no
long-range 1/r^2 force; without a long-range force there is contact-interaction
clustering only, which matter_clumping.py measured at 1.10x against the
observed ~19x. The baryon-asymmetry failure is separate and expected -- A_M = 0
exactly by symmetry, which is just the Sakharov conditions restated.

On the entropy gap, which deserves care rather than a verdict
---------------------------------------------------------------
gravity_arrow.py found matter holding a CONSTANT 16.6% behind a receding
ceiling. The universe sits at S/S_max ~ 1e104/1e122, about 1e-18 of its
ceiling, and falling further behind.

These are different notions of ceiling -- the model's is matter on a fixed
substrate, the universe's is the de Sitter horizon bound -- so this is not a
clean refutation. But the discrepancy is structural rather than a units
problem. The model equilibrates to a fixed fractional lag; the real universe
diverges from its bound. A mechanism producing a constant fraction is the wrong
SHAPE for an observable 18 orders of magnitude from saturation, and that would
remain true after any rescaling.

What the exercise establishes
-------------------------------
Most of the model is not falsifiable against data, and saying so is the honest
result. Eight observables are excluded by construction because the model has no
clock, no length scale and no stress-energy. That is not a gap to be filled by
better measurement; it follows from the model being pre-geometric.

Of what remains, the one pass is the one quantity that depends on the growth
law alone and not on the geometry. Everything requiring a geometry fails, and
fails through the same bottleneck the rest of the project already located.
"""

import numpy as np


# ----------------------------------------------------------- model numbers

# shear_damping.py, local rule with weight 2^-g:
#   spread * sqrt(N) is constant, so spread = c / sqrt(N)
SHEAR_N = np.array([2500, 5000, 10000, 20000, 40000], float)
SHEAR_SPREAD = np.array([0.028, 0.019, 0.015, 0.013, 0.007])

OBS_SHEAR = 4.7e-11          # Saadeh et al. 2016, 95% CL
OBS_SHEAR_CRUDE = 1e-3       # Ellis & van Elst eq. (59) + CBR isotropy


def shear_crossing():
    """At what N does the local damping rule meet the observed bound?"""
    c = float(np.mean(SHEAR_SPREAD * np.sqrt(SHEAR_N)))
    resid = SHEAR_SPREAD * np.sqrt(SHEAR_N)
    return c, float(resid.std() / resid.mean()), (c / OBS_SHEAR) ** 2


def line(name, obs, model, verdict, note=""):
    print("  %-24s %-18s %-20s %-11s %s"
          % (name, obs, model, verdict, note))


def main():
    print("=" * 100)
    print("1. COMPARABLE: dimensionless ratios and scaling exponents")
    print("=" * 100)
    print("  %-24s %-18s %-20s %-11s %s"
          % ("quantity", "observed", "model", "verdict", "source module"))
    print("  " + "-" * 96)

    c, spread_cv, n_cross = shear_crossing()

    line("effective dimension", "3 (to ~50 um)", "1.00 / no selection",
         "FAIL", "scale_factor, dimension_attractor")
    line("shear sigma/H", "< 4.7e-11", "1.4 / sqrt(N)",
         "PASS*", "shear_damping (see below)")
    line("force law exponent", "1/r^2", "1/r^2 exactly",
         "PASS", "self_amplification (not independent)")
    line("clustering vs random", "~19x at 1 Mpc", "1.10x (0.905 sep.)",
         "FAIL", "matter_clumping")
    line("density contrast", "1e-5 at recomb.", "8-15%, non-decaying",
         "FAIL", "six_criteria")
    line("baryon asymmetry", "6.1e-10", "exactly 0",
         "FAIL", "sakharov, conjugation")
    print()

    print("  * The shear result is the one genuine pass, and it needs its")
    print("    scaling stated rather than a single number.")
    print("    spread * sqrt(N) = %.3f, constant to %.1f%% across the runs,"
          % (c, 100 * spread_cv))
    print("    so spread = %.2f / sqrt(N) and the observed bound 4.7e-11 is"
          % c)
    print("    reached at N = %.1e nodes." % n_cross)
    print("    For scale: the observable universe holds ~1e80 baryons and")
    print("    ~1e185 Planck volumes, so that threshold is met with room to")
    print("    spare. Against the cruder 1e-3 bound it is reached at N = %.1e."
          % ((c / OBS_SHEAR_CRUDE) ** 2))
    print()

    print("=" * 100)
    print("2. NOT COMPARABLE: no counterpart exists in the model")
    print("=" * 100)
    print("  %-26s %-22s %s" % ("quantity", "observed", "why not comparable"))
    print("  " + "-" * 96)
    for q, o, why in [
        ("H0", "67-73 km/s/Mpc", "has units of 1/time; the model has no clock"),
        ("age of the universe", "13.8 Gyr", "same"),
        ("w0, wa", "w0>-1, wa<0 (DESI)", "w_eff not reparameterisation-invariant"),
        ("spectral index n_s", "0.965 +- 0.004", "requires an inflaton and a"),
        ("", "", "  primordial power spectrum"),
        ("curvature Omega_k", "|Omega_k| < 0.002", "no metric, so no curvature scalar"),
        ("Omega_m, Omega_L", "0.31, 0.69", "no stress-energy budget"),
        ("entropy S", "~1e104 k_B", "different ceiling: the model's is matter"),
        ("", "", "  on a fixed substrate, the universe's is"),
        ("", "", "  the de Sitter horizon bound"),
    ]:
        print("  %-26s %-22s %s" % (q, o, why))
    print()

    print("=" * 100)
    print("3. THE ENTROPY GAP, WHICH IS WORTH ONE LINE OF ITS OWN")
    print("=" * 100)
    print("  gravity_arrow.py found matter holding a CONSTANT 16.6% behind a")
    print("  receding ceiling. The universe sits at S/S_max ~ 1e104/1e122,")
    print("  i.e. about 1e-18 of its ceiling, and the fraction is falling.")
    print()
    print("  Those are different notions of ceiling, so this is not a clean")
    print("  refutation. But the STRUCTURE differs in a way that is not a")
    print("  units problem: the model equilibrates to a fixed fractional lag,")
    print("  while the real universe falls further behind its bound over time.")
    print("  A mechanism that produces a constant fraction is the wrong shape")
    print("  for an observable that is 18 orders of magnitude from saturation.")
    print()

    print("=" * 100)
    print("4. SCORE")
    print("=" * 100)
    print("  comparable quantities   : 6")
    print("  passes                  : 2  (shear, and the force law which is")
    print("                               not independent of assuming d=3)")
    print("  fails                   : 4  (dimension, clustering, density")
    print("                               contrast, baryon asymmetry)")
    print("  not comparable          : 8")
    print()
    print("  Every failure traces to the same root already measured: nothing")
    print("  selects a dimension, so there is no 3D substrate for a long-range")
    print("  force, and without a long-range force there is no clustering and")
    print("  no structure. The single pass is the one quantity that depends on")
    print("  the growth law alone and not on the geometry.")
    print()


if __name__ == "__main__":
    main()
