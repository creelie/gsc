# Least Loss of Passive Phase Shifters

Deep Bhattacharjee (corresponding author), Sanjeevan Singha Roy and Ushashi
Bhattacharya.

The paper proves that a passive causal system without zeros in the upper
half-plane cannot shift every frequency of a band [a, b] by the same phase
θ > 0 without loss. Its mean attenuation over the band is at least

    θ · (K(κ) − E(κ)) / κ   nepers,   κ = (b − a)/(b + a),

with K, E the complete elliptic integrals, and an explicit system attains the
bound. Delays are free and lossless resonances raise the bound. For a phase
θ + τω on the band the same explicit system, or its inversion ω ↦ ab/ω for a
lead, has the least attenuation at every frequency of the band at once. With
zeros in the upper half-plane allowed, one explicit system with a single zero
has the least attenuation at every frequency of the band at once; for a right
angle it is algebraic, its mean attenuation is (π/2)(√b − √a)/(√b + √a)
nepers, and no passive system with this phase passes more than the fraction
a/b of the power at the geometric centre √(ab) (at least 3.01 dB of loss over
an octave). The paper also checks the formulas of the earlier paper by
S. N. Thakur and D. Bhattacharjee, "Phase shift and infinitesimal wave energy
loss equations" (J. Phys. Chem. Biophys. 13 (2023) 365), and corrects them.

Release 1.2.0 has a new title, adds Ushashi Bhattacharya as third author and
adds Section 12, on a phase θ + τω with a delay τ > 0, with several channels
and with several bands. Proved there: without zeros the least mean loss in
closed form for every slope, and a positive floor below which it never falls;
with zeros, the optimal system exists and is unique for every phase, with a
test for optimality; its loss falls exponentially with the delay, and the
exponent lies between −πab/(b − a) and −a/e; a trace bound for several
channels and the absence of a least attenuation in the Loewner order; a
compensation criterion for several bands. Not known: the optimal system with
zeros and a delay and its exact exponent, nonaffine phases without zeros where
the velocity changes sign, the least trace for an indefinite weight, the sharp
multiband constant, joint bounds with reflection, and a formal proof of the
analysis.

## Layout

- `LeastLossOfPassivePhaseShifters.pdf`: the compiled paper (REVTeX, two columns, 40 pages).
- `paper/`: the LaTeX source (`main.tex`, `sections/`, `figures/`).
- `paper/anc/checks.py`: high-precision re-evaluation of every identity and
  numerical value of the paper (needs `mpmath`, and `scipy` for
  `paper/anc/zero_lp.py`; 421 checks).
- `paper/anc/make_figures.py`, `paper/anc/make_schematics.py`,
  `paper/anc/make_figures_zeros.py`,
  `paper/anc/make_figures_pointwise.py`, `paper/anc/make_figures_beyond.py`: the figures (needs `numpy`, `scipy`,
  `matplotlib`, `mpmath`).
- `lean/PhaseShiftCost.lean`: a Lean 4 certificate of the finite and
  algebraic parts (87 theorems, no `sorry`, core library only); the axiom
  report is `lean/axioms.txt`.

## Running the checks

    cd paper && python3 anc/checks.py        # 421 checks passed, 0 failed
    cd lean && lean PhaseShiftCost.lean      # Lean 4.34.0, a few seconds
    cd paper && pdflatex -jobname=output main.tex   # three passes

## What is and is not machine-checked

Every theorem has a proof in the paper; no statement rests on the programs.
The Lean file checks the dimensional analysis and arithmetic of the audit,
the coefficients of the series of the constant and their monotonicity for
every n, enclosures of the constant to nine decimals for six band ratios,
and the polynomial identities used in the proofs, among them those behind the
closed form of the optimal system at a right angle and the signs and
comparisons of the pointwise bounds without zeros, and the constants of
Section 12. Of its 87 theorems, 34
depend on no axiom and the rest on the standard axioms of Lean
(`propext`, `Classical.choice`, `Quot.sound`). The analytic part of the paper
(the sum rule, the reconstruction of Section 9, the convexity and decay
arguments of Section 12, boundary values of Cauchy integrals, the M. Riesz theorem,
Hardy spaces, the conformal map and extremal length used for systems with
zeros, and the elliptic-integral identities) is not formalised.

## Archive

| Release | Tag | Zenodo DOI |
| --- | --- | --- |
| 1.0.0 | `v1.0.0` | [10.5281/zenodo.23761044](https://doi.org/10.5281/zenodo.23761044) |
| 1.1.0 | `v1.1.0` | [10.5281/zenodo.23764158](https://doi.org/10.5281/zenodo.23764158) |
| 1.2.0 | `v1.2.0` | [10.5281/zenodo.23768958](https://doi.org/10.5281/zenodo.23768958) |

Concept DOI of all versions: [10.5281/zenodo.23761043](https://doi.org/10.5281/zenodo.23761043).

## Contact

GitHub-related contact: egspl@proton.me (GitHub: creelie). Corresponding
author: d.bhattacharjee@erl-forschung.de; also itsdeep@live.com.
