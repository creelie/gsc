# The Exact Price of a Phase Shift

Deep Bhattacharjee (corresponding author), Soumendra Nath Thakur and Ushashi
Bhattacharya.

A passive causal system without zeros in the upper half-plane cannot shift
every frequency of a band [a, b] by the same phase θ > 0 without loss. The
paper finds the least loss exactly. Its mean attenuation over the band is at
least

    θ · (K(κ) − E(κ)) / κ   nepers,   κ = (b − a)/(b + a),

with K, E the complete elliptic integrals, and an explicit system attains the
bound at every frequency of the band at once. Delays are free and lossless
resonances raise the bound. On any finite union B of bands the least mean
loss for the phase θ + τω is θ·Λ_B for an explicit constant Λ_B, attained by
a system with one lossless resonance in each gap, placed where the test
function of the bands vanishes; for two bands [a₁, b₁] ∪ [a₂, b₂] the
resonance is at c² = (b₁²b₂² − a₁²a₂²)/(b₁² + b₂² − a₁² − a₂²). One device
for several bands always costs more than separate devices. With zeros
allowed, one explicit system with a single zero has the least attenuation at
every frequency of the band; for a right angle its mean attenuation is
(π/2)(√b − √a)/(√b + √a) nepers and at least 3.01 dB is lost at the centre of
an octave. With a delay and zeros the optimum exists, is unique, is the limit
of finite convex problems, and its loss decays exponentially in the delay.

The paper also checks and corrects the formulas of its earlier version,
S. N. Thakur and D. Bhattacharjee, "Phase shift and infinitesimal wave energy
loss equations", Preprints.org, doi:10.20944/preprints202309.1831.v1 (2023),
published in J. Phys. Chem. Biophys. 13 (2023) 365. A later version is the
Hexagon preprint 2610.00657v1.

Still open: the optimum with zeros and a delay in closed form and its exact
decay exponent, its behaviour at the ends of the band, nonaffine phases
without zeros where the velocity changes sign, a strong lead and nonaffine
phases on several bands, the least trace for an indefinite weight that does
not commute with the slope, joint bounds with reflection, and a formal proof
of the analysis.

## Layout

- `TheExactPriceOfAPhaseShift.pdf`: the compiled paper (REVTeX, two columns, 46 pages).
- `paper/`: the LaTeX source (`main.tex`, `sections/`, `figures/`).
- `paper/anc/checks.py`: high-precision re-evaluation of every identity and
  numerical value of the paper (needs `mpmath`, and `scipy` for
  `paper/anc/zero_lp.py`; 456 checks).
- `paper/anc/make_figures.py`, `paper/anc/make_schematics.py`,
  `paper/anc/make_figures_zeros.py`,
  `paper/anc/make_figures_pointwise.py`, `paper/anc/make_figures_beyond.py`,
  `paper/anc/make_figures_bands.py`: the figures (needs `numpy`, `scipy`,
  `matplotlib`, `mpmath`).
- `lean/PhaseShiftCost.lean`: a Lean 4 certificate of the finite and
  algebraic parts (95 theorems, no `sorry`, core library only); the axiom
  report is `lean/axioms.txt`.

## Running the checks

    cd paper && python3 anc/checks.py        # 456 checks passed, 0 failed
    cd lean && lean PhaseShiftCost.lean      # Lean 4.34.0, a few seconds
    cd paper && pdflatex -jobname=output main.tex   # three passes

## What is and is not machine-checked

Every theorem has a proof in the paper; no statement rests on the programs.
The Lean file checks the dimensional analysis and arithmetic of the audit,
the coefficients of the series of the constant and their monotonicity for
every n, enclosures of the constant to nine decimals for six band ratios,
and the polynomial identities used in the proofs, among them those behind the
closed form of the optimal system at a right angle and the signs and
comparisons of the pointwise bounds without zeros, the constants of
Section 12 and the exact rational bound of the six-tap device. Of its 95
theorems, 40
depend on no axiom and the rest on the standard axioms of Lean
(`propext`, `Classical.choice`, `Quot.sound`). The analytic part of the paper
(the sum rule, the reconstruction of Section 9, the convexity, construction and decay
arguments of Section 12, boundary values of Cauchy integrals, the M. Riesz theorem,
Hardy spaces, the conformal map and extremal length used for systems with
zeros, and the elliptic-integral identities) is not formalised.

## Archive

| Release | Tag | Zenodo DOI |
| --- | --- | --- |
| 1.0.0 | `v1.0.0` | [10.5281/zenodo.23761044](https://doi.org/10.5281/zenodo.23761044) |
| 1.1.0 | `v1.1.0` | [10.5281/zenodo.23764158](https://doi.org/10.5281/zenodo.23764158) |
| 1.2.0 | `v1.2.0` | [10.5281/zenodo.23768958](https://doi.org/10.5281/zenodo.23768958) |
| 1.3.0 | `v1.3.0` | [10.5281/zenodo.23772796](https://doi.org/10.5281/zenodo.23772796) |
| 2.0.0 | `v2.0.0` | [10.5281/zenodo.23800767](https://doi.org/10.5281/zenodo.23800767) |

Concept DOI of all versions: [10.5281/zenodo.23761043](https://doi.org/10.5281/zenodo.23761043).

## Contact

GitHub-related contact: egspl@proton.me (GitHub: creelie). Corresponding
author: d.bhattacharjee@erl-forschung.de; also itsdeep@live.com.
