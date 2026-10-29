# The Energy Cost of a Phase Shift

Deep Bhattacharjee (corresponding author) and Sanjeevan Singha Roy.

The paper proves that a passive causal system without zeros in the upper
half-plane cannot shift every frequency of a band [a, b] by the same phase
θ > 0 without loss. Its mean attenuation over the band is at least

    θ · (K(κ) − E(κ)) / κ   nepers,   κ = (b − a)/(b + a),

with K, E the complete elliptic integrals, and an explicit system attains the
bound. Delays are free and lossless resonances raise the bound. For a phase
θ + τω on the band the same explicit system, or its inversion ω ↦ ab/ω for a
lead, has the least attenuation at every frequency of the band at once: every
unit of spectral measure off the band raises the attenuation at every point
of the band. This gives the sharp worst-frequency loss (1.273 θΛ for an
octave) and the sharp plain mean for a lead, θ(E − κ'²K)/κ. With zeros
in the upper half-plane allowed, one explicit system with a single zero has
the least attenuation at every frequency of the band at once; for a right
angle it is algebraic, its mean attenuation is (π/2)(√b − √a)/(√b + √a)
nepers, and no passive system with this phase passes more than the fraction a/b of the power at
the geometric centre √(ab) (at least 3.01 dB of loss over an octave). The paper also checks
the formulas of the earlier paper by S. N. Thakur and D. Bhattacharjee,
"Phase shift and infinitesimal wave energy loss equations" (J. Phys. Chem.
Biophys. 13 (2023) 365), and corrects them.

Release 1.1.0 adds Section 9 (the pointwise bounds without zeros; the positive
reconstruction behind it was found with OpenAI's ChatGPT) and corrects three
points of 1.0.0: a missing factor 2 in the proof of the Gauss lemma, the
values of zε and zμ on the imaginary axis, and the reading of the energy
bound, which concerns the transmitted channel (the power may be absorbed,
reflected or routed to other ports).

## Layout

- `EnergyCostOfAPhaseShift.pdf`: the compiled paper (REVTeX, two columns, 34 pages).
- `paper/`: the LaTeX source (`main.tex`, `sections/`, `figures/`).
- `paper/anc/checks.py`: high-precision re-evaluation of every identity and
  numerical value of the paper (needs `mpmath`; 244 checks, about six minutes).
- `paper/anc/make_figures.py`, `paper/anc/make_schematics.py`,
  `paper/anc/make_figures_zeros.py`,
  `paper/anc/make_figures_pointwise.py`: the figures (needs `numpy`, `scipy`,
  `matplotlib`, `mpmath`).
- `lean/PhaseShiftCost.lean`: a Lean 4 certificate of the finite and
  algebraic parts (75 theorems, no `sorry`, core library only); the axiom
  report is `lean/axioms.txt`.

## Running the checks

    cd paper && python3 anc/checks.py        # 244 checks passed, 0 failed
    cd lean && lean PhaseShiftCost.lean      # Lean 4.34.0, a few seconds
    cd paper && pdflatex -jobname=output main.tex   # three passes

## What is and is not machine-checked

Every theorem has a proof in the paper; no statement rests on the programs.
The Lean file checks the dimensional analysis and arithmetic of the audit,
the coefficients of the series of the constant and their monotonicity for
every n, enclosures of the constant to nine decimals for six band ratios,
and the polynomial identities used in the proofs, among them those behind the
closed form of the optimal system at a right angle and the signs and
comparisons of the pointwise bounds without zeros. Of its 75 theorems, 30
depend on no axiom and the rest on the standard axioms of Lean
(`propext`, `Classical.choice`, `Quot.sound`). The analytic part of the paper
(the sum rule, the reconstruction of Section 9, boundary values of Cauchy integrals, the M. Riesz theorem,
Hardy spaces, the conformal map and extremal length used for systems with
zeros, and the elliptic-integral identities) is not formalised.

## Archive

| Release | Tag | Zenodo DOI |
| --- | --- | --- |
| 1.0.0 | `v1.0.0` | [10.5281/zenodo.23761044](https://doi.org/10.5281/zenodo.23761044) |
| 1.1.0 | `v1.1.0` | [10.5281/zenodo.23764158](https://doi.org/10.5281/zenodo.23764158) |

Concept DOI of all versions: [10.5281/zenodo.23761043](https://doi.org/10.5281/zenodo.23761043).

## Contact

GitHub-related contact: egspl@proton.me (GitHub: creelie). Corresponding
author: d.bhattacharjee@erl-forschung.de; also itsdeep@live.com.
