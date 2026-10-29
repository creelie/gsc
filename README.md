# The Energy Cost of a Phase Shift

Deep Bhattacharjee (corresponding author) and Sanjeevan Singha Roy.

The paper proves that a passive causal system without zeros in the upper
half-plane cannot shift every frequency of a band [a, b] by the same phase
θ > 0 without loss. Its mean attenuation over the band is at least

    θ · (K(κ) − E(κ)) / κ   nepers,   κ = (b − a)/(b + a),

with K, E the complete elliptic integrals, and an explicit system attains the
bound. Delays are free and lossless resonances raise the bound. With zeros
in the upper half-plane allowed, one explicit system with a single zero has
the least attenuation at every frequency of the band at once; for a right
angle it is algebraic, its mean attenuation is (π/2)(√b − √a)/(√b + √a)
nepers, and no passive system with this phase passes more than the fraction a/b of the power at
the geometric centre √(ab) (at least 3.01 dB of loss over an octave). The paper also checks
the formulas of the earlier paper by S. N. Thakur and D. Bhattacharjee,
"Phase shift and infinitesimal wave energy loss equations" (J. Phys. Chem.
Biophys. 13 (2023) 365), and corrects them.

## Layout

- `EnergyCostOfAPhaseShift.pdf`: the compiled paper (REVTeX, two columns, 30 pages).
- `paper/`: the LaTeX source (`main.tex`, `sections/`, `figures/`).
- `paper/anc/checks.py`: high-precision re-evaluation of every identity and
  numerical value of the paper (needs `mpmath`; 161 checks, about three minutes).
- `paper/anc/make_figures.py`, `paper/anc/make_schematics.py`,
  `paper/anc/make_figures_zeros.py`: the figures (needs `numpy`, `scipy`,
  `matplotlib`, `mpmath`).
- `lean/PhaseShiftCost.lean`: a Lean 4 certificate of the finite and
  algebraic parts (70 theorems, no `sorry`, core library only); the axiom
  report is `lean/axioms.txt`.

## Running the checks

    cd paper && python3 anc/checks.py        # 161 checks passed, 0 failed
    cd lean && lean PhaseShiftCost.lean      # Lean 4.34.0, a few seconds
    cd paper && pdflatex -jobname=output main.tex   # three passes

## What is and is not machine-checked

Every theorem has a proof in the paper; no statement rests on the programs.
The Lean file checks the dimensional analysis and arithmetic of the audit,
the coefficients of the series of the constant and their monotonicity for
every n, enclosures of the constant to nine decimals for six band ratios,
and the polynomial identities used in the proofs, among them those behind the
closed form of the optimal system at a right angle. Of its 70 theorems, 30
depend on no axiom and the rest on the standard axioms of Lean
(`propext`, `Classical.choice`, `Quot.sound`). The analytic part of the paper
(the sum rule, boundary values of Cauchy integrals, the M. Riesz theorem,
Hardy spaces, the conformal map and extremal length used for systems with
zeros, and the elliptic-integral identities) is not formalised.

## Contact

GitHub-related contact: egspl@proton.me (GitHub: creelie). Corresponding
author: d.bhattacharjee@erl-forschung.de; also itsdeep@live.com.
