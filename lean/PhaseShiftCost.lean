/-
PhaseShiftCost.lean

A machine check, by the kernel of Lean 4 (core library only, no Mathlib),
of the finite and algebraic parts of "The Energy Cost of a Phase Shift"
by D. Bhattacharjee and S. Singha Roy.

What is checked here, and what is not.

  * Checked: the arithmetic and the dimensional analysis of the audit of
    the earlier paper (Section 2 of the paper); the coefficients of the
    power series of Lambda and their monotonicity, for every n; certified
    enclosures of Lambda/pi for six band ratios, from partial sums of that
    series with a proved tail bound; and the polynomial identities used in
    the proofs (moments of the weight, the Landen step, the two forms of
    tau_star, the upper bound for Lambda, the value of W at 0, the inversion
    of the band, the plain mean of the lead extremal, the threshold of the
    one-zero corollary, the constants pi^2/720 and pi/24, the closed form
    of the optimum at a right angle, and the signs, comparisons and
    linear coefficients of the pointwise bounds without zeros).

  * Not checked: the analytic theorems themselves (the sum rule, the
    identification of Lambda with an elliptic integral and with its power
    series, the extremal system, the one-zero proposition, the
    reconstruction of Section 9).  These rest on
    Cauchy integrals, the M. Riesz theorem and complete elliptic integrals,
    which are not formalised here.  Where a theorem below is the finite or
    algebraic step of an analytic argument, its comment names that step.

Axioms.  The theorems checked by evaluation in the kernel (`decide`) and
by elementary lemmas of the core library depend on no axiom.  Those that
use the tactics `grind` or `omega` (the polynomial identities and the
monotonicity of the series coefficients for every n) depend on Lean's three
standard axioms propext, Classical.choice and Quot.sound.  No theorem uses
`sorry`, `native_decide` or an axiom of its own.  The output of `#print axioms`
at the end of the file lists the dependencies of each theorem.
-/

namespace PhaseShiftCost

/-! ## Section 1. Dimensional analysis of the energy formulas

A physical dimension is a vector of exponents of mass, length and time. -/

structure Dim where
  M : Int
  L : Int
  T : Int
deriving DecidableEq, Repr

def Dim.mul (x y : Dim) : Dim := ⟨x.M + y.M, x.L + y.L, x.T + y.T⟩
def Dim.inv (x : Dim) : Dim := ⟨-x.M, -x.L, -x.T⟩

def energy : Dim := ⟨1, 2, -2⟩      -- joule
def action : Dim := ⟨1, 2, -1⟩      -- joule second
def planck : Dim := action          -- h
def freq   : Dim := ⟨0, 0, -1⟩      -- hertz
def time   : Dim := ⟨0, 0, 1⟩       -- second
def number : Dim := ⟨0, 0, 0⟩       -- dimensionless

/-- `h f` is an energy. -/
theorem hf_is_energy : planck.mul freq = energy := by decide

/-- `ΔE = h f Δt` is an action, not an energy. -/
theorem hfdt_is_action : (planck.mul freq).mul time = action ∧
    (planck.mul freq).mul time ≠ energy := by decide

/-- `ΔE = (2π h f₁/360) T(deg)` is an action (T(deg) a time, 2π/360 a number). -/
theorem second_formula_is_action :
    ((number.mul planck).mul freq).mul time = action := by decide

/-- `ΔE = (2π h/360) T(deg)/Δt` is an action. -/
theorem third_formula_is_action :
    (planck.mul time).mul time.inv = action := by decide

/-- `f₁ T(deg)` is a pure number. -/
theorem f_times_Tdeg_is_number : freq.mul time = number := by decide

/-! ## Section 2. The numbers of the audit

Constants are scaled to integers.  Planck's constant is
h = 662607015 × 10^(-42) J s (exact, SI 2019). -/

/-- The delay for one degree at 5 MHz is 1/(360 · 5·10^6) s = 555.55... ps:
    555.5 ps ≤ T(deg) < 555.6 ps, and 9 T(deg) = 5000 ps exactly. -/
theorem delay_one_degree :
    5555 * (360 * 5000000) ≤ 10 ^ 13 ∧ 10 ^ 13 < 5556 * (360 * 5000000) ∧
    9 * 10 ^ 12 = 5000 * (360 * 5000000) := by decide

/-- The product evaluated in the earlier paper:
    (6.626e-34)(4.986e6)(555e-12) = 18335665980 × 10^(-46) J s = 1.834e-36. -/
theorem hfdt_value :
    6626 * 4986 * 555 = 18335665980 ∧
    18335000000 ≤ 6626 * 4986 * 555 ∧ 6626 * 4986 * 555 < 18345000000 := by decide

/-- The stated value 1.848e-27 exceeds the correct product 1.8336e-36 by a
    factor between 1.0 × 10^9 and 1.01 × 10^9. -/
theorem factor_ten_to_nine :
    10 ^ 9 * 18335665980 ≤ 1848 * 10 ^ 16 ∧
    1848 * 10 ^ 16 < 101 * 10 ^ 7 * 18335665980 := by decide

/-- 1.848e-27 J would be 55.8% of the quantum E₁ = 3.313e-27 J. -/
theorem fifty_six_percent :
    5575 * 3313 ≤ 1848 * 10000 ∧ 1848 * 10000 < 5585 * 3313 := by decide

/-- E₁ = h f₁ at f₁ = 5 MHz: 662607015 × 5 = 3313035075, i.e.
    E₁ = 3.313035075 × 10^(-27) J. -/
theorem E1_value : 662607015 * 5 = 3313035075 := by decide

/-- f₂ = f₁ - f₁/360 = 4986111.11 Hz. -/
theorem f2_value :
    498611111 * 360 ≤ 5000000 * 359 * 100 ∧
    5000000 * 359 * 100 < 498611112 * 360 := by decide

/-- E₂ = h f₂ = (359/360) E₁ lies in [3.3035, 3.3045) × 10^(-27) J,
    so it is 3.304 × 10^(-27) J and not 3.313 × 10^(-27) J. -/
theorem E2_value :
    33035 * 10 ^ 5 * 360 ≤ 3313035075 * 359 ∧
    3313035075 * 359 < 33045 * 10 ^ 5 * 360 := by decide

/-- E₁ - E₂ = E₁/360 lies in [9.195, 9.205) × 10^(-30) J. -/
theorem E1_minus_E2 :
    9195000 * 360 ≤ 3313035075 ∧ 3313035075 < 9205000 * 360 := by decide

/-- The frequency shift of a phase slipping 1° per period at 5 MHz is
    f₁/360 = 13888.89 Hz. -/
theorem slip_frequency :
    1388888 * 360 ≤ 5000000 * 100 ∧ 5000000 * 100 < 1388889 * 360 := by decide

/-! ## Section 3. The series of Lambda

With c_n = C(2n,n)/4^n the paper proves
  Lambda = (pi/2) Σ_{n≥1} d_n κ^(2n-1),   d_n = c_n² · 2n/(2n-1).
Here d_n is a rational number with numerator C(2n,n)² · 2n and
denominator 16^n (2n-1). -/

def fact : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fact n

/-- The central binomial coefficient C(2n, n). -/
def cb (n : Nat) : Nat := fact (2 * n) / (fact n * fact n)

/-- A nonnegative rational number as a pair (numerator, denominator). -/
structure Q where
  num : Nat
  den : Nat

def Q.add (x y : Q) : Q := ⟨x.num * y.den + y.num * x.den, x.den * y.den⟩
def Q.le (x y : Q) : Bool := Nat.ble (x.num * y.den) (y.num * x.den)
def Q.lt (x y : Q) : Bool := Nat.blt (x.num * y.den) (y.num * x.den)
def Q.eq (x y : Q) : Bool := x.num * y.den == y.num * x.den

/-- c_n = C(2n,n)/4^n = Π_{k<n} (2k+1)/(2k+2), as numerator cn n and
    denominator cd n (not reduced).  The recurrence
    c_{n+1}/c_n = (2n+1)/(2n+2) holds by definition. -/
def cn : Nat → Nat
  | 0 => 1
  | n + 1 => cn n * (2 * n + 1)

def cd : Nat → Nat
  | 0 => 1
  | n + 1 => cd n * (2 * n + 2)

theorem cn_pos : ∀ n, 0 < cn n
  | 0 => by decide
  | n + 1 => Nat.mul_pos (cn_pos n) (Nat.succ_pos _)

theorem cd_pos : ∀ n, 0 < cd n
  | 0 => by decide
  | n + 1 => Nat.mul_pos (cd_pos n) (Nat.succ_pos _)

/-- d_n = c_n² · 2n/(2n-1), numerator and denominator. -/
def dnum (n : Nat) : Nat := cn n * cn n * (2 * n)
def dden (n : Nat) : Nat := cd n * cd n * (2 * n - 1)
def dcoef (n : Nat) : Q := ⟨dnum n, dden n⟩

/-- `cb` is the central binomial coefficient (the division in its
    definition is exact), and cn n / cd n = C(2n,n)/4^n, for n ≤ 60. -/
theorem cb_exact : (List.range 61).all (fun n =>
    cb n * (fact n * fact n) == fact (2 * n) && cn n * 4 ^ n == cb n * cd n) = true := by
  decide

/-- The first three coefficients are 1/2, 3/16, 15/128. -/
theorem first_coefficients :
    (dcoef 1).eq ⟨1, 2⟩ = true ∧ (dcoef 2).eq ⟨3, 16⟩ = true ∧
    (dcoef 3).eq ⟨15, 128⟩ = true := by decide

/-- After division by d_1 = 1/2 the expansion of Lambda / ((pi/4) κ) is
    1 + (3/8) κ² + (15/64) κ⁴ + ... . -/
theorem expansion_coefficients :
    (Q.eq ⟨3 * 2, 16⟩ ⟨3, 8⟩) = true ∧ (Q.eq ⟨15 * 2, 128⟩ ⟨15, 64⟩) = true := by
  decide

/-- The ratio d_{n+1}/d_n = (2n+1)(2n-1)/((2n+2)(2n)) is below 1: with
    m = 2n - 1 the inequality reads (m+2) m < (m+3)(m+1). -/
theorem ratio_lt_one (m : Nat) : (m + 2) * m < (m + 3) * (m + 1) := by
  have h1 : (m + 2) * m ≤ (m + 2) * (m + 1) :=
    Nat.mul_le_mul_left _ (Nat.le_succ m)
  have h2 : (m + 2) * (m + 1) < (m + 3) * (m + 1) :=
    Nat.mul_lt_mul_of_pos_right (Nat.lt_succ_self _) (Nat.succ_pos m)
  exact Nat.lt_of_le_of_lt h1 h2

/-- The coefficients decrease, for every n ≥ 1: d_{n+1} < d_n, i.e.
    dnum(n+1) dden(n) < dnum(n) dden(n+1) (stated for n = m + 1).  This is
    Proposition 7.2(ii) of the paper: d_n ≤ d_1 = 1/2. -/
theorem d_decreasing (m : Nat) :
    dnum (m + 2) * dden (m + 1) < dnum (m + 1) * dden (m + 2) := by
  have hP : 0 < cn (m+1) * cn (m+1) * (cd (m+1) * cd (m+1)) * ((2*m+3)*(2*m+4)) :=
    Nat.mul_pos (Nat.mul_pos (Nat.mul_pos (cn_pos _) (cn_pos _))
      (Nat.mul_pos (cd_pos _) (cd_pos _))) (Nat.mul_pos (Nat.succ_pos _) (Nat.succ_pos _))
  have e1 : dnum (m + 2) * dden (m + 1) =
      (cn (m+1) * cn (m+1) * (cd (m+1) * cd (m+1)) * ((2*m+3)*(2*m+4))) *
        ((2*m+3)*(2*m+1)) := by
    simp only [dnum, dden, cn, cd]
    have : 2 * (m + 1) - 1 = 2*m+1 := by omega
    rw [this]
    grind
  have e2 : dnum (m + 1) * dden (m + 2) =
      (cn (m+1) * cn (m+1) * (cd (m+1) * cd (m+1)) * ((2*m+3)*(2*m+4))) *
        ((2*m+2)*(2*m+4)) := by
    simp only [dnum, dden, cn, cd]
    have : 2 * (m + 2) - 1 = 2*m+3 := by omega
    rw [this]
    grind
  rw [e1, e2]
  apply Nat.mul_lt_mul_of_pos_left _ hP
  have h1 : (2*m+3)*(2*m+1) ≤ (2*m+3)*(2*m+2) := Nat.mul_le_mul_left _ (Nat.le_succ _)
  have h2 : (2*m+3)*(2*m+2) < (2*m+4)*(2*m+2) :=
    Nat.mul_lt_mul_of_pos_right (Nat.lt_succ_self _) (Nat.succ_pos _)
  have h3 : (2*m+4)*(2*m+2) = (2*m+2)*(2*m+4) := Nat.mul_comm _ _
  omega

theorem dden_pos (m : Nat) : 0 < dden (m + 1) := by
  unfold dden
  have : 2 * (m + 1) - 1 = 2*m+1 := by omega
  rw [this]
  exact Nat.mul_pos (Nat.mul_pos (cd_pos _) (cd_pos _)) (Nat.succ_pos _)

/-- Comparison of fractions by cross-multiplication is transitive when the
    middle denominator is positive. -/
theorem cross_trans {a b c d e f : Nat} (hd : 0 < d)
    (h1 : a * d ≤ c * b) (h2 : c * f ≤ e * d) : a * f ≤ e * b := by
  have s1 : a * d * f ≤ c * b * f := Nat.mul_le_mul_right _ h1
  have s2 : c * f * b ≤ e * d * b := Nat.mul_le_mul_right _ h2
  have s3 : (a * f) * d ≤ (e * b) * d := by
    have q1 : a * d * f = (a * f) * d := Nat.mul_right_comm _ _ _
    have q2 : c * b * f = c * f * b := Nat.mul_right_comm _ _ _
    have q3 : e * d * b = (e * b) * d := Nat.mul_right_comm _ _ _
    rw [q1, q2] at s1; rw [q3] at s2
    exact Nat.le_trans s1 s2
  exact Nat.le_of_mul_le_mul_right s3 hd

/-- d_{N+1+k} ≤ d_{N+1} for all N, k: the coefficient bound behind the tail
    estimate.  With it the tail Σ_{n>N} (1/2) d_n κ^(2n-1) is at most
    (1/2) d_{N+1} Σ_{n>N} κ^(2n-1) = (1/2) d_{N+1} κ^(2N+1)/(1-κ²); the
    summation of the geometric series is not formalised. -/
theorem d_le_first (N : Nat) :
    ∀ k, dnum (N + 1 + k) * dden (N + 1) ≤ dnum (N + 1) * dden (N + 1 + k)
  | 0 => Nat.le_refl _
  | k + 1 => by
    have ih := d_le_first N k
    have step : dnum (N + k + 2) * dden (N + k + 1) ≤
        dnum (N + k + 1) * dden (N + k + 2) := Nat.le_of_lt (d_decreasing (N + k))
    have e1 : N + 1 + (k + 1) = N + k + 2 := by omega
    have e2 : N + 1 + k = N + k + 1 := by omega
    rw [e1]; rw [e2] at ih
    exact cross_trans (dden_pos (N + k)) step ih

/-- The n-th term (1/2) d_n κ^(2n-1) of Lambda/pi at κ = p/q. -/
def term (p q n : Nat) : Q :=
  ⟨(dcoef n).num * p ^ (2 * n - 1), 2 * (dcoef n).den * q ^ (2 * n - 1)⟩

/-- The partial sum Σ_{n=1}^{N} (1/2) d_n κ^(2n-1). -/
def psum (p q : Nat) : Nat → Q
  | 0 => ⟨0, 1⟩
  | N + 1 => (psum p q N).add (term p q (N + 1))

/-- The tail bound (1/2) d_{N+1} κ^(2N+1) / (1 - κ²): since the d_n
    decrease, the tail Σ_{n>N} is at most the geometric series with first
    term (1/2) d_{N+1} κ^(2N+1) and ratio κ². -/
def tail (p q N : Nat) : Q :=
  ⟨(term p q (N + 1)).num * (q * q), (term p q (N + 1)).den * (q * q - p * p)⟩

/-- An enclosure lo ≤ S_N and S_N + tail ≤ hi; with the series of the
    paper it gives lo ≤ Lambda/pi ≤ hi. -/
def encloses (p q N : Nat) (lo hi : Q) : Bool :=
  Q.le lo (psum p q N) && Q.le ((psum p q N).add (tail p q N)) hi

/-- b/a = 1.01, κ = 1/201: Lambda/pi ∈ [0.001243792, 0.001243793]. -/
theorem lambda_1_01 : encloses 1 201 2 ⟨1243792, 10^9⟩ ⟨1243793, 10^9⟩ = true := by
  decide

/-- b/a = 1.1, κ = 1/21: Lambda/pi ∈ [0.011914899, 0.011914900]. -/
theorem lambda_1_1 : encloses 1 21 3 ⟨11914899, 10^9⟩ ⟨11914900, 10^9⟩ = true := by
  decide

/-- b/a = 11/9 (the band 4.5-5.5 MHz), κ = 1/10:
    Lambda/pi ∈ [0.025094340, 0.025094341]. -/
theorem lambda_5MHz : encloses 1 10 4 ⟨25094340, 10^9⟩ ⟨25094341, 10^9⟩ = true := by
  decide

/-- b/a = 1.5, κ = 1/5: Lambda/pi ∈ [0.050769314, 0.050769315]. -/
theorem lambda_1_5 : encloses 1 5 6 ⟨50769314, 10^9⟩ ⟨50769315, 10^9⟩ = true := by
  decide

/-- b/a = 2 (an octave), κ = 1/3: Lambda/pi ∈ [0.087068100, 0.087068101]. -/
theorem lambda_octave : encloses 1 3 9 ⟨87068100, 10^9⟩ ⟨87068101, 10^9⟩ = true := by
  decide

/-- b/a = 3, κ = 1/2: Lambda/pi ∈ [0.138966549, 0.138966550]. -/
theorem lambda_3 : encloses 1 2 14 ⟨138966549, 10^9⟩ ⟨138966550, 10^9⟩ = true := by
  decide

/-! ## Section 4. Identities of the test function and the moments

All variables are integers; an identity of rational functions is checked
after clearing denominators.  A = a², B = b², σ = t². -/

/-- ρ² = (z²-a²)/(z²-b²) = 1 + (b²-a²)/(z²-b²) (Lemma 5.1(i)), cleared:
    z² - a² = (z² - b²) + (b² - a²). -/
theorem rho_square (z a b : Int) :
    z*z - a*a = (z*z - b*b) + (b*b - a*a) := by grind

/-- The two forms of the weight p (equation (5.2)): with
    s = sqrt((σ-A)/(B-σ)), 1/s - s = (A+B-2σ)/sqrt((σ-A)(B-σ)), which after
    multiplication by sqrt((σ-A)(B-σ)) is (B - σ) - (σ - A) = A + B - 2σ. -/
theorem weight_two_forms (A B σ : Int) : (B - σ) - (σ - A) = A + B - 2*σ := by
  grind

/-- ∫ p/x = (1/2)((A+B) pi/sqrt(AB) - 2 pi) = pi (b-a)²/(2ab)
    (Lemma 5.3(iii)); the algebraic step is a² + b² - 2ab = (b-a)². -/
theorem moment_inverse (a b : Int) : a*a + b*b - 2*a*b = (b - a)*(b - a) := by
  grind

/-- ∫ x³ p = -(pi/8)(B-A)² (Lemma 5.3(iv)), from ∫σ/√ = pi(A+B)/2 and
    ∫σ²/√ = pi(3A²+2AB+3B²)/8:
    (1/2)[(A+B) pi (A+B)/2 - 2 pi (3A²+2AB+3B²)/8] = -(pi/8)(B-A)²,
    i.e. 2(A+B)² - (3A²+2AB+3B²) = -(B-A)². -/
theorem moment_cubic (A B : Int) :
    2*(A+B)*(A+B) - (3*A*A + 2*A*B + 3*B*B) = -((B - A)*(B - A)) := by grind

/-- The antisymmetry ϖ(A+B-σ) = -ϖ(σ) of Lemma 5.3(i): the numerator
    changes sign, A + B - 2(A+B-σ) = -(A+B-2σ), and the two factors of the
    denominator are exchanged, ((A+B-σ)-A)(B-(A+B-σ)) = (σ-A)(B-σ). -/
theorem varpi_antisymmetry (A B σ : Int) :
    A + B - 2*(A + B - σ) = -(A + B - 2*σ) ∧
    ((A + B - σ) - A) * (B - (A + B - σ)) = (σ - A) * (B - σ) := by
  constructor <;> grind

/-- W(0) = (ρ(0)-1)²/ρ(0) with ρ(0) = a/b equals (b-a)²/(ab); cleared:
    (a - b)² · (ab) = (b - a)² · (ab) and ((a/b) - 1)² / (a/b) =
    (a-b)²/(ab) amounts to (a - b)² b² · (ab) = (b-a)² · a b · b². -/
theorem W_at_zero (a b : Int) :
    (a - b)*(a - b) * (b*b) * (a*b) = (b - a)*(b - a) * (a*b) * (b*b) := by grind

/-! ## Section 5. The constant Lambda: Landen step and bounds

k' = a/b, κ = (b-a)/(b+a).  We write the identities in a, b. -/

/-- 1 + κ = 2b/(a+b): (b+a) + (b-a) = 2b. -/
theorem one_plus_kappa (a b : Int) : (b + a) + (b - a) = 2*b := by grind

/-- The Landen step of Proposition 7.1: b(1+κ)(1+k')² = 2(a+b), i.e.
    b · 2b · (a+b)² = 2(a+b) · (a+b) · b² after clearing (a+b) b². -/
theorem landen_K_coefficient (a b : Int) :
    b * (2*b) * ((a + b)*(a + b)) = 2*(a + b) * ((a + b) * (b*b)) := by grind

/-- The Landen step, E part: 2b(1+k') = 2(a+b), i.e. 2b(b+a) = 2(a+b) b. -/
theorem landen_E_coefficient (a b : Int) : 2*b*(b + a) = 2*(a + b)*b := by grind

/-- (1+k'²) + 2k' = (1+k')², used to combine (1+k'²)K(k) with -2E(k) =
    -2(1+k')E(κ) + 2k'K(k); in a, b: (b²+a²) + 2ab = (a+b)². -/
theorem one_plus_kprime_square (a b : Int) :
    (b*b + a*a) + 2*a*b = (a + b)*(a + b) := by grind

/-- 2(b-a) = 2κ(a+b), the divisor in Lambda = J/(2(b-a)). -/
theorem two_b_minus_a (a b : Int) : 2*(b - a)*(b + a) = 2*(b - a)*(a + b) := by grind

/-- The upper bound κ/(1-κ²) = (b²-a²)/(4ab) of Proposition 7.2(ii):
    (b+a)² - (b-a)² = 4ab. -/
theorem upper_bound_form (a b : Int) :
    (b + a)*(b + a) - (b - a)*(b - a) = 4*a*b := by grind

/-- The two forms of τ⋆ in (6.1): 2θK(k)/(πb) = 4θK(κ)/(π(a+b)), from
    K(k) = (1+κ)K(κ) and (1+κ)/b = 2/(a+b); cleared: 2 · 2b · (a+b) =
    4 · b · (a+b). -/
theorem tau_star_forms (a b : Int) : 2*(2*b)*(a + b) = 4*b*(a + b) := by grind

/-- Legendre's relation (Lemma A.1): kk'^2 dL/dk, expanded with the four
    derivative formulas, equals (k² + k'² - 1)(E K' - K E'), which vanishes
    since k² + k'² = 1.  Here E, K are at k and E1, K1 at k'. -/
theorem legendre_derivative (k k1 E K E1 K1 : Int) :
    k1*k1*(E - K)*K1 - E*(E1 - k*k*K1) + k*k*K*(K1 - E1)
      + E1*(E - k1*k1*K) - K1*(E - k1*k1*K) + K*(E1 - k*k*K1)
      = (k*k + k1*k1 - 1) * (E*K1 - K*E1) := by grind

/-- Landen's transformation for E (Proposition A.3): the coefficient of
    K(κ) is ((1-κ)² + 2κ(1-κ) - 2(1-κ)(1+κ))/(1+κ) = -(1-κ). -/
theorem landen_E_bracket (q : Int) :
    (1 - q)*(1 - q) + 2*q*(1 - q) - 2*(1 - q)*(1 + q) = -((1 - q)*(1 + q)) := by
  grind

/-- Landen's transformation for E: with k = 2√κ/(1+κ) and k' = (1-κ)/(1+κ),
    (dk/dκ)/(k k'²) = (1+κ)/(2κ(1-κ)).  Writing s = √κ (so κ = s²),
    dk/dκ = (1-s²)/(s(1+s²)²) and k k'² = 2s(1-s²)²/(1+s²)³; the claim,
    cleared of denominators, is the identity below. -/
theorem landen_dk (s : Int) :
    (1 - s*s) * ((1 + s*s)*(1 + s*s)*(1 + s*s)) * (2*(s*s)*(1 - s*s)) =
      (1 + s*s) * (s * ((1 + s*s)*(1 + s*s))) * (2*s*((1 - s*s)*(1 - s*s))) := by
  grind

/-- k' = (1-κ)/(1+κ) gives k'²(1+κ) = (1-κ)²/(1+κ), and 2/(1+κ) = 1 + k':
    (1+κ) + (1-κ) = 2. -/
theorem landen_kprime (q : Int) : (1 + q) + (1 - q) = 2 := by grind

/-- The mean of the extremal attenuation (Appendix C): with
    τ⋆ pi/4 = θK/(2b), τ⋆ pi m²/2 = θ b E, I0 = K/b, I2 = bE,
    τ⋆(pi/4)(a²+b²-2m²) = (θ/2)((a²+b²) I0 - 2 I2); the coefficients are
    (a²+b²)/(2b) = (a²+b²)·(1/b)/2 and b = 2b/2. -/
theorem extremal_mean (a b θ K E : Int) :
    θ*K*(a*a + b*b) - 2*θ*b*E*b = θ*((a*a + b*b)*K - 2*(b*E)*b) := by grind

/-- The narrow band: κ = Δω/(2ω̄), so (pi/4)κ = (pi/8) Δω/ω̄; and one degree
    is pi/180 rad, so 2θΛ ≈ 2 (pi/180)(pi/8) Δf/f = (pi²/720) Δf/f:
    the rational coefficients 1/4 · 1/2 = 1/8 and 2 · 1/180 · 1/8 = 1/720. -/
theorem narrow_band_constants : 1 * 8 = 4 * 2 ∧ 2 * 720 = 180 * 8 := by decide

/-! ## Section 6. Positivity of the coefficients and of J at the series level

For the record: every coefficient d_n is positive, which the paper uses
to conclude that Lambda increases with κ. -/

theorem coefficients_positive : ((List.range 60).map (· + 1)).all
    (fun n => Nat.blt 0 (dcoef n).num && Nat.blt 0 (dcoef n).den) = true := by
  decide

/-! ## Section 7. The lead branch and the zeros -/

/-- The inversion x = ab/ω of Theorem 6.6: r(ab/ω) = ab r(ω)/ω², from
    ((ab)² - a²ω²)(b²ω² - (ab)²) = a²b² (ω²-a²)(b²-ω²). -/
theorem inversion_r (a b w : Int) :
    ((a*b)*(a*b) - a*a*(w*w)) * (b*b*(w*w) - (a*b)*(a*b)) =
      a*a*(b*b) * ((w*w - a*a) * (b*b - w*w)) := by grind

/-- The weight after inversion: a² + b² - 2(ab/ω)² times ω² is
    (a²+b²)ω² - 2a²b², which is ω² times the numerator of p*. -/
theorem inversion_weight (a b w : Int) :
    (a*a + b*b)*(w*w) - 2*((a*b)*(a*b)) = (w*w)*(a*a + b*b) - 2*(a*a)*(b*b) := by
  grind

/-- The plain mean of the lead extremal (Remark 6.8):
    -|θ|Λ + (pi/4)(b-a)τ⋆ with Λ = (K-E)/κ and (pi/4)(b-a)τ⋆ = |θ|κK gives
    |θ|(κK - (K-E)/κ) = |θ|(E - (1-κ²)K)/κ; cleared by κ:
    κ²K - (K - E) = E - (1 - κ²)K. -/
theorem lead_mean (κ K E : Int) : κ*κ*K - (K - E) = E - (1 - κ*κ)*K := by grind

/-- (pi/4)(b-a)τ⋆ = |θ|κK: with τ⋆ = 4|θ|K/(pi(a+b)) one gets
    (pi/4)(b-a) · 4K/(pi(a+b)) = (b-a)K/(a+b) = κK; cleared:
    (b-a) · 4 · (b+a) = 4 · (b-a) · (a+b). -/
theorem lead_atom_cost (a b : Int) :
    (b - a) * 4 * (b + a) = 4 * (b - a) * (a + b) := by grind

/-- Corollary 10.3: (a/b)(pi-θ) > θ iff θ < a pi/(a+b); the two sides differ
    by the identity a(pi - θ) - bθ = a pi - (a+b)θ. -/
theorem threshold_identity (a b pi θ : Int) :
    a*(pi - θ) - b*θ = a*pi - (a + b)*θ := by grind

/-- The asymptotics (7.6): arctan u = u - u³/3 + ..., ∫x p = 0 and
    ∫x³p = -(pi/8)(B-A)², so δ(η) ~ (1/3)(pi/8)(B-A)²/((b-a)η³) =
    (pi/24)(B-A)²/((b-a)η³): the coefficient (1/3)(1/8) = 1/24. -/
theorem delta_coefficient : 1 * 24 = 3 * 8 := by decide

/-- The one-zero phase: g = θ - 2 arctan(ω/η) is nonnegative on the band
    iff 2 arctan(b/η) ≤ θ, i.e. η ≥ b cot(θ/2); the extremal choice
    η_θ = b cot(θ/2) solves 2 arctan(b/η) = θ.  The finite content is the
    half-angle identity tan(θ/2) = b/η_θ, recorded here as
    b = η_θ tan(θ/2) ⟺ η_θ = b cot(θ/2), i.e. b · cot = η and tan · cot = 1
    imply b = η · tan. -/
theorem half_angle (b η t c : Int) (h1 : η = b*c) (h2 : t*c = 1) :
    b = η*t := by
  subst h1
  grind

/-! ### The least loss with zeros (Section 8): the right angle -/

/-- `thm:rightangle`: R² + ((a+b)z)² = (z²+ab)² for R² = (z²-a²)(z²-b²), so the
    optimal system has modulus one off the band. -/
theorem right_modulus (a b z : Int) :
    (z*z - a*a)*(z*z - b*b) + (a + b)*(a + b)*(z*z) = (z*z + a*b)*(z*z + a*b) := by grind

/-- Q(z) = (z+a)(z+b) = z² + ab + (a+b)z. -/
theorem right_Q (a b z : Int) : (z + a)*(z + b) = z*z + a*b + (a + b)*z := by grind

/-- The two factorisations of R²: Q(z)(z-a)(z-b) = (z²-a²)(z²-b²). -/
theorem right_R_square (a b z : Int) :
    (z + a)*(z + b)*((z - a)*(z - b)) = (z*z - a*a)*(z*z - b*b) := by grind

/-- The real part of the cross-multiplication in the proof of `thm:rightangle`:
    R² + (a+b)zQ = (z²+ab)Q when R² = Q(z-a)(z-b). -/
theorem right_cross (a b z : Int) :
    (z + a)*(z + b)*((z - a)*(z - b)) + (a + b)*z*((z + a)*(z + b)) =
      (z*z + a*b)*((z + a)*(z + b)) := by grind

/-- On the band: ((a+b)ω)² - r² = (ω²+ab)² with r² = (ω²-a²)(b²-ω²). -/
theorem band_modulus (a b w : Int) :
    ((a + b)*w)*((a + b)*w) - (w*w - a*a)*(b*b - w*w) = (w*w + a*b)*(w*w + a*b) := by grind

/-- (a+b)ω - (ω²+ab) = (ω-a)(b-ω) and (a+b)ω + (ω²+ab) = (ω+a)(ω+b), whose
    product is r². -/
theorem band_factors (a b w : Int) :
    (a + b)*w - (w*w + a*b) = (w - a)*(b - w) ∧
    (a + b)*w + (w*w + a*b) = (w + a)*(w + b) ∧
    (w - a)*(b - w)*((w + a)*(w + b)) = (w*w - a*a)*(b*b - w*w) := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- The second form of (8.4): if r² = PQ then (P+r)²Q = (Q+r)²P, that is
    (P+r)/(Q+r) = √(P/Q). -/
theorem tanh_form (P Q r : Int) (h : r*r = P*Q) :
    (P + r)*(P + r)*Q = (Q + r)*(Q + r)*P := by
  have e : (P + r)*(P + r)*Q - (Q + r)*(Q + r)*P = (P - Q)*(P*Q - r*r) := by grind
  have : (P + r)*(P + r)*Q - (Q + r)*(Q + r)*P = 0 := by rw [e, h]; grind
  grind

/-- The footnote to the first form of (8.4): with c = (a+b)ω, m = ω²+ab and
    c² - m² = r², (Q+r)m = (Q-r)(c+r) for Q = m + c. -/
theorem first_form (c m r : Int) (h : c*c - m*m = r*r) :
    (m + c + r)*m = (m + c - r)*(c + r) := by
  have e : (m + c - r)*(c + r) - (m + c + r)*m = (c*c - m*m) - r*r := by grind
  grind

/-- c² - m² = r²: (a+b)²ω² - (ω²+ab)² = (ω²-a²)(b²-ω²). -/
theorem first_form_hyp (a b w : Int) :
    ((a + b)*w)*((a + b)*w) - (w*w + a*b)*(w*w + a*b) = (w*w - a*a)*(b*b - w*w) := by grind

/-- The symmetry ω ↦ ab/ω (multiplied out by ω²): P and Q both acquire the
    factor ab, so P/Q and α♯ are invariant. -/
theorem right_symmetry (a b w : Int) :
    (a*b - a*w)*(b*w - a*b) = a*b*((w - a)*(b - w)) ∧
    (a*b + a*w)*(b*w + a*b) = a*b*((w + a)*(w + b)) := by
  refine ⟨?_, ?_⟩ <;> grind

/-- At the geometric centre ω² = ab: r = ω(b-a), and e^{2α} = b/a. -/
theorem right_centre (a b w : Int) (h : w*w = a*b) :
    (w*w - a*a)*(b*b - w*w) = (w*(b - a))*(w*(b - a)) ∧
    ((a + b)*w + w*(b - a))*((a + b)*w + w*(b - a))*a = (w*w + a*b)*(w*w + a*b)*b := by
  constructor
  · have e : (w*(b - a))*(w*(b - a)) = (w*w)*((b - a)*(b - a)) := by grind
    rw [e, h]; grind
  · have hw : w*w = a*b := h
    have e1 : ((a + b)*w + w*(b - a)) = 2*b*w := by grind
    rw [e1, hw]
    have e2 : 2*b*w*(2*b*w)*a = 4*(b*b)*a*(w*w) := by grind
    rw [e2, hw]; grind

/-- `eq:lambdaright`: (π/2)(√b-√a)²/(b-a) = (π/2)(√b-√a)/(√b+√a), from
    (q-p)²(q+p) = (q-p)(q²-p²). -/
theorem right_mean (p q : Int) : (q - p)*(q - p)*(q + p) = (q - p)*(q*q - p*p) := by grind

/-- For an octave: (√2-1)² = 3-2√2 in ℤ[√2], (x+y√2)² = (x²+2y²) + 2xy√2. -/
theorem octave_right : (-1)*(-1) + 2*(1*1) = (3 : Int) ∧ 2*(-1)*1 = (-2 : Int) := by decide

/-- Lemma 9.1: for t off the band, R(t) and t² - x² have the same sign, so
    the kernel k_t is positive; below the band both factors are negative. -/
theorem kernel_sign (R d : Int) :
    (R < 0 → d < 0 → 0 < R*d) ∧ (0 < R → 0 < d → 0 < R*d) :=
  ⟨fun h1 h2 => Int.mul_pos_of_neg_of_neg h1 h2, fun h1 h2 => Int.mul_pos h1 h2⟩

/-- Corollary 9.3: for 0 < t < a < x, 2x/(x² - t²) exceeds 2/x because,
    cleared by x(x² - t²) > 0, 2x² - 2(x² - t²) = 2t² > 0. -/
theorem low_kernel_gap (x t : Int) :
    2*(x*x) - 2*(x*x - t*t) = 2*(t*t) ∧ 0 ≤ t*t := by
  refine ⟨by grind, ?_⟩
  rcases Int.le_total 0 t with h | h
  · exact Int.mul_nonneg h h
  · exact Int.mul_nonneg_of_nonpos_of_nonpos h h

/-- Corollary 9.7(iv): ∫ r/x = (pi/2)(√B - √A)²/2 with A = a², B = b²;
    the step 2c - 2√(AB) = (√B - √A)², c = (A + B)/2, with p = √A, q = √B. -/
theorem r_over_x_square (p q : Int) : p*p + q*q - 2*(p*q) = (q - p)*(q - p) := by grind

/-- Theorems 9.4 and 9.6: at the threshold the compensated linear
    coefficient θτ∘ + τ + C vanishes, C = -θτ∘ - τ, and for a lead
    -sτ∘ + τ + (sτ∘ - τ) = 0. -/
theorem compensated_linear (θτ τ s : Int) :
    θτ + τ + (-θτ - τ) = 0 ∧ -s + τ + (s - τ) = 0 := by
  refine ⟨?_, ?_⟩ <;> omega

/-- Theorem 9.8: the two mean costs ϑΛ and (pi - ϑ)Λ₋ are equal exactly
    when ϑ(Λ + Λ₋) = pi Λ₋, which defines ϑ_c. -/
theorem branch_threshold (pi t L M : Int) :
    t*L - (pi - t)*M = t*(L + M) - pi*M := by grind

end PhaseShiftCost

open PhaseShiftCost in
section
#print axioms hf_is_energy
#print axioms hfdt_is_action
#print axioms second_formula_is_action
#print axioms third_formula_is_action
#print axioms f_times_Tdeg_is_number
#print axioms delay_one_degree
#print axioms hfdt_value
#print axioms factor_ten_to_nine
#print axioms fifty_six_percent
#print axioms E1_value
#print axioms f2_value
#print axioms E2_value
#print axioms E1_minus_E2
#print axioms slip_frequency
#print axioms cn_pos
#print axioms cd_pos
#print axioms dden_pos
#print axioms cb_exact
#print axioms first_coefficients
#print axioms expansion_coefficients
#print axioms ratio_lt_one
#print axioms d_decreasing
#print axioms cross_trans
#print axioms d_le_first
#print axioms lambda_1_01
#print axioms lambda_1_1
#print axioms lambda_5MHz
#print axioms lambda_1_5
#print axioms lambda_octave
#print axioms lambda_3
#print axioms rho_square
#print axioms weight_two_forms
#print axioms moment_inverse
#print axioms moment_cubic
#print axioms varpi_antisymmetry
#print axioms W_at_zero
#print axioms one_plus_kappa
#print axioms landen_K_coefficient
#print axioms landen_E_coefficient
#print axioms one_plus_kprime_square
#print axioms two_b_minus_a
#print axioms upper_bound_form
#print axioms tau_star_forms
#print axioms legendre_derivative
#print axioms landen_E_bracket
#print axioms landen_dk
#print axioms landen_kprime
#print axioms extremal_mean
#print axioms narrow_band_constants
#print axioms coefficients_positive
#print axioms inversion_r
#print axioms inversion_weight
#print axioms lead_mean
#print axioms lead_atom_cost
#print axioms threshold_identity
#print axioms delta_coefficient
#print axioms half_angle
#print axioms right_modulus
#print axioms right_Q
#print axioms right_R_square
#print axioms right_cross
#print axioms band_modulus
#print axioms band_factors
#print axioms tanh_form
#print axioms first_form
#print axioms first_form_hyp
#print axioms right_symmetry
#print axioms right_centre
#print axioms right_mean
#print axioms octave_right
#print axioms kernel_sign
#print axioms low_kernel_gap
#print axioms r_over_x_square
#print axioms compensated_linear
#print axioms branch_threshold
end
