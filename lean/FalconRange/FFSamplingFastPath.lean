module
public import Mathlib
public import FalconRange.FixedPoint
public import FalconRange.ErrorPropagation
@[expose] public section

/-!
# Fast-path range proofs for recursive masked ffSampling

The C backend uses these theorems only for the recursive LDL and iXadj
call sites. Root complex multiplication and the root self-adjoint division
remain on the signed-192-bit fallback.
-/

namespace FalconRange

/-- The current signer accepts only keys with this squared norm bound. -/
def acceptedFGNorm : ℤ := 16822

def maxStandardDegree : ℤ := 1024

theorem root_gram_numeric_budget :
    maxStandardDegree * acceptedFGNorm < pow2 25 := by
  norm_num [maxStandardDegree, acceptedFGNorm, pow2]

/-- Cauchy--Schwarz bound for an unnormalised Fourier evaluation with
unit-magnitude (or smaller) twiddles. -/
theorem weighted_fourier_sq_le_card_mul_energy
    {ι : Type*} [DecidableEq ι] (s : Finset ι)
    (f : ι → ℝ) (w : ι → ℂ)
    (hw : ∀ i ∈ s, ‖w i‖ ≤ 1) :
    ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * w i‖ ^ 2 ≤
      (s.card : ℝ) * ∑ i ∈ s, (f i) ^ 2 := by
  have hnorm : ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * w i‖ ≤
      ∑ i ∈ s, |f i| := by
    calc
      ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * w i‖ ≤
          ∑ i ∈ s, ‖((f i : ℝ) : ℂ) * w i‖ := norm_sum_le _ _
      _ ≤ ∑ i ∈ s, |f i| := by
        apply Finset.sum_le_sum
        intro i hi
        calc
          ‖((f i : ℝ) : ℂ) * w i‖ = |f i| * ‖w i‖ := by
            rw [norm_mul, Complex.norm_real, Real.norm_eq_abs]
          _ ≤ |f i| * 1 := mul_le_mul_of_nonneg_left (hw i hi) (abs_nonneg _)
          _ = |f i| := by ring
  have hsum0 : 0 ≤ ∑ i ∈ s, |f i| := Finset.sum_nonneg fun _ _ => abs_nonneg _
  have hsq : ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * w i‖ ^ 2 ≤
      (∑ i ∈ s, |f i|) ^ 2 := by
    nlinarith [norm_nonneg (∑ i ∈ s, ((f i : ℝ) : ℂ) * w i)]
  calc
    ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * w i‖ ^ 2 ≤
        (∑ i ∈ s, |f i|) ^ 2 := hsq
    _ ≤ (s.card : ℝ) * ∑ i ∈ s, |f i| ^ 2 :=
      sq_sum_le_card_mul_sum_sq
    _ = (s.card : ℝ) * ∑ i ∈ s, (f i) ^ 2 := by
      congr 1
      apply Finset.sum_congr rfl
      intro i hi
      simp [sq_abs]

/-- The two-polynomial Gram value is below 2^25 for both standard degrees,
using only the key norm check already performed by the signer. -/
theorem accepted_root_gram_lt_pow2_25
    {ι : Type*} [DecidableEq ι] (s : Finset ι)
    (f g : ι → ℝ) (wf wg : ι → ℂ)
    (hwf : ∀ i ∈ s, ‖wf i‖ ≤ 1)
    (hwg : ∀ i ∈ s, ‖wg i‖ ≤ 1)
    (hcard : s.card ≤ 1024)
    (henergy : (∑ i ∈ s, (f i) ^ 2) +
      (∑ i ∈ s, (g i) ^ 2) ≤ 16822) :
    ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * wf i‖ ^ 2 +
      ‖∑ i ∈ s, ((g i : ℝ) : ℂ) * wg i‖ ^ 2 <
      (2 : ℝ) ^ 25 := by
  have hf := weighted_fourier_sq_le_card_mul_energy s f wf hwf
  have hg := weighted_fourier_sq_le_card_mul_energy s g wg hwg
  have hc : (s.card : ℝ) ≤ 1024 := by exact_mod_cast hcard
  calc
    ‖∑ i ∈ s, ((f i : ℝ) : ℂ) * wf i‖ ^ 2 +
        ‖∑ i ∈ s, ((g i : ℝ) : ℂ) * wg i‖ ^ 2 ≤
        (s.card : ℝ) * ((∑ i ∈ s, (f i) ^ 2) +
          (∑ i ∈ s, (g i) ^ 2)) := by nlinarith
    _ ≤ 1024 * 16822 := by
      have he0 : 0 ≤ (∑ i ∈ s, (f i) ^ 2) +
          (∑ i ∈ s, (g i) ^ 2) := by positivity
      nlinarith
    _ < (2 : ℝ) ^ 25 := by norm_num

/-- Arithmetic-mean child of a positive split stays within its parent
magnitude budget. -/
theorem arithmetic_mean_preserves_budget {x y B : ℝ}
    (hx0 : 0 ≤ x) (hy0 : 0 ≤ y) (hx : x ≤ B) (hy : y ≤ B) :
    0 ≤ (x + y) / 2 ∧ (x + y) / 2 ≤ B := by
  constructor <;> nlinarith

/-- Harmonic child of a positive split stays within the larger parent. -/
theorem harmonic_child_preserves_budget {x y B : ℝ}
    (hx0 : 0 < x) (hy0 : 0 < y) (hx : x ≤ B) (hy : y ≤ B) :
    0 < (2 * x * y) / (x + y) ∧ (2 * x * y) / (x + y) ≤ B := by
  have hxy : 0 < x + y := by linarith
  constructor
  · positivity
  · apply (div_le_iff₀ hxy).2
    by_cases h : x ≤ y
    · nlinarith
    · have hyx : y ≤ x := le_of_not_ge h
      nlinarith

/-- The non-root LDL correlation ratio has magnitude strictly below one. -/
theorem correlation_ratio_abs_lt_one {x y : ℝ}
    (hx : 0 < x) (hy : 0 < y) :
    |(x - y) / (x + y)| < 1 := by
  rw [abs_lt]
  constructor
  · apply (lt_div_iff₀ (by linarith : 0 < x + y)).2
    nlinarith
  · apply (div_lt_iff₀ (by linarith : 0 < x + y)).2
    nlinarith

/-- SamplerZ writes z=floor(t)+zi with -18≤zi≤19; the residual is therefore
strictly smaller than 20 independently of the sampled trace. -/
theorem samplerZ_residual_lt_20 {base zi : ℤ} {r : ℝ}
    (hr0 : 0 ≤ r) (hr1 : r < 1) (hzi0 : -18 ≤ zi) (hzi1 : zi ≤ 19) :
    |(((base : ℝ) + r) - ((base + zi : ℤ) : ℝ))| < 20 := by
  have hzi0r : (-18 : ℝ) ≤ (zi : ℝ) := by exact_mod_cast hzi0
  have hzi1r : (zi : ℝ) ≤ 19 := by exact_mod_cast hzi1
  rw [abs_lt]
  norm_num [Int.cast_add]
  constructor <;> nlinarith

/-- A length-k linear transform of residuals bounded by 20 has magnitude at
most 20*k. This is the split/merge envelope used by iXadj. -/
theorem residual_sum_le_twenty_card {ι : Type*} [DecidableEq ι]
    (s : Finset ι) (v : ι → ℝ)
    (hv : ∀ i ∈ s, |v i| ≤ 20) :
    |∑ i ∈ s, v i| ≤ 20 * s.card := by
  calc
    |∑ i ∈ s, v i| ≤ ∑ i ∈ s, |v i| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i ∈ s, (20 : ℝ) := Finset.sum_le_sum hv
    _ = 20 * s.card := by simp [mul_comm]

/-- Ten levels of a conservative doubling recurrence, with one ulp local
error per level, stay below one Q32.64 real unit. -/
theorem depth10_rounding_error_lt_one :
    errorRec 2 ((2 : ℝ) ^ (-64 : ℤ)) 10 < 1 := by
  rw [errorRec_two]
  norm_num

/-- A generic bridge from the exact analytic value to the fixed-point value.
It is instantiated below with the depth-ten error bound. -/
theorem abs_approx_lt_of_abs_error {x xhat B eps C : ℝ}
    (hx : |x| ≤ B) (herr : |xhat - x| < eps) (hbudget : B + eps < C) :
    |xhat| < C := by
  calc
    |xhat| = |x + (xhat - x)| := by ring_nf
    _ ≤ |x| + |xhat - x| := abs_add_le _ _
    _ < B + eps := add_lt_add_of_le_of_lt hx herr
    _ < C := hbudget

/-- The standard-key Gram budget retains more than one full real unit of
rounding margin at the 2^25 boundary. -/
theorem root_gram_with_rounding_lt_pow2_25 {x xhat : ℝ}
    (hx : |x| ≤ 1024 * 16822) (herr : |xhat - x| < 1) :
    |xhat| < (2 : ℝ) ^ 25 := by
  apply abs_approx_lt_of_abs_error hx herr
  norm_num

/-- A merged length-1024 SamplerZ residual retains a full-unit rounding
margin below 2^15. -/
theorem residual_with_rounding_lt_pow2_15 {x xhat : ℝ}
    (hx : |x| ≤ 20 * 1024) (herr : |xhat - x| < 1) :
    |xhat| < (2 : ℝ) ^ 15 := by
  apply abs_approx_lt_of_abs_error hx herr
  norm_num

/-- Rounding can move the strict correlation ratio across one, but never
outside the conservative magnitude-two envelope. -/
theorem rounded_correlation_lt_two {x xhat : ℝ}
    (hx : |x| < 1) (herr : |xhat - x| < 1) :
    |xhat| < 2 := by
  calc
    |xhat| = |x + (xhat - x)| := by ring_nf
    _ ≤ |x| + |xhat - x| := abs_add_le _ _
    _ < 1 + 1 := add_lt_add hx herr
    _ = 2 := by norm_num

/-- The two terms forming w in iXadj fit the conservative 2^26 envelope. -/
theorem ixadj_w_lt_pow2_26 {a mua : ℝ}
    (ha : |a| < (2 : ℝ) ^ 15) (hmua : |mua| < (2 : ℝ) ^ 25) :
    |a + mua| < (2 : ℝ) ^ 26 := by
  calc
    |a + mua| ≤ |a| + |mua| := abs_add_le _ _
    _ < (2 : ℝ) ^ 15 + (2 : ℝ) ^ 25 := add_lt_add ha hmua
    _ < (2 : ℝ) ^ 26 := by norm_num

/-- Numerical envelopes used by the C call-site map. -/
theorem ldl_theta_real_budget : 4 * 1024 * 16822 < (2 : ℤ) ^ 27 := by norm_num

theorem ixadj_mu_real_budget : 20 * (1024 : ℤ) ^ 2 < (2 : ℤ) ^ 25 := by norm_num

/-- LDL public theta times a diagonal: raw exponents 66+89=155. -/
theorem ldl_theta_a0_fits_s160 {theta a0 : ℤ}
    (htheta : |theta| < pow2 66) (ha0 : |a0| < pow2 89) :
    InSigned 160 (theta * a0) := by
  apply product_fits_s160_of_magnitude htheta ha0
  norm_num

/-- The rounded LDL coefficient is below two in real magnitude and the
off-diagonal Gram value is below 2^25. Their Q32.64 representatives therefore
have raw exponents 65 and 89; this derives the Schur-product bound instead of
assuming that bound as a premise. -/
theorem ldl_schur_product_fits_s160 {e b1 : ℤ}
    (he : |e| < pow2 65) (hb1 : |b1| < pow2 89) :
    InSigned 160 (e * b1) := by
  apply product_fits_s160_of_magnitude he hb1
  norm_num

/-- iXadj theta has raw exponent 66, and its rounded LDL correlation operand
has raw exponent 65. Thus the exact product before the public half has raw
exponent 131. -/
theorem ixadj_theta_b_fits_s160 {theta b : ℤ}
    (htheta : |theta| < pow2 66) (hb : |b| < pow2 65) :
    InSigned 160 (theta * b) := by
  apply product_fits_s160_of_magnitude htheta hb
  norm_num

/-- The public mu table is below 2^10 and the residual transform below 2^15
in real units: raw exponents 74+79=153. -/
theorem ixadj_mu_a_fits_s160 {mu a : ℤ}
    (hmu : |mu| < pow2 74) (ha : |a| < pow2 79) :
    InSigned 160 (mu * a) := by
  apply product_fits_s160_of_magnitude hmu ha
  norm_num

/-- Final iXadj multiply: conservative real bounds 4 and 2^26 become raw
exponents 66 and 90, for a 156-bit exact product. -/
theorem ixadj_final_fits_s160 {bt w : ℤ}
    (hbt : |bt| < pow2 66) (hw : |w| < pow2 90) :
    InSigned 160 (bt * w) := by
  apply product_fits_s160_of_magnitude hbt hw
  norm_num

end FalconRange
