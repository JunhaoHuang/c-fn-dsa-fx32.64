module
public import Mathlib
public import FalconRange.FixedPoint
public import FalconRange.RangeBudgets
@[expose] public section

/-!
# Range obligations for the masked dynamic sampler backend
This module records the integer bounds that decide whether a C call site may
use the 160-bit ring multiply. Recursive LDL/iXadj sites are enabled only
through the theorems in FFSamplingFastPath; root complex products stay on the
192-bit fallback. No empirical maximum is promoted to a premise.

All exponents below refer to raw signed integers. A real bound 2^m in Q32.64
therefore has raw exponent m+64.
-/

namespace FalconRange

def InUQ0_64 (x : ℤ) : Prop := 0 ≤ x ∧ x < pow2 64

/-- Arithmetic halving cannot take a Q32.64 representative out of range. -/
theorem q32_64_half_fits {x : ℤ} (hx : InQ32_64 x) :
    InQ32_64 (x / 2) := by
  rcases hx with ⟨hlo, hhi⟩
  constructor <;>
    norm_num [InQ32_64, InSigned, pow2] at hlo hhi ⊢ <;> omega

/-- A genuine unsigned Q0.64 public twiddle times any Q32.64 value has an
exact raw product in signed 160 bits. This is the theorem attached to the
x_SIN and IRNJ fast call sites in mask_xpoly.c. -/
theorem q0_64_twiddle_product_fits_s160 {x c : ℤ}
    (hx : InQ32_64 x) (hc : InUQ0_64 c) :
    InSigned 160 (x * c) := by
  have hax : |x| ≤ pow2 95 := q32_64_abs_le hx
  have hcabs : |c| < pow2 64 := by
    rw [abs_of_nonneg hc.1]
    exact hc.2
  by_cases hz : x = 0
  · simp [hz, InSigned, pow2]
  · have hxpos : 0 < |x| := abs_pos.mpr hz
    have hprod : |x * c| < pow2 159 := by
      rw [abs_mul]
      calc
        |x| * |c| < |x| * pow2 64 :=
          mul_lt_mul_of_pos_left hcabs hxpos
        _ ≤ pow2 95 * pow2 64 :=
          mul_le_mul_of_nonneg_right hax (pow2_pos _).le
        _ = pow2 159 := by norm_num [pow2, pow_add]
    rw [abs_lt] at hprod
    exact ⟨le_of_lt hprod.1, hprod.2⟩

/-- Denominator/reciprocal products in the normalized Newton division use
at most raw exponents 64 and 66. -/
theorem div_newton_den_product_fits_s160 {b x : ℤ}
    (hb : |b| < pow2 64) (hx : |x| < pow2 66) :
    InSigned 160 (b * x) := by
  apply product_fits_s160_of_magnitude hb hx
  norm_num

/-- The second multiplication of a division Newton round uses two values
with raw exponent at most 66. -/
theorem div_newton_update_product_fits_s160 {x y : ℤ}
    (hx : |x| < pow2 66) (hy : |y| < pow2 66) :
    InSigned 160 (x * y) := by
  apply product_fits_s160_of_magnitude hx hy
  norm_num

/-- The documented square-root input interval implies a raw exponent below
79 before the public twelve-bit down-scaling. -/
theorem sqrt_input_raw_lt_2pow79 {a : ℤ}
    (hlo : 6327 * pow2 64 < a)
    (hhi : a < 16823 * pow2 64) :
    |a| < pow2 79 := by
  have ha0 : 0 < a := by
    have : 0 < 6327 * pow2 64 := mul_pos (by norm_num) (pow2_pos 64)
    omega
  rw [abs_of_pos ha0]
  calc
    a < 16823 * pow2 64 := hhi
    _ < pow2 79 := by norm_num [pow2]

/-- After the public twelve-bit scale, the square-root argument has raw
exponent below 67. -/
theorem sqrt_scaled_raw_lt_2pow67 {a : ℤ}
    (ha0 : 0 ≤ a) (ha : a < pow2 79) :
    |a / pow2 12| < pow2 67 := by
  rw [abs_of_nonneg (Int.ediv_nonneg ha0 (pow2_pos 12).le)]
  norm_num [pow2] at ha ⊢
  omega

/-- Local product budget used by the eight inverse-square-root rounds. -/
theorem sqrt_newton_product_fits_s160 {x y : ℤ}
    (hx : |x| < pow2 67) (hy : |y| < pow2 67) :
    InSigned 160 (x * y) := by
  apply product_fits_s160_of_magnitude hx hy
  norm_num

/-- Multiplication by the public inverse-sigma and inverse-q constants is a
Q0.64 fast call whenever the secret operand remains Q32.64. -/
theorem sampler_public_scale_fits_s160 {x c : ℤ}
    (hx : InQ32_64 x) (hc : InUQ0_64 c) :
    InSigned 160 (x * c) :=
  q0_64_twiddle_product_fits_s160 hx hc

/-- One linear FFT layer raises a real-magnitude exponent by at most one;
public Q0.64 multiplication itself does not raise it. -/
def fftMagnitudeExponent (base depth : ℕ) : ℕ := base + depth

@[simp] theorem fftMagnitudeExponent_zero (base : ℕ) :
    fftMagnitudeExponent base 0 = base := by
  simp [fftMagnitudeExponent]

theorem fftMagnitudeExponent_succ (base depth : ℕ) :
    fftMagnitudeExponent base (depth + 1) =
      fftMagnitudeExponent base depth + 1 := by
  simp [fftMagnitudeExponent, Nat.add_assoc]

/-- A nine-layer FN-DSA-512 FFT remains eligible for Q0.64 fast twiddles when
its real input exponent is at most 22. -/
theorem fft512_public_twiddle_budget {base : ℕ} (hbase : base ≤ 22) :
    fftMagnitudeExponent base 9 ≤ 31 := by
  simp [fftMagnitudeExponent]
  omega

/-- A ten-layer FN-DSA-1024 FFT remains eligible for Q0.64 fast twiddles when
its real input exponent is at most 21. -/
theorem fft1024_public_twiddle_budget {base : ℕ} (hbase : base ≤ 21) :
    fftMagnitudeExponent base 10 ≤ 31 := by
  simp [fftMagnitudeExponent]
  omega

/-- Executable C call-site policy. The product exponent records a strict
upper bound on the exact, unshifted raw product, not merely an output bound. -/
inductive MaskedMulSite
  | publicQ0_64
  | divisionNewton
  | sqrtNewton
  | ldlThetaA0
  | ldlSchur
  | ixadjThetaB
  | ixadjMuA
  | ixadjFinal
  | rootComplex
  | rootSelfAdjDiv
  deriving DecidableEq

def MaskedMulSite.fastEnabled : MaskedMulSite → Bool
  | .publicQ0_64 | .divisionNewton | .sqrtNewton
  | .ldlThetaA0 | .ldlSchur | .ixadjThetaB | .ixadjMuA
  | .ixadjFinal => true
  | .rootComplex | .rootSelfAdjDiv => false

def MaskedMulSite.rawProductExponent : MaskedMulSite → ℕ
  | .publicQ0_64 => 159
  | .divisionNewton => 132
  | .sqrtNewton => 134
  | .ldlThetaA0 => 155
  | .ldlSchur => 154
  | .ixadjThetaB => 131
  | .ixadjMuA => 153
  | .ixadjFinal => 156
  | .rootComplex | .rootSelfAdjDiv => 0

theorem enabled_site_raw_product_exponent_lt_160 (site : MaskedMulSite)
    (h : site.fastEnabled = true) :
    site.rawProductExponent < 160 := by
  cases site <;> simp [MaskedMulSite.fastEnabled,
    MaskedMulSite.rawProductExponent] at h ⊢

/-- Every arbitrary Q32.64-by-Q32.64 call remains functionally safe through
the 192-bit fallback, including unclosed root sites and general polynomial APIs. -/
theorem dynamic_fallback_covers_q32_products {x y : ℤ}
    (hx : InQ32_64 x) (hy : InQ32_64 y) :
    InSigned 192 (x * y) :=
  q32_64_product_fits_s192 hx hy

end FalconRange
