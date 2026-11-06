module
public import Mathlib
@[expose] public section

/-!
# Exact fixed-point model used by `c-fn-dsa-fx32.64`

An `rx` value is a signed 96-bit integer scaled by `2^64`.  The masked
backend stores arithmetic shares in `ZMod (2^160)`; the 192-bit fallback
multiplies sign-extended 96-bit representatives before applying the same
rounding rule as `rx_mul`.
-/

namespace FalconRange

def pow2 (n : ℕ) : ℤ := (2 : ℤ) ^ n

def InSigned (w : ℕ) (x : ℤ) : Prop :=
  -(pow2 (w - 1)) ≤ x ∧ x < pow2 (w - 1)

def InQ32_64 (x : ℤ) : Prop := InSigned 96 x

def roundShift (n : ℕ) (x : ℤ) : ℤ :=
  (x + pow2 (n - 1)) / pow2 n

def roundShift64 (x : ℤ) : ℤ := roundShift 64 x

/-- Real value denoted by a signed Q32.64 raw integer. -/
noncomputable def decodeQ32_64 (x : ℤ) : ℝ := (x : ℝ) / (2 : ℝ) ^ 64

@[simp] theorem pow2_pos (n : ℕ) : 0 < pow2 n := by
  simp [pow2]

theorem q32_64_fits_s160 {x : ℤ} (hx : InQ32_64 x) : InSigned 160 x := by
  rcases hx with ⟨hx0, hx1⟩
  constructor <;> norm_num [InQ32_64, InSigned, pow2] at hx0 hx1 ⊢ <;> omega

theorem q32_64_abs_le {x : ℤ} (hx : InQ32_64 x) : |x| ≤ pow2 95 := by
  rcases hx with ⟨hx0, hx1⟩
  rw [abs_le]
  constructor
  · exact hx0
  · exact le_of_lt hx1

theorem q32_64_product_fits_s192 {x y : ℤ}
    (hx : InQ32_64 x) (hy : InQ32_64 y) : InSigned 192 (x * y) := by
  have hax := q32_64_abs_le hx
  have hay := q32_64_abs_le hy
  have habs : |x * y| ≤ pow2 190 := by
    rw [abs_mul]
    calc
      |x| * |y| ≤ pow2 95 * pow2 95 := mul_le_mul hax hay (abs_nonneg _) (pow2_pos _).le
      _ = pow2 190 := by norm_num [pow2, pow_add]
  rw [abs_le] at habs
  rcases habs with ⟨hlo, hhi⟩
  constructor <;> norm_num [InSigned, pow2] at hlo hhi ⊢ <;> omega

/-- A sufficient fast-path condition: strict magnitude exponents whose sum is
at most 159 make the exact raw product fit in signed 160 bits. -/
theorem product_fits_s160_of_magnitude {x y : ℤ} {mx my : ℕ}
    (hx : |x| < pow2 mx) (hy : |y| < pow2 my) (hsum : mx + my ≤ 159) :
    InSigned 160 (x * y) := by
  have habs : |x * y| < pow2 159 := by
    rw [abs_mul]
    have hp : pow2 (mx + my) ≤ pow2 159 := by
      simp only [pow2]
      exact pow_le_pow_right₀ (by norm_num) hsum
    calc
      |x| * |y| < pow2 mx * pow2 my := by
        gcongr
      _ = pow2 (mx + my) := by simp [pow2, pow_add]
      _ ≤ pow2 159 := hp
  rw [abs_lt] at habs
  simpa [InSigned] using And.intro (le_of_lt habs.1) habs.2

/-- `rx_mul` adds `2^63` and then arithmetically shifts by 64.  Therefore the
integer error is at most half an ulp and exact ties go toward `+∞`. -/
theorem roundShift64_error (x : ℤ) :
    |pow2 64 * roundShift64 x - x| ≤ pow2 63 := by
  norm_num [roundShift64, roundShift, pow2, abs_le]
  omega

/-- Universal local precision theorem for `rx_mul`: independently of operand
magnitudes, the decoded rounded product differs from exact real multiplication
by at most `2^-65`.  Range hypotheses are needed separately to rule out
storage wraparound. -/
theorem qmul_real_error (x y : ℤ) :
    |decodeQ32_64 (roundShift64 (x * y)) -
      decodeQ32_64 x * decodeQ32_64 y| ≤ (2 : ℝ) ^ (-65 : ℤ) := by
  have h := roundShift64_error (x * y)
  have hR : |((pow2 64 * roundShift64 (x * y) - x * y : ℤ) : ℝ)| ≤
      ((pow2 63 : ℤ) : ℝ) := by exact_mod_cast h
  norm_num [pow2] at hR
  rw [show decodeQ32_64 (roundShift64 (x * y)) -
      decodeQ32_64 x * decodeQ32_64 y =
      (((roundShift64 (x * y) : ℤ) : ℝ) * (2 : ℝ) ^ 64 -
        ((x * y : ℤ) : ℝ)) / ((2 : ℝ) ^ 64) ^ 2 by
        simp [decodeQ32_64]
        ring]
  rw [abs_div, abs_of_pos (by positivity : 0 < ((2 : ℝ) ^ 64) ^ 2)]
  norm_num
  apply (div_le_iff₀ (by norm_num : (0 : ℝ) <
    340282366920938463463374607431768211456)).2
  norm_num [mul_comm] at hR ⊢
  exact hR

theorem roundShift64_tie_pos (k : ℤ) :
    roundShift64 (pow2 64 * k + pow2 63) = k + 1 := by
  norm_num [roundShift64, roundShift, pow2]
  rw [show 18446744073709551616 * k + 9223372036854775808 +
      9223372036854775808 = 18446744073709551616 * (k + 1) by ring]
  simp

theorem roundShift64_tie_neg (k : ℤ) :
    roundShift64 (pow2 64 * k - pow2 63) = k := by
  norm_num [roundShift64, roundShift, pow2]

/-- If the exact real product stays in Q32.64 range with a half-ulp margin,
the rounded result is again a valid signed 96-bit representative. -/
theorem rounded_product_fits_q32_64 {p : ℤ}
    (hlo : -(pow2 159) + pow2 63 ≤ p)
    (hhi : p < pow2 159 - pow2 63) :
    InQ32_64 (roundShift64 p) := by
  have herr := roundShift64_error p
  rw [abs_le] at herr
  constructor <;> norm_num [InQ32_64, InSigned, pow2, roundShift64, roundShift] at hlo hhi herr ⊢ <;> omega

end FalconRange
