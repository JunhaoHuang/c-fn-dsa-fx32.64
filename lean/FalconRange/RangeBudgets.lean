module
public import Mathlib
public import FalconRange.FixedPoint
@[expose] public section

/-!
# Falcon-512 magnitude budgets

The numerical statements below are exact integer/rational consequences of
the four key-generation filters used by the fixed-point Falcon paper:
`γ_fg = 255`, `α_h = 4`, `γ_FG = 3500`, and `γ_root = 24`.

No sampled key and no floating-point oracle appears in these proofs.
-/

namespace FalconRange

def q : ℤ := 12289
def n : ℤ := 512
def gammaFGSmall : ℤ := 255
def alphaHybrid : ℤ := 4
def gammaFGLarge : ℤ := 3500
def gammaRoot : ℤ := 24

theorem gram00_lt_2pow17 : 2 * gammaFGSmall ^ 2 < pow2 17 := by
  norm_num [gammaFGSmall, pow2]

theorem gram01_lt_2pow21 :
    2 * gammaFGSmall * gammaFGLarge < pow2 21 := by
  norm_num [gammaFGSmall, gammaFGLarge, pow2]

theorem root_l10_lt_2pow5 : gammaRoot < pow2 5 := by
  norm_num [gammaRoot, pow2]

theorem inner_diagonal_lt_2pow18 : alphaHybrid ^ 2 * q < pow2 18 := by
  norm_num [alphaHybrid, q, pow2]

theorem leaf_diagonal_lt_2pow15 :
    (117 : ℤ) ^ 2 * q < 100 ^ 2 * pow2 15 := by
  norm_num [q, pow2]

theorem centered_mod_q_lt_2pow13 : (q - 1) / 2 < pow2 13 := by
  norm_num [q, pow2]

theorem hash_point_lt_2pow14 : q - 1 < pow2 14 := by
  norm_num [q, pow2]

/-- A rational envelope for the Lemma-13 drift.  The premises are analytic
facts supplied by that lemma (`τ < 15`, `σ < 166`, `√n < 23`, `√2 < 3/2`,
and `√q > 110`).  Everything after those premises is exact arithmetic. -/
theorem drift_lt_222000
    {tau sigma sqrtN sqrtTwo sqrtQ : ℝ}
    (htau : tau < 15) (hsigma : sigma < 166)
    (hsqrtN : sqrtN < 23) (hsqrtTwo : sqrtTwo < 3 / 2)
    (hsqrtQ : 110 < sqrtQ)
    (hpos : 0 ≤ tau ∧ 0 ≤ sigma ∧ 0 ≤ sqrtN ∧ 0 ≤ sqrtTwo) :
    2 * ((2 * tau * sigma * 4) / sqrtQ) * sqrtN *
      (24 + sqrtTwo + 1) < 222000 := by
  rcases hpos with ⟨ht0, hs0, hn0, h20⟩
  have hq0 : 0 < sqrtQ := lt_trans (by norm_num) hsqrtQ
  have hA : (2 * tau * sigma * 4) / sqrtQ < 182 := by
    apply (div_lt_iff₀ hq0).2
    have hnum : 2 * tau * sigma * 4 < 2 * 15 * 166 * 4 := by
      nlinarith
    nlinarith
  have hA0 : 0 ≤ (2 * tau * sigma * 4) / sqrtQ := by positivity
  have hAN : ((2 * tau * sigma * 4) / sqrtQ) * sqrtN < 182 * 23 := by
    calc
      ((2 * tau * sigma * 4) / sqrtQ) * sqrtN ≤ 182 * sqrtN :=
        mul_le_mul_of_nonneg_right hA.le hn0
      _ < 182 * 23 := mul_lt_mul_of_pos_left hsqrtN (by norm_num)
  have hAN0 : 0 ≤ ((2 * tau * sigma * 4) / sqrtQ) * sqrtN := mul_nonneg hA0 hn0
  have hC : 24 + sqrtTwo + 1 < (53 : ℝ) / 2 := by linarith
  have hC0 : 0 ≤ 24 + sqrtTwo + 1 := by nlinarith
  have hANC : ((2 * tau * sigma * 4) / sqrtQ) * sqrtN *
      (24 + sqrtTwo + 1) < (182 * 23) * ((53 : ℝ) / 2) := by
    calc
      ((2 * tau * sigma * 4) / sqrtQ) * sqrtN * (24 + sqrtTwo + 1) ≤
          (182 * 23) * (24 + sqrtTwo + 1) := mul_le_mul_of_nonneg_right hAN.le hC0
      _ < (182 * 23) * ((53 : ℝ) / 2) := mul_lt_mul_of_pos_left hC (by norm_num)
  calc
    2 * ((2 * tau * sigma * 4) / sqrtQ) * sqrtN * (24 + sqrtTwo + 1) =
        2 * ((((2 * tau * sigma * 4) / sqrtQ) * sqrtN) * (24 + sqrtTwo + 1)) := by ring
    _ < 2 * ((182 * 23) * ((53 : ℝ) / 2)) := mul_lt_mul_of_pos_left hANC (by norm_num)
    _ < 222000 := by norm_num

theorem ffsampling_tweak_lt_2pow18
    {drift : ℝ} (_h0 : 0 ≤ drift) (hd : drift < 222000) :
    (512 : ℝ) / 2 + drift < (2 : ℝ) ^ 18 := by
  norm_num
  linarith

theorem ffsampling_standard_lt_2pow21
    {drift : ℝ} (_h0 : 0 ≤ drift) (hd : drift < 222000) :
    (512 : ℝ) * 3500 * (12288 / 12289) + drift < (2 : ℝ) ^ 21 := by
  norm_num
  nlinarith

/-- Raw-Q32.64 fast-path criterion.  Real magnitude tags `m₁,m₂` become raw
integer exponents `m₁+64,m₂+64`; therefore `m₁+m₂ ≤ 31` suffices for the
unshifted product to fit in signed 160 bits. -/
theorem q32_64_fast_path {x y : ℤ} {m1 m2 : ℕ}
    (hx : |x| < pow2 (m1 + 64)) (hy : |y| < pow2 (m2 + 64))
    (hsum : m1 + m2 ≤ 31) : InSigned 160 (x * y) := by
  apply product_fits_s160_of_magnitude hx hy
  omega

end FalconRange
