module
public import Mathlib
public import FalconRange.FixedPoint
@[expose] public section

/-!
# Exact range obligations for masked Falcon NTT and CRT

These theorems justify the integer widths used by `mask_mq.c`. They contain
no sampled maximum, proof placeholder, or custom axiom. Security simulation arguments
remain in `security/NTT_CRT_PROOFS.md`; this module proves the arithmetic and
range obligations on which those arguments rely.
-/

namespace FalconRange

def nttQ : ℤ := 12289
def nttQ2 : ℤ := 18433
def nttM : ℤ := 226523137
def nttH : ℤ := 113261568

/-- Falcon's CRT modulus is exactly the product of its two NTT primes. -/
theorem nttM_eq : nttM = nttQ * nttQ2 := by
  norm_num [nttM, nttQ, nttQ2]

/-- A product of canonical residues is a valid input to 32-bit Montgomery
REDC: it is below q*R and therefore below 2^64. -/
theorem montgomery_product_bound {q a b : ℤ}
    (hq0 : 0 < q) (hqR : q < pow2 32)
    (ha : 0 ≤ a ∧ a < q) (hb : 0 ≤ b ∧ b < q) :
    0 ≤ a * b ∧ a * b < q * pow2 32 := by
  constructor
  · exact mul_nonneg ha.1 hb.1
  · calc
      a * b < q * q := by
        nlinarith [mul_nonneg ha.1 hb.1]
      _ < q * pow2 32 := by
        exact mul_lt_mul_of_pos_left hqR hq0

/-- Both Falcon primes satisfy the public Montgomery precondition. -/
theorem falcon_primes_lt_R :
    0 < nttQ ∧ nttQ < pow2 32 ∧ 0 < nttQ2 ∧ nttQ2 < pow2 32 := by
  norm_num [nttQ, nttQ2, pow2]

/-- A butterfly sum of two canonical q' residues fits in uint32_t; q' is the
larger prime, so this also covers q. -/
theorem butterfly_sum_fits_u32 {a b : ℤ}
    (ha : 0 ≤ a ∧ a < nttQ2) (hb : 0 ≤ b ∧ b < nttQ2) :
    0 ≤ a + b ∧ a + b < pow2 32 := by
  norm_num [nttQ2, pow2] at ha hb ⊢
  omega

/-- Per-share CRT reconstruction is canonical modulo M. -/
theorem crt_component_range {a u : ℤ}
    (ha : 0 ≤ a ∧ a < nttQ) (hu : 0 ≤ u ∧ u < nttQ2) :
    0 ≤ a + nttQ * u ∧ a + nttQ * u < nttM := by
  norm_num [nttQ, nttQ2, nttM] at ha hu ⊢
  omega

/-- Adding two canonical CRT shares cannot wrap a 32-bit arithmetic sharing. -/
theorem crt_share_sum_range {c0 c1 : ℤ}
    (h0 : 0 ≤ c0 ∧ c0 < nttM)
    (h1 : 0 ≤ c1 ∧ c1 < nttM) :
    0 ≤ c0 + c1 ∧ c0 + c1 < 2 * nttM := by
  constructor <;> omega

/-- The C implementation's stronger bound: the unreduced share sum is below
2^29, so signed-bit tests after subtracting M are unambiguous. -/
theorem twice_crt_modulus_lt_pow2_29 : 2 * nttM < pow2 29 := by
  norm_num [nttM, pow2]

def reduceM (s : ℤ) : ℤ :=
  if nttM ≤ s then s - nttM else s

def centerM (s : ℤ) : ℤ :=
  let x := reduceM s
  if nttH < x then x - nttM else x

/-- One conditional subtraction maps the proven input interval into [0,M). -/
theorem reduceM_range {s : ℤ}
    (hs : 0 ≤ s ∧ s < 2 * nttM) :
    0 ≤ reduceM s ∧ reduceM s < nttM := by
  unfold reduceM
  split_ifs with h
  · norm_num [nttM] at hs h ⊢
    omega
  · norm_num [nttM] at hs h ⊢
    omega

/-- The second masked comparison/subtraction produces exactly the signed
Falcon CRT interval [-H,H]. -/
theorem centerM_range {s : ℤ}
    (hs : 0 ≤ s ∧ s < 2 * nttM) :
    -nttH ≤ centerM s ∧ centerM s ≤ nttH := by
  have hx := reduceM_range hs
  unfold centerM
  dsimp
  split_ifs with h
  · norm_num [nttM, nttH] at hx h ⊢
    omega
  · norm_num [nttM, nttH] at hx h ⊢
    omega

/-- Constants used by the C sign-bit comparisons are below 2^31. -/
theorem center_subtractions_fit_s32 :
    nttM < pow2 31 ∧ nttH + 1 < pow2 31 := by
  norm_num [nttM, nttH, pow2]

end FalconRange
