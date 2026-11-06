module
public import Mathlib
public import FalconRange.FixedPoint
public import FalconRange.MaskedNTTCRT
@[expose] public section

/-!
# Range obligations for the integrated masked signer

These theorems cover the arithmetic connection layer added by
`sign_masked.c`: signed `mi32` modulus switching, q/q-prime centering,
64-bit norm/encoding accumulators, and the caller-owned arena. They contain
no sampled premise, `sorry`, or custom axiom. They are functional range
proofs, not probing-security proofs.
-/

namespace FalconRange

/-- One limb of an arithmetic sharing is below 2^32. -/
def InU32 (x : ℤ) : Prop := 0 ≤ x ∧ x < pow2 32

/-- Two canonical uint32 shares have a carry bit, never a larger carry. -/
theorem two_u32_sum_range {a0 a1 : ℤ}
    (h0 : InU32 a0) (h1 : InU32 a1) :
    0 ≤ a0 + a1 ∧ a0 + a1 < 2 * pow2 32 := by
  unfold InU32 at h0 h1
  constructor <;> omega

/-- Algebraic justification for the C modulus-switch correction. If
`low + carry*2^32` is the unsigned share sum and `sign` is the low-word sign
bit, subtracting `(carry+sign)*(2^32 mod p)` modulo p represents the signed
32-bit value `low-sign*2^32`. -/
theorem signed_share_mod_correction {p a0 a1 low carry sign : ℤ}
    (_hp : 0 < p)
    (hdecomp : a0 + a1 = low + carry * pow2 32) :
    Int.ModEq p (low - sign * pow2 32)
      (a0 + a1 - (carry + sign) * pow2 32) := by
  rw [hdecomp]
  ring_nf
  exact Int.ModEq.rfl

/-- The signed representative selected from a 32-bit low word is in the
ordinary int32 interval. -/
theorem signed_low_word_range {low sign : ℤ}
    (hlow : 0 ≤ low ∧ low < pow2 32)
    (hsign0 : sign = 0 → low < pow2 31)
    (hsign1 : sign = 1 → pow2 31 ≤ low)
    (hsign : sign = 0 ∨ sign = 1) :
    -(pow2 31) ≤ low - sign * pow2 32 ∧
      low - sign * pow2 32 < pow2 31 := by
  rcases hsign with rfl | rfl
  · norm_num [pow2] at hsign0 ⊢
    omega
  · norm_num [pow2] at hlow hsign1 ⊢
    omega

def centerQ (x : ℤ) : ℤ :=
  if 6144 < x then x - nttQ else x

def centerQ2 (x : ℤ) : ℤ :=
  if 9216 < x then x - nttQ2 else x

/-- Joint centering modulo q produces Falcon's exact symmetric interval. -/
theorem centerQ_range {x : ℤ} (hx : 0 ≤ x ∧ x < nttQ) :
    -6144 ≤ centerQ x ∧ centerQ x ≤ 6144 := by
  unfold centerQ
  split_ifs <;> norm_num [nttQ] at hx ⊢ <;> omega

/-- Joint centering modulo q-prime produces its exact symmetric interval. -/
theorem centerQ2_range {x : ℤ} (hx : 0 ≤ x ∧ x < nttQ2) :
    -9216 ≤ centerQ2 x ∧ centerQ2 x ≤ 9216 := by
  unfold centerQ2
  split_ifs <;> norm_num [nttQ2] at hx ⊢ <;> omega

/-- Squaring one accepted coefficient cannot overflow a 32-bit unsigned
integer; the implementation nevertheless accumulates in a masked 160-bit
value and compares against a public 32-bit bound. -/
theorem accepted_coefficient_square_bound {x : ℤ} (hx : |x| ≤ 840) :
    0 ≤ x * x ∧ x * x ≤ 705600 := by
  have h := abs_le.mp hx
  constructor <;> nlinarith [sq_nonneg x]
/-- Even the conservative sum of 2048 maximal coefficient squares fits in
an unsigned 64-bit accumulator. -/
theorem signature_sqnorm_fits_u64 :
    2 * 1024 * 705600 < pow2 64 := by
  norm_num [pow2]

/-- `comp_encode` uses at most 15 bits for an accepted coefficient. -/
theorem compressed_coefficient_bits_bound {x : ℤ}
    (hx0 : 0 ≤ x) (hx : x ≤ 840) :
    9 + x / 128 ≤ 15 := by
  omega

/-- The conservative 1024-coefficient compressed-length sum fits in the
same 64-bit Boolean accumulator used by the C gadget. -/
theorem compressed_length_fits_u64 :
    1024 * 15 < pow2 64 := by
  norm_num [pow2]

/-- Exact byte formula used by `fndsa_masked_tmp_size()`. -/
def maskedArenaSize (n sigSize : ℕ) : ℕ :=
  216 * n + 63 + sigSize

/-- Exact standard-parameter arena sizes, including alignment allowance and
staged signature storage. -/
theorem masked_arena_sizes :
    maskedArenaSize 512 666 = 111321 ∧
    maskedArenaSize 1024 1280 = 222527 := by
  norm_num [maskedArenaSize]

/-- The largest arena is below 218 KiB and far below `size_t` on a 32-bit
platform, so every intermediate size calculation is non-wrapping. -/
theorem masked_arena_1024_size_bound :
    maskedArenaSize 1024 1280 ≤ 218 * 1024 ∧
    maskedArenaSize 1024 1280 < pow2 32 := by
  norm_num [maskedArenaSize, pow2]

/-- The v2 key formula is likewise non-wrapping for the largest standard
parameter (verification-key size is 1793 bytes at logn=10). -/
theorem masked_key_v2_size_bound :
    72 + 1793 + 16 * 1024 = 18249 ∧
    72 + 1793 + 16 * 1024 < pow2 32 := by
  norm_num [pow2]

end FalconRange
