module
public import Mathlib
public import FalconRange.FixedPoint
public import FalconRange.MaskedRing
public import FalconRange.DynamicBackend
public import FalconRange.FFSamplingFastPath
@[expose] public section

/-!
# Sharewise public-constant Q32.64 multiplication

The C gadget multiplies both arithmetic shares by the same public integer,
then performs one masked add-half/right-shift conversion.  These theorems
cover the linear ring step, the 160-bit period removed by the public shift,
and every fast public-constant range class used by the signer.
-/

namespace FalconRange

def Shares160.mulPublic (a : Shares160) (c : R160) : Shares160 where
  s0 := a.s0 * c
  s1 := a.s1 * c

/-- Sharewise multiplication by a public constant reconstructs to the ring
product; no secret-secret cross-products are required. -/
theorem reconstruct_mulPublic (a : Shares160) (c : R160) :
    (a.mulPublic c).reconstruct = a.reconstruct * c := by
  simp [Shares160.mulPublic, Shares160.reconstruct]
  ring

/-- Replacing an exact product by its representative modulo 2^160 changes the
rounded, shifted integer only by a multiple of 2^96.  Thus the low signed
96-bit Q32.64 payload is identical whenever the exact product satisfies the
fast-path range theorem. -/
theorem roundShift64_add_pow2_160 (p k : ℤ) :
    roundShift64 (p + k * pow2 160) =
      roundShift64 p + k * pow2 96 := by
  rw [roundShift64, roundShift]
  have hpow : pow2 160 = pow2 96 * pow2 64 := by
    norm_num [pow2]
  rw [show p + k * pow2 160 + pow2 63 =
      (p + pow2 63) + (k * pow2 96) * pow2 64 by
        rw [hpow]
        ring]
  rw [Int.add_mul_ediv_right _ _ (ne_of_gt (pow2_pos 64))]
  rfl

/-- The wide public-constant path is sufficient for arbitrary Q32.64 input
and public Q32.64 constants. -/
theorem public_wide_product_fits_s192 {x c : ℤ}
    (hx : InQ32_64 x) (hc : InQ32_64 c) :
    InSigned 192 (x * c) :=
  q32_64_product_fits_s192 hx hc

/-- The normalized inverse/division seed multiplies |ub|<2^64 by the exact
32/17 table constant whose raw magnitude is below 2^66. -/
theorem inverse_seed_public_fits_s160 {ub c : ℤ}
    (hub : |ub| < pow2 64) (hc : |c| < pow2 66) :
    InSigned 160 (ub * c) := by
  apply product_fits_s160_of_magnitude hub hc
  norm_num

/-- BerExp reduction products are bounded by raw exponents 72 and 65. -/
theorem ber_exp_public_fits_s160 {x c : ℤ}
    (hx : |x| < pow2 72) (hc : |c| < pow2 65) :
    InSigned 160 (x * c) := by
  apply product_fits_s160_of_magnitude hx hc
  norm_num

/-- isigma is Q0.64 and both standard sigma-min constants are below two. -/
theorem sampler_sigma_min_public_fits_s160 {x c : ℤ}
    (hx : |x| < pow2 64) (hc : |c| < pow2 66) :
    InSigned 160 (x * c) := by
  apply product_fits_s160_of_magnitude hx hc
  norm_num

/-- BaseSampler's z0^2 term is below 2^74 raw and the public inverse scale is
Q0.64, leaving a large signed-160 margin. -/
theorem sampler_inv_2sqr_public_fits_s160 {x c : ℤ}
    (hx : |x| < pow2 74) (hc : |c| < pow2 64) :
    InSigned 160 (x * c) := by
  apply product_fits_s160_of_magnitude hx hc
  norm_num

/-- Exact raw encodings used in mask160.c/mask_sampler.c. -/
def raw96 (lo mid hi : ℤ) : ℤ :=
  lo + mid * pow2 32 + hi * pow2 64

def c32_17_raw : ℤ := raw96 3789677026 3789677025 1
def invLog2_raw : ℤ := raw96 3090145655 1901360722 1
def log2_raw : ℤ := raw96 3520035244 2977044471 0
def sigmaMin9_raw : ℤ := raw96 4154935963 1193286641 1
def sigmaMin10_raw : ℤ := raw96 209075216 1281104281 1
def inv2SqrSigma0_raw : ℤ := raw96 125273110 647960451 0

theorem c32_17_raw_bound : |c32_17_raw| < pow2 66 := by
  norm_num [c32_17_raw, raw96, pow2]

theorem invLog2_raw_bound : |invLog2_raw| < pow2 65 := by
  norm_num [invLog2_raw, raw96, pow2]

theorem log2_raw_bound : |log2_raw| < pow2 64 := by
  norm_num [log2_raw, raw96, pow2]

theorem sigmaMin9_raw_bound : |sigmaMin9_raw| < pow2 66 := by
  norm_num [sigmaMin9_raw, raw96, pow2]

theorem sigmaMin10_raw_bound : |sigmaMin10_raw| < pow2 66 := by
  norm_num [sigmaMin10_raw, raw96, pow2]

theorem inv2SqrSigma0_raw_bound : |inv2SqrSigma0_raw| < pow2 64 := by
  norm_num [inv2SqrSigma0_raw, raw96, pow2]

/-- C call-site policy for public-constant products.  Disabled sites always
use the 192-bit public fallback. -/
inductive PublicMulSite
  | splitTwiddle
  | treeLdlTheta
  | treeIxadjTheta
  | treeIxadjMu
  | invSigmaOrInvQ
  | inverseOrDivisionSeed
  | berExpReduction
  | samplerSigmaMin
  | samplerInv2SqrSigma0
  | fftSqrt2OrRNJ
  | generalLdlOrIxadj
  deriving DecidableEq

def PublicMulSite.fastEnabled : PublicMulSite -> Bool
  | .splitTwiddle | .treeLdlTheta | .treeIxadjTheta | .treeIxadjMu
  | .invSigmaOrInvQ | .inverseOrDivisionSeed | .berExpReduction
  | .samplerSigmaMin | .samplerInv2SqrSigma0 => true
  | .fftSqrt2OrRNJ | .generalLdlOrIxadj => false

def PublicMulSite.rawProductExponent : PublicMulSite -> Nat
  | .splitTwiddle | .invSigmaOrInvQ => 159
  | .treeLdlTheta => 155
  | .treeIxadjTheta => 131
  | .treeIxadjMu => 153
  | .inverseOrDivisionSeed => 130
  | .berExpReduction => 137
  | .samplerSigmaMin => 130
  | .samplerInv2SqrSigma0 => 138
  | .fftSqrt2OrRNJ | .generalLdlOrIxadj => 0

theorem enabled_public_site_raw_product_exponent_lt_160
    (site : PublicMulSite) (h : site.fastEnabled = true) :
    site.rawProductExponent < 160 := by
  cases site <;> simp [PublicMulSite.fastEnabled,
    PublicMulSite.rawProductExponent] at h ⊢

end FalconRange
