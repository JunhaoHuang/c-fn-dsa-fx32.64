module
public import Mathlib
@[expose] public section

/-!
# Two-share arithmetic in the 160-bit ring

This is the algebraic model of `mrx160`.  It proves functional correctness
of refresh and the first-order multiplication layout independently of the
randomness generator and side-channel assumptions.
-/

namespace FalconRange

abbrev R160 := ZMod (2 ^ 160)

structure Shares160 where
  s0 : R160
  s1 : R160

def Shares160.reconstruct (x : Shares160) : R160 := x.s0 + x.s1

def Shares160.refresh (x : Shares160) (r : R160) : Shares160 where
  s0 := x.s0 + r
  s1 := x.s1 - r

/-- Two-share ISW-style multiplication.  The cross-products are assigned to
the second output share and `r` refreshes the result. -/
def Shares160.mul (a b : Shares160) (r : R160) : Shares160 where
  s0 := a.s0 * b.s0 + r
  s1 := a.s0 * b.s1 + a.s1 * b.s0 + a.s1 * b.s1 - r

theorem reconstruct_refresh (x : Shares160) (r : R160) :
    (x.refresh r).reconstruct = x.reconstruct := by
  simp [Shares160.refresh, Shares160.reconstruct]

theorem reconstruct_mul (a b : Shares160) (r : R160) :
    (a.mul b r).reconstruct = a.reconstruct * b.reconstruct := by
  simp [Shares160.mul, Shares160.reconstruct]
  ring

/-- Refresh does not change a product's reconstructed value. -/
theorem reconstruct_mul_refresh (a b : Shares160) (rMul rRefresh : R160) :
    ((a.mul b rMul).refresh rRefresh).reconstruct =
      a.reconstruct * b.reconstruct := by
  rw [reconstruct_refresh, reconstruct_mul]

/-- Add and subtract are linear and therefore may run sharewise. -/
theorem reconstruct_add (a b : Shares160) :
    ({ s0 := a.s0 + b.s0, s1 := a.s1 + b.s1 } : Shares160).reconstruct =
      a.reconstruct + b.reconstruct := by
  change (a.s0 + b.s0) + (a.s1 + b.s1) = (a.s0 + a.s1) + (b.s0 + b.s1)
  ring


end FalconRange
