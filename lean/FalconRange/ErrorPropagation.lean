module
public import Mathlib
@[expose] public section

/-!
# Deterministic absolute-error propagation

These lemmas replace "observed RMSE plus a margin" with universal bounds.
They are intentionally independent of a particular Falcon call graph: the
algorithm-specific proof must instantiate `A`, `B`, `ea`, `eb`, and the local
rounding error at every node.
-/

namespace FalconRange

theorem add_abs_error {a ahat b bhat ea eb : ℝ}
    (ha : |ahat - a| ≤ ea) (hb : |bhat - b| ≤ eb) :
    |(ahat + bhat) - (a + b)| ≤ ea + eb := by
  calc
    |(ahat + bhat) - (a + b)| = |(ahat - a) + (bhat - b)| := by ring_nf
    _ ≤ |ahat - a| + |bhat - b| := abs_add_le _ _
    _ ≤ ea + eb := add_le_add ha hb

theorem sub_abs_error {a ahat b bhat ea eb : ℝ}
    (ha : |ahat - a| ≤ ea) (hb : |bhat - b| ≤ eb) :
    |(ahat - bhat) - (a - b)| ≤ ea + eb := by
  calc
    |(ahat - bhat) - (a - b)| = |(ahat - a) - (bhat - b)| := by ring_nf
    _ ≤ |ahat - a| + |bhat - b| := abs_sub _ _
    _ ≤ ea + eb := add_le_add ha hb

/-- Absolute error for a rounded multiplication.  `A` and `B` are universal
magnitude bounds on the exact operands; `ρ` is the local rounding error. -/
theorem mul_abs_error {a ahat b bhat ea eb A B ρ r : ℝ}
    (ha : |ahat - a| ≤ ea) (hb : |bhat - b| ≤ eb)
    (haMag : |a| ≤ A) (hbMag : |b| ≤ B)
    (hr : |r - ahat * bhat| ≤ ρ)
    (hea : 0 ≤ ea) (_heb : 0 ≤ eb) :
    |r - a * b| ≤ ρ + A * eb + B * ea + ea * eb := by
  have hahat : |ahat| ≤ A + ea := by
    calc
      |ahat| = |a + (ahat - a)| := by ring_nf
      _ ≤ |a| + |ahat - a| := abs_add_le _ _
      _ ≤ A + ea := add_le_add haMag ha
  have hA0 : 0 ≤ A := (abs_nonneg a).trans haMag
  have hB0 : 0 ≤ B := (abs_nonneg b).trans hbMag
  have hdecomp : ahat * bhat - a * b =
      ahat * (bhat - b) + b * (ahat - a) := by ring
  calc
    |r - a * b| = |(r - ahat * bhat) + (ahat * bhat - a * b)| := by ring_nf
    _ ≤ |r - ahat * bhat| + |ahat * bhat - a * b| := abs_add_le _ _
    _ ≤ ρ + |ahat * (bhat - b) + b * (ahat - a)| := by
      rw [← hdecomp]
      exact add_le_add hr (le_refl _)
    _ ≤ ρ + (|ahat| * |bhat - b| + |b| * |ahat - a|) := by
      have htriangle :
          |ahat * (bhat - b) + b * (ahat - a)| ≤
            |ahat| * |bhat - b| + |b| * |ahat - a| := by
        calc
          |ahat * (bhat - b) + b * (ahat - a)| ≤
              |ahat * (bhat - b)| + |b * (ahat - a)| := abs_add_le _ _
          _ = |ahat| * |bhat - b| + |b| * |ahat - a| := by rw [abs_mul, abs_mul]
      simpa [add_comm] using add_le_add_left htriangle ρ
    _ ≤ ρ + ((A + ea) * eb + B * ea) := by
      have hproducts := add_le_add
        (mul_le_mul hahat hb (abs_nonneg _) (add_nonneg hA0 hea))
        (mul_le_mul hbMag ha (abs_nonneg _) hB0)
      simpa [add_comm] using add_le_add_left hproducts ρ
    _ = ρ + A * eb + B * ea + ea * eb := by ring

/-- One radix-2 butterfly with exact unit-magnitude twiddle. -/
theorem butterfly_error {a ahat b bhat w what out ea eb ew ρ B : ℝ}
    (ha : |ahat - a| ≤ ea) (hb : |bhat - b| ≤ eb)
    (hw : |what - w| ≤ ew) (hwMag : |w| ≤ 1)
    (hbMag : |b| ≤ B) (hr : |out - (ahat + what * bhat)| ≤ ρ)
    (_hea : 0 ≤ ea) (heb : 0 ≤ eb) (hew : 0 ≤ ew) :
    |out - (a + w * b)| ≤ ρ + ea + eb + B * ew + ew * eb := by
  have hmul := mul_abs_error hw hb hwMag hbMag (r := what * bhat) (ρ := 0)
    (by simp) hew heb
  have hadd := add_abs_error ha hmul
  calc
    |out - (a + w * b)| ≤ |out - (ahat + what * bhat)| +
        |(ahat + what * bhat) - (a + w * b)| := by
      calc
        |out - (a + w * b)| = |(out - (ahat + what * bhat)) +
            ((ahat + what * bhat) - (a + w * b))| := by ring_nf
        _ ≤ _ := abs_add_le _ _
    _ ≤ ρ + (ea + (0 + 1 * eb + B * ew + ew * eb)) := add_le_add hr hadd
    _ = ρ + ea + eb + B * ew + ew * eb := by ring

/-- Linear recurrence used to close a depth-`d` FFT error proof. -/
def errorRec (α β : ℝ) : ℕ → ℝ
  | 0 => 0
  | n + 1 => α * errorRec α β n + β

theorem errorRec_nonneg {α β : ℝ} (hα : 0 ≤ α) (hβ : 0 ≤ β) :
    ∀ n, 0 ≤ errorRec α β n := by
  intro n
  induction n with
  | zero => simp [errorRec]
  | succ n ih => simp [errorRec]; positivity

theorem errorRec_two (β : ℝ) (n : ℕ) :
    errorRec 2 β n = ((2 : ℝ) ^ n - 1) * β := by
  induction n with
  | zero => simp [errorRec]
  | succ n ih => simp [errorRec, ih, pow_succ]; ring

end FalconRange
