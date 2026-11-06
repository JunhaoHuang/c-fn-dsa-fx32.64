# Lean proofs for the fixed-point Falcon range analysis

This directory turns the reusable, deterministic part of the fixed-point
analysis into kernel-checked Lean theorems. It deliberately does **not** turn
a finite-precision experiment into a universal claim.

## What is proved

- every signed Q32.64 representative embeds exactly in signed 160 bits;
- the exact product of two signed Q32.64 representatives fits in signed
  192 bits, justifying the wide fallback;
- a generic magnitude criterion for the 160-bit fast multiplication path;
- the exact `rx_mul` rounding rule (add `2^63`, arithmetic shift by 64), its
  half-ulp error, and its tie-toward-`+∞` behavior;
- a decoded Q32.64 multiplication error bound of `2^-65`;
- reconstruction correctness for 160-bit two-share refresh and multiplication;
- exact Falcon-512 magnitude budgets implied by the four filtered-key bounds;
- deterministic add, subtract, multiply, FFT-butterfly, recurrence-error, and
  error-propagation lemmas;
- arithmetic halving, Q0.64 public-twiddle products, division/square-root local
  product budgets, and 512/1024 public FFT-depth budgets;
- Montgomery product/REDC input bounds, q/q' butterfly uint32 bounds, per-share
  CRT bounds, the `2*M < 2^29` invariant, and the exact centered interval;
- signed `mi32` modulus-switch correction, exact q/q-prime centered ranges,
  accepted-coefficient square and compressed-bit bounds, and exact v2 key/
  512/1024 arena sizes;
- the standard 512/1024 Gram-root bound, positive LDL mean/harmonic children,
  `|L10| < 1`, SamplerZ residual and depth-ten rounding margins;
- sharewise public-constant reconstruction, 160-bit rounding-period equivalence,
  exact table bounds and public fast/wide call-site policy;
- an executable per-call policy enabling recursive LDL/iXadj fast products
  while retaining the 192-bit fallback for unproved root products.

The C-to-Lean review map is in [CALL_SITES.md](CALL_SITES.md). Recursive
ffSampling theorems are in `FalconRange/FFSamplingFastPath.lean`, executable
range policy in `FalconRange/DynamicBackend.lean`, and public-constant
obligations in `FalconRange/PublicConstant.lean`. Modular transform/CRT and
real-signer connection obligations are in `MaskedNTTCRT.lean` and
`IntegratedSigning.lean`.

## What remains an explicit proof obligation

The paper's empirical RMSE values for the full
`Gram → ffLDL → σ_i → ffSampling` pipeline are measurements over finitely many
keys. To replace them by a theorem, each call site must provide:

1. an exact interval for every operand and divisor;
2. the local rounding error of the concrete fixed-point operation;
3. a recursive invariant connecting split/merge and ffLDL child bounds;
4. a SamplerZ stability or distributional argument near integer and Bernoulli
   decision boundaries.

`ErrorPropagation.lean` supplies universal composition rules, while
`FFSamplingFastPath.lean` closes the width obligations for recursive LDL and
iXadj without assuming observed ceilings. This is not a proof that finite
precision leaves every SamplerZ/BerExp decision or output distribution
unchanged for all real-valued executions. Root complex and root-division
products also remain on the exact 192-bit fallback.

## Build

```sh
lake update
lake build
```

The project pins Lean and Mathlib to `v4.35.0`. The 2026-11-05 build completed
successfully with 9057 jobs, including `FalconRange.DynamicBackend`,
`FalconRange.FFSamplingFastPath`, `FalconRange.PublicConstant`, `FalconRange.MaskedNTTCRT`, and
`FalconRange.IntegratedSigning` in the default target. Project sources,
excluding downloaded dependencies under `.lake`, contain no `sorry` or custom
`axiom` declarations.
