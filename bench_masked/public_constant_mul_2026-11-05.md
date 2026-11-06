# Q32.64 public-constant masked multiplication benchmark

Date: 2026-11-05
Host: AMD Ryzen 9 9950X3D 16-Core Processor
Compiler: clang, `-O2`; masking modules use `-fno-lto`
Counter: invariant TSC/reference ticks, **not PMU core cycles**
Command: `./speed_masked -tsc`

The public paths retain the old mask-RNG endpoint: public-fast consumes 285
core words and discards 5; public-wide consumes 621 and discards 288. The
benchmark therefore includes the compatibility cost rather than reporting an
idealized no-padding number.

## Direct gadget comparison

| Operation | Generic | Public constant | Reduction |
| --- | ---: | ---: | ---: |
| signed-160 Q32.64 multiply | 1,965.96 | 1,738.06 | 11.59% |
| signed-192 Q32.64 multiply | 5,723.30 | 4,533.63 | 20.79% |

The dedicated symbols do not call `ring_mul()` or `ring_mul192()`. Results,
half-tie behavior, in-place operation and final RNG structs are checked against
the corresponding generic gadgets by `test_masked`.

## Integrated measurements

The comparison baseline is the preceding 160-bit ffSampling report
`ffsampling_fast160_2026-11-05.md`. Since these are separate benchmark runs,
the following whole-function percentages include ordinary host variation; they
are supporting performance evidence, not a hard correctness criterion.

| Operation | Previous 160-bit baseline | Public-constant build | Reduction |
| --- | ---: | ---: | ---: |
| left tree, FN-DSA-512 | 140,133,646 | 136,242,017 | 2.78% |
| right tree, FN-DSA-512 | 146,920,250 | 143,068,224 | 2.62% |
| seeded sign, FN-DSA-512 | 352,691,418 | 339,184,430 | 3.83% |
| random sign, FN-DSA-512 | 353,524,328 | 339,887,910 | 3.86% |
| left tree, FN-DSA-1024 | 302,203,828 | 293,376,874 | 2.92% |
| right tree, FN-DSA-1024 | 316,066,727 | 306,128,352 | 3.14% |
| seeded sign, FN-DSA-1024 | 761,967,224 | 730,999,140 | 4.06% |
| random sign, FN-DSA-1024 | 751,302,020 | 725,978,675 | 3.37% |

No complete-signature measurement regressed by 1%; all four decreased. The
test-only fast-product oracle still observed maxima `145/143/129/143/135`
bits for the five recursive sites, below the signed-160 limit. These observed
values are regression data only; Lean range proofs do not assume them.

A Cortex-M4 object build is checked separately, but no DWT number is claimed
here because this host has no Cortex-M4 DWT hardware. The fixed demo LFSR
remains **DEMO ONLY** and these timings do not support a production leakage
claim.
