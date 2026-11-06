# Masked floating-point versus masked Q32.64 performance

## Scope

This comparison uses `burak231708/masked_FALCON` commit
`4c103fb774adfb41edc406d8d52283b9b74418c9` and this repository's portable
160-bit Q32.64 masking backend. Both configurations use first-order masking
with two shares and were compiled by Clang 18.1.3 with `-O2 -fno-lto`.
Measurements were taken in one process pinned to logical CPU 2 on an AMD Ryzen
9 9950X3D. The reported unit is invariant-TSC reference ticks, not hardware
core cycles. Each complete run uses medians of 101 samples after 20 warm-ups;
the table reports the median of three complete runs.
The microbenchmark measures steady-state call throughput rather than a single
dependency-chained operation latency.

The benchmark source is [`compare_float_fixed.c`](compare_float_fixed.c). The
original primitive results are preserved in
[`float_vs_fixed_tsc.txt`](float_vs_fixed_tsc.txt); the new corresponding
sampler-gadget results, inputs, ratio ranges, and metadata are in
[`float_vs_fixed_sampler_gadgets_tsc.txt`](float_vs_fixed_sampler_gadgets_tsc.txt).
The benchmark calls masked_FALCON's public API; it copies none of that
repository's source code.

## Comparable numerical operations

| Numerical operation | masked floating point | masked Q32.64 | Result |
|---|---:|---:|---:|
| addition | 25,206.60 | 7.85 | fixed is about 3,214x faster |
| multiplication | 15,755.20 | 1,198.84 (160-bit fast path) | fixed is 13.14x faster |
| multiplication | 15,755.20 | 3,844.20 (192-bit fallback) | fixed is 4.10x faster |

These rows compare the same mathematical operations on exactly representable
positive inputs (`5 + 3` and `5 * 3`), and both reconstructed results are
validated before timing. They do not imply identical domains: IEEE-754 has a
dynamic exponent, while Q32.64 has a fixed range. In particular, masked
floating-point addition must perform secret exponent comparison/alignment,
normalization, and rounding; Q32.64 addition is a linear, share-wise multi-limb
addition. The very large addition ratio is therefore a representation-level
advantage, not a like-for-like comparison of the same internal circuit.

The fast fixed multiplier may be used only where range analysis proves that
the pre-shift product fits signed 160 bits. The 192-bit fallback is the
appropriate conservative comparison when that proof is unavailable.

## Representation support operations

| Operation | masked floating point | masked Q32.64 | Interpretation |
|---|---:|---:|---|
| create sharing | 53.47 (64-bit Boolean) | 24.08 (160-bit arithmetic) | fixed demo LFSR makes this 2.22x faster despite greater width |
| refresh | 58.44 (64-bit Boolean) | 25.80 (160-bit arithmetic) | fixed is 2.27x faster; RNG cost is included |
| A2B | 1,894.58 (16-bit exponent) | 926.22 (160-bit value) | fixed is 2.05x faster, but widths and circuits differ |
| B2A | 113.73 (16-bit exponent) | 174.15 (160-bit value) | floating is 1.53x faster, but converts only one tenth as many bits |

The B2A row is not evidence that the floating representation has a faster
full-value conversion: masked_FALCON converts a 16-bit exponent there, while
the Q32.64 backend converts all 160 bits. Per converted bit, the fixed Psi-B2A
is substantially cheaper, although per-bit normalization is only a rough
diagnostic and not a cryptographic cost model.

## Corresponding floor and sampling gadgets

The fixed-point implementations are now measured in the same process and with
the same compiler as their masked floating-point counterparts. Values are TSC
reference ticks per operation; the ratio is the median of the three complete-run
ratios.

| Operation | masked floating point | masked Q32.64 | float/fixed |
|---|---:|---:|---:|
| floor | 17,591.30 | 1,196.69 | 14.70x |
| inverse | 1,423,042.00 | 30,710.60 | 46.36x |
| ApproxExp / `expm_p63` | 468,227.00 | 13,235.40 | 35.36x |
| BerExp | 558,484.00 | 19,066.20 | 29.29x |
| SamplerZ, isigma already available | 3,660,633.00 | 26,660.00 | 137.25x |
| SamplerZ, including masked inverse | 3,660,633.00 | 57,534.00 | 63.63x |

The last two rows need different interpretations. masked_FALCON's `SamplerZ`
accepts `sigma` and computes its masked inverse internally. The fixed FN-DSA
sampler API accepts `isigma`. Thus **63.63x is the like-interface pipeline
comparison**; 137.25x isolates the sampling work when `isigma` is already
available from the tree. The three-run full-pipeline ratio ranged from 62.07x
to 63.88x.

Floor used 5.5, inverse used 5.0, and both exponential paths used `x=0.5` and
`ccs=0.5`. BerExp used the same 64-bit uniform word on both sides and the
accept decisions were checked equal. The fixed implementation's lazy byte
comparison may consume fewer bytes, whereas the floating API accepts an
already-masked full word; that structural difference is part of the measured
implementations.

These are gadget-level results only. At the time of this older capture, the
top-level API still used the compatibility bridge. That status is historical:
the v2 real masked signing path is now connected and measured separately in
`integrated_signing_summary.md`; this table remains a gadget comparison.

## Validity and security limits

- The floating `SamplerZ(mu, sigma)` includes inverse; fixed
  `mrx160_sampler_z(mu, isigma)` does not. Both the isolated and full-pipeline
  figures are reported so this interface difference cannot inflate the main
  claim.

- masked_FALCON does not build with GCC 13 under its own CMake `-Werror`
  settings: `B2Aext128` triggers `stringop-overflow` diagnostics at
  `gadgets.c:573-574`. The comparison uses Clang for both implementations.
- UBSan stops masked_FALCON at `gadgets.c:161`, where `SecAdd` shifts a
  `uint64_t` by 64. The optimized Clang build produces correct results for the
  benchmark's validation cases, but the timings are not measurements of an
  undefined-behavior-free reference.
- An AddressSanitizer build does not report an explicit memory violation before
  validation, but the same masked floating-point `5 * 3` test reconstructs the
  wrong result. This sanitizer-sensitive behavior is additional evidence that
  the reference must be repaired before publication-quality cycle claims.
- masked_FALCON's `rand()`/`rand64()` and this project's fixed LFSR-127 are both
  demo-only non-cryptographic generators. Replacing either with production
  randomness can materially change performance.
- Both implementations expose data-dependent sampling control flow and do not
  establish a production side-channel security claim. Equal masking order
  does not by itself establish equivalent leakage resistance.

## Conclusion

On this host, the fixed-point representation is faster for every corresponding
masked gadget measured. The strongest comparable sampler claim is **about
63.6x for inverse plus SamplerZ**, with a 62.1x to 63.9x three-run range.
Individually, floor is about 14.7x faster, inverse 46.4x, `expm_p63` 35.4x,
and BerExp 29.3x. When `isigma` is already available, the fixed SamplerZ alone
is about 137.3x faster.

The performance gap comes primarily from avoiding secret IEEE-754 exponent
alignment, normalization, repacking, and repeated floating gadgets. It should
not be generalized into a full signing claim until masked FFT, LDL, dynamic
tree, modular expansion, and top-level signing are integrated and measured.
The masked_FALCON UB/sanitizer issues and both demo-only random generators also
prevent treating these TSC measurements as publication-quality security or
cycle claims without further repair and hardware-cycle experiments.
