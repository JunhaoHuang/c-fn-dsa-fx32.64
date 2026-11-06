# Boolean SNI hardening overhead (2026-11-04)

The post-fix run is `mask_boolean_sni_tsc_2026-11-04.txt`, captured with
Ubuntu Clang 18.1.3, `-O2 -fno-lto` for the masking modules, and
`taskset -c 2 ./speed_masked -tsc`. Values are TSC/reference ticks, not PMU
core cycles. The comparison baseline is `mask_integrated_signing_tsc.txt`;
therefore ratios are indicative single-run deltas, not a statistical result.

| Component | Previous | Input-isolated shared AND | Ratio | Increase |
|---|---:|---:|---:|---:|
| A2B-160 | 859.57 | 1,388.04 | 1.615x | 61.5% |
| Q32.64 multiply (160-bit) | 1,303.33 | 1,961.66 | 1.505x | 50.5% |
| floor | 1,027.27 | 1,468.88 | 1.430x | 43.0% |
| SamplerZ-512 | 34,228 | 44,720 | 1.307x | 30.7% |
| SamplerZ-1024 | 34,056 | 44,634 | 1.311x | 31.1% |
| seeded masked sign-512 | 285,805,606 | 365,916,240 | 1.280x | 28.0% |
| seeded masked sign-1024 | 613,793,825 | 789,379,681 | 1.286x | 28.6% |

The increase is expected: each 32-bit nonlinear AND now consumes independent
input-refresh masks `rx,ry` in addition to the output mask `r`, and composite
sampler outputs receive explicit refreshes. Signature bytes remain unchanged
because the mask RNG is separate from the signature SHAKE stream.
