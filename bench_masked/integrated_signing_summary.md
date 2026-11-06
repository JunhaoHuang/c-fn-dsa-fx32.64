# Integrated masked signing benchmark (TSC reference)

Recorded on 2026-10-31 with the optimized portable-C build using
`./speed_masked -tsc`. These are TSC/reference ticks, **not hardware core
cycles**; they are suitable for same-host bottleneck comparisons only. The full
raw capture is in `mask_integrated_signing_tsc.txt`.

| Operation | FN-DSA-512 | FN-DSA-1024 |
| --- | ---: | ---: |
| Plain seeded sign | 1,127,116 | 2,388,177 |
| Masked seeded sign | 285,805,606 | 613,793,825 |
| Masked/plain ratio | 253.6x | 257.0x |
| Masked basis expansion | 207,862 | 411,338 |
| Masked final vectors (`s1`, `s0`) | 1,983,762 | 3,998,871 |
| Masked aggregate range/norm/length acceptance | 5,376,333 | 10,743,077 |
| Masked left ffSampling | 114,565,889 | 246,291,487 |
| Masked right ffSampling | 121,000,237 | 259,345,298 |
| Both dynamic trees / total masked sign | 82.4% | 82.4% |

The dominant bottleneck is the pair of recursive dynamic trees, not modular
NTT/CRT or final-vector generation. Within the tree, 192-bit wide masked
multiplication, A2B/B2A conversion, division/square-root iterations and
SamplerZ account for the bulk of work. The aggregate acceptance scan is the
largest new post-sampling stage (about 1.9%/1.8% of total); final-vector
construction is below 1%.

For scale, one masked q/q-prime Gram at logn 10 costs roughly 102k/103k ticks,
whereas one complete left/right tree costs about 246M/259M ticks. Consequently,
optimizing the dynamic fixed-point gadgets and recursive tree offers much more
leverage than further NTT micro-optimization.
