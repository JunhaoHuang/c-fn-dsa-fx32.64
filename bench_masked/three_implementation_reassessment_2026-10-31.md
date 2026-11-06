# Three-implementation masked-gadget reassessment

> Historical pre-hardening result. On 2026-11-04 the common Boolean AND was
> changed to independently refresh both input sharings and composite sampler
> outputs gained explicit refreshes. The fixed-point columns below must be
> remeasured with the same-process three-way harness before further citation.
> The post-fix internal reference run and overhead are recorded in
> `mask_boolean_sni_tsc_2026-11-04.txt` and
> `mask_boolean_sni_overhead_2026-11-04.md`.

## Result

The current Q32.64 backend remains substantially faster than Burak's masked
IEEE-754 implementation for every corresponding gadget tested. Against the
Chen--Chen floating implementation, the fixed 160-bit fast multiplier is
about 2.00 times faster, but the conservative 192-bit fallback is about
1.40 times slower. This makes range proofs that enable the 160-bit path a
performance-critical part of the design, not merely a proof artifact.

All implementations use first-order masking with two shares. The third-party
sources are pinned to Chen--Chen commit
`24be736ad5544f02f2c11f217de708740d988c8e` and Burak commit
`4c103fb774adfb41edc406d8d52283b9b74418c9`.

## Strict three-way comparison

The table reports the median of five complete runs. Every complete run uses
101 timing samples after 20 warm-ups, and reports the sample median. The inputs
are exactly representable values 5 and 3, and reconstruction is checked before
timing.

| Operation | Current fixed Q32.64 | Chen--Chen float | Burak float |
|---|---:|---:|---:|
| addition | 7.80 | 6,611.25 | 25,234.55 |
| multiply, proved 160-bit path | 1,681.30 | 3,362.60 | 15,781.00 |
| multiply, 192-bit fallback | 4,714.23 | 3,362.60 | 15,781.00 |

The units are invariant-TSC reference ticks per operation, not PMU hardware
cycles. The measurements were made in one process, pinned to logical CPU 2 on
an AMD Ryzen 9 9950X3D, with Clang 18.1.3 and
`-O2 -fno-lto -std=c99`.

Derived comparisons:

- Fixed addition is about 848 times faster than Chen--Chen and 3,235 times
  faster than Burak. This reflects the representation: fixed addition is a
  linear limb operation, while floating addition must align and normalize a
  secret exponent.
- The fixed 160-bit multiplier is 2.00 times faster than Chen--Chen and
  9.39 times faster than Burak.
- The fixed 192-bit fallback is 1.40 times slower than Chen--Chen, but remains
  3.35 times faster than Burak.
- Burak multiplication is about 4.69 times slower than Chen--Chen.

## Extended gadget comparison

Chen--Chen exposes masked floating addition, multiplication and their support
circuits, but not standalone floor, inverse, exponential, BerExp or SamplerZ
interfaces matching the other two repositories. These rows are therefore an
honest two-way comparison rather than fabricated three-way entries.

| Operation | Burak float | Current fixed | Float/fixed |
|---|---:|---:|---:|
| floor | 14,073.90 | 1,083.17 | 12.99x |
| inverse | 1,422,139.00 | 37,453.00 | 37.97x |
| ApproxExp / `expm_p63` | 470,076.00 | 17,733.20 | 26.51x |
| BerExp | 560,634.00 | 24,570.20 | 22.82x |
| SamplerZ, fixed `isigma` precomputed | 3,667,900.00 | 34,271.00 | 107.02x |
| SamplerZ, fixed inverse included | 3,667,900.00 | 71,595.00 | 51.23x |

The last row is the better like-interface comparison because Burak's SamplerZ
accepts `sigma` and performs the inverse internally, whereas the fixed sampler
normally accepts `isigma`.

## Current bottlenecks

Within one full current-backend run, the expensive scalar gadgets were square
root (47,157 ticks), inverse (37,032), SamplerZ (35,260), BerExp (about
25,300), division (24,304) and `expm_p63` (17,492). A 192-bit multiply cost
3,775 ticks versus 1,321 for the proved 160-bit path. A2B cost 882 ticks and
B2A 138 ticks.

At polynomial level, the dominant cost is recursive dynamic sampling:

| Backend operation | FN-DSA-512 | FN-DSA-1024 |
|---|---:|---:|
| masked FFT | 8,652,632 | 19,168,583 |
| masked LDL | 8,300,763 | 16,613,222 |
| ffSampling left | 117,455,747 | 252,966,764 |
| ffSampling right | 124,101,827 | 266,204,615 |

The right recursion is about 5--6% slower because it includes the additional
deepest-level normalization/division work. The dynamic tree is dominated by
the repeated fixed-point nonlinear gadgets, especially multiplication,
division, inverse/square-root preparation and A2B-dependent shifts.

The new Montgomery modular multiplication is not a leading bottleneck: it
costs about 42 ticks for 12289/18433 and 46 ticks for 226523137, far below a
single 160-bit fixed-point multiplication.

## Change from the older report

The earlier saved table measured the fixed fast and wide multipliers at
1,227.65 and 3,812.67 ticks. The new values are respectively about 37% and
24% higher. The primary reason is that the current benchmark includes the
random-first first-order repairs to ring multiplication and the latest
conversion path. The new costs should be used for planning; the old numbers
should not be cited for the repaired implementation.

## Limits

- The host denies PMU access (`perf_event_paranoid=4`), so these are same-host
  relative TSC measurements, not Cortex-M4 cycles or true core-cycle counts.
- The three random sources are demo/test generators and their costs are
  included. They are not equivalent production CSPRNGs.
- Representation, range and security contracts differ. Equal masking order
  does not imply identical leakage guarantees.
- Burak's reference retains previously identified undefined-behaviour and
  sanitizer sensitivity; its optimized benchmark validates the selected
  inputs but should not be treated as a repaired production reference.
- At the time of this three-gadget capture, the top-level API still used the
  compatibility bridge. That statement is historical; current end-to-end v2
  measurements are in `integrated_signing_summary.md`.

Full raw medians and the five strict three-way runs are saved in
`three_implementation_recheck_tsc.txt`.
