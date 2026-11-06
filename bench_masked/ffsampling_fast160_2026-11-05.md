# ffSampling 160-bit fast-path A/B benchmark

Date: 2026-11-05
Host: AMD Ryzen 9 9950X3D 16-Core Processor
Compiler: clang, `-O2`; masking modules use `-fno-lto`
Counter: invariant TSC/reference ticks, **not PMU core cycles**
Command: `./speed_masked -tsc`

Both builds use the real signing flags: left tree `NOOUT`, right tree `INIT0`.
The baseline is compiled with `FNDSA_MASK_TREE_FORCE_WIDE=1`; the optimized
build uses the Lean-approved 160-bit sites. Each reported operation is the
benchmark harness median after its normal warm-up/sampling schedule.
The optimized path discards the fixed 619-word consumption difference after
each fast multiply, so both builds leave the mask RNG in exactly the same
state. This compatibility cost is included in every optimized measurement.

| Operation | Forced 192-bit | 160-bit sites | Reduction |
| --- | ---: | ---: | ---: |
| left tree, FN-DSA-512 | 143,197,783 | 140,133,646 | 2.14% |
| right tree, FN-DSA-512 | 149,557,053 | 146,920,250 | 1.76% |
| both trees, FN-DSA-512 | 292,754,836 | 287,053,896 | 1.95% |
| seeded sign, FN-DSA-512 | 356,602,440 | 352,691,418 | 1.10% |
| random sign, FN-DSA-512 | 356,675,884 | 353,524,328 | 0.88% |
| left tree, FN-DSA-1024 | 307,782,906 | 302,203,828 | 1.81% |
| right tree, FN-DSA-1024 | 320,971,264 | 316,066,727 | 1.53% |
| both trees, FN-DSA-1024 | 628,754,170 | 618,270,555 | 1.67% |
| seeded sign, FN-DSA-1024 | 774,840,005 | 761,967,224 | 1.66% |
| random sign, FN-DSA-1024 | 766,853,701 | 751,302,020 | 2.03% |

The test-only exact 96x96-to-192 oracle observed maximum product bit lengths
`145/143/129/143/135` for LDL theta, LDL Schur, iXadj theta, iXadj mu, and
iXadj final respectively while running synthetic trees and real seeded
FN-DSA-512/1024 signing. These maxima are regression data, not Lean premises.

The signature SHAKE stream, output bytes, and internal mask-RNG state remain
unchanged. The 160-bit gadget consumes 290 fresh 32-bit words and then discards
the 619-word difference from the former 909-word wide path. No consumed word is
reused as a mask. The fixed LFSR remains DEMO ONLY.
