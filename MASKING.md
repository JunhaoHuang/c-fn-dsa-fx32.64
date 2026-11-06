# Experimental 160-bit Masking Backend

## Status and scope

This directory contains the first implementation milestone for a two-share,
first-order, 160-bit masking backend for the Q32.64 FN-DSA signer. It is a
research scaffold and is **DEMO ONLY**.

The v2 milestone implements the complete standard FN-DSA-512/1024 dynamic
signing data path. The signer expands separate q/q-prime shares and keeps the
basis, Gram/cross polynomials, FFT/LDL tree, samples, final `s0/s1`, squared
norm and compressed-length decision shared until the documented public
boundaries. It never reconstructs an ordinary private key or calls the
unmasked signer. Existing stable APIs and their behavior are unchanged.

| Component | Current status |
| --- | --- |
| Separate experimental API | Implemented in `fndsa_masked.h` |
| V2 seed-compressed masked key for logn 9/10 | Implemented; V1 rejected |
| Stored public key, 40-byte hedge token and corrections for f/g/F/G | Implemented |
| Two-share arithmetic values in Z/(2^160) | Implemented |
| Word-sliced A2B/B2A, refresh and 160/192-bit multiplication | Implemented |
| Floor, inverse, division, square root, expm, BerExp and SamplerZ | Implemented |
| Masked FFT, LDL and dynamic ffSampling | Implemented and wired to signing |
| Separate-share q/q-prime NTT/iNTT, Gram/cross and CRT centering | Implemented and wired |
| Final masked s0/s1, range, norm and compressed-length checks | Implemented |
| Bit-exact seeded FN-DSA-512/1024 signatures | Implemented by the masked core |
| Cortex-M4 assembly optimization and physical leakage validation | Not implemented |

The new API remains outside `fndsa.h` because the fixed LFSR and observable
compatibility control flow make this a reproducible research backend, not a
production side-channel-secure API.

## Masked-key encoding

For degree `n = 2^logn`, the object size is:

```
72 + FNDSA_VRFY_KEY_SIZE(logn) + 16*n
```

All multi-byte integers are little-endian.

| Offset | Length | Contents |
| --- | ---: | --- |
| 0 | 4 | Magic `FMK2` |
| 4 | 1 | Version 2 |
| 5 | 1 | `logn`, restricted to 9 or 10 |
| 6 | 1 | Share count, exactly 2 |
| 7 | 1 | Demo flags |
| 8 | 4 | Total encoded length |
| 12 | 4 | Reserved, zero |
| 16 | variable | Standard encoded FN-DSA verification key |
| next | 16 | Share-expansion seed |
| next | 40 | `SHAKE256(f || g)` hedging compatibility token |
| next | `16*n` | Four canonical correction polynomials, uint32 modulo M |

Here `M = 12289 * 18433 = 226523137`. SHAKE256 expands share 1 with
domain separation by format version, degree, and polynomial number. The
stored canonical correction is `share0 = secret - share1 mod M`. V1 objects
are rejected; no legacy-key fallback exists.

The 40-byte hedge token reproduces the ordinary signer's `"20261006"` seed
schedule. It is deliberately stored unmasked and is an explicit trusted
compatibility boundary; it is not claimed to be protected by masked SHAKE.

Import is trusted provisioning. It decodes `f`, `g`, and `F`, derives `G` and
the public key with the existing modular implementation, validates the key,
and writes the v2 object. Import is outside the masking boundary.

## Arithmetic behavior

`mrx160` holds two arithmetic shares, each as five little-endian 32-bit
limbs. Their sum is interpreted modulo `2^160`. `mb160` has the same
shape and reconstructs with XOR.

A2B follows the proven construction recalled by Schneider--Paglialonga--Oder--
Gueneysu (SPOG19): each arithmetic input share is first independently encoded
as a complete Boolean sharing, then the two encodings are added. The adder is a
fixed-topology, word-sliced Kogge--Stone network in the style of the secure
logarithmic adder used by CGTV15. The 160- and 192-bit variants operate on five
or six 32-bit words and have eight public carry-prefix stages instead of a
bit-serial 160/192-step chain. Every nonlinear word AND uses a fresh 32-bit ISW
mask, and the Boolean result receives a final fresh refresh. The implementation
never represents an arithmetic share as `(a_i, 0)`: that former shortcut made
the second Boolean output share's least-significant bit constant zero and put
the reconstructed secret bit directly in the other share.

B2A is exactly the `n = 2` specialization of BCZ18 Algorithm 1. Its `C_2`
instance expands the two inputs to three refreshed Boolean shares, evaluates
`Psi(x,y) = (x XOR y) - y`, and reduces both `C_1` leaves with XOR. The same
word-oriented routine covers `2^160` and `2^192`; multi-limb subtraction is
exact, and all conversion temporaries are explicitly erased.

The literature result instantiated here is **first-order SNI only** because
`FNDSA_MASK_SHARES` is 2. It assumes fresh independent randomness and the
software probing model with 32-bit words. It is not a proof for the fixed demo
LFSR, compiler-induced transition leakage, glitches, register spills, or the
compiled binary. The conversion functions and nonlinear word-AND boundary are
`noinline`, use compiler value barriers, and the masked module is built without
LTO, but generated assembly and physical leakage still require separate audit.
The regression suite checks exact reconstruction and verifies that every bit of
each conversion output share takes both values over repeated fresh masks; this
catches the old constant-LSB defect but is intentionally not presented as a
statistical or formal SNI proof.

Primary references: [SPOG19](https://doi.org/10.1007/978-3-030-17259-6_18),
[BCZ18](https://tches.iacr.org/index.php/TCHES/article/view/873), and the
[FALCON masked-floor gadget table](https://cic.iacr.org/p/1/4/9/pdf). Historical
performance comparisons remain in
[`bench_masked/masked_falcon_b2a_comparison.md`](bench_masked/masked_falcon_b2a_comparison.md)
and
[`bench_masked/float_vs_fixed_comparison.md`](bench_masked/float_vs_fixed_comparison.md).

The fast `mrx160_mul()` requires a call-site range proof that the true
pre-shift product fits in signed 160 bits. `mrx160_mul_wide()` is the fixed
fallback: it converts both inputs to Boolean shares, sign-extends the 96-bit
Q32.64 payloads to 192 bits, multiplies in the 192-bit shared ring, adds
`2^63`, performs the arithmetic shift by 64, truncates to the legacy
96-bit result, and converts back to 160-bit arithmetic shares. Both functions
match the actual legacy rule: add `2^63`, then arithmetic-right-shift by 64.

Public constants use dedicated `mrx160_mul_public_fast()` and
`mrx160_mul_public_wide()` paths. They multiply each arithmetic share by the
same sign-extended public constant, then perform one global masked rounding and
refresh; they never call `ring_mul()` or `ring_mul192()`. The fast and wide
cores consume 285 and 621 random words, then discard 5 and 288 fresh words so
the final mask-RNG states remain compatible with the old 290/909-word generic
calls. Split/merge, proved recursive tree sites, inverse/division seeds,
BerExp, SamplerZ and leaf scales use public-fast; ordinary FFT/LDL sites without
a signed-160 theorem use public-wide. The C-to-Lean map is in
[`lean/CALL_SITES.md`](lean/CALL_SITES.md).

The same sharing must not be supplied as both logical inputs of the fast
multiplier. `mrx160_sqr()` instead makes two independently refreshed sharings
of its input before invoking `mrx160_mul()`; the square operations in the
square-root iteration and SamplerZ use this dedicated gadget.

### Standalone fixed-point sampling gadgets

`mask_sampler.c` adds a clean-room fixed-point gadget layer. It takes only the
Boolean data-flow ideas (secret barrel shifting, masked comparison, and a
masked acceptance bit) from masked_FALCON commit
`4c103fb774adfb41edc406d8d52283b9b74418c9`; that repository does not state a
license, so no IEEE-754 encoding, source code, or LDL implementation is copied.

All sampler Boolean nonlinear nodes now call the common
`fndsa_mask_bool_and_words()` boundary. It independently refreshes both input
sharings, seeds the output before the first cross-product, and accumulates the
four products separately behind value barriers. OR/select/nonzero and 64-bit
select have explicit output refreshes; the byte comparator refreshes `eq` and
`lt` independently. The conditional MIMO-SNI fan-out map is in
[`security/SAMPLER_BOOLEAN_MIMO.md`](security/SAMPLER_BOOLEAN_MIMO.md).

`mrx160_floor()` extracts bits 64 through 95 after A2B and sign-extends the
integer result. `mrx160_inv_spec()` uses fixed-topology secret normalization,
the legacy `48/17 - 32*b/17` seed, and the same five Newton steps and rounding
order as `rx_inv_spec()`. `mrx160_expm_p63()` evaluates the original 13-term
FACCT polynomial with an exact masked unsigned high-half multiply; it does not
reuse the rounded Q32.64 multiply.

BerExp keeps the threshold Boolean-shared, saturates and applies the secret
shift with a barrel shifter, and performs the legacy most-significant-byte-first
lazy comparison. Only the per-byte equality needed to continue and the final
acceptance decision are opened at the explicit compatibility boundary.
SamplerZ keeps `mu`, `isigma`, the rejection exponent, and the returned integer
shared. Its BaseSampler/sign-bit bytes come from the independent signature
source callback, while all gadget masks come from `fndsa_mask_rng`.

The signer invokes these sampler gadgets through the masked dynamic tree.
SamplerZ/BerExp compatibility openings remain observable, so their integration
does not make the overall trace oblivious.

Masked modules are compiled with `-fno-lto`; nonlinear entry points are
`noinline` and contain compiler value barriers.

## Internal masked dynamic sampling backend

The research-only mask_xpoly.c and mask_tree.c modules now implement the
split FFT representation, complex multiplication, split/merge, self-adjoint
LDL, general masked division, square root, and recursive left/right dynamic
ffSampling. They mirror the operation order and special logn=1/2 branches of
the legacy fixed-point implementation.

The dynamic temporary-area query is fndsa_masked_dyn_tmp_size(). After
32-byte alignment, the area contains 2*n mrx160 work values and n mi32 sample
values: 45,056 bytes for FN-DSA-512 and 90,112 bytes for FN-DSA-1024, plus
32 bytes of alignment allowance. Recursion allocates no arrays proportional
to n on the C stack.

SamplerZ is invoked through its existing byte callback, so signature SHAKE
bytes and mask randomness remain separate. The only opened decisions remain
the documented BerExp prefix-continuation and SamplerZ acceptance boundaries.
The dynamic module itself never calls mrx160_open_rx().

Range selection remains call-site specific. Public unsigned Q0.64 x_SIN/IRNJ
products and the recursive-tree LDL/iXadj products use the 160-bit gadget.
Their exact raw-product bounds are proved in Lean from the existing key-norm
check, positive split/LDL invariants, and the SamplerZ residual bound. General
polynomial helpers, root complex products, RNJ constants, and root
self-adjoint-division products retain the exact 192-bit fallback. Test builds
also run an independent 96x96-to-192 range oracle. The proof-to-C map is in
lean/CALL_SITES.md.

This backend is now the dynamic-sampling engine used by the public experimental
masked signer. Its q/q-prime inputs come from the separate-share modular layer;
its samples feed masked final-vector generation and aggregate acceptance.

## Constant-time reduction audit

The masked production sources (`mask160.c`, `mask_sampler.c`,
`sign_masked.c`, `mask_mq.c`, `mask_xpoly.c`, and `mask_tree.c`) contain no
secret C integer division or remainder expressions. Persistent q/q' values and
all NTT butterflies/pointwise products use 32-bit Montgomery arithmetic with
R=2^32. Generic 64-bit reduction, uniform mask sampling and the small CRT
helper use Barrett reduction with precomputed `floor(2^64/q)` reciprocals. The
high half of the 64-by-64 product is assembled from four 32-bit products, so
the same source does not require native 128-bit integers. Other public 32-bit
moduli remain supported through a fixed 64-step reciprocal builder. Final
conditional subtraction uses masks/conditional moves, not a secret branch.
Generated x86-64 objects contain no `div` or `idiv` instruction.

The small-modulus fast path uses Raccoon-style signed Montgomery REDC.
For `q <= INT32_MAX`, canonical products and decoding need only a conditional
addition after REDC. Arbitrary-word encoding uses this path only when its
public `R2` constant satisfies the signed input bound; other encodings and
larger moduli retain the unsigned fallback.

NTT butterflies keep signed lazy values internally, removing per-butterfly
canonical corrections. At most ten stages give `|x| < 11*q <= 202763`;
products are explicitly 64-bit. Forward outputs are normalized once using
Montgomery one; inverse outputs need only a conditional addition. All API
boundaries and every secret-multiplication accumulator remain canonical.
Uniform random masks are sampled directly as Montgomery residues, avoiding
redundant encoding without changing their distribution or RNG draw count.
See `security/MONTGOMERY_PROOFS.md` and `security/NTT_CRT_PROOFS.md` for the
range arguments and masking constraints. These are not new machine-code
or end-to-end side-channel security claims.

Signed-path validation includes the full `test_masked` and `test_fndsa`
suites, Clang ASan/UBSan instrumentation of the changed arithmetic and test
translation units, and GCC compilation/testing of those same units linked
against unchanged support objects. Tests cover signed REDC input limits,
moduli on both sides of `INT32_MAX`, arbitrary-word imports, in-place products,
and patterned NTT/iNTT inputs for both primes at every supported degree.
Cortex-M4 and Windows runtime validation remain outstanding.

The focused driver `bench_masked/bench_signed_mont.c` reuses the existing
benchmark routines. `bench_masked/signed_montgomery_tsc.csv` records medians
of 11 paired runs on CPU 0 of an AMD Ryzen 9 9950X3D, Clang 18.1.3 `-O2`,
on 2026-11-04. The baseline is the working tree saved immediately before
this change, not the repository's clean HEAD. Units are TSC/reference ticks,
not core cycles; these are local microbenchmarks, not whole-signing speedups.
Multiplication takes 16.7–17.7% fewer ticks and masked NTT/iNTT takes
26.8–30.4% fewer ticks in this run; tiny add/sub/open differences are noise-scale.
To build and run the current focused driver from the repository root:

```sh
make test_masked
objects=$(find . -maxdepth 1 -name '*.o' ! -name 'test*' ! -name 'speed*' ! -name 'sign_masked_test.o' -print)
clang -O2 -I. bench_masked/bench_signed_mont.c $objects -o /tmp/bench_signed_mont
taskset -c 0 /tmp/bench_signed_mont
```

Canonical-key validation and reconstructed-coefficient range checks scan every
coefficient and aggregate an error bit; they no longer return at the first bad
secret coefficient. Centering and modular correction are branchless. Loops,
array indices, shift amounts, and modulus-selection branches in the arithmetic
gadgets depend only on public dimensions or public moduli.

The following variable-time behavior intentionally remains and limits the
security claim:

- Exact uniform sampling modulo a non-power-of-two modulus uses rejection in
  `fndsa_mask_rng_mod()` and masked-key SHAKE expansion. Its iteration count
  depends on fresh mask randomness, not the masked value. Removing it without
  bias requires a separately specified fixed-batch sampler and failure policy.
- RNG reservoir branches depend on the public gadget call schedule. A future
  secret-dependent gadget schedule would violate this assumption.
- Standalone BerExp deliberately opens the equal-prefix continuation bit and
  SamplerZ opens the final acceptance bit. Thus equal-prefix length, rejection
  count, and the enclosing retry schedule remain observable by design.
- The real masked signer opens only the aggregate key-valid bit, documented
  BerExp/SamplerZ control bits, the aggregate attempt retry bit, and the
  accepted `s1` vector that the signature publishes. Total retries and lazy
  comparison prefixes remain observable, so the signer is not trace-oblivious.
- Ordinary-key import and public-key/G derivation are trusted provisioning and
  remain outside the masking boundary.

## Demo randomness warning

The mask source is the LFSR-127 with polynomial
`x^127 + x^64 + 1`, initialized from the fixed bytes:

```
F0 E1 D2 C3 B4 A5 96 87 78 69 5A 4B 3C 2D 1E 0F
```

The RNG APIs expose a single continuous little-endian bit stream. A 64-bit
reservoir serves mixed 1-, 32-, and 64-bit requests without returning a bit
twice or discarding partial words. Word-sliced gadgets consume full independent
mask words; future bit-oriented gadgets must use fndsa_mask_rng_next_bit()
instead of generating a word and discarding 31 bits.

It is deterministic, predictable, and repeats across top-level calls. It is
suitable only for reproducible functional tests. It is not a cryptographic
PRNG and provides no production masking claim. The signature SHAKE stream is
separate, so the demo mask stream does not change nonce or signature sampling
bytes.

## Building and testing

Run:

```sh
make
./test_fndsa
./test_masked
./speed_masked
```

`speed_masked` measures the individual masking gadgets, including A2B/B2A,
The benchmark includes direct generic-fast/public-fast and
generic-wide/public-wide comparisons.
fixed-point arithmetic, sampler operations, FFT/LDL, and complete left/right
dynamic trees for FN-DSA-512/1024. The internal 192-bit helpers are measured
through `mrx160_mul_wide()`.

On Linux x86, the default mode opens a pinned `PERF_COUNT_HW_CPU_CYCLES`
event, so the kernel selects the correct PMU on both AMD and Intel. Run it
with root/CAP_PERFMON, or temporarily set `kernel.perf_event_paranoid=0`.
This avoids the Intel-only fixed-counter index used by the original benchmark.
`./speed_masked -tsc` remains the unprivileged fallback; those results are
explicitly labelled `TSC-ticks` and must not be reported as core cycles.

The current integrated full capture is stored in
`bench_masked/mask_integrated_signing_tsc.txt`, with a concise bottleneck
summary in `bench_masked/integrated_signing_summary.md`. Both are explicitly
labelled as TSC/reference baselines, not hardware core-cycle claims.

`test_masked` covers the LFSR/reservoir KAT; exact five-limb conversion;
per-bit output-share randomization regression for A2B/B2A; arithmetic; fast and
192-bit multiplication; floor, inverse, division and square-root differentials;
FACCT, BerExp and 100000 SamplerZ samples per parameter set; FFT/split/merge
coefficient differentials for logn 1 through 10; LDL differentials for logn 1
through 9; exact 512/1024 left/right dynamic tree samples, FFT outputs and
post-call SHAKE states; q/q-prime NTT/iNTT, pointwise, Gram/cross and CRT
centering for logn 1 through 10; real 512/1024 masked-key basis expansion;
modular arithmetic; separate coefficient-range, squared-norm and compressed-
length rejection; empty/raw, maximum-context, randomized-message and SHA-256
prehash signature cases; mask-RNG perturbation invariance; forced 27-attempt
cleanup; malformed v1/v2 keys and failure
zeroization; exact seeded signature compatibility; and verification with the
existing verifier.

The RNG-perturbation and forced-retry controls are compiled only into
`sign_masked_test.o` under `FNDSA_MASK_TEST_HOOKS`; `test_masked` links that
object instead of `sign_masked.o`. Production binaries contain neither the
control function nor the forced-rejection branch.

## Remaining production work

The functional end-to-end masked data path is connected. A production security
claim still requires replacing the fixed LFSR with independently seeded
cryptographic masks or physical entropy, protecting or redesigning the
unmasked hedging boundary, closing the FFT/LDL MIMO-SNI and compiled-machine
proof obligations, and validating physical leakage on each target. A
trace-oblivious profile would additionally need fixed-trace BerExp, SamplerZ
and attempt scheduling. The current implementation deliberately makes none of
those production claims.
