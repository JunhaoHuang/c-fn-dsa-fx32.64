# Masked multiplication call-site range policy

This is the review map between the C backend, `DynamicBackend.lean`, and
`FFSamplingFastPath.lean`. A call may use `mrx160_mul()` only when its exact,
pre-shift raw product is covered below. Every unlisted or disabled Q32.64 call
uses `mrx160_mul_wide()` and `dynamic_fallback_covers_q32_products`.

## Recursive ffSampling calls

The fast functions are deliberately separate from the general polynomial API:
`mxpoly_ldl_tree_fft()` and `mxpoly_mul_iXadj_tree_fft()` are called only by
`mask_tree.c`; `mxpoly_ldl_fft()` and `mxpoly_mul_iXadj_fft()` remain wide-safe
for arbitrary research/test inputs.

| C operation | Strict raw-product bound | C path | Lean theorem |
| --- | ---: | --- | --- |
| split/merge `x_SIN` and `IRNJ` | `< 2^159` | fast 160 | `q0_64_twiddle_product_fits_s160` |
| tree LDL `theta*a0` | `< 2^155` | fast 160 | `ldl_theta_a0_fits_s160` |
| tree LDL `e*b1` Schur term | `2^65 * 2^89 = 2^154` | fast 160 | `ldl_schur_product_fits_s160` |
| tree iXadj `theta*b` | `2^66 * 2^65 = 2^131` | fast 160 | `ixadj_theta_b_fits_s160` |
| tree iXadj `mu*a` | `< 2^153` | fast 160 | `ixadj_mu_a_fits_s160` |
| tree iXadj final `bt*w` | `< 2^156` | fast 160 | `ixadj_final_fits_s160` |
| root `mcomplex_mul()` | unclosed without extra F/G filters | wide 192 | disabled `rootComplex` site |
| root self-adjoint division result multiply | unclosed without root-L10 filter | wide 192 | disabled `rootSelfAdjDiv` site |

The recursive envelope uses only the existing standard-key check
`||f||^2+||g||^2 <= 16822`, degree `n <= 1024`, positive mean/harmonic LDL
children, `|L10| < 1`, and the deterministic SamplerZ residual `|t-z| < 20`.
It does not add the paper-specific gamma or alpha key filters. Exact C table
encodings of theta and mu are checked in `test_masked.c`.

`MaskedMulSite.fastEnabled` is executable policy metadata. Its enabled product
exponents are 159, 132, 134, 155, 154, 131, 153, and 156, all strictly below
160; the two root sites are false. Test builds additionally reconstruct the
96-bit operands, compute an exact 192-bit magnitude, and abort above 159 bits.
Those observed maxima are regression evidence only, never proof hypotheses.

A fast multiply consumes 290 fresh 32-bit RNG words versus 909 for the wide
path. `mmul_tree_fast()` discards the fixed 619-word difference, so the next
call sees exactly the same mask-RNG state as before the optimization.

## Public-constant Q32.64 products

`mrx160_mul_public_fast()` multiplies each arithmetic share by the same public
96-bit constant, then performs one shared add-half/A2B/shift/B2A/refresh.
`mrx160_mul_public_wide()` sign-extends the secret once to 192 bits and uses the
same sharewise public multiplication there. The former consumes 285 core mask
words plus 5 discarded words; the latter consumes 621 plus 288 discarded words,
so their final RNG states match the former generic 290/909-word calls.

| C call site | Path | Strict raw-product theorem |
| --- | --- | --- |
| split/merge `x_SIN`, `IRNJ` | public fast | `q0_64_twiddle_product_fits_s160` |
| recursive LDL theta and iXadj theta/mu | public fast + existing 619-word tree pad | `ldl_theta_a0_fits_s160`, `ixadj_theta_b_fits_s160`, `ixadj_mu_a_fits_s160` |
| leaf `RX_INV_SIGMA`, `RX_IQ` | public fast | `sampler_public_scale_fits_s160` |
| inverse/division `C32_17*ub` | public fast | `inverse_seed_public_fits_s160`, `c32_17_raw_bound` |
| BerExp `INV_LOG2`, `LOG2` | public fast | `ber_exp_public_fits_s160`, exact table bounds in `PublicConstant.lean` |
| SamplerZ `SIGMA_MIN*isigma`, `INV_2SQRSIGMA0*y` | public fast | `sampler_sigma_min_public_fits_s160`, `sampler_inv_2sqr_public_fits_s160` |
| ordinary FFT `SQRT2`, `RNJ` | public wide | `public_wide_product_fits_s192` |
| general LDL/iXadj theta/mu | public wide | `public_wide_product_fits_s192` |

`reconstruct_mulPublic` proves the sharewise ring step, and
`roundShift64_add_pow2_160` proves that a 160-bit representative changes the
post-shift value only by a multiple of `2^96`. `PublicMulSite.fastEnabled` is
the executable proof policy: every enabled exponent is strictly below 160;
unproved/general sites are disabled and fixed to the 192-bit public fallback.

## Masked modular-polynomial obligations

| C area | Integer invariant | Lean theorem |
| --- | --- | --- |
| `mmod32` product before REDC | canonical residues and `q < 2^32` imply `a*b < q*2^32` | `montgomery_product_bound` |
| q/q' NTT butterfly add/sub | both operands canonical; q' is the larger modulus | `butterfly_sum_fits_u32` |
| one per-share CRT component | `0 <= a+q*u < M` | `crt_component_range` |
| sum of two CRT shares | `0 <= c0+c1 < 2*M < 2^29` | `crt_share_sum_range`, `twice_crt_modulus_lt_pow2_29` |
| whole-secret center | two conditional subtractions produce `[-H,H]` | `reduceM_range`, `centerM_range` |

These theorems justify widths and sign tests only. The corresponding first-order
NI/SNI simulations are documented separately in `security/NTT_CRT_PROOFS.md`.

## Integrated signing connection-layer obligations

| C area | Integer invariant | Lean theorem |
| --- | --- | --- |
| `mi32_to_mmod32` | two uint32 shares sum below `2*2^32`; carry/sign correction represents the signed word modulo p | `two_u32_sum_range`, `signed_share_mod_correction`, `signed_low_word_range` |
| q/q-prime iNTT center | one whole-secret conditional center lands in the symmetric interval | `centerQ_range`, `centerQ2_range` |
| masked square/norm accumulation | every accepted coefficient square is at most 705600; 2048 such terms fit in 64 bits | `accepted_coefficient_square_bound`, `signature_sqnorm_fits_u64` |
| compressed-length accumulation | each accepted coefficient consumes at most 15 bits; 1024 terms fit in 64 bits | `compressed_coefficient_bits_bound`, `compressed_length_fits_u64` |
| v2 key and caller arena | exact 512/1024 formulas do not overflow 32-bit `size_t` | `masked_arena_sizes`, `masked_arena_1024_size_bound`, `masked_key_v2_size_bound` |

These theorems establish arithmetic correctness and range only. The aggregate
acceptance gadget and complete signing call graph remain subject to the
MIMO-SNI and trace qualifications in `security/COMPOSITION.md`.
