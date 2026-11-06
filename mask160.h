#ifndef FNDSA_MASK160_H__
#define FNDSA_MASK160_H__

/*
 * Portable reference gadgets for two-share, first-order masking.
 *
 * Persistent fixed-point values are arithmetic shares in Z/(2^160), held
 * as five little-endian 32-bit limbs per share. Their reconstructed payload
 * is a sign-extended Q32.64 value. Boolean shares reconstruct with XOR.
 *
 * WARNING: the bundled LFSR is deterministic and is DEMO ONLY.
 */

#include <stddef.h>
#include <stdint.h>

#define FNDSA_MASK160_LIMBS   5
#define FNDSA_MASK_SHARES     2

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS];
} mrx160;

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS];
} mb160;

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES];
} mbit;

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES][2];
} mb64;

typedef struct {
	uint64_t s[2];
	uint64_t reservoir;
	unsigned reservoir_len;
} fndsa_mask_rng;

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES];
} mmod32;

/* Compact arithmetic sharing of a signed integer modulo 2^32. */
typedef struct {
	uint32_t share[FNDSA_MASK_SHARES];
} mi32;

void fndsa_mask_rng_init_demo(fndsa_mask_rng *rng);
uint64_t fndsa_mask_rng_next_u64(fndsa_mask_rng *rng);
uint32_t fndsa_mask_rng_next_u32(fndsa_mask_rng *rng);
uint32_t fndsa_mask_rng_next_bit(fndsa_mask_rng *rng);
uint32_t fndsa_mask_rng_mod(fndsa_mask_rng *rng, uint32_t modulus);

uint32_t fndsa_mask_reduce_u64(uint64_t x, uint32_t modulus);

/*
 * Shared Boolean nonlinear boundary. The AND implementation independently
 * refreshes both input sharings before evaluating any cross-product, then
 * accumulates every product into an output initialized with fresh randomness.
 * The refresh helper also isolates MIMO outputs. These are internal research
 * interfaces, not stable FN-DSA APIs.
 */
void fndsa_mask_bool_and_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *x0, const uint32_t *x1,
	const uint32_t *y0, const uint32_t *y1,
	size_t n, fndsa_mask_rng *rng);
void fndsa_mask_bool_refresh_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *x0, const uint32_t *x1,
	size_t n, fndsa_mask_rng *rng);

void mrx160_share_rx(mrx160 *d, const uint32_t x[3], fndsa_mask_rng *rng);
void mrx160_open_rx(uint32_t d[3], const mrx160 *a);
void mrx160_copy(mrx160 *d, const mrx160 *a);
void mrx160_add(mrx160 *d, const mrx160 *a, const mrx160 *b);
void mrx160_sub(mrx160 *d, const mrx160 *a, const mrx160 *b);
void mrx160_neg(mrx160 *d, const mrx160 *a);
void mrx160_refresh(mrx160 *a, fndsa_mask_rng *rng);

void mrx160_floor(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng);
void mrx160_half(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng);
void mrx160_inv_spec(mrx160 *d, const mrx160 *a,
	fndsa_mask_rng *rng);
void mrx160_div_spec(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng);
void mrx160_sqrt_spec(mrx160 *d, const mrx160 *a,
	fndsa_mask_rng *rng);

void mrx160_a2b(mb160 *d, const mrx160 *a, fndsa_mask_rng *rng);
void mrx160_b2a(mrx160 *d, const mb160 *a, fndsa_mask_rng *rng);

/*
 * Q32.64 multiply: add 2^63 and arithmetic-right-shift by 64. The fast
 * function requires the true raw product to fit signed 160 bits.
 */
void mrx160_mul(mrx160 *d, const mrx160 *a, const mrx160 *b,
        fndsa_mask_rng *rng);
/*
 * Q32.64 square with the same rounding rule as mrx160_mul(). The input is
 * independently refreshed into two logical sharings before multiplication;
 * callers must use this function instead of mrx160_mul(d, a, a, rng).
 */
void mrx160_sqr(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng);
void mrx160_mul_wide(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng);
/*
 * Public-constant Q32.64 multiplication. The fast form requires the exact
 * signed raw product to fit in 160 bits; the wide form accepts any two
 * Q32.64 payloads. Both preserve the generic gadgets' random-word budgets
 * (290 and 909 words respectively), and support d == a.
 */
void mrx160_mul_public_fast(mrx160 *d, const mrx160 *a,
	const uint32_t c[3], fndsa_mask_rng *rng);
void mrx160_mul_public_wide(mrx160 *d, const mrx160 *a,
	const uint32_t c[3], fndsa_mask_rng *rng);
void mrx160_mul_public_u64_fast(mrx160 *d, const mrx160 *a,
	uint64_t c, fndsa_mask_rng *rng);
void mrx160_mul_u64_hi(mb64 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng);

uint32_t mbit_open_compat(const mbit *a);
uint64_t mb64_open_compat(const mb64 *a);

/*
 * Masked arithmetic modulo an odd public modulus. Internal shares are
 * canonical Montgomery residues with R = 2^32; share/open are the domain
 * boundaries. Unsupported even moduli produce zero.
 */
void mmod32_share(mmod32 *d, uint32_t x, uint32_t modulus,
	fndsa_mask_rng *rng);
uint32_t mmod32_open(const mmod32 *a, uint32_t modulus);
void mmod32_add(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus);
void mmod32_sub(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus);
void mmod32_mul(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus, fndsa_mask_rng *rng);
void mmod32_refresh(mmod32 *a, uint32_t modulus,
	fndsa_mask_rng *rng);

/*
 * Internal Montgomery-domain helpers for the masked polynomial backend.
 * The public multiplier is already Montgomery encoded. The share import
 * and export functions never reconstruct their two input shares; they only
 * change the representation of each share independently.
 */
uint32_t mmod32_encode_public(uint32_t x, uint32_t modulus);
void mmod32_set_residue_shares(mmod32 *d, uint32_t s0, uint32_t s1,
	uint32_t modulus, fndsa_mask_rng *rng);
void mmod32_decode_residue_shares(uint32_t d[FNDSA_MASK_SHARES],
	const mmod32 *a, uint32_t modulus);
void mmod32_mul_public_mont(mmod32 *d, const mmod32 *a,
	uint32_t b_mont, uint32_t modulus);

void mi32_share(mi32 *d, int32_t x, fndsa_mask_rng *rng);
int32_t mi32_open_compat(const mi32 *a);
void mi32_to_mrx160(mrx160 *d, const mi32 *a, fndsa_mask_rng *rng);
void mrx160_to_mi32(mi32 *d, const mrx160 *a, fndsa_mask_rng *rng);

/*
 * Convert a signed Z/(2^32) sharing into an odd public prime field without
 * opening the sign or the carry of the arithmetic-share sum. The result is
 * returned in the canonical Montgomery representation used by mmod32.
 */
void mi32_to_mmod32(mmod32 *d, const mi32 *a, uint32_t modulus,
	fndsa_mask_rng *rng);

/*
 * Reduce a non-negative two-share integer s with 0 <= s < 2*M and center
 * it into [-(M-1)/2, +(M-1)/2]. M is the Falcon CRT modulus q*q'.
 * This is a masked comparison/select operation; it does not center shares
 * separately and does not open the comparison bits.
 */
void mi32_center_mod_m(mi32 *d, const mi32 *a, fndsa_mask_rng *rng);

/* Jointly reduce 0 <= share0+share1 < 2*modulus and center the result. */
void mi32_center_mod_public(mi32 *d, const mi32 *a, uint32_t modulus,
	fndsa_mask_rng *rng);

/*
 * Masked vector predicates used at the signing boundary. The returned bit is
 * Boolean shared and must be opened only as the aggregate key-valid/retry
 * decision. No coefficient or partial norm is opened by these functions.
 */
void mi32_poly_sqnorm_le(mbit *accept, unsigned logn,
	const mi32 *a, const mi32 *b, uint32_t bound,
	fndsa_mask_rng *rng);
void mi32_signature_is_acceptable(mbit *accept, unsigned logn,
	const mi32 *s0, const mi32 *s1, size_t encoded_len,
	fndsa_mask_rng *rng);

#endif
