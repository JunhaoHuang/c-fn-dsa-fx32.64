/*
 * Masked Falcon NTT/iNTT and two-prime CRT backend.
 *
 * Every linear transform is written as two explicit share-major passes: one
 * complete pass over share 0 followed by one complete pass over share 1. The
 * transform never reconstructs or combines shares. Boundary refresh is a
 * separate final pass.
 *
 * Gram/cross contain secret-by-secret products and therefore deliberately use
 * mmod32_mul(), which evaluates all four share products under fresh masking.
 * Computing only matching-share products would omit the two cross terms.
 *
 * WARNING: the bundled mask RNG is DEMO ONLY.
 */

#include <stddef.h>
#include <stdint.h>

#include "inner.h"
#include "sign_inner.h"
#include "mask_mq.h"
#include "mask_mont32.h"

static int
mod_data(fndsa_mask_modulus mod, uint32_t *q,
	const uint16_t **gm, const uint16_t **igm)
{
	switch (mod) {
	case FNDSA_MASK_MOD_Q:
		*q = FNDSA_MASK_Q;
		*gm = mq_GM;
		*igm = mq_iGM;
		return 1;
	case FNDSA_MASK_MOD_Q2:
		*q = FNDSA_MASK_Q2;
		*gm = mq2_GM;
		*igm = mq2_iGM;
		return 1;
	default:
		*q = 0;
		*gm = NULL;
		*igm = NULL;
		return 0;
	}
}

#define MASK_Q_Q0I   UINT32_C(0xF7002FFF)
#define MASK_Q2_Q0I  UINT32_C(0xEBC047FF)

/* q/q'-specialized operations on exactly one Montgomery share. */
static inline uint32_t
mq_add_share(uint32_t a, uint32_t b, uint32_t q)
{
	return mask_mont32_add(a, b, q);
}

static inline uint32_t
mq_sub_share(uint32_t a, uint32_t b, uint32_t q)
{
	return mask_mont32_sub(a, b, q);
}

static inline uint32_t
mq_redc_share(int64_t t, uint32_t q, uint32_t q0i)
{
	return mask_mont32_cadd(mask_mont32_redc(t, q, q0i), q);
}

static inline uint32_t
mq_mul_public_share(uint32_t a, uint32_t b_mont,
	uint32_t q, uint32_t q0i)
{
	return mq_redc_share((int64_t)a * b_mont, q, q0i);
}

static inline uint32_t
mq_decode_share(uint32_t a, uint32_t q, uint32_t q0i)
{
	return mq_redc_share(a, q, q0i);
}

static int
valid_logn(unsigned logn)
{
	return logn >= 1 && logn <= 10;
}

void
mmpoly_refresh(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mmod32_refresh(&a[i], q, rng);
	}
}

int
mmpoly_ntt(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return 0;
	}
	(void)igm;
	size_t n = (size_t)1 << logn;
	uint32_t mont_one = mmod32_encode_public(1, q);
	uint32_t q0i = mod == FNDSA_MASK_MOD_Q
		? MASK_Q_Q0I : MASK_Q2_Q0I;

	/* Two independent linear transforms; no share is read by the other. */
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		size_t t = n;
		for (unsigned lm = 0; lm < logn; lm ++) {
			size_t m = (size_t)1 << lm;
			size_t ht = t >> 1;
			size_t j0 = 0;
			for (size_t i = 0; i < m; i ++) {
				uint32_t s = gm[i + m];
				for (size_t j = 0; j < ht; j ++) {
					size_t j1 = j0 + j;
					size_t j2 = j1 + ht;
					int32_t x1 = mask_mont32_signed_word(a[j1].share[sh]);
					int32_t x2 = mask_mont32_redc((int64_t)
						mask_mont32_signed_word(a[j2].share[sh]) * s,
						q, q0i);
					a[j1].share[sh] = (uint32_t)(x1 + x2);
					a[j2].share[sh] = (uint32_t)(x1 - x2);
				}
				j0 += t;
			}
			t = ht;
		}
		for (size_t i = 0; i < n; i ++) {
			a[i].share[sh] = mq_redc_share((int64_t)
				mask_mont32_signed_word(a[i].share[sh]) * mont_one,
				q, q0i);
		}
	}
	mmpoly_refresh(logn, a, mod, rng);
	return 1;
}

int
mmpoly_intt(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return 0;
	}
	(void)gm;
	size_t n = (size_t)1 << logn;
	uint32_t q0i = mod == FNDSA_MASK_MOD_Q
		? MASK_Q_Q0I : MASK_Q2_Q0I;
	uint32_t inv2 = mmod32_encode_public((q + 1) >> 1, q);

	/* As above, complete share 0 before touching share 1. */
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		size_t t = 1;
		for (unsigned lm = 0; lm < logn; lm ++) {
			size_t hm = (size_t)1 << (logn - 1 - lm);
			size_t dt = t << 1;
			size_t j0 = 0;
			for (size_t i = 0; i < hm; i ++) {
				uint32_t s = igm[i + hm];
				for (size_t j = 0; j < t; j ++) {
					size_t j1 = j0 + j;
					size_t j2 = j1 + t;
					int32_t x1 = mask_mont32_signed_word(a[j1].share[sh]);
					int32_t x2 = mask_mont32_signed_word(a[j2].share[sh]);
					a[j1].share[sh] = (uint32_t)mask_mont32_redc(
						(int64_t)(x1 + x2) * inv2, q, q0i);
					a[j2].share[sh] = (uint32_t)mask_mont32_redc(
						(int64_t)(x1 - x2) * s, q, q0i);
				}
				j0 += dt;
			}
			t = dt;
		}
		for (size_t i = 0; i < n; i ++) {
			a[i].share[sh] = mask_mont32_cadd(
				mask_mont32_signed_word(a[i].share[sh]), q);
		}
	}
	mmpoly_refresh(logn, a, mod, rng);
	return 1;
}

void
mmpoly_add(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b, fndsa_mask_modulus mod)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || b == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			d[i].share[sh] = mq_add_share(
				a[i].share[sh], b[i].share[sh], q);
		}
	}
}

void
mmpoly_sub(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b, fndsa_mask_modulus mod)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || b == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			d[i].share[sh] = mq_sub_share(
				a[i].share[sh], b[i].share[sh], q);
		}
	}
}

void
mmpoly_neg(unsigned logn, mmod32 *d,
	const mmod32 *a, fndsa_mask_modulus mod)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			d[i].share[sh] = mq_sub_share(
				0, a[i].share[sh], q);
		}
	}
}

void
mmpoly_mul_ntt(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || b == NULL
		|| rng == NULL || !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		/* Nonlinear: all four share products are required. */
		mmod32_mul(&d[i], &a[i], &b[i], q, rng);
	}
}

void
mmpoly_mul_public_ntt(unsigned logn, mmod32 *d,
	const mmod32 *a, const uint32_t *b_mont,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || b_mont == NULL
		|| rng == NULL || !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	uint32_t q0i = mod == FNDSA_MASK_MOD_Q
		? MASK_Q_Q0I : MASK_Q2_Q0I;
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			d[i].share[sh] = mq_mul_public_share(
				a[i].share[sh], b_mont[i], q, q0i);
		}
	}
	for (size_t i = 0; i < n; i ++) {
		mmod32_refresh(&d[i], q, rng);
	}
}

void
mmpoly_from_mi32(unsigned logn, mmod32 *d, const mi32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mi32_to_mmod32(&d[i], &a[i], q, rng);
	}
}

int
mmpoly_public_ntt(unsigned logn, uint32_t *d, const uint16_t *a,
	fndsa_mask_modulus mod)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return 0;
	}
	(void)igm;
	size_t n = (size_t)1 << logn;
	uint32_t q0i = mod == FNDSA_MASK_MOD_Q
		? MASK_Q_Q0I : MASK_Q2_Q0I;

	/* Reverse order permits a uint16 source to alias the uint32 output. */
	for (size_t i = n; i -- > 0;) {
		d[i] = mmod32_encode_public(a[i], q);
	}
	for (size_t t = n, m = 1; m < n; m <<= 1) {
		size_t ht = t >> 1;
		size_t j0 = 0;
		for (size_t i = 0; i < m; i ++) {
			uint32_t s = gm[i + m];
			for (size_t j = 0; j < ht; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + ht;
				int32_t x1 = mask_mont32_signed_word(d[j1]);
				int32_t x2 = mask_mont32_redc((int64_t)
					mask_mont32_signed_word(d[j2]) * s, q, q0i);
				d[j1] = (uint32_t)(x1 + x2);
				d[j2] = (uint32_t)(x1 - x2);
			}
			j0 += t;
		}
		t = ht;
	}
	uint32_t mont_one = mmod32_encode_public(1, q);
	for (size_t i = 0; i < n; i ++) {
		d[i] = mq_redc_share((int64_t)
			mask_mont32_signed_word(d[i]) * mont_one, q, q0i);
	}
	return 1;
}

void
mmpoly_mul_public_const(unsigned logn, mmod32 *d,
	const mmod32 *a, uint32_t c, fndsa_mask_modulus mod,
	fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	uint32_t cm = mmod32_encode_public(c, q);
	size_t n = (size_t)1 << logn;
	uint32_t q0i = mod == FNDSA_MASK_MOD_Q
		? MASK_Q_Q0I : MASK_Q2_Q0I;
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			d[i].share[sh] = mq_mul_public_share(a[i].share[sh],
				cm, q, q0i);
		}
	}
	for (size_t i = 0; i < n; i ++) {
		mmod32_refresh(&d[i], q, rng);
	}
}

void
mmpoly_public_sub_ntt(unsigned logn, mmod32 *d,
	const uint32_t *a_public_mont, const mmod32 *b,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a_public_mont == NULL
		|| b == NULL || rng == NULL || !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		d[i].share[0] = mq_sub_share(
			a_public_mont[i], b[i].share[0], q);
		d[i].share[1] = mq_sub_share(0, b[i].share[1], q);
		mmod32_refresh(&d[i], q, rng);
	}
}

int
mmpoly_intt_center(unsigned logn, mi32 *d, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || a == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm)
		|| !mmpoly_intt(logn, a, mod, rng))
	{
		return 0;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mi32 raw;
		raw.share[0] = mq_decode_share(a[i].share[0], q,
			mod == FNDSA_MASK_MOD_Q ? MASK_Q_Q0I : MASK_Q2_Q0I);
		raw.share[1] = mq_decode_share(a[i].share[1], q,
			mod == FNDSA_MASK_MOD_Q ? MASK_Q_Q0I : MASK_Q2_Q0I);
		mi32_center_mod_public(&d[i], &raw, q, rng);
		raw.share[0] = 0;
		raw.share[1] = 0;
	}
	return 1;
}

void
mmpoly_gram_ntt(unsigned logn, mmod32 *d,
	const mmod32 *f, const mmod32 *g,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || f == NULL || g == NULL
		|| rng == NULL || !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mmod32 x, y;
		/* Secret products cannot be split into matching-share products. */
		mmod32_mul(&x, &f[i], &f[n - 1 - i], q, rng);
		mmod32_mul(&y, &g[i], &g[n - 1 - i], q, rng);
		mmod32_add(&d[i], &x, &y, q);
		mmod32_refresh(&d[i], q, rng);
	}
}

void
mmpoly_cross_ntt(unsigned logn, mmod32 *d,
	const mmod32 *F, const mmod32 *f,
	const mmod32 *G, const mmod32 *g,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng)
{
	uint32_t q;
	const uint16_t *gm, *igm;
	if (!valid_logn(logn) || d == NULL || F == NULL || f == NULL
		|| G == NULL || g == NULL || rng == NULL
		|| !mod_data(mod, &q, &gm, &igm))
	{
		return;
	}
	(void)gm;
	(void)igm;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mmod32 x, y;
		/* Each gadget evaluates 00, 01, 10 and 11 under fresh masking. */
		mmod32_mul(&x, &F[i], &f[n - 1 - i], q, rng);
		mmod32_mul(&y, &G[i], &g[n - 1 - i], q, rng);
		mmod32_add(&d[i], &x, &y, q);
		mmod32_refresh(&d[i], q, rng);
	}
}

void
mmpoly_crt_center(unsigned logn, mi32 *d,
	const mmod32 *a_q, const mmod32 *a_q2,
	fndsa_mask_rng *rng)
{
	if (!valid_logn(logn) || d == NULL || a_q == NULL || a_q2 == NULL
		|| rng == NULL)
	{
		return;
	}
	size_t n = (size_t)1 << logn;

	/* CRT is a ring homomorphism: map each share independently first. */
	for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
		for (size_t i = 0; i < n; i ++) {
			uint32_t aq = mq_decode_share(
				a_q[i].share[sh], FNDSA_MASK_Q, MASK_Q_Q0I);
			uint32_t aq2 = mq_decode_share(
				a_q2[i].share[sh], FNDSA_MASK_Q2, MASK_Q2_Q0I);
			uint64_t z = (uint64_t)aq2 + FNDSA_MASK_Q2 - aq;
			uint64_t ge = (uint64_t)0
				- (uint64_t)(z >= FNDSA_MASK_Q2);
			uint32_t w = (uint32_t)(z
				- ((uint64_t)FNDSA_MASK_Q2 & ge));
			uint32_t u = fndsa_mask_reduce_u64(
				(uint64_t)w * UINT32_C(3), FNDSA_MASK_Q2);
			d[i].share[sh] = aq + FNDSA_MASK_Q * u;
		}
	}

	/* Centering is nonlinear: compare the reconstructed shared value once. */
	for (size_t i = 0; i < n; i ++) {
		mi32 raw = d[i];
		mi32_center_mod_m(&d[i], &raw, rng);
		raw.share[0] = 0;
		raw.share[1] = 0;
	}
}
