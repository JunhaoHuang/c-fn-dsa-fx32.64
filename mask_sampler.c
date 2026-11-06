/*
 * Standalone masked fixed-point sampler gadgets.
 *
 * This is a clean-room fixed-point design inspired by the Boolean data-flow
 * of masked_FALCON.  No IEEE-754 packing code is used or copied.
 */
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "mask_sampler.h"

int32_t fndsa_masked_gaussian0_helper(uint64_t lo, uint32_t hi);

#if defined(__GNUC__) || defined(__clang__)
#define MASK_SAMPLER_NOINLINE   __attribute__((noinline))
#else
#define MASK_SAMPLER_NOINLINE
#endif

static void
sampler_clear(void *buf, size_t len)
{
	volatile uint8_t *d = (volatile uint8_t *)buf;
	while (len -- > 0) {
		*d ++ = 0;
	}
}


static void
/* Put a public Q32.64 constant in share 0; share 1 is zero. */
mrx_public_raw(mrx160 *d, const uint32_t x[3])
{
	memset(d, 0, sizeof *d);
	d->share[0][0] = x[0];
	d->share[0][1] = x[1];
	d->share[0][2] = x[2];
	d->share[0][3] = d->share[0][4]
		= (uint32_t)(0 - (x[2] >> 31));
}

static void
mrx_public_u64(mrx160 *d, uint64_t x)
{
	uint32_t t[3] = { (uint32_t)x, (uint32_t)(x >> 32), 0 };
	mrx_public_raw(d, t);
}

static void
mb64_to_mrx(mrx160 *d, const mb64 *a, fndsa_mask_rng *rng)
{
	mb160 t;
	memset(&t, 0, sizeof t);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		t.share[j][0] = a->share[j][0];
		t.share[j][1] = a->share[j][1];
	}
	mrx160_b2a(d, &t, rng);
	sampler_clear(&t, sizeof t);
}

static void
mrx_to_mb64(mb64 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mb160 t;
	mrx160_a2b(&t, a, rng);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j][0] = t.share[j][0];
		d->share[j][1] = t.share[j][1];
	}
	sampler_clear(&t, sizeof t);
}

static void
mask_and_word(uint32_t *z0, uint32_t *z1,
	uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1,
	fndsa_mask_rng *rng)
{
	fndsa_mask_bool_and_words(
		z0, z1, &x0, &x1, &y0, &y1, 1, rng);
}

static void
mbit_set_public(mbit *d, uint32_t x, fndsa_mask_rng *rng)
{
	uint32_t r = fndsa_mask_rng_next_bit(rng);
	d->share[0] = r;
	d->share[1] = r ^ (x & 1);
}

static void
mbit_refresh(mbit *d, const mbit *a, fndsa_mask_rng *rng)
{
	uint32_t z0, z1;
	fndsa_mask_bool_refresh_words(&z0, &z1,
		&a->share[0], &a->share[1], 1, rng);
	d->share[0] = z0 & 1;
	d->share[1] = z1 & 1;
}

static void
mbit_not(mbit *d, const mbit *a)
{
	d->share[0] = (a->share[0] ^ 1) & 1;
	d->share[1] = a->share[1] & 1;
}

static void
mbit_and(mbit *d, const mbit *a, const mbit *b,
	fndsa_mask_rng *rng)
{
	uint32_t z0, z1;
	mask_and_word(&z0, &z1,
		a->share[0], a->share[1], b->share[0], b->share[1], rng);
	d->share[0] = z0 & 1;
	d->share[1] = z1 & 1;
}

static void
mbit_or(mbit *d, const mbit *a, const mbit *b,
	fndsa_mask_rng *rng)
{
	mbit t, o;
	mbit_and(&t, a, b, rng);
	o.share[0] = (a->share[0] ^ b->share[0] ^ t.share[0]) & 1;
	o.share[1] = (a->share[1] ^ b->share[1] ^ t.share[1]) & 1;
	mbit_refresh(d, &o, rng);
	sampler_clear(&t, sizeof t);
	sampler_clear(&o, sizeof o);
}

static void
mbit_select(mbit *d, const mbit *a, const mbit *b,
	const mbit *ctl, fndsa_mask_rng *rng)
{
	mbit x, t, o;
	x.share[0] = a->share[0] ^ b->share[0];
	x.share[1] = a->share[1] ^ b->share[1];
	mbit_and(&t, &x, ctl, rng);
	o.share[0] = (a->share[0] ^ t.share[0]) & 1;
	o.share[1] = (a->share[1] ^ t.share[1]) & 1;
	mbit_refresh(d, &o, rng);
	sampler_clear(&x, sizeof x);
	sampler_clear(&t, sizeof t);
	sampler_clear(&o, sizeof o);
}

static void
mbit_nonzero32(mbit *d, uint32_t x0, uint32_t x1,
	fndsa_mask_rng *rng)
{
	static const unsigned SH[] = { 16, 8, 4, 2, 1 };
	for (size_t i = 0; i < sizeof SH / sizeof SH[0]; i ++) {
		uint32_t y0 = x0 >> SH[i];
		uint32_t y1 = x1 >> SH[i];
		uint32_t z0, z1;
		mask_and_word(&z0, &z1, x0, x1, y0, y1, rng);
		x0 ^= y0 ^ z0;
		x1 ^= y1 ^ z1;
	}
	mbit t = { { x0 & 1, x1 & 1 } };
	mbit_refresh(d, &t, rng);
	sampler_clear(&t, sizeof t);
}

static void
mb64_rshift_public(mb64 *d, const mb64 *a, unsigned n)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		uint64_t x = (uint64_t)a->share[j][0]
			| ((uint64_t)a->share[j][1] << 32);
		x >>= n;
		d->share[j][0] = (uint32_t)x;
		d->share[j][1] = (uint32_t)(x >> 32);
	}
}

static void
mb64_refresh(mb64 *d, const mb64 *a, fndsa_mask_rng *rng)
{
	fndsa_mask_bool_refresh_words(d->share[0], d->share[1],
		a->share[0], a->share[1], 2, rng);
}

static void
mb64_select(mb64 *d, const mb64 *a, const mb64 *b,
	const mbit *ctl, fndsa_mask_rng *rng)
{
	mb64 t, o;
	uint32_t c0 = (uint32_t)0 - (ctl->share[0] & 1);
	uint32_t c1 = (uint32_t)0 - (ctl->share[1] & 1);
	for (size_t i = 0; i < 2; i ++) {
		uint32_t x0 = a->share[0][i] ^ b->share[0][i];
		uint32_t x1 = a->share[1][i] ^ b->share[1][i];
		mask_and_word(&t.share[0][i], &t.share[1][i],
			x0, x1, c0, c1, rng);
		o.share[0][i] = a->share[0][i] ^ t.share[0][i];
		o.share[1][i] = a->share[1][i] ^ t.share[1][i];
	}
	mb64_refresh(d, &o, rng);
	sampler_clear(&t, sizeof t);
	sampler_clear(&o, sizeof o);
}

static void
/* Six masked-select stages implement a secret shift in the range 0..63. */
mb64_secret_rshift(mb64 *d, const mb64 *a, const mbit shift[6],
	fndsa_mask_rng *rng)
{
	mb64 x = *a;
	for (unsigned i = 0; i < 6; i ++) {
		mb64 y, z;
		mb64_rshift_public(&y, &x, 1u << i);
		mb64_select(&z, &x, &y, &shift[i], rng);
		x = z;
		sampler_clear(&y, sizeof y);
		sampler_clear(&z, sizeof z);
	}
	*d = x;
	sampler_clear(&x, sizeof x);
}

static void
/* Compare one secret threshold byte with public w, most-significant bit first. */
mb64_compare_byte(mbit *eq, mbit *lt, const mb64 *a,
	unsigned byte_index, uint8_t w, fndsa_mask_rng *rng)
{
	mbit e, l;
	mbit_set_public(&e, 1, rng);
	mbit_set_public(&l, 0, rng);
	unsigned word = byte_index >> 2;
	unsigned off = (byte_index & 3) << 3;
	for (int i = 7; i >= 0; i --) {
		mbit b, diff, ndiff, take, ne, nl;
		b.share[0] = (a->share[0][word] >> (off + (unsigned)i)) & 1;
		b.share[1] = (a->share[1][word] >> (off + (unsigned)i)) & 1;
		diff = b;
		diff.share[0] ^= (w >> i) & 1;
		mbit_not(&ndiff, &diff);
		mbit_and(&take, &e, &diff, rng);
		if (((w >> i) & 1) == 0) {
			mbit_or(&nl, &l, &take, rng);
			l = nl;
		}
		mbit_and(&ne, &e, &ndiff, rng);
		e = ne;
		sampler_clear(&b, sizeof b);
		sampler_clear(&diff, sizeof diff);
		sampler_clear(&ndiff, sizeof ndiff);
		sampler_clear(&take, sizeof take);
		sampler_clear(&ne, sizeof ne);
		sampler_clear(&nl, sizeof nl);
	}
	mbit_refresh(eq, &e, rng);
	mbit_refresh(lt, &l, rng);
	sampler_clear(&e, sizeof e);
	sampler_clear(&l, sizeof l);
}

static void
mrx_half_positive(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mb160 x, y;
	mrx160_a2b(&x, a, rng);
	memset(&y, 0, sizeof y);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			y.share[j][i] = x.share[j][i] >> 1;
			if (i + 1 < FNDSA_MASK160_LIMBS) {
				y.share[j][i] |= x.share[j][i + 1] << 31;
			}
		}
	}
	mrx160_b2a(d, &y, rng);
	mrx160_refresh(d, rng);
	sampler_clear(&x, sizeof x);
	sampler_clear(&y, sizeof y);
}

MASK_SAMPLER_NOINLINE void
/*
 * Horner-evaluate the legacy 13-coefficient FACCT polynomial. Each multiply
 * is an unrounded unsigned high-half product, matching the original order.
 */
mrx160_expm_p63(mb64 *d, const mrx160 *x, const mrx160 *ccs,
	fndsa_mask_rng *rng)
{
	static const uint64_t C[] = {
		UINT64_C(0x00000004741183A3),
		UINT64_C(0x00000036548CFC06),
		UINT64_C(0x0000024FDCBF140A),
		UINT64_C(0x0000171D939DE045),
		UINT64_C(0x0000D00CF58F6F84),
		UINT64_C(0x000680681CF796E3),
		UINT64_C(0x002D82D8305B0FEA),
		UINT64_C(0x011111110E066FD0),
		UINT64_C(0x0555555555070F00),
		UINT64_C(0x155555555581FF00),
		UINT64_C(0x400000000002B400),
		UINT64_C(0x7FFFFFFFFFFF4800),
		UINT64_C(0x8000000000000000)
	};
	mrx160 y, p, c;
	/* All operands are unsigned Q0.64 words; every exact product fits
	   128 bits, so mrx160_mul_u64_hi() cannot wrap its 160-bit ring. */
	mb64 bp;
	mrx_public_u64(&y, C[0]);
	for (size_t i = 1; i < sizeof C / sizeof C[0]; i ++) {
		mrx160_mul_u64_hi(&bp, &y, x, rng);
		mb64_to_mrx(&p, &bp, rng);
		mrx_public_u64(&c, C[i]);
		mrx160_sub(&y, &c, &p);
	}
	mrx160_mul_u64_hi(d, &y, ccs, rng);
	sampler_clear(&y, sizeof y);
	sampler_clear(&p, sizeof p);
	sampler_clear(&c, sizeof c);
	sampler_clear(&bp, sizeof bp);
}

MASK_SAMPLER_NOINLINE void
/*
 * Reduce x modulo log(2), build a secret 64-bit threshold, saturate the
 * secret shift to 63, then compare signature-random bytes lazily.
 */
mrx160_ber_exp(mbit *accept, fndsa_mask_sampler_source *src,
	const mrx160 *x, const mrx160 *ccs, fndsa_mask_rng *rng)
{
	static const uint32_t INV_LOG2[3] = {
		UINT32_C(3090145655), UINT32_C(1901360722), 1
	};
	static const uint32_t LOG2[3] = {
		UINT32_C(3520035244), UINT32_C(2977044471), 0
	};
	mrx160 tx, s, t, r, yz, twoz, ulp;
	mb160 sb;
	mb64 y, z, zs;
	mbit ge, sh[6];

	if (src == NULL || src->next_u8 == NULL) {
		mbit_set_public(accept, 0, rng);
		return;
	}
	/* BerExp is called with x in the sampler rejection range; all
	   reduction products are below 2^136, hence fit signed 160 bits. */
	mrx160_mul_public_fast(&tx, x, INV_LOG2, rng);
	mrx160_floor(&s, &tx, rng);
	mrx160_mul_public_fast(&t, &s, LOG2, rng);
	mrx160_sub(&r, x, &t);
	mrx160_expm_p63(&y, &r, ccs, rng);

	/* z <- (2*y - 1) mod 2^64. */
	mb64_to_mrx(&yz, &y, rng);
	mrx160_add(&twoz, &yz, &yz);
	mrx_public_u64(&ulp, 1);
	mrx160_sub(&twoz, &twoz, &ulp);
	mrx_to_mb64(&z, &twoz, rng);

	/* Saturate the secret shift to 63, exactly as the legacy code does. */
	mrx160_a2b(&sb, &s, rng);
	mbit_nonzero32(&ge, sb.share[0][2] & UINT32_C(0xFFFFFFC0),
		sb.share[1][2] & UINT32_C(0xFFFFFFC0), rng);
	for (unsigned i = 0; i < 6; i ++) {
		mbit bit, one;
		bit.share[0] = (sb.share[0][2] >> i) & 1;
		bit.share[1] = (sb.share[1][2] >> i) & 1;
		mbit_set_public(&one, 1, rng);
		mbit_select(&sh[i], &bit, &one, &ge, rng);
		sampler_clear(&bit, sizeof bit);
		sampler_clear(&one, sizeof one);
	}
	mb64_secret_rshift(&zs, &z, sh, rng);

	/* Legacy-compatible lazy comparison, most-significant byte first. */
	for (int i = 7; i >= 0; i --) {
		uint8_t w = src->next_u8(src->ctx);
		mbit eq, lt;
		mb64_compare_byte(&eq, &lt, &zs, (unsigned)i, w, rng);
		uint32_t same = mbit_open_compat(&eq);
		if (same == 0) {
			*accept = lt;
			sampler_clear(&lt, sizeof lt);
			sampler_clear(&eq, sizeof eq);
			goto cleanup;
		}
		sampler_clear(&eq, sizeof eq);
		sampler_clear(&lt, sizeof lt);
	}
	mbit_set_public(accept, 0, rng);

cleanup:
	sampler_clear(&tx, sizeof tx);
	sampler_clear(&s, sizeof s);
	sampler_clear(&t, sizeof t);
	sampler_clear(&r, sizeof r);
	sampler_clear(&yz, sizeof yz);
	sampler_clear(&twoz, sizeof twoz);
	sampler_clear(&ulp, sizeof ulp);
	sampler_clear(&sb, sizeof sb);
	sampler_clear(&y, sizeof y);
	sampler_clear(&z, sizeof z);
	sampler_clear(&zs, sizeof zs);
	sampler_clear(&ge, sizeof ge);
	sampler_clear(sh, sizeof sh);
}

static uint64_t
/* Assemble signature-PRNG bytes in the legacy little-endian order. */
source_next_u64(fndsa_mask_sampler_source *src)
{
	uint64_t x = 0;
	for (unsigned i = 0; i < 8; i ++) {
		x |= (uint64_t)src->next_u8(src->ctx) << (8 * i);
	}
	return x;
}

static uint32_t
source_next_u16(fndsa_mask_sampler_source *src)
{
	uint32_t x = src->next_u8(src->ctx);
	return x | ((uint32_t)src->next_u8(src->ctx) << 8);
}

MASK_SAMPLER_NOINLINE void
/*
 * Masked SamplerZ: split mu, derive dss/ccs, draw BaseSampler candidates,
 * and apply BerExp rejection without opening candidate values.
 */
mrx160_sampler_z(mrx160 *zout, fndsa_mask_sampler_source *src,
	const mrx160 *mu, const mrx160 *isigma, fndsa_mask_rng *rng)
{
	static const uint32_t INV_2SQRSIGMA0[3] = {
		UINT32_C(125273110), UINT32_C(647960451), 0
	};
	static const uint32_t SIGMA_MIN_9[3] = {
		UINT32_C(4154935963), UINT32_C(1193286641), 1
	};
	static const uint32_t SIGMA_MIN_10[3] = {
		UINT32_C(209075216), UINT32_C(1281104281), 1
	};
	mrx160 s, r, dss, ccs;

	if (src == NULL || src->next_u8 == NULL
		|| (src->logn != 9 && src->logn != 10))
	{
		memset(zout, 0, sizeof *zout);
		return;
	}
	/* Decompose the center once: mu = floor(mu) + fractional remainder. */
	mrx160_floor(&s, mu, rng);
	mrx160_sub(&r, mu, &s);
	/* Production isigma is in (0,1); these setup products are < 2^128. */
	mrx160_sqr(&dss, isigma, rng);
	mrx_half_positive(&dss, &dss, rng);
	mrx160_mul_public_fast(&ccs, isigma,
		src->logn == 9 ? SIGMA_MIN_9 : SIGMA_MIN_10, rng);

	/*
	 * Candidate bytes come only from src; gadget masks come only from rng.
	 * This separation preserves the legacy signature-random byte sequence.
	 */
	for (;;) {
		uint64_t lo = source_next_u64(src);
		uint32_t hi = source_next_u16(src);
		int32_t z0 = fndsa_masked_gaussian0_helper(lo, hi);
		int32_t b = (int32_t)lo & 1;
		int32_t zi = b + ((b << 1) - 1) * z0;
		uint32_t zr[3] = { 0, 0, (uint32_t)zi };
		uint32_t z0r[3] = { 0, 0, (uint32_t)(z0 * z0) };
		mrx160 zm, x, y, t;
		mbit accept;
		mrx160_share_rx(&zm, zr, rng);
		mrx160_sub(&x, &zm, &r);
		/* BaseSampler candidates are bounded by its 79-bit table.  For
		   FN-DSA-512/1024 the squared/reweighted rejection terms remain
		   far below signed 2^159 before the Q32.64 shift. */
		mrx160_sqr(&x, &x, rng);
		mrx160_mul(&x, &x, &dss, rng);
		mrx160_share_rx(&y, z0r, rng);
		sampler_clear(&lo, sizeof lo);
		sampler_clear(&hi, sizeof hi);
		sampler_clear(&z0, sizeof z0);
		sampler_clear(&b, sizeof b);
		sampler_clear(&zi, sizeof zi);
		sampler_clear(zr, sizeof zr);
		sampler_clear(z0r, sizeof z0r);
		mrx160_mul_public_fast(&t, &y, INV_2SQRSIGMA0, rng);
		mrx160_sub(&x, &x, &t);
		mrx160_ber_exp(&accept, src, &x, &ccs, rng);
		/* Observable rejection-loop boundary retained for compatibility. */
		uint32_t ok = mbit_open_compat(&accept);
		if (ok != 0) {
			mrx160_add(zout, &s, &zm);
			mrx160_refresh(zout, rng);
			sampler_clear(&zm, sizeof zm);
			sampler_clear(&x, sizeof x);
			sampler_clear(&y, sizeof y);
			sampler_clear(&t, sizeof t);
			sampler_clear(&accept, sizeof accept);
			break;
		}
		sampler_clear(&zm, sizeof zm);
		sampler_clear(&x, sizeof x);
		sampler_clear(&y, sizeof y);
		sampler_clear(&t, sizeof t);
		sampler_clear(&accept, sizeof accept);
	}

	sampler_clear(&s, sizeof s);
	sampler_clear(&r, sizeof r);
	sampler_clear(&dss, sizeof dss);
	sampler_clear(&ccs, sizeof ccs);
}
