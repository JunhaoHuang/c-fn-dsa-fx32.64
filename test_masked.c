/*
 * Tests for the experimental 160-bit masking reference backend.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "fndsa.h"
#include "fndsa_masked.h"
#include "mask160.h"
#include "mask_mont32.h"
#include "mask_mq.h"
#include "mask_sampler.h"
#include "mask_xpoly.h"
#include "mask_tree.h"
#include "sign_inner.h"

extern const rx fndsa_X_PLUS_1_SQUARED_OVER_X[];
extern const rx fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[];

static void
fail(const char *what)
{
	fprintf(stderr, "test_masked: %s failed\n", what);
	exit(EXIT_FAILURE);
}

static uint64_t
prng64(uint64_t *s)
{
	uint64_t x = *s;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	*s = x;
	return x;
}

static int
eq_rx(const uint32_t a[3], const uint32_t b[3])
{
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

static unsigned
raw96_abs_bits(const uint32_t a[3])
{
	uint32_t m[3];
	uint32_t neg = a[2] >> 31;
	uint32_t mask = (uint32_t)0 - neg;
	uint64_t cc = neg;
	for (size_t i = 0; i < 3; i ++) {
		uint64_t z = (uint64_t)(a[i] ^ mask) + cc;
		m[i] = (uint32_t)z;
		cc = z >> 32;
	}
	for (size_t i = 3; i > 0; i --) {
		if (m[i - 1] != 0) {
			return (unsigned)((i - 1) * 32 + 32
				- (size_t)__builtin_clz(m[i - 1]));
		}
	}
	return 0;
}

static int
all_zero(const uint8_t *buf, size_t len)
{
	uint32_t x = 0;
	for (size_t i = 0; i < len; i ++) {
		x |= buf[i];
	}
	return x == 0;
}

static void
test_lfsr(void)
{
	static const uint64_t KAT[8] = {
		UINT64_C(0x1E3C5A7896B4D2F1),
		UINT64_C(0x3355FF98AACC6602),
		UINT64_C(0x5AD34BC078F169E6),
		UINT64_C(0xD30D68B1A47A1FC9),
		UINT64_C(0x13BC46E3B916EC5F),
		UINT64_C(0x81625CA43AD9E72D),
		UINT64_C(0x25BC348F079E16E5),
		UINT64_C(0x49BCD0567A8FE390)
	};
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	for (size_t i = 0; i < sizeof KAT / sizeof KAT[0]; i ++) {
		if (fndsa_mask_rng_next_u64(&rng) != KAT[i]) {
			fail("LFSR KAT");
		}
	}

	/*
	 * Mixed-width requests must consume one continuous stream without
	 * discarding or reusing bits.
	 */
	fndsa_mask_rng ref, mix;
	fndsa_mask_rng_init_demo(&ref);
	fndsa_mask_rng_init_demo(&mix);
	uint64_t w0 = fndsa_mask_rng_next_u64(&ref);
	uint64_t w1 = fndsa_mask_rng_next_u64(&ref);
	uint64_t w2 = fndsa_mask_rng_next_u64(&ref);
	for (unsigned i = 0; i < 7; i ++) {
		if (fndsa_mask_rng_next_bit(&mix) != ((w0 >> i) & 1)) {
			fail("LFSR bit reservoir");
		}
	}
	if (fndsa_mask_rng_next_u32(&mix) != (uint32_t)(w0 >> 7)
		|| fndsa_mask_rng_next_u64(&mix)
			!= ((w0 >> 39) | (w1 << 25))
		|| fndsa_mask_rng_next_bit(&mix) != ((w1 >> 39) & 1)
		|| fndsa_mask_rng_next_u32(&mix)
			!= ((uint32_t)(w1 >> 40) | (uint32_t)(w2 << 24)))
	{
		fail("LFSR mixed-width reservoir");
	}
}

/*
 * Structural regression for the single Boolean nonlinear boundary. This checks
 * correctness with independent, correlated and aliased input sharings; it is
 * not a statistical leakage test or a substitute for the SNI argument.
 */
static void
test_mask_bool_nonlinear(void)
{
	fndsa_mask_rng rng;
	uint64_t rs = UINT64_C(0xB713D4E90A5C862F);
	uint32_t changed = 0;
	fndsa_mask_rng_init_demo(&rng);
	for (unsigned k = 0; k < 2000; k ++) {
		uint32_t x0[5], x1[5], y0[5], y1[5];
		uint32_t z0[5], z1[5], q0[5], q1[5];
		uint32_t io0[5], io1[5];
		for (size_t i = 0; i < 5; i ++) {
			uint32_t x = (uint32_t)prng64(&rs);
			uint32_t y = (uint32_t)prng64(&rs);
			x0[i] = (uint32_t)prng64(&rs);
			x1[i] = x ^ x0[i];
			y0[i] = (uint32_t)prng64(&rs);
			y1[i] = y ^ y0[i];
		}
		fndsa_mask_bool_and_words(
			z0, z1, x0, x1, y0, y1, 5, &rng);
		for (size_t i = 0; i < 5; i ++) {
			if ((z0[i] ^ z1[i])
				!= ((x0[i] ^ x1[i]) & (y0[i] ^ y1[i])))
			{
				fail("masked Boolean AND differential");
			}
		}
		fndsa_mask_bool_refresh_words(q0, q1, z0, z1, 5, &rng);
		for (size_t i = 0; i < 5; i ++) {
			if ((q0[i] ^ q1[i]) != (z0[i] ^ z1[i])) {
				fail("masked Boolean refresh differential");
			}
			changed |= q0[i] ^ z0[i];
			io0[i] = x0[i];
			io1[i] = x1[i];
		}
		fndsa_mask_bool_and_words(
			io0, io1, io0, io1, io0, io1, 5, &rng);
		for (size_t i = 0; i < 5; i ++) {
			if ((io0[i] ^ io1[i]) != (x0[i] ^ x1[i])) {
				fail("masked Boolean AND aliased sharing");
			}
		}
		for (size_t i = 0; i < 5; i ++) {
			y0[i] = x0[i] >> 1;
			y1[i] = x1[i] >> 1;
		}
		fndsa_mask_bool_and_words(
			z0, z1, x0, x1, y0, y1, 5, &rng);
		for (size_t i = 0; i < 5; i ++) {
			uint32_t x = x0[i] ^ x1[i];
			if ((z0[i] ^ z1[i]) != (x & (x >> 1))) {
				fail("masked Boolean AND correlated fan-out");
			}
		}
	}
	if (changed == 0) {
		fail("masked Boolean refresh changes shares");
	}
}
/*
 * Structural regression for the conversion boundary.  This is not a leakage
 * proof; it catches the former A2B shortcut whose second Boolean output share
 * had a constant zero least-significant bit (and therefore put the secret LSB
 * directly in the other output share).
 */
static void
test_conversion_share_distribution(void)
{
	static const uint32_t AX[3] = {
		UINT32_C(0xC7A13E59), UINT32_C(0x6B2084FD),
		UINT32_C(0x12345678)
	};
	static const uint32_t BV[FNDSA_MASK160_LIMBS] = {
		UINT32_C(0x91D43A6F), UINT32_C(0x02BE7851),
		UINT32_C(0xF0C35A97), UINT32_C(0x4E1268BD),
		UINT32_C(0xA5793CE0)
	};
	fndsa_mask_rng share_rng, rng;
	mrx160 a, ar;
	mb160 b;
	uint32_t a_one[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS] = { { 0 } };
	uint32_t a_zero[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS] = { { 0 } };
	uint32_t b_one[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS] = { { 0 } };
	uint32_t b_zero[FNDSA_MASK_SHARES][FNDSA_MASK160_LIMBS] = { { 0 } };
	uint64_t rs = UINT64_C(0x69EF4B82C1357AD0);

	fndsa_mask_rng_init_demo(&share_rng);
	mrx160_share_rx(&a, AX, &share_rng);
	fndsa_mask_rng_init_demo(&rng);
	for (unsigned k = 0; k < 512; k ++) {
		mrx160_a2b(&b, &a, &rng);
		for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
			for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
				a_one[j][i] |= b.share[j][i];
				a_zero[j][i] |= ~b.share[j][i];
				uint32_t ex = i < 3 ? AX[i] : 0;
				if ((b.share[0][i] ^ b.share[1][i]) != ex) {
					fail("A2B share-distribution correctness");
				}
			}
		}
	}
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			if (a_one[j][i] != UINT32_MAX
				|| a_zero[j][i] != UINT32_MAX)
			{
				fail("A2B output share not fully randomized");
			}
		}
	}

	for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
		b.share[0][i] = (uint32_t)prng64(&rs);
		b.share[1][i] = b.share[0][i] ^ BV[i];
	}
	fndsa_mask_rng_init_demo(&rng);
	for (unsigned k = 0; k < 512; k ++) {
		uint64_t cc = 0;
		mrx160_b2a(&ar, &b, &rng);
		for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
			for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
				b_one[j][i] |= ar.share[j][i];
				b_zero[j][i] |= ~ar.share[j][i];
			}
		}
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			uint64_t z = (uint64_t)ar.share[0][i]
				+ ar.share[1][i] + cc;
			if ((uint32_t)z != BV[i]) {
				fail("B2A share-distribution correctness");
			}
			cc = z >> 32;
		}
	}
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			if (b_one[j][i] != UINT32_MAX
				|| b_zero[j][i] != UINT32_MAX)
			{
				fail("B2A output share not fully randomized");
			}
		}
	}
}

static void
test_mask160_edges(fndsa_mask_rng *rng)
{
	static const rx edge[] = {
		{ 0, 0, 0 },
		{ 1, 0, 0 },
		{ UINT32_MAX, UINT32_MAX, UINT32_MAX },
		{ 0, UINT32_C(0x80000000), 0 },
		{ 0, UINT32_C(0x80000000), UINT32_MAX },
		{ 0, 0, 1 },
		{ 0, 0, UINT32_MAX },
		{ UINT32_MAX, UINT32_MAX, UINT32_C(0x7FFFFFFF) },
		{ 0, 0, UINT32_C(0x80000000) }
	};
	for (size_t i = 0; i < sizeof edge / sizeof edge[0]; i ++) {
		mrx160 a, t;
		mb160 b;
		rx got;
		mrx160_share_rx(&a, edge[i], rng);
		mrx160_open_rx(got, &a);
		if (!eq_rx(got, edge[i])) {
			fail("edge share/open");
		}
		mrx160_a2b(&b, &a, rng);
		mrx160_b2a(&t, &b, rng);
		mrx160_open_rx(got, &t);
		if (!eq_rx(got, edge[i])) {
			fail("edge A2B/B2A");
		}
	}
	for (size_t i = 0; i < sizeof edge / sizeof edge[0]; i ++) {
		for (size_t j = 0; j < sizeof edge / sizeof edge[0]; j ++) {
			mrx160 a, b, c;
			rx ref, got;
			mrx160_share_rx(&a, edge[i], rng);
			mrx160_share_rx(&b, edge[j], rng);
			mrx160_mul_wide(&c, &a, &b, rng);
			rx_mul(ref, edge[i], edge[j]);
			mrx160_open_rx(got, &c);
			if (!eq_rx(got, ref)) {
				fail("192-bit Q32.64 multiply fallback");
			}
			mb160 full;
			mrx160_a2b(&full, &c, rng);
			uint32_t sw = (uint32_t)(0 - (ref[2] >> 31));
			for (size_t k = 0; k < FNDSA_MASK160_LIMBS; k ++) {
				uint32_t x = full.share[0][k] ^ full.share[1][k];
				uint32_t ex = k < 3 ? ref[k] : sw;
				if (x != ex) {
					fail("Q32.64 sign extension");
				}
			}
		}
	}
}

static void
test_public_raw(mrx160 *d, const uint32_t x[3])
{
	uint32_t sw = (uint32_t)(0 - (x[2] >> 31));
	memset(d, 0, sizeof *d);
	d->share[0][0] = x[0];
	d->share[0][1] = x[1];
	d->share[0][2] = x[2];
	d->share[0][3] = sw;
	d->share[0][4] = sw;
}

static void
test_public_constant_mul(void)
{
	static const rx constants[] = {
		{ 0, 0, 0 },
		{ 1, 0, 0 },
		{ UINT32_MAX, UINT32_MAX, UINT32_MAX },
		{ 0, UINT32_C(0x80000000), 0 },
		{ 0, UINT32_C(0x80000000), UINT32_MAX },
		{ UINT32_MAX, UINT32_MAX, 0 },
		{ 0, 0, 1 },
		{ 0, 0, UINT32_MAX }
	};
	fndsa_mask_rng setup;
	uint64_t rs = UINT64_C(0xD91B5E7A24C6038F);
	fndsa_mask_rng_init_demo(&setup);

	for (unsigned k = 0; k < 1000; k ++) {
		rx a, ref, got0, got1;
		a[0] = (uint32_t)prng64(&rs);
		a[1] = (uint32_t)prng64(&rs);
		a[2] = (uint32_t)((int32_t)(prng64(&rs) % 63) - 31);
		const uint32_t *c = constants[k
			% (sizeof constants / sizeof constants[0])];
		mrx160 ma, mc, dg, dp, alias;
		fndsa_mask_rng rg, rp;
		mrx160_share_rx(&ma, a, &setup);
		test_public_raw(&mc, c);
		rg = setup;
		rp = setup;
		mrx160_mul(&dg, &ma, &mc, &rg);
		mrx160_mul_public_fast(&dp, &ma, c, &rp);
		rx_mul(ref, a, c);
		mrx160_open_rx(got0, &dg);
		mrx160_open_rx(got1, &dp);
		if (!eq_rx(got0, ref) || !eq_rx(got1, ref)) {
			fail("public constant fast differential");
		}
		if (memcmp(&rg, &rp, sizeof rg) != 0) {
			fail("public constant fast RNG budget");
		}
		alias = ma;
		rp = setup;
		mrx160_mul_public_fast(&alias, &alias, c, &rp);
		mrx160_open_rx(got1, &alias);
		if (!eq_rx(got1, ref)) {
			fail("public constant fast alias");
		}
	}

	for (unsigned k = 0; k < 1000; k ++) {
		rx a, c, ref, got0, got1;
		for (size_t i = 0; i < 3; i ++) {
			a[i] = (uint32_t)prng64(&rs);
			c[i] = (uint32_t)prng64(&rs);
		}
		mrx160 ma, mc, dg, dp, alias;
		fndsa_mask_rng rg, rp;
		mrx160_share_rx(&ma, a, &setup);
		test_public_raw(&mc, c);
		rg = setup;
		rp = setup;
		mrx160_mul_wide(&dg, &ma, &mc, &rg);
		mrx160_mul_public_wide(&dp, &ma, c, &rp);
		rx_mul(ref, a, c);
		mrx160_open_rx(got0, &dg);
		mrx160_open_rx(got1, &dp);
		if (!eq_rx(got0, ref) || !eq_rx(got1, ref)) {
			fail("public constant wide differential");
		}
		if (memcmp(&rg, &rp, sizeof rg) != 0) {
			fail("public constant wide RNG budget");
		}
		alias = ma;
		rp = setup;
		mrx160_mul_public_wide(&alias, &alias, c, &rp);
		mrx160_open_rx(got1, &alias);
		if (!eq_rx(got1, ref)) {
			fail("public constant wide alias");
		}
	}

	/* Exact half ties: round after adding 2^63, hence toward +infinity. */
	{
		static const rx ulp = { 1, 0, 0 };
		static const rx pos_half = { 0, UINT32_C(0x80000000), 0 };
		static const rx neg_half = {
			0, UINT32_C(0x80000000), UINT32_MAX
		};
		mrx160 a, d;
		rx got;
		mrx160_share_rx(&a, pos_half, &setup);
		mrx160_mul_public_fast(&d, &a, ulp, &setup);
		mrx160_open_rx(got, &d);
		if (got[0] != 1 || got[1] != 0 || got[2] != 0) {
			fail("public constant positive half tie");
		}
		mrx160_share_rx(&a, neg_half, &setup);
		mrx160_mul_public_fast(&d, &a, ulp, &setup);
		mrx160_open_rx(got, &d);
		if (got[0] != 0 || got[1] != 0 || got[2] != 0) {
			fail("public constant negative half tie");
		}
	}
}

static void
test_mask160_arithmetic(void)
{
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	uint64_t rs = UINT64_C(0xA7C318D29B4605EF);
	test_conversion_share_distribution();
	test_mask160_edges(&rng);
	test_public_constant_mul();

	/* Exercise all five limbs, including carries across every word boundary. */
	for (unsigned k = 0; k < 2000; k ++) {
		mrx160 ma, mt;
		mb160 ba;
		uint32_t refa[FNDSA_MASK160_LIMBS];
		uint32_t refb[FNDSA_MASK160_LIMBS];
		uint64_t cc = 0;
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			ma.share[0][i] = (uint32_t)prng64(&rs);
			ma.share[1][i] = (uint32_t)prng64(&rs);
			uint64_t z = (uint64_t)ma.share[0][i]
				+ ma.share[1][i] + cc;
			refa[i] = (uint32_t)z;
			cc = z >> 32;
		}
		mrx160_a2b(&ba, &ma, &rng);
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			if ((ba.share[0][i] ^ ba.share[1][i]) != refa[i]) {
				fail("full-width A2B");
			}
			ba.share[0][i] = (uint32_t)prng64(&rs);
			ba.share[1][i] = (uint32_t)prng64(&rs);
			refb[i] = ba.share[0][i] ^ ba.share[1][i];
		}
		mrx160_b2a(&mt, &ba, &rng);
		cc = 0;
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			uint64_t z = (uint64_t)mt.share[0][i]
				+ mt.share[1][i] + cc;
			if ((uint32_t)z != refb[i]) {
				fail("full-width B2A");
			}
			cc = z >> 32;
		}
	}

	for (unsigned k = 0; k < 400; k ++) {
		rx a, b, ref, got;
		a[0] = (uint32_t)prng64(&rs);
		a[1] = (uint32_t)prng64(&rs);
		a[2] = (uint32_t)((int32_t)(prng64(&rs) % 63) - 31);
		b[0] = (uint32_t)prng64(&rs);
		b[1] = (uint32_t)prng64(&rs);
		b[2] = (uint32_t)((int32_t)(prng64(&rs) % 63) - 31);

		mrx160 ma, mb, mc, mt;
		mb160 ba;
		mrx160_share_rx(&ma, a, &rng);
		mrx160_share_rx(&mb, b, &rng);
		mrx160_open_rx(got, &ma);
		if (!eq_rx(a, got)) {
			fail("share/open");
		}

		mrx160_add(&mc, &ma, &mb);
		rx_add(ref, a, b);
		mrx160_open_rx(got, &mc);
		if (!eq_rx(ref, got)) {
			fail("add");
		}

		mrx160_sub(&mc, &ma, &mb);
		rx_sub(ref, a, b);
		mrx160_open_rx(got, &mc);
		if (!eq_rx(ref, got)) {
			fail("sub");
		}

		mrx160_neg(&mc, &ma);
		rx_neg(ref, a);
		mrx160_open_rx(got, &mc);
		if (!eq_rx(ref, got)) {
			fail("neg");
		}

		mrx160_a2b(&ba, &ma, &rng);
		mrx160_b2a(&mt, &ba, &rng);
		mrx160_open_rx(got, &mt);
		if (!eq_rx(a, got)) {
			fail("A2B/B2A");
		}

		mrx160_mul(&mc, &ma, &mb, &rng);
		rx_mul(ref, a, b);
		mrx160_open_rx(got, &mc);
		if (!eq_rx(ref, got)) {
			fail("Q32.64 multiply");
		}

                        mrx160_sqr(&mc, &ma, &rng);
                        rx_sqr(ref, a);
                        mrx160_open_rx(got, &mc);
                        if (!eq_rx(ref, got)) {
                                fail("Q32.64 square");
                        }

                        /* The dedicated square supports an in-place destination. */
                        mrx160_sqr(&ma, &ma, &rng);
                        mrx160_open_rx(got, &ma);
                        if (!eq_rx(ref, got)) {
                                fail("Q32.64 square alias");
                }

		uint32_t before[3], after[3];
		mrx160_open_rx(before, &ma);
		mrx160_refresh(&ma, &rng);
		mrx160_open_rx(after, &ma);
		if (!eq_rx(before, after)) {
			fail("refresh");
		}
	}
}

static void
test_barrett_reduce(void)
{
	static const uint32_t MODS[] = {
		2, 3, 4, 5, 12289, 18433, 65536, 65537,
		226523137, UINT32_C(0x80000000), UINT32_MAX
	};
	static const uint64_t EDGES[] = {
		0, 1, UINT32_MAX, UINT64_C(0x100000000),
		UINT64_C(0x7FFFFFFFFFFFFFFF), UINT64_C(0x8000000000000000),
		UINT64_MAX
	};
	uint64_t rs = UINT64_C(0x6C8E9CF570932BD5);
	if (fndsa_mask_reduce_u64(UINT64_MAX, 0) != 0
		|| fndsa_mask_reduce_u64(UINT64_MAX, 1) != 0)
	{
		fail("Barrett invalid modulus");
	}
	for (size_t j = 0; j < sizeof MODS / sizeof MODS[0]; j ++) {
		uint32_t q = MODS[j];
		for (size_t i = 0; i < sizeof EDGES / sizeof EDGES[0]; i ++) {
			if (fndsa_mask_reduce_u64(EDGES[i], q) != EDGES[i] % q) {
				fail("Barrett edge reduction");
			}
		}
		for (unsigned i = 0; i < 5000; i ++) {
			uint64_t x = prng64(&rs);
			if (fndsa_mask_reduce_u64(x, q) != x % q) {
				fail("Barrett random reduction");
			}
		}
	}
}

static void
plain_ntt(unsigned logn, uint16_t *a, fndsa_mask_modulus mod)
{
	if (mod == FNDSA_MASK_MOD_Q) {
		mqpoly_ext_to_int(logn, a);
		mqpoly_int_to_ntt(logn, a);
	} else {
		mq2poly_unsigned_to_mod(logn, a);
		mq2poly_NTT(logn, a);
	}
}

static void
plain_intt(unsigned logn, uint16_t *a, fndsa_mask_modulus mod)
{
	if (mod == FNDSA_MASK_MOD_Q) {
		mqpoly_ntt_to_int(logn, a);
		mqpoly_int_to_ext(logn, a);
	} else {
		mq2poly_iNTT(logn, a);
	}
}

static void
test_masked_ntt_one(unsigned logn, fndsa_mask_modulus mod)
{
	uint32_t q = (uint32_t)mod;
	size_t n = (size_t)1 << logn;
	mmod32 *a = malloc(n * sizeof *a);
	mmod32 *b = malloc(n * sizeof *b);
	mmod32 *d = malloc(n * sizeof *d);
	uint16_t *ra = malloc(n * sizeof *ra);
	uint16_t *rb = malloc(n * sizeof *rb);
	uint16_t *orig = malloc(n * sizeof *orig);
	uint32_t *pub = malloc(n * sizeof *pub);
	if (a == NULL || b == NULL || d == NULL || ra == NULL || rb == NULL
		|| orig == NULL || pub == NULL)
	{
		fail("masked NTT allocation");
	}
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	uint64_t rs = UINT64_C(0xA729C85E13D4B60F)
		^ ((uint64_t)logn << 32) ^ q;
	for (size_t i = 0; i < n; i ++) {
		orig[i] = (uint16_t)(prng64(&rs) % q);
		ra[i] = orig[i];
		rb[i] = (uint16_t)(prng64(&rs) % q);
		mmod32_share(&a[i], ra[i], q, &rng);
		mmod32_share(&b[i], rb[i], q, &rng);
	}
	plain_ntt(logn, ra, mod);
	plain_ntt(logn, rb, mod);
	if (!mmpoly_public_ntt(logn, pub, orig, mod)) {
		fail("public NTT accepted parameters");
	}
	for (size_t i = 0; i < n; i ++) {
		mmod32 t = { { pub[i], 0 } };
		if (mmod32_open(&t, q) != (uint32_t)ra[i] % q) {
			fprintf(stderr, "public NTT mismatch logn=%u mod=%u i=%zu got=%u ref=%u\\n",
				logn, q, i, mmod32_open(&t, q), (uint32_t)ra[i] % q);
			fail("public NTT differential");
		}
	}
	if (!mmpoly_ntt(logn, a, mod, &rng)
		|| !mmpoly_ntt(logn, b, mod, &rng))
	{
		fail("masked NTT accepted parameters");
	}
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&a[i], q) != (uint32_t)ra[i] % q
			|| mmod32_open(&b[i], q) != (uint32_t)rb[i] % q)
		{
			fail("masked NTT differential");
		}
	}

	/* Secret and public pointwise products have the same reconstruction. */
	mmpoly_mul_ntt(logn, d, a, b, mod, &rng);
	for (size_t i = 0; i < n; i ++) {
		uint32_t ref = (uint32_t)(((uint64_t)(ra[i] % q)
			* (rb[i] % q)) % q);
		if (mmod32_open(&d[i], q) != ref) {
			fail("masked NTT pointwise multiplication");
		}
		pub[i] = mmod32_encode_public(rb[i], q);
	}
	mmpoly_mul_public_ntt(logn, d, a, pub, mod, &rng);
	for (size_t i = 0; i < n; i ++) {
		uint32_t ref = (uint32_t)(((uint64_t)(ra[i] % q)
			* (rb[i] % q)) % q);
		if (mmod32_open(&d[i], q) != ref) {
			fail("masked public NTT multiplication");
		}
	}

	mmpoly_gram_ntt(logn, d, a, b, mod, &rng);
	for (size_t i = 0; i < n; i ++) {
		size_t j = n - 1 - i;
		uint32_t ref = (uint32_t)(((uint64_t)(ra[i] % q)
			* (ra[j] % q) + (uint64_t)(rb[i] % q)
			* (rb[j] % q)) % q);
		if (mmod32_open(&d[i], q) != ref) {
			fail("masked NTT Gram");
		}
	}
	mmpoly_cross_ntt(logn, d, a, b, b, a, mod, &rng);
	for (size_t i = 0; i < n; i ++) {
		size_t j = n - 1 - i;
		uint32_t ref = (uint32_t)(((uint64_t)(ra[i] % q)
			* (rb[j] % q) + (uint64_t)(rb[i] % q)
			* (ra[j] % q)) % q);
		if (mmod32_open(&d[i], q) != ref) {
			fail("masked NTT cross product");
		}
	}

	/* Linear operations are checked before the inverse transform. */
	mmpoly_add(logn, d, a, b, mod);
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&d[i], q)
			!= ((uint32_t)(ra[i] % q) + (rb[i] % q)) % q)
		{
			fail("masked polynomial add");
		}
	}
	mmpoly_sub(logn, d, a, b, mod);
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&d[i], q)
			!= ((uint32_t)(ra[i] % q) + q - (rb[i] % q)) % q)
		{
			fail("masked polynomial subtract");
		}
	}
	mmpoly_neg(logn, d, a, mod);
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&d[i], q)
			!= (q - (uint32_t)(ra[i] % q)) % q)
		{
			fail("masked polynomial negation");
		}
	}

	plain_intt(logn, ra, mod);
	if (!mmpoly_intt(logn, a, mod, &rng)) {
		fail("masked iNTT accepted parameters");
	}
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&a[i], q) != orig[i]
			|| (uint32_t)ra[i] % q != orig[i])
		{
			fail("masked iNTT round trip");
		}
	}

	free(pub);
	free(orig);
	free(rb);
	free(ra);
	free(d);
	free(b);
	free(a);
}

static void
test_signed_ntt_edges(unsigned logn, fndsa_mask_modulus mod)
{
	uint32_t modulus = (uint32_t)mod;
	size_t count = (size_t)1 << logn;
	mmod32 values[1024], scaled[1024];
	uint16_t original[1024], reference[1024], public_input[1024];
	uint32_t *public_ntt = malloc(count * sizeof *public_ntt);
	if (public_ntt == NULL) {
		fail("signed NTT allocation");
	}
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	for (unsigned pattern = 0; pattern < 4; pattern ++) {
		for (size_t index = 0; index < count; index ++) {
			for (size_t share = 0; share < FNDSA_MASK_SHARES; share ++) {
				uint32_t residue = pattern == 0 ? 0 : modulus - 1;
				if (pattern == 2) {
					residue = ((index + share) & 1) ? modulus - 1 : 0;
				} else if (pattern == 3) {
					residue = (uint32_t)(index * 7919 + share * 11) % modulus;
				}
				values[index].share[share] = residue;
			}
			original[index] = (uint16_t)mmod32_open(&values[index], modulus);
			reference[index] = original[index];
			public_input[index] = (uint16_t)(original[index]
				+ ((UINT16_MAX - original[index]) / modulus) * modulus);
		}
		plain_ntt(logn, reference, mod);
		memcpy(public_ntt, public_input, count * sizeof *public_input);
		if (!mmpoly_public_ntt(logn, public_ntt,
			(const uint16_t *)public_ntt, mod)
			|| !mmpoly_ntt(logn, values, mod, &rng))
		{
			fail("signed NTT edge parameters");
		}
		for (size_t index = 0; index < count; index ++) {
			uint32_t expected = reference[index] % modulus;
			if (values[index].share[0] >= modulus
				|| values[index].share[1] >= modulus
				|| mmod32_open(&values[index], modulus) != expected
				|| public_ntt[index] != mmod32_encode_public(expected, modulus))
			{
				fail("signed NTT edges, canonicality or public alias");
			}
		}
		memcpy(scaled, values, count * sizeof *values);
		mmpoly_mul_public_const(logn, scaled, scaled, UINT32_MAX, mod, &rng);
		for (size_t index = 0; index < count; index ++) {
			if (scaled[index].share[0] >= modulus
				|| scaled[index].share[1] >= modulus
				|| mmod32_open(&scaled[index], modulus)
					!= (uint64_t)(reference[index] % modulus) * UINT32_MAX % modulus)
			{
				fail("signed public constant multiply alias");
			}
		}
		if (!mmpoly_intt(logn, values, mod, &rng)) {
			fail("signed inverse NTT edge parameters");
		}
		for (size_t index = 0; index < count; index ++) {
			if (values[index].share[0] >= modulus
				|| values[index].share[1] >= modulus
				|| mmod32_open(&values[index], modulus) != original[index])
			{
				fail("signed inverse NTT edges or canonicality");
			}
		}
	}
	free(public_ntt);
}

static void
test_masked_ntt_crt(void)
{
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	mmod32 dummy = { { 1, 2 } };
	if (mmpoly_ntt(0, &dummy, FNDSA_MASK_MOD_Q, &rng)
		|| mmpoly_ntt(11, &dummy, FNDSA_MASK_MOD_Q, &rng)
		|| mmpoly_ntt(1, &dummy, (fndsa_mask_modulus)17, &rng))
	{
		fail("masked NTT public parameter rejection");
	}
	for (unsigned logn = 1; logn <= 10; logn ++) {
		test_masked_ntt_one(logn, FNDSA_MASK_MOD_Q);
		test_masked_ntt_one(logn, FNDSA_MASK_MOD_Q2);
		test_signed_ntt_edges(logn, FNDSA_MASK_MOD_Q);
		test_signed_ntt_edges(logn, FNDSA_MASK_MOD_Q2);
	}

	enum { CRT_LOGN = 6, CRT_N = 1 << CRT_LOGN };
	mmod32 aq[CRT_N], aq2[CRT_N];
	mi32 out[CRT_N];
	static const uint32_t edge[] = {
		0, 1, FNDSA_MASK_QQ2_H, FNDSA_MASK_QQ2_H + 1,
		FNDSA_MASK_QQ2 - 2, FNDSA_MASK_QQ2 - 1
	};
	uint32_t x[CRT_N];
	uint64_t rs = UINT64_C(0xD0B561A7C84E932F);
	for (size_t i = 0; i < CRT_N; i ++) {
		x[i] = i < sizeof edge / sizeof edge[0]
			? edge[i] : (uint32_t)(prng64(&rs) % FNDSA_MASK_QQ2);
		mmod32_share(&aq[i], x[i] % FNDSA_MASK_Q,
			FNDSA_MASK_Q, &rng);
		mmod32_share(&aq2[i], x[i] % FNDSA_MASK_Q2,
			FNDSA_MASK_Q2, &rng);
	}
	mmpoly_crt_center(CRT_LOGN, out, aq, aq2, &rng);
	for (size_t i = 0; i < CRT_N; i ++) {
		int32_t ref = x[i] > FNDSA_MASK_QQ2_H
			? (int32_t)(x[i] - FNDSA_MASK_QQ2) : (int32_t)x[i];
		if (mi32_open_compat(&out[i]) != ref) {
			fail("masked CRT and centering differential");
		}
	}
}

static void
test_masked_acceptance(void)
{
	enum { LOGN = 9, N = 1 << LOGN };
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	mi32 *s0 = calloc(N, sizeof *s0);
	mi32 *s1 = calloc(N, sizeof *s1);
	if (s0 == NULL || s1 == NULL) {
		fail("masked acceptance allocation");
	}
	mbit accept;
	size_t payload_len = FNDSA_SIGNATURE_SIZE(LOGN) - 41;
	for (size_t i = 0; i < N; i ++) {
		mi32_share(&s0[i], 0, &rng);
		mi32_share(&s1[i], 0, &rng);
	}
	mi32_signature_is_acceptable(&accept, LOGN, s0, s1,
		payload_len, &rng);
	if (mbit_open_compat(&accept) != 1) {
		fail("masked acceptance zero vector");
	}

	mi32_share(&s0[0], 841, &rng);
	mi32_signature_is_acceptable(&accept, LOGN, s0, s1,
		payload_len, &rng);
	if (mbit_open_compat(&accept) != 0) {
		fail("masked coefficient range rejection");
	}
	mi32_share(&s0[0], 0, &rng);

	for (size_t i = 0; i < N; i ++) {
		mi32_share(&s0[i], 300, &rng);
	}
	mi32_signature_is_acceptable(&accept, LOGN, s0, s1,
		payload_len, &rng);
	if (mbit_open_compat(&accept) != 0) {
		fail("masked norm rejection");
	}
	for (size_t i = 0; i < N; i ++) {
		mi32_share(&s0[i], 0, &rng);
		mi32_share(&s1[i], 128, &rng);
	}
	mi32_signature_is_acceptable(&accept, LOGN, s0, s1,
		payload_len, &rng);
	if (mbit_open_compat(&accept) != 0) {
		fail("masked compressed-length rejection");
	}

	mi32 pair[2];
	mi32_share(&pair[0], 3, &rng);
	mi32_share(&pair[1], 4, &rng);
	mi32_poly_sqnorm_le(&accept, 1, pair, NULL, 25, &rng);
	if (mbit_open_compat(&accept) != 1) {
		fail("masked norm exact bound");
	}
	mi32_poly_sqnorm_le(&accept, 1, pair, NULL, 24, &rng);
	if (mbit_open_compat(&accept) != 0) {
		fail("masked norm boundary rejection");
	}
	free(s1);
	free(s0);
}

static void
check_signed_redc(int64_t product, uint32_t modulus, uint32_t q0i)
{
	int32_t reduced = mask_mont32_redc(product, modulus, q0i);
	if ((int64_t)reduced <= -(int64_t)modulus || (int64_t)reduced >= modulus) {
		fail("signed REDC range");
	}
	int64_t expected = product % modulus;
	if (expected < 0) {
		expected += modulus;
	}
	uint32_t canonical = mask_mont32_cadd(reduced, modulus);
	if (canonical >= modulus
		|| ((uint64_t)canonical << 32) % modulus != (uint64_t)expected)
	{
		fail("signed REDC congruence");
	}
}

static void
test_montgomery_boundaries(uint32_t modulus, fndsa_mask_rng *rng)
{
	uint32_t edges[] = { 0, 1, modulus / 2, modulus - 1, modulus, UINT32_MAX };
	for (size_t left = 0; left < sizeof edges / sizeof edges[0]; left ++) {
		uint32_t encoded = mmod32_encode_public(edges[left], modulus);
		if (encoded != ((uint64_t)edges[left] << 32) % modulus) {
			fail("Montgomery encode boundary");
		}
		for (size_t right = 0; right < sizeof edges / sizeof edges[0]; right ++) {
			mmod32 first = { { edges[left] % modulus, edges[right] % modulus } };
			mmod32 second = { { edges[right] % modulus, edges[left] % modulus } };
			uint32_t opened = mmod32_open(&first, modulus);
			uint32_t expected = (uint32_t)((uint64_t)opened * opened % modulus);
			mmod32_mul(&first, &first, &second, modulus, rng);
			if (first.share[0] >= modulus || first.share[1] >= modulus
				|| mmod32_open(&first, modulus) != expected)
			{
				fail("Montgomery boundary multiply alias");
			}
			mmod32_mul(&second, &first, &second, modulus, rng);
			if (second.share[0] >= modulus || second.share[1] >= modulus
				|| mmod32_open(&second, modulus)
					!= (uint64_t)expected * opened % modulus)
			{
				fail("Montgomery boundary right multiply alias");
			}
			mmod32_set_residue_shares(&first, edges[left], edges[right], modulus, rng);
			if (first.share[0] >= modulus || first.share[1] >= modulus
				|| mmod32_open(&first, modulus)
					!= ((uint64_t)edges[left] + edges[right]) % modulus)
			{
				fail("Montgomery noncanonical residue import");
			}
			if (modulus <= INT32_MAX) {
				uint32_t first_residue = edges[left] % modulus;
				uint32_t second_residue = edges[right] % modulus;
				if (mask_mont32_add(first_residue, second_residue, modulus)
						!= ((uint64_t)first_residue + second_residue) % modulus
					|| mask_mont32_sub(first_residue, second_residue, modulus)
						!= ((uint64_t)first_residue + modulus - second_residue) % modulus)
				{
					fail("signed modular add/sub boundaries");
				}
			}
		}
	}
	if (modulus <= INT32_MAX) {
		uint32_t inverse = 1;
		for (unsigned round = 0; round < 5; round ++) {
			inverse *= 2 - modulus * inverse;
		}
		uint32_t q0i = (uint32_t)0 - inverse;
		int64_t limit = (int64_t)modulus * INT64_C(2147483648);
		int64_t products[] = {
			-limit + 1, -limit + 2, -INT64_C(4294967296), -1, 0, 1,
			INT64_C(4294967296), limit - 2, limit - 1,
			(int64_t)(modulus - 1) * (modulus - 1)
		};
		for (size_t index = 0; index < sizeof products / sizeof products[0]; index ++) {
			check_signed_redc(products[index], modulus, q0i);
		}
		uint64_t state = UINT64_C(0x65E9B71A23740C8D);
		for (unsigned iteration = 0; iteration < 10000; iteration ++) {
			int64_t product = (int64_t)(prng64(&state) % (uint64_t)(2 * limit - 1))
				- (limit - 1);
			check_signed_redc(product, modulus, q0i);
		}
	}
}

static void
test_modular(void)
{
	test_barrett_reduce();
	static const uint32_t MODS[] = {
		3, 5, 19, 29, 65537, 12289, 18433, 226523137,
		UINT32_C(1073741823), UINT32_C(2147483645), UINT32_C(2147483647),
		UINT32_C(2147483649), UINT32_C(2147483659), UINT32_C(4294967291),
		UINT32_MAX
	};
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	{
		static const int32_t tv[] = {
			INT32_MIN, -18434, -12290, -1, 0, 1, 12288, 18432, INT32_MAX
		};
		for (size_t ti = 0; ti < sizeof tv / sizeof tv[0]; ti ++) {
			for (size_t qm = 0; qm < 2; qm ++) {
				uint32_t q = qm == 0 ? FNDSA_MASK_Q : FNDSA_MASK_Q2;
				mi32 a;
				mmod32 b;
				mi32_share(&a, tv[ti], &rng);
				mi32_to_mmod32(&b, &a, q, &rng);
				int64_t rr = tv[ti];
				while (rr < 0) rr += q;
				while (rr >= q) rr -= q;
				if (mmod32_open(&b, q) != (uint32_t)rr) {
					fprintf(stderr, "mi32 mod: x=%d q=%u got=%u ref=%u\\n",
						tv[ti], q, mmod32_open(&b, q), (uint32_t)rr);
					fail("mi32 to masked modulus");
				}
			}
		}
	}
	uint64_t rs = UINT64_C(0xC9E48B7215AF603D);
	mmod32 invalid;
	mmod32_share(&invalid, UINT32_MAX, 2, &rng);
	if (invalid.share[0] != 0 || invalid.share[1] != 0
		|| mmod32_open(&invalid, 2) != 0)
	{
		fail("Montgomery even modulus accepted");
	}

	for (size_t j = 0; j < sizeof MODS / sizeof MODS[0]; j ++) {
		uint32_t q = MODS[j];
		test_montgomery_boundaries(q, &rng);
		for (unsigned k = 0; k < 1000; k ++) {
			uint32_t x = (uint32_t)(prng64(&rs) % q);
			uint32_t y = (uint32_t)(prng64(&rs) % q);
			mmod32 a, b, c;
			mmod32_share(&a, x, q, &rng);
			mmod32_share(&b, y, q, &rng);
			uint32_t raw = (uint32_t)prng64(&rs);
			mmod32_share(&c, raw, q, &rng);
			if (a.share[0] >= q || a.share[1] >= q
				|| b.share[0] >= q || b.share[1] >= q
				|| c.share[0] >= q || c.share[1] >= q)
			{
				fail("Montgomery share canonicality");
			}
			if (mmod32_open(&c, q) != raw % q) {
				fail("modular non-canonical share");
			}
			mmod32_add(&c, &a, &b, q);
			if (c.share[0] >= q || c.share[1] >= q) {
				fail("Montgomery add canonicality");
			}
			if (mmod32_open(&c, q) != ((uint64_t)x + y) % q) {
				fail("modular add");
			}
			mmod32_sub(&c, &a, &b, q);
			if (c.share[0] >= q || c.share[1] >= q) {
				fail("Montgomery subtract canonicality");
			}
			if (mmod32_open(&c, q) != ((uint64_t)x + q - y) % q) {
				fail("modular sub");
			}
			mmod32_mul(&c, &a, &b, q, &rng);
			if (c.share[0] >= q || c.share[1] >= q) {
				fail("Montgomery multiply canonicality");
			}
			if (mmod32_open(&c, q)
				!= (uint32_t)(((uint64_t)x * y) % q))
			{
				fail("modular multiply");
			}
			uint32_t z = mmod32_open(&a, q);
			mmod32_refresh(&a, q, &rng);
			if (a.share[0] >= q || a.share[1] >= q) {
				fail("Montgomery refresh canonicality");
			}
			if (mmod32_open(&a, q) != z) {
				fail("modular refresh");
			}
		}
	}
}

static uint64_t
ref_expm_p63(uint64_t x, uint64_t ccs)
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
	uint64_t y = C[0];
	for (size_t i = 1; i < sizeof C / sizeof C[0]; i ++) {
		unsigned __int128 p = (unsigned __int128)y * x;
		y = C[i] - (uint64_t)(p >> 64);
	}
	return (uint64_t)(((unsigned __int128)y * ccs) >> 64);
}

typedef struct {
	uint8_t buf[16];
	size_t ptr;
} scripted_source;

static uint8_t
scripted_next_u8(void *ctx)
{
	scripted_source *s = ctx;
	return s->buf[s->ptr ++];
}

static uint8_t
shake_source_next_u8(void *ctx)
{
	uint8_t x;
	shake_extract((shake_context *)ctx, &x, 1);
	return x;
}

static uint64_t
ref_ber_threshold(const rx x, uint64_t ccs, uint32_t *sp)
{
	static const rx INV_LOG2 = {
		UINT32_C(3090145655), UINT32_C(1901360722), 1
	};
	static const rx LOG2 = {
		UINT32_C(3520035244), UINT32_C(2977044471), 0
	};
	rx t, r;
	rx_mul(t, x, INV_LOG2);
	int32_t si = rx_floor(t);
	rx_of32(t, si);
	rx_mul(t, t, LOG2);
	rx_sub(r, x, t);
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;
	s &= 63;
	if (sp != NULL) {
		*sp = s;
	}
	return ((ref_expm_p63(rx_ufrac(r), ccs) << 1) - 1) >> s;
}

static void
test_masked_floor_inverse_expm(void)
{
	static const rx FLOOR_EDGE[] = {
		{ 0, 0, 0 },
		{ 1, 0, 0 },
		{ UINT32_MAX, UINT32_MAX, UINT32_MAX },
		{ 0, UINT32_C(0x80000000), 0 },
		{ 0, UINT32_C(0x80000000), UINT32_MAX },
		{ UINT32_MAX, UINT32_MAX, UINT32_C(0x7FFFFFFF) },
		{ 0, 0, UINT32_C(0x80000000) }
	};
	static const rx INV_EDGE[] = {
		{ UINT32_MAX, 3, 0 },
		{ 0, 4, 0 },
		{ UINT32_MAX, UINT32_C(0x7FFFFFFF), 0 },
		{ 0, UINT32_C(0x80000000), 0 },
		{ UINT32_MAX, UINT32_MAX, 0 },
		{ 0, 0, 1 },
		{ 1, 0, 1 },
		{ UINT32_MAX, UINT32_MAX, 1 },
		{ 0, 0, 2 },
		{ UINT32_MAX, UINT32_MAX, UINT32_C(0x00007FFF) }
	};
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	for (size_t i = 0; i < sizeof INV_EDGE / sizeof INV_EDGE[0]; i ++) {
		mrx160 a, d;
		rx got, ref;
		mrx160_share_rx(&a, INV_EDGE[i], &rng);
		mrx160_inv_spec(&d, &a, &rng);
		mrx160_open_rx(got, &d);
		rx_inv_spec(ref, INV_EDGE[i]);
		if (!eq_rx(got, ref)) {
			fail("masked inverse edge");
		}
	}
	fndsa_mask_rng_init_demo(&rng);
	for (size_t i = 0; i < sizeof FLOOR_EDGE / sizeof FLOOR_EDGE[0]; i ++) {
		mrx160 a, d;
		rx got, ref = { 0, 0, FLOOR_EDGE[i][2] };
		mrx160_share_rx(&a, FLOOR_EDGE[i], &rng);
		mrx160_floor(&d, &a, &rng);
		mrx160_open_rx(got, &d);
		if (!eq_rx(got, ref)) {
			fail("masked floor edge");
		}
	}

	uint64_t rs = UINT64_C(0xD86A3C19570BE42F);
	for (unsigned i = 0; i < 256; i ++) {
		rx a, ref, got;
		a[0] = (uint32_t)prng64(&rs);
		a[1] = (uint32_t)prng64(&rs);
		a[2] = (uint32_t)(1 + (prng64(&rs) % 32767));
		mrx160 ma, md;
		mrx160_share_rx(&ma, a, &rng);
		mrx160_inv_spec(&md, &ma, &rng);
		mrx160_open_rx(got, &md);
		rx_inv_spec(ref, a);
		if (!eq_rx(got, ref)) {
			fail("masked inverse differential");
		}
	}

	for (unsigned i = 0; i < 1000; i ++) {
		uint64_t x = prng64(&rs) % UINT64_C(0xB17217F7D1CF79AC);
		uint64_t ccs = prng64(&rs);
		rx ax = { (uint32_t)x, (uint32_t)(x >> 32), 0 };
		rx ac = { (uint32_t)ccs, (uint32_t)(ccs >> 32), 0 };
		mrx160 mx, mc;
		mb64 d;
		mrx160_share_rx(&mx, ax, &rng);
		mrx160_share_rx(&mc, ac, &rng);
		mrx160_mul_u64_hi(&d, &mx, &mc, &rng);
		uint64_t href = (uint64_t)(((unsigned __int128)x * ccs) >> 64);
		if (mb64_open_compat(&d) != href) {
			fail("masked unsigned high64 differential");
		}
		mrx160_expm_p63(&d, &mx, &mc, &rng);
		if (mb64_open_compat(&d) != ref_expm_p63(x, ccs)) {
			fail("masked expm_p63 differential");
		}
	}
}

static void
test_masked_ber_exp(void)
{
	static const rx LOG2 = {
		UINT32_C(3520035244), UINT32_C(2977044471), 0
	};
	static const uint64_t CCS = UINT64_C(0xC000000000000000);
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	rx cases[3];
	memset(cases, 0, sizeof cases);
	cases[0][1] = UINT32_C(0x20000000);
	rx k;
	rx_of32(k, 63);
	rx_mul(cases[1], k, LOG2);
	cases[1][1] += UINT32_C(0x10000000);
	rx_of32(k, 65);
	rx_mul(cases[2], k, LOG2);

	for (size_t ci = 0; ci < 3; ci ++) {
		uint32_t shift;
		uint64_t z = ref_ber_threshold(cases[ci], CCS, &shift);
		if ((ci == 0 && shift != 0) || (ci == 1 && shift != 63)
			|| (ci == 2 && shift != 63))
		{
			fail("BerExp shift coverage");
		}
		for (unsigned prefix = 0; prefix <= 8; prefix ++) {
			scripted_source ss;
			memset(&ss, 0, sizeof ss);
			for (unsigned i = 0; i < 8; i ++) {
				ss.buf[i] = (uint8_t)(z >> (56 - 8 * i));
			}
			if (prefix < 8) {
				uint8_t bz = ss.buf[prefix];
				ss.buf[prefix] = bz == 0 ? 1 : (uint8_t)(bz - 1);
			}
			int ref = 0;
			for (unsigned i = 0; i < 8; i ++) {
				uint8_t bz = (uint8_t)(z >> (56 - 8 * i));
				uint8_t w = ss.buf[i];
				if (w != bz) {
					ref = w < bz;
					break;
				}
			}
			mrx160 mx, mc;
			rx ccs = { (uint32_t)CCS, (uint32_t)(CCS >> 32), 0 };
			mrx160_share_rx(&mx, cases[ci], &rng);
			mrx160_share_rx(&mc, ccs, &rng);
			ss.ptr = 0;
			fndsa_mask_sampler_source src = {
				&ss, scripted_next_u8, 9
			};
			mbit accept;
			mrx160_ber_exp(&accept, &src, &mx, &mc, &rng);
			if ((int)mbit_open_compat(&accept) != ref
				|| ss.ptr != (prefix < 8 ? prefix + 1 : 8))
			{
				fail("masked BerExp lazy differential");
			}
		}
	}
}

#ifndef MASK_SAMPLER_DIFF_ITERS
#define MASK_SAMPLER_DIFF_ITERS   100000
#endif

static void
test_masked_sampler_z(unsigned logn)
{
	uint8_t seed[48];
	for (size_t i = 0; i < sizeof seed; i ++) {
		seed[i] = (uint8_t)(0x5D + 13 * i + logn);
	}
	shake_context pc_ref, pc_mask;
	shake_init(&pc_ref, 256);
	shake_inject(&pc_ref, seed, sizeof seed);
	shake_flip(&pc_ref);
	pc_mask = pc_ref;
	sampler_state ss;
	memset(&ss, 0, sizeof ss);
	ss.pc = pc_ref;
	ss.logn = logn;
	fndsa_mask_sampler_source src = {
		&pc_mask, shake_source_next_u8, logn
	};
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	uint64_t rs = UINT64_C(0x8E4F71C2A93D650B) ^ logn;
	for (unsigned i = 0; i < MASK_SAMPLER_DIFF_ITERS; i ++) {
		uint64_t mf = prng64(&rs);
		int32_t mi = (int32_t)(prng64(&rs) % 201) - 100;
		rx mu = { (uint32_t)mf, (uint32_t)(mf >> 32), (uint32_t)mi };
		uint64_t isf = UINT64_C(0x9000000000000000)
			| (prng64(&rs) & UINT64_C(0x1FFFFFFFFFFFFFFF));
		rx isigma = { (uint32_t)isf, (uint32_t)(isf >> 32), 0 };
		mrx160 mmu, mis, mz;
		mrx160_share_rx(&mmu, mu, &rng);
		mrx160_share_rx(&mis, isigma, &rng);
		int32_t ref = x_sampler_next(&ss, mu, isigma);
		mrx160_sampler_z(&mz, &src, &mmu, &mis, &rng);
		rx got;
		mrx160_open_rx(got, &mz);
		if (got[0] != 0 || got[1] != 0 || (int32_t)got[2] != ref) {
			fail("masked SamplerZ differential");
		}
	}
	uint8_t a[32], b[32];
	shake_extract(&ss.pc, a, sizeof a);
	shake_extract(&pc_mask, b, sizeof b);
	if (memcmp(a, b, sizeof a) != 0) {
		fail("masked SamplerZ PRNG state");
	}
}

static void
test_masked_div_sqrt_mi32(void)
{
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	uint64_t rs = UINT64_C(0x17D4A6E93CB8520F);
	for (unsigned i = 0; i < 2000; i ++) {
		rx a, b, ref, got;
		a[0] = (uint32_t)prng64(&rs);
		a[1] = (uint32_t)prng64(&rs);
		a[2] = (uint32_t)((int32_t)(prng64(&rs) % 2001) - 1000);
		b[0] = (uint32_t)prng64(&rs);
		b[1] = (uint32_t)prng64(&rs);
		b[2] = (uint32_t)(1 + prng64(&rs) % 1000);
		mrx160 ma, mb, md;
		mrx160_share_rx(&ma, a, &rng);
		mrx160_share_rx(&mb, b, &rng);
		mrx160_div_spec(&md, &ma, &mb, &rng);
		mrx160_open_rx(got, &md);
		rx_div_spec(ref, a, b);
		if (!eq_rx(got, ref)) {
			fail("masked general division differential");
		}

		rx s = {
			(uint32_t)prng64(&rs), (uint32_t)prng64(&rs),
			(uint32_t)(6328 + prng64(&rs) % (16822 - 6328 + 1))
		};
		mrx160_share_rx(&ma, s, &rng);
		mrx160_sqrt_spec(&md, &ma, &rng);
		mrx160_open_rx(got, &md);
		rx_sqrt_spec(ref, s);
		if (!eq_rx(got, ref)) {
			fail("masked square-root differential");
		}

		int32_t z = (int32_t)prng64(&rs);
		mi32 mz, mz2;
		mi32_share(&mz, z, &rng);
		mi32_to_mrx160(&ma, &mz, &rng);
		mrx160_to_mi32(&mz2, &ma, &rng);
		if (mi32_open_compat(&mz2) != z) {
			fail("mi32 conversion round trip");
		}
	}
}

static void
check_masked_poly_value(const char *what, unsigned logn, size_t i,
	const mrx160 *a, const rx ref)
{
	rx got;
	mrx160_open_rx(got, a);
	if (!eq_rx(got, ref)) {
		fprintf(stderr, "%s: logn=%u index=%zu\n", what, logn, i);
		fail("masked polynomial differential");
	}
}


static void
test_xpoly_constant_ranges(void)
{
	for (size_t i = 1; i < 512; i ++) {
		if (raw96_abs_bits(fndsa_X_PLUS_1_SQUARED_OVER_X[i]) > 66
			|| raw96_abs_bits(
				fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[i]) > 74)
		{
			fail("masked tree public constant range");
		}
	}
}
static void
test_masked_xpoly(void)
{
	test_xpoly_constant_ranges();
	uint64_t rs = UINT64_C(0xC54A19E872B30D6F);
	for (unsigned logn = 1; logn <= 10; logn ++) {
		size_t n = (size_t)1 << logn;
		size_t hn = n >> 1;
		int32_t *v = malloc(n * sizeof *v);
		rx *u = malloc(n * sizeof *u);
		rx *u1 = malloc(hn * sizeof *u1);
		rx *um = malloc(n * sizeof *um);
		mi32 *mv = malloc(n * sizeof *mv);
		mi32 *scratch = malloc((n + hn) * sizeof *scratch);
		mrx160 *m = malloc(n * sizeof *m);
		mrx160 *m1 = malloc(hn * sizeof *m1);
		mrx160 *mm = malloc(n * sizeof *mm);
		if (v == NULL || u == NULL || u1 == NULL || um == NULL
			|| mv == NULL || scratch == NULL || m == NULL
			|| m1 == NULL || mm == NULL)
		{
			fail("masked xpoly allocation");
		}
		fndsa_mask_rng rng;
		fndsa_mask_rng_init_demo(&rng);
		for (size_t i = 0; i < n; i ++) {
			v[i] = (int32_t)(prng64(&rs) % 2001) - 1000;
			mi32_share(&mv[i], v[i], &rng);
		}
		xpoly_FFT_int32(logn, u, v);
		mxpoly_fft_i32(logn, m, mv, scratch, &rng);
		for (size_t i = 0; i < n; i ++) {
			check_masked_poly_value("FFT", logn, i, &m[i], u[i]);
		}
		xpoly_split_fft(logn, u, u1);
		mxpoly_split_fft(logn, m, m1, &rng);
		for (size_t i = 0; i < hn; i ++) {
			check_masked_poly_value("split-0", logn, i, &m[i], u[i]);
			check_masked_poly_value("split-1", logn, i, &m1[i], u1[i]);
		}
		xpoly_merge_fft(logn, um, u, u1);
		mxpoly_merge_fft(logn, mm, m, m1, &rng);
		for (size_t i = 0; i < n; i ++) {
			check_masked_poly_value("merge", logn, i, &mm[i], um[i]);
		}
		free(mm); free(m1); free(m); free(scratch); free(mv);
		free(um); free(u1); free(u); free(v);
	}

	for (unsigned logn = 1; logn <= 9; logn ++) {
		size_t hn = (size_t)1 << (logn - 1);
		rx *a = malloc(hn * sizeof *a);
		rx *b = malloc(hn * sizeof *b);
		rx *d = malloc(hn * sizeof *d);
		rx *e = malloc(hn * sizeof *e);
		mrx160 *ma = malloc(hn * sizeof *ma);
		mrx160 *mb = malloc(hn * sizeof *mb);
		mrx160 *md = malloc(hn * sizeof *md);
		mrx160 *me = malloc(hn * sizeof *me);
		mrx160 *mdt = malloc(hn * sizeof *mdt);
		mrx160 *met = malloc(hn * sizeof *met);
		if (a == NULL || b == NULL || d == NULL || e == NULL
			|| ma == NULL || mb == NULL || md == NULL || me == NULL
			|| mdt == NULL || met == NULL)
		{
			fail("masked LDL allocation");
		}
		fndsa_mask_rng rng;
		fndsa_mask_rng_init_demo(&rng);
		for (size_t i = 0; i < hn; i ++) {
			rx_of32(a[i], 100 + (int32_t)(prng64(&rs) % 100));
			rx_of32(b[i], (int32_t)(prng64(&rs) % 11) - 5);
			mrx160_share_rx(&ma[i], a[i], &rng);
			mrx160_share_rx(&mb[i], b[i], &rng);
		}
		fndsa_mask_rng rng_wide = rng;
		fndsa_mask_rng rng_tree = rng;
		xpoly_LDL_fft(logn, d, e, a, b);
		mxpoly_ldl_fft(logn, md, me, ma, mb, &rng_wide);
		mxpoly_ldl_tree_fft(logn, mdt, met, ma, mb, &rng_tree);
		if (memcmp(&rng_wide, &rng_tree, sizeof rng_wide) != 0) {
			fail("tree fast/wide mask RNG consumption");
		}
		for (size_t i = 0; i < hn; i ++) {
			check_masked_poly_value("LDL-d", logn, i, &md[i], d[i]);
			check_masked_poly_value("LDL-e", logn, i, &me[i], e[i]);
			check_masked_poly_value("LDL-tree-d", logn, i,
				&mdt[i], d[i]);
			check_masked_poly_value("LDL-tree-e", logn, i,
				&met[i], e[i]);
		}
		free(met); free(mdt); free(me); free(md); free(mb); free(ma);
		free(e); free(d); free(b); free(a);
	}
}

static int
test_dynamic_tree_one(unsigned logn, int right)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t mlen = fndsa_masked_dyn_tmp_size(logn);
	void *mtmp = malloc(mlen);
	rx *utmp = calloc(2 * n, sizeof *utmp);
	int16_t *udz = calloc(n, sizeof *udz);
	if (mtmp == NULL || utmp == NULL || udz == NULL) {
		fail("masked dynamic tree allocation");
	}
	uint8_t seed[40];
	for (size_t i = 0; i < sizeof seed; i ++) {
		seed[i] = (uint8_t)(3 + 7 * i + logn + (unsigned)right);
	}
	sampler_state us;
	memset(&us, 0, sizeof us);
	shake_init(&us.pc, 256);
	shake_inject(&us.pc, seed, sizeof seed);
	shake_flip(&us.pc);
	us.dz = udz;
	us.logn = logn;

	shake_context mpc;
	shake_init(&mpc, 256);
	shake_inject(&mpc, seed, sizeof seed);
	shake_flip(&mpc);
	fndsa_mask_sampler_source src = {
		&mpc, shake_source_next_u8, logn
	};
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	fndsa_mask_tree_state ms;
	if (!fndsa_mask_tree_bind(&ms, logn, &src, &rng, mtmp, mlen)) {
		fail("masked dynamic tree bind");
	}
	mrx160 *mw = fndsa_mask_tree_work(&ms);
	mi32 *mdz = fndsa_mask_tree_samples(&ms);
	rx a, t;
	for (size_t i = 0; i < n; i ++) {
		rx_of32(t, (int32_t)(i % 5) - 2);
		rx_copy(utmp[i], t);
		mrx160_share_rx(&mw[i], t, &rng);
	}
	rx_of32(a, 10000);
	for (size_t i = 0; i < hn; i ++) {
		rx_copy(utmp[n + i], a);
		mrx160_share_rx(&mw[n + i], a, &rng);
	}
	for (size_t i = 0; i < n; i ++) {
		mi32_share(&mdz[i], 0, &rng);
	}

	if (right) {
		x_ffsamp_fft_right(&us, 0, utmp);
		fndsa_mask_ffsampling_right(&ms, 0);
	} else {
		x_ffsamp_fft_left(&us, 0, utmp);
		fndsa_mask_ffsampling_left(&ms, 0);
	}
	for (size_t i = 0; i < n; i ++) {
		if (mi32_open_compat(&mdz[i]) != udz[i]) {
			fail("masked dynamic tree sample");
		}
		check_masked_poly_value(right ? "right-tree" : "left-tree",
			logn, i, &mw[i], utmp[i]);
	}
	uint8_t aout[32], bout[32];
	shake_extract(&us.pc, aout, sizeof aout);
	shake_extract(&mpc, bout, sizeof bout);
	if (memcmp(aout, bout, sizeof aout) != 0) {
		fail("masked dynamic tree SHAKE state");
	}
	fndsa_mask_tree_clear(&ms);
	if (!all_zero(mtmp, mlen)) {
		fail("masked dynamic tree clear");
	}
	free(udz);
	free(utmp);
	free(mtmp);
	return 1;
}

static void
test_masked_dynamic_tree(void)
{
	fndsa_mask_range_stats_reset();
	if (fndsa_masked_dyn_tmp_size(8) != 0
		|| fndsa_masked_dyn_tmp_size(11) != 0
		|| fndsa_masked_dyn_tmp_size(9)
			!= 2 * 512 * sizeof(mrx160) + 512 * sizeof(mi32) + 32
		|| fndsa_masked_dyn_tmp_size(10)
			!= 2 * 1024 * sizeof(mrx160) + 1024 * sizeof(mi32) + 32)
	{
		fail("masked dynamic workspace size");
	}
	{
		uint8_t tiny[1];
		fndsa_mask_tree_state st;
		fndsa_mask_sampler_source src = {
			NULL, shake_source_next_u8, 9
		};
		fndsa_mask_rng rng;
		memset(&st, 0xA5, sizeof st);
		fndsa_mask_rng_init_demo(&rng);
		if (fndsa_mask_tree_bind(&st, 9, &src, &rng,
			tiny, sizeof tiny) || !all_zero((const uint8_t *)&st, sizeof st))
		{
			fail("masked dynamic short workspace");
		}
	}
	(void)test_dynamic_tree_one(9, 0);
	(void)test_dynamic_tree_one(9, 1);
	(void)test_dynamic_tree_one(10, 0);
	(void)test_dynamic_tree_one(10, 1);
}


static void
check_fast_range_stats(void)
{
	fndsa_mask_range_stats stats;
	fndsa_mask_range_stats_get(&stats);
	for (unsigned i = 0; i < FNDSA_MASK_RANGE_SITE_COUNT; i ++) {
		if (stats.calls[i] == 0 || stats.max_product_bits[i] > 159) {
			fail("masked dynamic fast-multiply range instrumentation");
		}
	}
	printf("[fast-product bits %u/%u/%u/%u/%u]",
		stats.max_product_bits[FNDSA_MASK_RANGE_LDL_THETA_A0],
		stats.max_product_bits[FNDSA_MASK_RANGE_LDL_SCHUR],
		stats.max_product_bits[FNDSA_MASK_RANGE_IXADJ_THETA_B],
		stats.max_product_bits[FNDSA_MASK_RANGE_IXADJ_MU_A],
		stats.max_product_bits[FNDSA_MASK_RANGE_IXADJ_FINAL]);
}
static uint32_t
signed_to_mod(int32_t x, uint32_t q)
{
	return x < 0 ? (uint32_t)(x + (int32_t)q) : (uint32_t)x;
}

static void
test_masked_key_mq_expansion(unsigned logn,
	const uint8_t *sk, size_t sk_len,
	const uint8_t *mk, size_t mk_len,
	const uint8_t *vk, size_t vk_len)
{
	(void)sk_len;
	size_t n = (size_t)1 << logn;
	mmod32 *p = calloc(10 * n, sizeof *p);
	if (p == NULL) {
		fail("masked key MQ expansion allocation");
	}
	fndsa_masked_basis_mq basis = {
		p + 0 * n, p + 1 * n, p + 2 * n, p + 3 * n,
		p + 4 * n, p + 5 * n, p + 6 * n, p + 7 * n
	};
	mmod32 *t0 = p + 8 * n;
	mmod32 *t1 = p + 9 * n;
	const uint8_t *got_vk;
	size_t got_vk_len;
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	if (fndsa_masked_key_expand_mq_internal(mk, mk_len - 1,
		&basis, &got_vk, &got_vk_len, &rng))
	{
		fail("truncated masked key MQ expansion");
	}
	if (!fndsa_masked_key_expand_mq_internal(mk, mk_len,
		&basis, &got_vk, &got_vk_len, &rng)
		|| got_vk_len != vk_len || memcmp(got_vk, vk, vk_len) != 0)
	{
		fail("masked key MQ expansion");
	}

	int8_t f[1024], g[1024], F[1024];
	size_t off = 1;
	unsigned nbits = fg_nbits(logn);
	off += trim_i8_decode(logn, sk + off, f, nbits);
	off += trim_i8_decode(logn, sk + off, g, nbits);
	off += trim_i8_decode(logn, sk + off, F, 8);
	if (off + 64 != FNDSA_SIGN_KEY_SIZE(logn)) {
		fail("masked key MQ test decode");
	}
	for (size_t i = 0; i < n; i ++) {
		int32_t sv[3] = { f[i], g[i], F[i] };
		mmod32 *qv[3] = { basis.f_q, basis.g_q, basis.F_q };
		mmod32 *q2v[3] = { basis.f_q2, basis.g_q2, basis.F_q2 };
		for (size_t j = 0; j < 3; j ++) {
			if (mmod32_open(&qv[j][i], FNDSA_MASK_Q)
				!= signed_to_mod(sv[j], FNDSA_MASK_Q)
				|| mmod32_open(&q2v[j][i], FNDSA_MASK_Q2)
				!= signed_to_mod(sv[j], FNDSA_MASK_Q2))
			{
				fail("masked key MQ coefficient expansion");
			}
		}
		uint32_t gq = mmod32_open(&basis.G_q[i], FNDSA_MASK_Q);
		uint32_t gq2 = mmod32_open(&basis.G_q2[i], FNDSA_MASK_Q2);
		int32_t sgq = gq > (FNDSA_MASK_Q >> 1)
			? (int32_t)(gq - FNDSA_MASK_Q) : (int32_t)gq;
		int32_t sgq2 = gq2 > (FNDSA_MASK_Q2 >> 1)
			? (int32_t)(gq2 - FNDSA_MASK_Q2) : (int32_t)gq2;
		if (sgq != sgq2 || sgq < -127 || sgq > 127) {
			fail("masked key G expansion");
		}
	}

	/* Stored G must satisfy f*G = g*F modulo Falcon q. */
	if (!mmpoly_ntt(logn, basis.f_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.g_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.F_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.G_q, FNDSA_MASK_MOD_Q, &rng))
	{
		fail("masked key MQ NTT");
	}
	mmpoly_mul_ntt(logn, t0, basis.f_q, basis.G_q,
		FNDSA_MASK_MOD_Q, &rng);
	mmpoly_mul_ntt(logn, t1, basis.g_q, basis.F_q,
		FNDSA_MASK_MOD_Q, &rng);
	mmpoly_sub(logn, t0, t0, t1, FNDSA_MASK_MOD_Q);
	for (size_t i = 0; i < n; i ++) {
		if (mmod32_open(&t0[i], FNDSA_MASK_Q) != 0) {
			fail("masked key NTRU basis relation");
		}
	}

	memset(f, 0, sizeof f);
	memset(g, 0, sizeof g);
	memset(F, 0, sizeof F);
	memset(p, 0, 10 * n * sizeof *p);
	free(p);
}

static void
check_seeded_sign_case(unsigned logn,
	const uint8_t *sk, size_t sk_len,
	const uint8_t *mk, size_t mk_len,
	const uint8_t *vk, size_t vk_len,
	const void *ctx, size_t ctx_len, const char *id,
	const void *hv, size_t hv_len,
	const uint8_t seed[40],
	uint8_t *sig0, uint8_t *sig1,
	void *tmp0, size_t tmp0_len, void *tmp1, size_t tmp1_len)
{
	size_t sig_len = FNDSA_SIGNATURE_SIZE(logn);
	size_t r0 = fndsa_sign_seeded_temp(sk, sk_len,
		ctx, ctx_len, id, hv, hv_len, seed, 40,
		sig0, sig_len, tmp0, tmp0_len);
	size_t r1 = fndsa_sign_masked_seeded_temp(mk, mk_len,
		ctx, ctx_len, id, hv, hv_len, seed, 40,
		sig1, sig_len, tmp1, tmp1_len);
	if (r0 != sig_len || r1 != sig_len || memcmp(sig0, sig1, sig_len) != 0) {
		fail("bit-exact masked signing case");
	}
	if (!fndsa_verify(sig1, sig_len, vk, vk_len,
		ctx, ctx_len, id, hv, hv_len))
	{
		fail("masked signing case verification");
	}
}

static void
test_sign_compat(unsigned logn)
{
	size_t sk_len = FNDSA_SIGN_KEY_SIZE(logn);
	size_t vk_len = FNDSA_VRFY_KEY_SIZE(logn);
	size_t sig_len = FNDSA_SIGNATURE_SIZE(logn);
	size_t mk_len = fndsa_masked_key_size(logn);
	size_t tmp0_len = FNDSA_SIGN_TMP_SIZE(logn);
	size_t tmp1_len = fndsa_masked_tmp_size(logn);
	uint8_t *sk = malloc(sk_len);
	uint8_t *vk = malloc(vk_len);
	uint8_t *mk = malloc(mk_len);
	uint8_t *tmp0 = malloc(tmp0_len);
	uint8_t *tmp1 = malloc(tmp1_len);
	uint8_t *sig0 = malloc(sig_len);
	uint8_t *sig1 = malloc(sig_len);
	if (sk == NULL || vk == NULL || mk == NULL || tmp0 == NULL
		|| tmp1 == NULL || sig0 == NULL || sig1 == NULL)
	{
		fail("allocation");
	}

	uint8_t kseed[48], sseed[40];
	for (size_t i = 0; i < sizeof kseed; i ++) {
		kseed[i] = (uint8_t)(i + 17 * logn);
	}
	for (size_t i = 0; i < sizeof sseed; i ++) {
		sseed[i] = (uint8_t)(0xA5 ^ i ^ logn);
	}
	fndsa_keygen_seeded(logn, kseed, sizeof kseed, sk, vk);
	if (fndsa_masked_key_size(8) != 0
		|| fndsa_masked_key_size(11) != 0
		|| fndsa_masked_tmp_size(8) != 0
		|| fndsa_masked_tmp_size(11) != 0)
	{
		fail("non-standard parameter rejection");
	}
	sk[0] ^= UINT8_C(0x10);
	if (fndsa_masked_key_import_demo(sk, sk_len, mk, mk_len)) {
		fail("private-key tag rejection");
	}
	sk[0] ^= UINT8_C(0x10);
	if (fndsa_masked_key_import_demo(sk, sk_len - 1, mk, mk_len)
		|| fndsa_masked_key_import_demo(sk, sk_len, mk, mk_len - 1))
	{
		fail("masked key import length rejection");
	}
	if (!fndsa_masked_key_import_demo(sk, sk_len, mk, mk_len)) {
		fail("masked key import");
	}
	uint32_t encoded_len = (uint32_t)mk[8]
		| ((uint32_t)mk[9] << 8)
		| ((uint32_t)mk[10] << 16)
		| ((uint32_t)mk[11] << 24);
	if (mk_len != (size_t)72 + vk_len + ((size_t)16 << logn)) {
		fail("masked key v2 size");
	}
	if (memcmp(mk, "FMK2", 4) != 0 || mk[4] != 2 || mk[5] != logn
		|| mk[6] != FNDSA_MASK_SHARES || mk[7] != 1
		|| encoded_len != mk_len
		|| mk[12] != 0 || mk[13] != 0 || mk[14] != 0 || mk[15] != 0
		|| memcmp(mk + 16, vk, vk_len) != 0)
	{
		fail("masked key serialization");
	}
	test_masked_key_mq_expansion(logn,
		sk, sk_len, mk, mk_len, vk, vk_len);
	if (fndsa_sign_masked_seeded_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, NULL, 0, NULL, 0) != sig_len)
	{
		fail("signature size query");
	}

	size_t r0 = fndsa_sign_seeded_temp(sk, sk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig0, sig_len, tmp0, tmp0_len);
	size_t r1 = fndsa_sign_masked_seeded_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len);
	if (r0 != sig_len || r1 != sig_len || memcmp(sig0, sig1, sig_len) != 0) {
		fail("bit-exact masked signing core");
	}
	if (!fndsa_verify(sig1, sig_len, vk, vk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7))
	{
		fail("masked signature verification");

	}
	fndsa_masked_test_configure(19, 0);
	r1 = fndsa_sign_masked_seeded_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len);
	fndsa_masked_test_configure(0, 0);
	if (r1 != sig_len || memcmp(sig0, sig1, sig_len) != 0) {
		fail("mask RNG independence");
	}
	uint8_t max_ctx[255], random_msg[37], prehash[32];
	for (size_t i = 0; i < sizeof max_ctx; i ++) {
		max_ctx[i] = (uint8_t)(3 * i + logn);
	}
	for (size_t i = 0; i < sizeof random_msg; i ++) {
		random_msg[i] = (uint8_t)(0xC7 ^ (11 * i) ^ logn);
	}
	for (size_t i = 0; i < sizeof prehash; i ++) {
		prehash[i] = (uint8_t)(0x5A + 7 * i + logn);
	}
	check_seeded_sign_case(logn, sk, sk_len, mk, mk_len, vk, vk_len,
		NULL, 0, FNDSA_HASH_ID_RAW, NULL, 0, sseed,
		sig0, sig1, tmp0, tmp0_len, tmp1, tmp1_len);
	check_seeded_sign_case(logn, sk, sk_len, mk, mk_len, vk, vk_len,
		max_ctx, sizeof max_ctx, FNDSA_HASH_ID_RAW,
		random_msg, sizeof random_msg, sseed,
		sig0, sig1, tmp0, tmp0_len, tmp1, tmp1_len);
	check_seeded_sign_case(logn, sk, sk_len, mk, mk_len, vk, vk_len,
		"ph", 2, FNDSA_HASH_ID_SHA256, prehash, sizeof prehash, sseed,
		sig0, sig1, tmp0, tmp0_len, tmp1, tmp1_len);

	r1 = fndsa_sign_masked_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sig1, sig_len, tmp1, tmp1_len);
	if (r1 != sig_len || !fndsa_verify(sig1, sig_len, vk, vk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7))
	{
		fail("unseeded masked signing core");
	}

	memset(sig1, 0xA5, sig_len);
	if (fndsa_sign_masked_seeded_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len - 1) != 0
		|| !all_zero(sig1, sig_len))
	{
		fail("short temporary area rejection");
	}

	memset(sig1, 0xA5, sig_len);
	if (fndsa_sign_masked_seeded_temp(mk, mk_len - 1,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len) != 0
		|| !all_zero(sig1, sig_len))
	{
		fail("truncated masked key rejection");
	}

	{
		uint8_t saved[5];
		memcpy(saved, mk, sizeof saved);
		memcpy(mk, "FMK1", 4);
		mk[4] = 1;
		memset(sig1, 0xA5, sig_len);
		if (fndsa_sign_masked_seeded_temp(mk, mk_len,
			"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
			sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len) != 0
			|| !all_zero(sig1, sig_len))
		{
			fail("masked key v1 rejection");
		}
		memcpy(mk, saved, sizeof saved);
	}
	static const size_t header_fields[] = { 0, 4, 5, 6, 7, 12 };
	for (size_t j = 0;
		j < sizeof header_fields / sizeof header_fields[0]; j ++)
	{
		size_t off = header_fields[j];
		mk[off] ^= 1;
		memset(sig1, 0xA5, sig_len);
		if (fndsa_sign_masked_seeded_temp(mk, mk_len,
			"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
			sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len) != 0
			|| !all_zero(sig1, sig_len))
		{
			fail("masked key header rejection");
		}
		mk[off] ^= 1;
	}

	if (logn == 9) {
		fndsa_masked_test_configure(0, 27);
		memset(sig1, 0xA5, sig_len);
		r1 = fndsa_sign_masked_seeded_temp(mk, mk_len,
			"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
			sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len);
		fndsa_masked_test_configure(0, 0);
		if (r1 != 0 || !all_zero(sig1, sig_len)) {
			fail("forced 27-attempt failure cleanup");
		}
	}

	size_t corr_off = 16 + vk_len + 16 + 40;
	uint8_t corr_save[4];
	memcpy(corr_save, mk + corr_off, sizeof corr_save);
	memset(mk + corr_off, 0xFF, sizeof corr_save);
	memset(sig1, 0xA5, sig_len);
	if (fndsa_sign_masked_seeded_temp(mk, mk_len,
		"mask-test", 9, FNDSA_HASH_ID_RAW, "message", 7,
		sseed, sizeof sseed, sig1, sig_len, tmp1, tmp1_len) != 0
		|| !all_zero(sig1, sig_len))
	{
		fail("non-canonical correction rejection");
	}
	memcpy(mk + corr_off, corr_save, sizeof corr_save);
	fndsa_masked_key_clear(mk, mk_len);
	for (size_t i = 0; i < mk_len; i ++) {
		if (mk[i] != 0) {
			fail("masked key clear");
		}
	}

	free(sk);
	free(vk);
	free(mk);
	free(tmp0);
	free(tmp1);
	free(sig0);
	free(sig1);
}

int
main(void)
{
	printf("Test masked LFSR: ");
	fflush(stdout);
	test_lfsr();
	printf("done.\nTest masked Boolean nonlinear boundary: ");
	fflush(stdout);
	test_mask_bool_nonlinear();
	printf("done.\nTest masked 160-bit arithmetic: ");
	fflush(stdout);
	test_mask160_arithmetic();
	printf("done.\nTest masked fixed-point sampler gadgets: ");
	fflush(stdout);
	test_masked_floor_inverse_expm();
	test_masked_ber_exp();
	printf("[SamplerZ-512]");
	fflush(stdout);
	test_masked_sampler_z(9);
	printf("[SamplerZ-1024]");
	fflush(stdout);
	test_masked_sampler_z(10);
	printf("done.\nTest masked division, sqrt and mi32: ");
	fflush(stdout);
	test_masked_div_sqrt_mi32();
	printf("done.\nTest masked FFT, split/merge and LDL: ");
	fflush(stdout);
	test_masked_xpoly();
	printf("done.\nTest masked dynamic ffSampling: [L512][R512][L1024][R1024]");
	fflush(stdout);
	test_masked_dynamic_tree();
	printf(" done.\nTest masked NTT, iNTT, Gram and CRT: ");
	fflush(stdout);
	test_masked_ntt_crt();
	printf("done.\nTest masked range, norm and length acceptance: ");
	fflush(stdout);
	test_masked_acceptance();
	printf("done.\nTest masked modular arithmetic: ");
	fflush(stdout);
	test_modular();
	printf("done.\nTest masked signing core: [512]");
	fflush(stdout);
	test_sign_compat(9);
	printf("[1024]");
	fflush(stdout);
	test_sign_compat(10);
	check_fast_range_stats();
	printf(" done.\n");
	return 0;
}
