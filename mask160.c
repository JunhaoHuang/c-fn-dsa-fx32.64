/*
 * Portable 160-bit masking reference implementation.
 *
 * WARNING: DEMO ONLY.  fndsa_mask_rng_init_demo() deliberately initializes
 * a fixed, predictable LFSR-127 stream.  Replace the randomness provider
 * before making any side-channel security claim.
 */

#include <stddef.h>
#include <string.h>
#if defined FNDSA_MASK_RANGE_CHECKS
#include <stdlib.h>
#endif

#include "mask160.h"
#include "mask_mont32.h"

#if FNDSA_MASK_SHARES != 2
#error "mask160.c implements only the proven two-share (first-order) gadgets"
#endif

#if defined(__GNUC__) || defined(__clang__)
#define MASK_NOINLINE   __attribute__((noinline))
#else
#define MASK_NOINLINE
#endif

static inline void
/* Prevent compiler algebra from recombining values across share boundaries. */
value_barrier_u32(uint32_t *x)
{
#if defined(__GNUC__) || defined(__clang__)
	__asm__ volatile("" : "+r" (*x));
#else
	volatile uint32_t y = *x;
	*x = y;
#endif
}


static void mask_clear(void *buf, size_t len);

/* x^127 + x^64 + 1, stepped 64 times. */
static uint64_t
lfsr_step_u64(fndsa_mask_rng *rng)
{
	uint64_t x = ((rng->s[1] << 1) | (rng->s[0] >> 63))
		^ (rng->s[1] >> 62);
	rng->s[1] = (x ^ rng->s[0]) & UINT64_C(0x7FFFFFFFFFFFFFFF);
	rng->s[0] = x;
	return x;
}

/*
 * Return the next 64 bits of the common random bit stream. If a narrower
 * request left bits in the reservoir, then those bits are consumed first;
 * no generated bit is returned twice or silently discarded.
 */
MASK_NOINLINE uint64_t
fndsa_mask_rng_next_u64(fndsa_mask_rng *rng)
{
	unsigned k = rng->reservoir_len;
	if (k == 0) {
		return lfsr_step_u64(rng);
	}
	if (k == 64) {
		uint64_t x = rng->reservoir;
		rng->reservoir = 0;
		rng->reservoir_len = 0;
		return x;
	}

	uint64_t lo = rng->reservoir;
	uint64_t hi = lfsr_step_u64(rng);
	uint64_t x = lo | (hi << k);
	rng->reservoir = hi >> (64 - k);
	rng->reservoir_len = k;
	return x;
}

void
fndsa_mask_rng_init_demo(fndsa_mask_rng *rng)
{
	/* F0 E1 D2 C3 B4 A5 96 87 78 69 5A 4B 3C 2D 1E 0F */
	rng->s[0] = UINT64_C(0x8796A5B4C3D2E1F0);
	rng->s[1] = UINT64_C(0x0F1E2D3C4B5A6978);
	rng->reservoir = 0;
	rng->reservoir_len = 0;
}

uint32_t
fndsa_mask_rng_next_u32(fndsa_mask_rng *rng)
{
	unsigned k = rng->reservoir_len;
	if (k == 0) {
		rng->reservoir = lfsr_step_u64(rng);
		rng->reservoir_len = 64;
	} else if (k < 32) {
		uint32_t lo = (uint32_t)rng->reservoir;
		uint64_t hi = lfsr_step_u64(rng);
		uint32_t x = lo | (uint32_t)(hi << k);
		rng->reservoir = hi >> (32 - k);
		rng->reservoir_len = 32 + k;
		return x;
	}
	uint32_t x = (uint32_t)rng->reservoir;
	rng->reservoir >>= 32;
	rng->reservoir_len -= 32;
	return x;
}

uint32_t
fndsa_mask_rng_next_bit(fndsa_mask_rng *rng)
{
	if (rng->reservoir_len == 0) {
		rng->reservoir = lfsr_step_u64(rng);
		rng->reservoir_len = 64;
	}
	uint32_t x = (uint32_t)rng->reservoir & 1;
	rng->reservoir >>= 1;
	rng->reservoir_len --;
	return x;
}

/* Return the high 64 bits of x*y using only 32x32->64 products. */
static uint64_t
mul64_hi(uint64_t x, uint64_t y)
{
	uint64_t x0 = (uint32_t)x;
	uint64_t x1 = x >> 32;
	uint64_t y0 = (uint32_t)y;
	uint64_t y1 = y >> 32;
	uint64_t w0 = x0 * y0;
	uint64_t t = x1 * y0 + (w0 >> 32);
	uint64_t w1 = (uint32_t)t;
	uint64_t w2 = t >> 32;
	w1 += x0 * y1;
	return x1 * y1 + w2 + (w1 >> 32);
}

/*
 * floor(2^64/m). The production moduli use compile-time constants. The
 * fallback builds the reciprocal with a fixed 64-step restoring division of
 * 2^64-1; the final correction covers the power-of-two divisors of 2^64.
 * The modulus is public, and no secret-dependent branch is used here.
 */
static uint64_t
barrett_reciprocal(uint32_t m)
{
	switch (m) {
	case UINT32_C(12289):
		return UINT64_C(0x00055538E425E9E0);
	case UINT32_C(18433):
		return UINT64_C(0x00038E2C3F62ACDA);
	case UINT32_C(226523137):
		return UINT64_C(0x00000012F5DC32CD);
	default:
		break;
	}

	uint64_t q = 0;
	uint64_t r = 0;
	for (unsigned i = 0; i < 64; i ++) {
		r = (r << 1) | 1;
		uint64_t t = r - m;
		uint64_t ge = (uint64_t)0 - (uint64_t)(r >= m);
		r = (t & ge) | (r & ~ge);
		q = (q << 1) | (ge & 1);
	}
	q += (uint64_t)((m & (m - 1)) == 0);
	return q;
}

static uint32_t
barrett_reduce_u64(uint64_t x, uint32_t m, uint64_t mu)
{
	uint64_t q = mul64_hi(x, mu);
	uint64_t r = x - q * m;
	uint64_t t = r - m;
	uint64_t ge = (uint64_t)0 - (uint64_t)(r >= m);
	return (uint32_t)((t & ge) | (r & ~ge));
}

/*
 * Portable 32-bit Montgomery arithmetic, following the REDC structure used
 * by Raccoon. Values held by mmod32 are canonical Montgomery residues x*R
 * modulo q, with R = 2^32. The modulus is public and must be odd.
 */

typedef struct {
	uint32_t q;
	uint32_t q0i;  /* -q^(-1) mod 2^32 */
	uint32_t r2;   /* 2^64 mod q */
} mont32_params;

static uint32_t
mont32_add(uint32_t a, uint32_t b, uint32_t q)
{
	if (q <= INT32_MAX) {
		return mask_mont32_add(a, b, q);
	}
	uint64_t z = (uint64_t)a + b;
	uint64_t ge = (uint64_t)0 - (uint64_t)(z >= q);
	z -= (uint64_t)q & ge;
	return (uint32_t)z;
}

static uint32_t
mont32_sub(uint32_t a, uint32_t b, uint32_t q)
{
	if (q <= INT32_MAX) {
		return mask_mont32_sub(a, b, q);
	}
	uint64_t z = (uint64_t)a + q - b;
	uint64_t ge = (uint64_t)0 - (uint64_t)(z >= q);
	z -= (uint64_t)q & ge;
	return (uint32_t)z;
}

/*
 * Return t/R modulo q. For all callers t < q*R; hence the unreduced result is
 * below 2*q and one branch-free conditional subtraction is sufficient.
 * Splitting t into high and low words avoids a 64-bit overflow in t+m*q.
 */
static uint32_t
mont32_redc(uint64_t t, const mont32_params *mp)
{
	uint32_t m = (uint32_t)t * mp->q0i;
	uint64_t u = (t >> 32)
		+ (((uint64_t)(uint32_t)t + (uint64_t)m * mp->q) >> 32);
	uint64_t ge = (uint64_t)0 - (uint64_t)(u >= mp->q);
	u -= (uint64_t)mp->q & ge;
	return (uint32_t)u;
}

static uint32_t
mont32_mul(uint32_t a, uint32_t b, const mont32_params *mp)
{
	if (mp->q <= INT32_MAX) {
		return mask_mont32_cadd(mask_mont32_redc(
			(int64_t)a * b, mp->q, mp->q0i), mp->q);
	}
	return mont32_redc((uint64_t)a * b, mp);
}

static uint32_t
mont32_encode(uint32_t x, const mont32_params *mp)
{
	/* x may be any 32-bit word; x*r2 is still strictly below q*R. */
	if (mp->q <= INT32_MAX && mp->r2 <= (mp->q >> 1)) {
		return mask_mont32_cadd(mask_mont32_redc(
			(int64_t)x * mp->r2, mp->q, mp->q0i), mp->q);
	}
	return mont32_redc((uint64_t)x * mp->r2, mp);
}

static uint32_t
mont32_decode(uint32_t x, const mont32_params *mp)
{
	if (mp->q <= INT32_MAX) {
		return mask_mont32_cadd(
			mask_mont32_redc(x, mp->q, mp->q0i), mp->q);
	}
	return mont32_redc(x, mp);
}

/*
 * Resolve public Montgomery constants. The FN-DSA q/q' and CRT modulus use
 * precomputed constants, as in Raccoon. The fixed-topology fallback handles
 * test/research odd moduli without division, remainder, or secret branches.
 */
static int
mont32_get_params(mont32_params *mp, uint32_t q)
{
	if (q < 3 || (q & 1) == 0) {
		memset(mp, 0, sizeof *mp);
		return 0;
	}
	mp->q = q;
	switch (q) {
	case UINT32_C(12289):
		mp->q0i = UINT32_C(0xF7002FFF);
		mp->r2 = UINT32_C(0x00001620);
		return 1;
	case UINT32_C(18433):
		mp->q0i = UINT32_C(0xEBC047FF);
		mp->r2 = UINT32_C(0x00000326);
		return 1;
	case UINT32_C(226523137):
		mp->q0i = UINT32_C(0xD54077FF);
		mp->r2 = UINT32_C(0x02D3B533);
		return 1;
	case UINT32_C(2147483659):
		mp->q0i = UINT32_C(0xC5D1745D);
		mp->r2 = UINT32_C(0x000001E4);
		return 1;
	case UINT32_C(4294967291):
		mp->q0i = UINT32_C(0xCCCCCCCD);
		mp->r2 = UINT32_C(0x00000019);
		return 1;
	default:
		break;
	}

	/* Newton inversion modulo 2^32: five steps grow 1 correct bit to 32. */
	uint32_t qi = 1;
	for (unsigned i = 0; i < 5; i ++) {
		qi *= UINT32_C(2) - q * qi;
	}
	mp->q0i = (uint32_t)0 - qi;

	/* r2 = 2^64 mod q, with a fixed public 64-step doubling schedule. */
	uint32_t r2 = 1;
	for (unsigned i = 0; i < 64; i ++) {
		r2 = mont32_add(r2, r2, q);
	}
	mp->r2 = r2;
	return 1;
}

static void
mmod32_refresh_mont(mmod32 *a, uint32_t modulus,
	fndsa_mask_rng *rng)
{
	uint32_t r = fndsa_mask_rng_mod(rng, modulus);
	a->share[0] = mont32_add(a->share[0], r, modulus);
	a->share[1] = mont32_sub(a->share[1], r, modulus);
}

uint32_t
fndsa_mask_reduce_u64(uint64_t x, uint32_t modulus)
{
	if (modulus < 2) {
		return 0;
	}
	return barrett_reduce_u64(x, modulus,
		barrett_reciprocal(modulus));
}

uint32_t
/* Uniform rejection sampler for masks modulo a public modulus. */
fndsa_mask_rng_mod(fndsa_mask_rng *rng, uint32_t modulus)
{
	if (modulus < 2) {
		return 0;
	}
	uint64_t mu = barrett_reciprocal(modulus);
	uint64_t lim = (mu >> 32) * modulus;
	uint32_t x;
	do {
		x = fndsa_mask_rng_next_u32(rng);
	} while ((uint64_t)x >= lim);
	return barrett_reduce_u64(x, modulus, mu);
}

static void
random160(uint32_t d[5], fndsa_mask_rng *rng)
{
	for (size_t i = 0; i < 5; i ++) {
		d[i] = fndsa_mask_rng_next_u32(rng);
	}
}

static void
add5(uint32_t d[5], const uint32_t a[5], const uint32_t b[5])
{
	uint64_t cc = 0;
	for (size_t i = 0; i < 5; i ++) {
		uint64_t z = (uint64_t)a[i] + b[i] + cc;
		d[i] = (uint32_t)z;
		cc = z >> 32;
	}
}

static void
sub5(uint32_t d[5], const uint32_t a[5], const uint32_t b[5])
{
	uint64_t cc = 0;
	for (size_t i = 0; i < 5; i ++) {
		uint64_t ai = a[i];
		uint64_t y = (uint64_t)b[i] + cc;
		d[i] = (uint32_t)(ai - y);
		cc = ai < y;
	}
}

static void
neg5(uint32_t d[5], const uint32_t a[5])
{
	uint32_t z[5] = { 0, 0, 0, 0, 0 };
	sub5(d, z, a);
}

static void
mul5(uint32_t d[5], const uint32_t a[5], const uint32_t b[5])
{
	memset(d, 0, 5 * sizeof *d);
	for (size_t i = 0; i < 5; i ++) {
		uint64_t cc = 0;
		for (size_t j = 0; j < 5 - i; j ++) {
			uint64_t z = (uint64_t)a[i] * b[j]
				+ d[i + j] + cc;
			d[i + j] = (uint32_t)z;
			cc = z >> 32;
		}
	}
}

/*
 * First-order Boolean AND on n parallel 32-bit lanes. Two independent
 * input refreshes cut alias/fan-out edges before the four cross-products are
 * evaluated. The output accumulator is seeded with a third fresh word before
 * any product is added; no RHS combines multiple cross-products.
 */
MASK_NOINLINE void
fndsa_mask_bool_and_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *x0, const uint32_t *x1,
	const uint32_t *y0, const uint32_t *y1,
	size_t n, fndsa_mask_rng *rng)
{
	for (size_t i = 0; i < n; i ++) {
		uint32_t rx = fndsa_mask_rng_next_u32(rng);
		uint32_t ry = fndsa_mask_rng_next_u32(rng);
		uint32_t r = fndsa_mask_rng_next_u32(rng);
		uint32_t ax0 = x0[i] ^ rx;
		uint32_t ax1 = x1[i] ^ rx;
		uint32_t ay0 = y0[i] ^ ry;
		uint32_t ay1 = y1[i] ^ ry;
		uint32_t a = r;
		uint32_t b = r;
		uint32_t t;
		value_barrier_u32(&ax0);
		value_barrier_u32(&ax1);
		value_barrier_u32(&ay0);
		value_barrier_u32(&ay1);
		value_barrier_u32(&a);
		value_barrier_u32(&b);
		t = ax0 & ay0;
		value_barrier_u32(&t);
		a ^= t;
		value_barrier_u32(&a);
		t = ax0 & ay1;
		value_barrier_u32(&t);
		b ^= t;
		value_barrier_u32(&b);
		t = ax1 & ay0;
		value_barrier_u32(&t);
		b ^= t;
		value_barrier_u32(&b);
		t = ax1 & ay1;
		value_barrier_u32(&t);
		b ^= t;
		value_barrier_u32(&b);
		z0[i] = a;
		z1[i] = b;
	}
}

/* Refresh a Boolean sharing with one fresh random word per 32-bit slice. */
MASK_NOINLINE void
fndsa_mask_bool_refresh_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *x0, const uint32_t *x1,
	size_t n, fndsa_mask_rng *rng)
{
	for (size_t i = 0; i < n; i ++) {
		uint32_t r = fndsa_mask_rng_next_u32(rng);
		uint32_t a = x0[i] ^ r;
		uint32_t b = x1[i] ^ r;
		value_barrier_u32(&a);
		value_barrier_u32(&b);
		z0[i] = a;
		z1[i] = b;
	}
}

/* Public logical left shift of a Boolean-shared multi-word value. */
static void
lshift_words(uint32_t *d, const uint32_t *a, size_t n, unsigned shift)
{
	size_t ws = shift >> 5;
	unsigned bs = shift & 31;
	for (size_t i = 0; i < n; i ++) {
		uint32_t x = 0;
		if (i >= ws) {
			x = a[i - ws] << bs;
			if (bs != 0 && i > ws) {
				x |= a[i - ws - 1] >> (32 - bs);
			}
		}
		d[i] = x;
	}
}

/*
 * Boolean masked addition with a word-sliced Kogge-Stone prefix network.
 * The topology depends only on the public width. For 160 and 192 bits the
 * carry depth is eight steps, instead of 160 or 192 serial steps.
 */
static MASK_NOINLINE void
bool_add_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *x0, const uint32_t *x1,
	const uint32_t *y0, const uint32_t *y1,
	size_t n, fndsa_mask_rng *rng)
{
	uint32_t po[2][6], p[2][6], g[2][6];
	uint32_t ps[2][6], gs[2][6], t[2][6], u[2][6];

	for (size_t i = 0; i < n; i ++) {
		po[0][i] = p[0][i] = x0[i] ^ y0[i];
		po[1][i] = p[1][i] = x1[i] ^ y1[i];
	}
	fndsa_mask_bool_and_words(g[0], g[1], x0, x1, y0, y1, n, rng);

	for (unsigned shift = 1; shift < 32 * n; shift <<= 1) {
		lshift_words(ps[0], p[0], n, shift);
		lshift_words(ps[1], p[1], n, shift);
		lshift_words(gs[0], g[0], n, shift);
		lshift_words(gs[1], g[1], n, shift);
		fndsa_mask_bool_and_words(t[0], t[1],
			p[0], p[1], gs[0], gs[1], n, rng);
		fndsa_mask_bool_and_words(u[0], u[1],
			p[0], p[1], ps[0], ps[1], n, rng);
		for (size_t i = 0; i < n; i ++) {
			/*
			 * Generate and propagate are disjoint, therefore XOR
			 * implements the prefix OR without another nonlinear
			 * gadget.
			 */
			g[0][i] ^= t[0][i];
			g[1][i] ^= t[1][i];
			p[0][i] = u[0][i];
			p[1][i] = u[1][i];
		}
	}

	lshift_words(gs[0], g[0], n, 1);
	lshift_words(gs[1], g[1], n, 1);
	for (size_t i = 0; i < n; i ++) {
		z0[i] = po[0][i] ^ gs[0][i];
		z1[i] = po[1][i] ^ gs[1][i];
	}
}

/*
 * Two-share, word-sliced A2B used for both the 160- and 192-bit paths.
 *
 * This is the n = 2 specialization of the proven A2B construction recalled
 * by Schneider et al. (SPOG19): independently encode every arithmetic input
 * share as a full Boolean sharing, add those encodings with a t-SNI masked
 * adder, then refresh the Boolean output.  In particular, neither operand is
 * ever represented as (a_i, 0); that unsafe shortcut exposes the sum LSB.
 */
static MASK_NOINLINE void
arith_to_bool_words(uint32_t *z0, uint32_t *z1,
	const uint32_t *a0, const uint32_t *a1,
	size_t n, fndsa_mask_rng *rng)
{
	uint32_t x0[6], x1[6], y0[6], y1[6], t0[6], t1[6];
	for (size_t i = 0; i < n; i ++) {
		uint32_t rx = fndsa_mask_rng_next_u32(rng);
		uint32_t ry = fndsa_mask_rng_next_u32(rng);
		x0[i] = a0[i] ^ rx;
		x1[i] = rx;
		y0[i] = a1[i] ^ ry;
		y1[i] = ry;
		value_barrier_u32(&x0[i]);
		value_barrier_u32(&x1[i]);
		value_barrier_u32(&y0[i]);
		value_barrier_u32(&y1[i]);
	}
	bool_add_words(t0, t1, x0, x1, y0, y1, n, rng);
	fndsa_mask_bool_refresh_words(z0, z1, t0, t1, n, rng);
	mask_clear(x0, sizeof x0);
	mask_clear(x1, sizeof x1);
	mask_clear(y0, sizeof y0);
	mask_clear(y1, sizeof y1);
	mask_clear(t0, sizeof t0);
	mask_clear(t1, sizeof t1);
}

void
mrx160_share_rx(mrx160 *d, const uint32_t x[3], fndsa_mask_rng *rng)
{
	uint32_t v[5];
	v[0] = x[0];
	v[1] = x[1];
	v[2] = x[2];
	v[3] = v[4] = (uint32_t)(0 - (x[2] >> 31));
	random160(d->share[1], rng);
	sub5(d->share[0], v, d->share[1]);
}

void
mrx160_open_rx(uint32_t d[3], const mrx160 *a)
{
	uint32_t v[5];
	add5(v, a->share[0], a->share[1]);
	d[0] = v[0];
	d[1] = v[1];
	d[2] = v[2];
}

void
mrx160_copy(mrx160 *d, const mrx160 *a)
{
	memcpy(d, a, sizeof *d);
}

void
mrx160_add(mrx160 *d, const mrx160 *a, const mrx160 *b)
{
	for (size_t j = 0; j < 2; j ++) {
		add5(d->share[j], a->share[j], b->share[j]);
	}
}

void
mrx160_sub(mrx160 *d, const mrx160 *a, const mrx160 *b)
{
	for (size_t j = 0; j < 2; j ++) {
		sub5(d->share[j], a->share[j], b->share[j]);
	}
}

void
mrx160_neg(mrx160 *d, const mrx160 *a)
{
	/* Swap-and-negate is alias safe through temporaries. */
	mrx160 t;
	neg5(t.share[0], a->share[0]);
	neg5(t.share[1], a->share[1]);
	*d = t;
}

void
mrx160_refresh(mrx160 *a, fndsa_mask_rng *rng)
{
	uint32_t r[5], t0[5], t1[5];
	random160(r, rng);
	add5(t0, a->share[0], r);
	sub5(t1, a->share[1], r);
	memcpy(a->share[0], t0, sizeof t0);
	memcpy(a->share[1], t1, sizeof t1);
}

MASK_NOINLINE void
/* Convert an arithmetic sharing without ever reconstructing an operand. */
mrx160_a2b(mb160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	arith_to_bool_words(d->share[0], d->share[1],
		a->share[0], a->share[1], 5, rng);
}

/*
 * BCZ18 Algorithm 1 specialized to n = 2 shares (the recursive C_2 case).
 * RefreshMasks expands (x0,x1) to the three Boolean shares y0,y1,y2; the
 * two C_1 leaves then reduce to XOR.  Psi(x,y) = (x XOR y) - y is evaluated
 * on five or six 32-bit words, exactly modulo 2^160 or 2^192.
 *
 * BCZ18 proves this construction (n-1)-SNI; with FNDSA_MASK_SHARES = 2 the
 * instantiated claim is first-order SNI, assuming fresh ideal randomness and
 * the software probing model.  The fixed demo LFSR does not meet that premise.
 */
static MASK_NOINLINE void
psi_words(uint32_t *d, const uint32_t *x, const uint32_t *y, size_t n)
{
	uint64_t cc = 0;
	for (size_t i = 0; i < n; i ++) {
		uint64_t xi = x[i] ^ y[i];
		uint64_t yi = (uint64_t)y[i] + cc;
		d[i] = (uint32_t)(xi - yi);
		cc = xi < yi;
		value_barrier_u32(&d[i]);
	}
}

static MASK_NOINLINE void
psi_b2a_words(uint32_t *d0, uint32_t *d1,
	const uint32_t *x0, const uint32_t *x1, size_t n,
	fndsa_mask_rng *rng)
{
	uint32_t r1[6], r2[6], y0[6], y1[6], y2[6];
	uint32_t p[6], z0[6], z1[6];
	for (size_t i = 0; i < n; i ++) {
		r1[i] = fndsa_mask_rng_next_u32(rng);
		r2[i] = fndsa_mask_rng_next_u32(rng);
		y0[i] = x0[i] ^ r1[i] ^ r2[i];
		y1[i] = x1[i] ^ r1[i];
		y2[i] = r2[i];
		value_barrier_u32(&y0[i]);
		value_barrier_u32(&y1[i]);
		value_barrier_u32(&y2[i]);
	}
	psi_words(p, y0, y1, n);
	for (size_t i = 0; i < n; i ++) {
		z0[i] = y0[i] ^ p[i];
		value_barrier_u32(&z0[i]);
	}
	psi_words(z1, y0, y2, n);
	for (size_t i = 0; i < n; i ++) {
		d0[i] = y1[i] ^ y2[i];
		d1[i] = z0[i] ^ z1[i];
		value_barrier_u32(&d0[i]);
		value_barrier_u32(&d1[i]);
	}
	mask_clear(r1, sizeof r1);
	mask_clear(r2, sizeof r2);
	mask_clear(y0, sizeof y0);
	mask_clear(y1, sizeof y1);
	mask_clear(y2, sizeof y2);
	mask_clear(p, sizeof p);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
}

MASK_NOINLINE void
mrx160_b2a(mrx160 *d, const mb160 *a, fndsa_mask_rng *rng)
{
	psi_b2a_words(d->share[0], d->share[1],
		a->share[0], a->share[1], 5, rng);
}

static MASK_NOINLINE void
/*
 * First-order SNI product in Z/(2^160).
 *
 * Each cross-product is safe to expose in the probing model because it uses
 * at most one share of each input. The accumulators are initialized with r
 * and -r before the first product is added, so every partial sum is a uniform
 * ring element. Do not combine cross-products before injecting r.
 */
ring_mul(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	uint32_t p[5], r[5], t0[5], t1[5];

	random160(r, rng);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&r[i]);
	}
	memcpy(t0, r, sizeof t0);
	neg5(t1, r);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&t0[i]);
		value_barrier_u32(&t1[i]);
	}

	mul5(p, a->share[0], b->share[0]);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p[i]);
	}
	add5(t0, t0, p);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&t0[i]);
	}

	mul5(p, a->share[0], b->share[1]);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p[i]);
	}
	add5(t1, t1, p);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&t1[i]);
	}
	mul5(p, a->share[1], b->share[0]);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p[i]);
	}
	add5(t1, t1, p);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&t1[i]);
	}
	mul5(p, a->share[1], b->share[1]);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p[i]);
	}
	add5(t1, t1, p);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&t1[i]);
	}

	memcpy(d->share[0], t0, sizeof t0);
	memcpy(d->share[1], t1, sizeof t1);
	mask_clear(p, sizeof p);
	mask_clear(r, sizeof r);
	mask_clear(t0, sizeof t0);
	mask_clear(t1, sizeof t1);
}

static void
bool_asr64(mb160 *d, const mb160 *a)
{
	memset(d, 0, sizeof *d);
	for (size_t j = 0; j < 2; j ++) {
		d->share[j][0] = a->share[j][2];
		d->share[j][1] = a->share[j][3];
		d->share[j][2] = a->share[j][4];
		uint32_t sw = (uint32_t)(0 - (a->share[j][4] >> 31));
		d->share[j][3] = sw;
		d->share[j][4] = sw;
	}
}

MASK_NOINLINE void
mrx160_mul(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	mrx160 p;
	mb160 bp, bs;
	ring_mul(&p, a, b, rng);

	/* Add the Q32.64 half-ulp before the exact arithmetic shift. */
	uint32_t c[5] = { 0, UINT32_C(0x80000000), 0, 0, 0 };
	add5(p.share[0], p.share[0], c);
	mrx160_a2b(&bp, &p, rng);
	bool_asr64(&bs, &bp);
	mrx160_b2a(d, &bs, rng);
	mrx160_refresh(d, rng);
}

#define PUBLIC_FAST_PAD_WORDS   5u
#define PUBLIC_WIDE_PAD_WORDS   288u

#if defined FNDSA_MASK_RANGE_CHECKS
static void
public_fast_range_check(const mrx160 *a, const uint32_t c[3])
{
	uint32_t ar[3], am[3], cm[3], p[6] = { 0 };
	uint32_t na, nc;
	uint64_t cc;

	mrx160_open_rx(ar, a);
	na = ar[2] >> 31;
	nc = c[2] >> 31;
	cc = na;
	for (size_t i = 0; i < 3; i ++) {
		uint64_t z = (uint64_t)(ar[i] ^ ((uint32_t)0 - na)) + cc;
		am[i] = (uint32_t)z;
		cc = z >> 32;
	}
	cc = nc;
	for (size_t i = 0; i < 3; i ++) {
		uint64_t z = (uint64_t)(c[i] ^ ((uint32_t)0 - nc)) + cc;
		cm[i] = (uint32_t)z;
		cc = z >> 32;
	}
	for (size_t i = 0; i < 3; i ++) {
		cc = 0;
		for (size_t j = 0; j < 3; j ++) {
			uint64_t z = (uint64_t)am[i] * cm[j]
				+ p[i + j] + cc;
			p[i + j] = (uint32_t)z;
			cc = z >> 32;
		}
		p[i + 3] = (uint32_t)cc;
	}
	/* Strict magnitude < 2^159: words 5 and 160..191 are zero,
	   and bit 159 (word 4 bit 31) is zero. */
	if (p[5] != 0 || (p[4] >> 31) != 0) {
		abort();
	}
	mask_clear(ar, sizeof ar);
	mask_clear(am, sizeof am);
	mask_clear(cm, sizeof cm);
	mask_clear(p, sizeof p);
}
#else
#define public_fast_range_check(a, c)   ((void)0)
#endif

static void
public_constant_words(uint32_t *d, size_t n, const uint32_t c[3])
{
	uint32_t sw = (uint32_t)(0 - (c[2] >> 31));
	d[0] = c[0];
	d[1] = c[1];
	d[2] = c[2];
	for (size_t i = 3; i < n; i ++) {
		d[i] = sw;
	}
}

static void
mask_rng_discard_words(fndsa_mask_rng *rng, unsigned count)
{
	for (unsigned i = 0; i < count; i ++) {
		(void)fndsa_mask_rng_next_u32(rng);
	}
}

MASK_NOINLINE void
mrx160_mul_public_fast(mrx160 *d, const mrx160 *a,
	const uint32_t c[3], fndsa_mask_rng *rng)
{
	static const uint32_t half[5] = {
		0, UINT32_C(0x80000000), 0, 0, 0
	};
	uint32_t c5[5];
	mrx160 p;
	mb160 bp, bs;

	public_fast_range_check(a, c);
	public_constant_words(c5, 5, c);
	mul5(p.share[0], a->share[0], c5);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p.share[0][i]);
	}
	mul5(p.share[1], a->share[1], c5);
	for (size_t i = 0; i < 5; i ++) {
		value_barrier_u32(&p.share[1][i]);
	}
	add5(p.share[0], p.share[0], half);
	mrx160_a2b(&bp, &p, rng);
	bool_asr64(&bs, &bp);
	mrx160_b2a(d, &bs, rng);
	mrx160_refresh(d, rng);

	/* Match mrx160_mul(): 285 core words plus five discarded words. */
	mask_rng_discard_words(rng, PUBLIC_FAST_PAD_WORDS);
	mask_clear(c5, sizeof c5);
	mask_clear(&p, sizeof p);
	mask_clear(&bp, sizeof bp);
	mask_clear(&bs, sizeof bs);
}

MASK_NOINLINE void
mrx160_mul_public_u64_fast(mrx160 *d, const mrx160 *a, uint64_t c,
	fndsa_mask_rng *rng)
{
	uint32_t raw[3] = { (uint32_t)c, (uint32_t)(c >> 32), 0 };
	mrx160_mul_public_fast(d, a, raw, rng);
	mask_clear(raw, sizeof raw);
}

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES][6];
} ma192;

typedef struct {
	uint32_t share[FNDSA_MASK_SHARES][6];
} mb192;

static void
random192(uint32_t d[6], fndsa_mask_rng *rng)
{
	for (size_t i = 0; i < 6; i ++) {
		d[i] = fndsa_mask_rng_next_u32(rng);
	}
}

static void
add6(uint32_t d[6], const uint32_t a[6], const uint32_t b[6])
{
	uint64_t cc = 0;
	for (size_t i = 0; i < 6; i ++) {
		uint64_t z = (uint64_t)a[i] + b[i] + cc;
		d[i] = (uint32_t)z;
		cc = z >> 32;
	}
}

static void
sub6(uint32_t d[6], const uint32_t a[6], const uint32_t b[6])
{
	uint64_t cc = 0;
	for (size_t i = 0; i < 6; i ++) {
		uint64_t ai = a[i];
		uint64_t y = (uint64_t)b[i] + cc;
		d[i] = (uint32_t)(ai - y);
		cc = ai < y;
	}
}


static void
mul6(uint32_t d[6], const uint32_t a[6], const uint32_t b[6])
{
	memset(d, 0, 6 * sizeof *d);
	for (size_t i = 0; i < 6; i ++) {
		uint64_t cc = 0;
		for (size_t j = 0; j < 6 - i; j ++) {
			uint64_t z = (uint64_t)a[i] * b[j]
				+ d[i + j] + cc;
			d[i + j] = (uint32_t)z;
			cc = z >> 32;
		}
	}
}


static void
bool_sign_extend_96(mb192 *d, const mb160 *a)
{
	memset(d, 0, sizeof *d);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		memcpy(d->share[j], a->share[j], 3 * sizeof(uint32_t));
		uint32_t sw = (uint32_t)(0 - (a->share[j][2] >> 31));
		d->share[j][3] = sw;
		d->share[j][4] = sw;
		d->share[j][5] = sw;
	}
}

static void
a2b192(mb192 *d, const ma192 *a, fndsa_mask_rng *rng)
{
	arith_to_bool_words(d->share[0], d->share[1],
		a->share[0], a->share[1], 6, rng);
}

static void
b2a192(ma192 *d, const mb192 *a, fndsa_mask_rng *rng)
{
	psi_b2a_words(d->share[0], d->share[1],
		a->share[0], a->share[1], 6, rng);
}

static MASK_NOINLINE void
/* 192-bit instance of the same random-first SNI multiplication gadget. */
ring_mul192(ma192 *d, const ma192 *a, const ma192 *b,
	fndsa_mask_rng *rng)
{
	uint32_t p[6], r[6], t0[6], t1[6];
	static const uint32_t zero[6] = { 0, 0, 0, 0, 0, 0 };

	random192(r, rng);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&r[i]);
	}
	memcpy(t0, r, sizeof t0);
	sub6(t1, zero, r);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&t0[i]);
		value_barrier_u32(&t1[i]);
	}

	mul6(p, a->share[0], b->share[0]);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p[i]);
	}
	add6(t0, t0, p);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&t0[i]);
	}

	mul6(p, a->share[0], b->share[1]);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p[i]);
	}
	add6(t1, t1, p);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&t1[i]);
	}
	mul6(p, a->share[1], b->share[0]);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p[i]);
	}
	add6(t1, t1, p);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&t1[i]);
	}
	mul6(p, a->share[1], b->share[1]);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p[i]);
	}
	add6(t1, t1, p);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&t1[i]);
	}

	memcpy(d->share[0], t0, sizeof t0);
	memcpy(d->share[1], t1, sizeof t1);
	mask_clear(p, sizeof p);
	mask_clear(r, sizeof r);
	mask_clear(t0, sizeof t0);
	mask_clear(t1, sizeof t1);
}

static void
bool_asr64_192_to_160(mb160 *d, const mb192 *a)
{
	memset(d, 0, sizeof *d);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j][0] = a->share[j][2];
		d->share[j][1] = a->share[j][3];
		d->share[j][2] = a->share[j][4];
		uint32_t sw = (uint32_t)(0 - (a->share[j][4] >> 31));
		d->share[j][3] = sw;
		d->share[j][4] = sw;
	}
}

static void
mask_clear(void *buf, size_t len)
{
	volatile uint8_t *d = (volatile uint8_t *)buf;
	while (len -- > 0) {
		*d ++ = 0;
	}
}

MASK_NOINLINE void
/*
 * Exact fallback: convert to Boolean form before sign-extending the 96-bit
 * payload. Sign-extending each arithmetic share independently is invalid.
 */
mrx160_mul_wide(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	mb160 ab, bb, out;
	mb192 awb, bwb, pb;
	ma192 awa, bwa, p;

	mrx160_a2b(&ab, a, rng);
	mrx160_a2b(&bb, b, rng);
	bool_sign_extend_96(&awb, &ab);
	bool_sign_extend_96(&bwb, &bb);
	b2a192(&awa, &awb, rng);
	b2a192(&bwa, &bwb, rng);
	ring_mul192(&p, &awa, &bwa, rng);

	/* Existing rx_mul adds 2^63 before its arithmetic right shift. */
	static const uint32_t half[6] = {
		0, UINT32_C(0x80000000), 0, 0, 0, 0
	};
	add6(p.share[0], p.share[0], half);
	a2b192(&pb, &p, rng);
	bool_asr64_192_to_160(&out, &pb);
	mrx160_b2a(d, &out, rng);
	mrx160_refresh(d, rng);

	mask_clear(&ab, sizeof ab);
	mask_clear(&bb, sizeof bb);
	mask_clear(&out, sizeof out);
	mask_clear(&awb, sizeof awb);
	mask_clear(&bwb, sizeof bwb);
	mask_clear(&pb, sizeof pb);
	mask_clear(&awa, sizeof awa);
	mask_clear(&bwa, sizeof bwa);
	mask_clear(&p, sizeof p);
}

MASK_NOINLINE void
mrx160_mul_public_wide(mrx160 *d, const mrx160 *a,
	const uint32_t c[3], fndsa_mask_rng *rng)
{
	static const uint32_t half[6] = {
		0, UINT32_C(0x80000000), 0, 0, 0, 0
	};
	uint32_t c6[6];
	mb160 ab, out;
	mb192 awb, pb;
	ma192 awa, p;

	mrx160_a2b(&ab, a, rng);
	bool_sign_extend_96(&awb, &ab);
	b2a192(&awa, &awb, rng);
	public_constant_words(c6, 6, c);
	mul6(p.share[0], awa.share[0], c6);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p.share[0][i]);
	}
	mul6(p.share[1], awa.share[1], c6);
	for (size_t i = 0; i < 6; i ++) {
		value_barrier_u32(&p.share[1][i]);
	}
	add6(p.share[0], p.share[0], half);
	a2b192(&pb, &p, rng);
	bool_asr64_192_to_160(&out, &pb);
	mrx160_b2a(d, &out, rng);
	mrx160_refresh(d, rng);

	/* Match mrx160_mul_wide(): 621 core words plus 288 discarded words. */
	mask_rng_discard_words(rng, PUBLIC_WIDE_PAD_WORDS);
	mask_clear(c6, sizeof c6);
	mask_clear(&ab, sizeof ab);
	mask_clear(&out, sizeof out);
	mask_clear(&awb, sizeof awb);
	mask_clear(&pb, sizeof pb);
	mask_clear(&awa, sizeof awa);
	mask_clear(&p, sizeof p);
}

/*
 * Square a shared Q32.64 value without aliasing both logical inputs of the
 * multiplication gadget. A direct ring_mul(a, a) would expose cross-products
 * that depend on both shares of the same input. The two independent refreshes
 * below cut both fan-out edges before the regular two-input SNI multiply.
 */
MASK_NOINLINE void
mrx160_sqr(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
        mrx160 x, y;
        mrx160_copy(&x, a);
        mrx160_copy(&y, a);
        mrx160_refresh(&x, rng);
        mrx160_refresh(&y, rng);
        mrx160_mul(d, &x, &y, rng);
        mask_clear(&x, sizeof x);
        mask_clear(&y, sizeof y);
}

/* Extract whether a Boolean-shared 32-bit word is non-zero. */
static void
mbit_nonzero_word(mbit *d, uint32_t x0, uint32_t x1,
	fndsa_mask_rng *rng)
{
	static const unsigned SH[] = { 16, 8, 4, 2, 1 };
	for (size_t i = 0; i < sizeof SH / sizeof SH[0]; i ++) {
		uint32_t y0 = x0 >> SH[i];
		uint32_t y1 = x1 >> SH[i];
		uint32_t z0, z1;
		fndsa_mask_bool_and_words(&z0, &z1, &x0, &x1, &y0, &y1, 1, rng);
		x0 ^= y0 ^ z0;
		x1 ^= y1 ^ z1;
	}
	d->share[0] = x0 & 1;
	d->share[1] = x1 & 1;
}

static void
mbit_not(mbit *d, const mbit *a)
{
	d->share[0] = (a->share[0] ^ 1) & 1;
	d->share[1] = a->share[1] & 1;
}

static void
mb160_lshift1(mb160 *d, const mb160 *a)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			uint32_t x = a->share[j][i] << 1;
			if (i != 0) {
				x |= a->share[j][i - 1] >> 31;
			}
			d->share[j][i] = x;
		}
	}
}

static void
mb160_rshift1(mb160 *d, const mb160 *a)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			uint32_t x = a->share[j][i] >> 1;
			if (i + 1 < FNDSA_MASK160_LIMBS) {
				x |= a->share[j][i + 1] << 31;
			}
			d->share[j][i] = x;
		}
	}
}

/* d <- ctl ? b : a, where ctl is a Boolean-shared bit. */
static void
mb160_select(mb160 *d, const mb160 *a, const mb160 *b,
	const mbit *ctl, fndsa_mask_rng *rng)
{
	uint32_t x0[5], x1[5], m0[5], m1[5], t0[5], t1[5];
	uint32_t c0 = (uint32_t)0 - (ctl->share[0] & 1);
	uint32_t c1 = (uint32_t)0 - (ctl->share[1] & 1);
	for (size_t i = 0; i < 5; i ++) {
		x0[i] = a->share[0][i] ^ b->share[0][i];
		x1[i] = a->share[1][i] ^ b->share[1][i];
		m0[i] = c0;
		m1[i] = c1;
	}
	fndsa_mask_bool_and_words(t0, t1, x0, x1, m0, m1, 5, rng);
	for (size_t i = 0; i < 5; i ++) {
		d->share[0][i] = a->share[0][i] ^ t0[i];
		d->share[1][i] = a->share[1][i] ^ t1[i];
	}
	mask_clear(x0, sizeof x0);
	mask_clear(x1, sizeof x1);
	mask_clear(m0, sizeof m0);
	mask_clear(m1, sizeof m1);
	mask_clear(t0, sizeof t0);
	mask_clear(t1, sizeof t1);
}

static void
mrx160_public_raw(mrx160 *d, const uint32_t x[3])
{
	memset(d, 0, sizeof *d);
	d->share[0][0] = x[0];
	d->share[0][1] = x[1];
	d->share[0][2] = x[2];
	d->share[0][3] = d->share[0][4]
		= (uint32_t)(0 - (x[2] >> 31));
}

MASK_NOINLINE void
/* Drop fractional bits 0..63 and sign-extend the retained integer. */
mrx160_floor(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mb160 x, y;
	mrx160_a2b(&x, a, rng);
	memset(&y, 0, sizeof y);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		y.share[j][2] = x.share[j][2];
		uint32_t sw = (uint32_t)(0 - (x.share[j][2] >> 31));
		y.share[j][3] = sw;
		y.share[j][4] = sw;
	}
	mrx160_b2a(d, &y, rng);
	mrx160_refresh(d, rng);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
}

MASK_NOINLINE void
/*
 * Secretly normalize the denominator into [1/2,1), apply the same power-of-
 * two scaling to numerator 1, then run the legacy five Newton iterations.
 */
mrx160_inv_spec(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	static const uint32_t C32_17[3] = {
		UINT32_C(3789677026), UINT32_C(3789677025), 1
	};
	static const uint32_t C48_17[3] = {
		UINT32_C(3537031891), UINT32_C(3537031890), 2
	};
	static const uint32_t ONE[3] = { 0, 0, 1 };
	mb160 ba, bb, ca, cb;
	mrx160 ua, ub, c48, one, x, y, t;

	mrx160_a2b(&bb, a, rng);
	memset(&ba, 0, sizeof ba);
	ba.share[0][2] = 1;

	/* Secret right normalization for inputs >= 1. */
	for (unsigned i = 0; i < 32; i ++) {
		mbit ctl;
		mbit_nonzero_word(&ctl, bb.share[0][2], bb.share[1][2], rng);
		mb160_rshift1(&ca, &ba);
		mb160_rshift1(&cb, &bb);
		mb160_select(&ba, &ba, &ca, &ctl, rng);
		mb160_select(&bb, &bb, &cb, &ctl, rng);
		mask_clear(&ctl, sizeof ctl);
	}

	/* Secret left normalization for inputs below 1/2. */
	for (unsigned i = 0; i < 64; i ++) {
		mbit top, ctl;
		top.share[0] = (bb.share[0][1] >> 31) & 1;
		top.share[1] = (bb.share[1][1] >> 31) & 1;
		mbit_not(&ctl, &top);
		mb160_lshift1(&ca, &ba);
		mb160_lshift1(&cb, &bb);
		mb160_select(&ba, &ba, &ca, &ctl, rng);
		mb160_select(&bb, &bb, &cb, &ctl, rng);
		mask_clear(&top, sizeof top);
		mask_clear(&ctl, sizeof ctl);
	}

	mrx160_b2a(&ua, &ba, rng);
	mrx160_b2a(&ub, &bb, rng);
	mrx160_public_raw(&c48, C48_17);
	mrx160_public_raw(&one, ONE);

	/*
	 * Range ledger for the fast 160-bit products below: after
	 * normalization, 1/2 <= ub < 1.  The affine seed and all five
	 * Newton iterates stay in (1,2.83), so their raw Q32.64 products
	 * are below 2^132.  For the last ua*x product, the API precondition
	 * 1/a < 2^31 implies the true raw product is below signed 2^159.
	 */
	mrx160_mul_public_fast(&t, &ub, C32_17, rng);
	mrx160_sub(&x, &c48, &t);
	for (unsigned i = 0; i < 5; i ++) {
		mrx160_mul(&t, &x, &ub, rng);
		mrx160_sub(&y, &one, &t);
		mrx160_mul(&t, &y, &x, rng);
		mrx160_add(&x, &x, &t);
	}
	mrx160_mul(d, &ua, &x, rng);
	mrx160_refresh(d, rng);

	mask_clear(&ba, sizeof ba);
	mask_clear(&bb, sizeof bb);
	mask_clear(&ca, sizeof ca);
	mask_clear(&cb, sizeof cb);
	mask_clear(&ua, sizeof ua);
	mask_clear(&ub, sizeof ub);
	mask_clear(&c48, sizeof c48);
	mask_clear(&one, sizeof one);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
	mask_clear(&t, sizeof t);
}

/* Public shifts of Boolean-shared 160-bit words. */
static void
mb160_lshift_public(mb160 *d, const mb160 *a, unsigned n)
{
	unsigned ws = n >> 5;
	unsigned bs = n & 31;
	memset(d, 0, sizeof *d);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			uint32_t x = 0;
			if (i >= ws) {
				x = a->share[j][i - ws] << bs;
				if (bs != 0 && i > ws) {
					x |= a->share[j][i - ws - 1] >> (32 - bs);
				}
			}
			d->share[j][i] = x;
		}
	}
}

static void
mb160_rshift_public(mb160 *d, const mb160 *a, unsigned n,
	int arithmetic)
{
	unsigned ws = n >> 5;
	unsigned bs = n & 31;
	memset(d, 0, sizeof *d);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		uint32_t fill = arithmetic
			? (uint32_t)(0 - (a->share[j][4] >> 31)) : 0;
		for (size_t i = 0; i < FNDSA_MASK160_LIMBS; i ++) {
			size_t si = i + ws;
			uint32_t lo = si < FNDSA_MASK160_LIMBS
				? a->share[j][si] : fill;
			uint32_t hi = si + 1 < FNDSA_MASK160_LIMBS
				? a->share[j][si + 1] : fill;
			uint32_t x = lo >> bs;
			if (bs != 0) {
				x |= hi << (32 - bs);
			}
			d->share[j][i] = x;
		}
	}
}

static void
mrx160_half_internal(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mb160 x, y;
	mrx160_a2b(&x, a, rng);
	mb160_rshift_public(&y, &x, 1, 1);
	mrx160_b2a(d, &y, rng);
	mrx160_refresh(d, rng);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
}

/*
 * Divide a signed Q32.64 numerator by a positive secret denominator.
 * Boolean barrel normalization is followed by the legacy five Newton rounds.
 */
MASK_NOINLINE void
mrx160_div_spec(mrx160 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	static const uint32_t C32_17[3] = {
		UINT32_C(3789677026), UINT32_C(3789677025), 1
	};
	static const uint32_t C48_17[3] = {
		UINT32_C(3537031891), UINT32_C(3537031890), 2
	};
	static const uint32_t ONE[3] = { 0, 0, 1 };
	static const unsigned RSH[] = { 16, 8, 4, 2, 1 };
	static const unsigned LSH[] = { 32, 16, 8, 4, 2, 1 };
	mb160 ba, bb, ca, cb;
	mrx160 ua, ub, c48, one, x, y, t;

	mrx160_a2b(&ba, a, rng);
	mrx160_a2b(&bb, b, rng);

	/* If b >= 1, shift both operands right by bit_length(floor(b)). */
	for (size_t i = 0; i < sizeof RSH / sizeof RSH[0]; i ++) {
		unsigned s = RSH[i];
		mbit ctl;
		mbit_nonzero_word(&ctl, bb.share[0][2] >> s,
			bb.share[1][2] >> s, rng);
		mb160_rshift_public(&ca, &ba, s, 1);
		mb160_rshift_public(&cb, &bb, s, 0);
		mb160_select(&ba, &ba, &ca, &ctl, rng);
		mb160_select(&bb, &bb, &cb, &ctl, rng);
		mask_clear(&ctl, sizeof ctl);
	}
	{
		mbit ctl;
		mbit_nonzero_word(&ctl, bb.share[0][2], bb.share[1][2], rng);
		mb160_rshift_public(&ca, &ba, 1, 1);
		mb160_rshift_public(&cb, &bb, 1, 0);
		mb160_select(&ba, &ba, &ca, &ctl, rng);
		mb160_select(&bb, &bb, &cb, &ctl, rng);
		mask_clear(&ctl, sizeof ctl);
	}

	/* If b < 1/2, left-normalize both operands with six public stages. */
	for (size_t i = 0; i < sizeof LSH / sizeof LSH[0]; i ++) {
		unsigned s = LSH[i];
		unsigned r = 32 - s;
		mbit nz, ctl;
		mbit_nonzero_word(&nz, bb.share[0][1] >> r,
			bb.share[1][1] >> r, rng);
		mbit_not(&ctl, &nz);
		mb160_lshift_public(&ca, &ba, s);
		mb160_lshift_public(&cb, &bb, s);
		mb160_select(&ba, &ba, &ca, &ctl, rng);
		mb160_select(&bb, &bb, &cb, &ctl, rng);
		mask_clear(&nz, sizeof nz);
		mask_clear(&ctl, sizeof ctl);
	}

	mrx160_b2a(&ua, &ba, rng);
	mrx160_b2a(&ub, &bb, rng);
	mrx160_public_raw(&c48, C48_17);
	mrx160_public_raw(&one, ONE);
	mrx160_mul_public_fast(&t, &ub, C32_17, rng);
	mrx160_sub(&x, &c48, &t);
	for (unsigned i = 0; i < 5; i ++) {
		mrx160_mul(&t, &x, &ub, rng);
		mrx160_sub(&y, &one, &t);
		mrx160_mul(&t, &y, &x, rng);
		mrx160_add(&x, &x, &t);
	}
	mrx160_mul_wide(d, &ua, &x, rng);
	mrx160_refresh(d, rng);

	mask_clear(&ba, sizeof ba);
	mask_clear(&bb, sizeof bb);
	mask_clear(&ca, sizeof ca);
	mask_clear(&cb, sizeof cb);
	mask_clear(&ua, sizeof ua);
	mask_clear(&ub, sizeof ub);
	mask_clear(&c48, sizeof c48);
	mask_clear(&one, sizeof one);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
	mask_clear(&t, sizeof t);
}

/* Eight-round masked inverse-square-root path, bit-exact with rx_sqrt_spec(). */
MASK_NOINLINE void
mrx160_sqrt_spec(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	static const uint32_t THREE[3] = { 0, 0, 3 };
	const uint64_t x0 = UINT64_C(11805916207174113034);
	uint32_t xr[3] = { (uint32_t)x0, (uint32_t)(x0 >> 32), 0 };
	mb160 ab, bs;
	mrx160 b, x, y, three, xs;

	mrx160_a2b(&ab, a, rng);
	mb160_rshift_public(&bs, &ab, 12, 0);
	mrx160_b2a(&b, &bs, rng);
	mrx160_public_raw(&x, xr);
	mrx160_public_raw(&three, THREE);
	for (unsigned i = 0; i < 8; i ++) {
		mrx160_sqr(&y, &x, rng);
		mrx160_mul(&y, &b, &y, rng);
		mrx160_sub(&y, &three, &y);
		mrx160_half_internal(&y, &y, rng);
		mrx160_mul(&x, &y, &x, rng);
	}
	mrx160_a2b(&ab, &x, rng);
	mb160_rshift_public(&bs, &ab, 6, 0);
	mrx160_b2a(&xs, &bs, rng);
	mrx160_mul_wide(d, a, &xs, rng);
	mrx160_refresh(d, rng);

	mask_clear(xr, sizeof xr);
	mask_clear(&ab, sizeof ab);
	mask_clear(&bs, sizeof bs);
	mask_clear(&b, sizeof b);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
	mask_clear(&three, sizeof three);
	mask_clear(&xs, sizeof xs);
}

MASK_NOINLINE void
/* Return exact product bits 64..127 without the Q32.64 rounding step. */
mrx160_mul_u64_hi(mb64 *d, const mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	mrx160 p;
	mb160 bp;
	ring_mul(&p, a, b, rng);
	mrx160_a2b(&bp, &p, rng);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j][0] = bp.share[j][2];
		d->share[j][1] = bp.share[j][3];
	}
	mask_clear(&p, sizeof p);
	mask_clear(&bp, sizeof bp);
}

uint32_t
mbit_open_compat(const mbit *a)
{
	return (a->share[0] ^ a->share[1]) & 1;
}

uint64_t
mb64_open_compat(const mb64 *a)
{
	uint32_t lo = a->share[0][0] ^ a->share[1][0];
	uint32_t hi = a->share[0][1] ^ a->share[1][1];
	return (uint64_t)lo | ((uint64_t)hi << 32);
}

void
mmod32_share(mmod32 *d, uint32_t x, uint32_t modulus,
	fndsa_mask_rng *rng)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	uint32_t r = fndsa_mask_rng_mod(rng, modulus);
	uint32_t xr = mont32_encode(x, &mp);
	d->share[1] = r;
	d->share[0] = mont32_sub(xr, r, modulus);
}

uint32_t
mmod32_open(const mmod32 *a, uint32_t modulus)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		return 0;
	}
	return mont32_decode(
		mont32_add(a->share[0], a->share[1], modulus), &mp);
}

void
mmod32_add(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus)
{
	if (modulus < 3 || (modulus & 1) == 0) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	for (size_t i = 0; i < 2; i ++) {
		d->share[i] = mont32_add(a->share[i], b->share[i], modulus);
	}
}

void
mmod32_sub(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus)
{
	if (modulus < 3 || (modulus & 1) == 0) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	for (size_t i = 0; i < 2; i ++) {
		d->share[i] = mont32_sub(a->share[i], b->share[i], modulus);
	}
}

MASK_NOINLINE void
/*
 * First-order SNI multiplication in the Montgomery domain for a public odd
 * modulus. It embeds a 32-bit analogue of Raccoon's Montgomery REDC in this
 * implementation's random-first masked accumulator.
 *
 * A Montgomery product uses at most one share of each input. Both
 * accumulators are masked before their first product, and each cross-product
 * is reduced and accumulated separately. No intermediate contains several
 * unmasked cross-products from different share pairs.
 */
mmod32_mul(mmod32 *d, const mmod32 *a, const mmod32 *b,
	uint32_t modulus, fndsa_mask_rng *rng)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	uint32_t r = fndsa_mask_rng_mod(rng, modulus);
	uint32_t z0 = r;
	uint32_t z1 = mont32_sub(0, r, modulus);
	uint32_t p;
	value_barrier_u32(&r);
	value_barrier_u32(&z0);
	value_barrier_u32(&z1);

	p = mont32_mul(a->share[0], b->share[0], &mp);
	value_barrier_u32(&p);
	z0 = mont32_add(z0, p, modulus);
	value_barrier_u32(&z0);

	p = mont32_mul(a->share[0], b->share[1], &mp);
	value_barrier_u32(&p);
	z1 = mont32_add(z1, p, modulus);
	value_barrier_u32(&z1);

	p = mont32_mul(a->share[1], b->share[0], &mp);
	value_barrier_u32(&p);
	z1 = mont32_add(z1, p, modulus);
	value_barrier_u32(&z1);

	p = mont32_mul(a->share[1], b->share[1], &mp);
	value_barrier_u32(&p);
	z1 = mont32_add(z1, p, modulus);
	value_barrier_u32(&z1);

	d->share[0] = z0;
	d->share[1] = z1;
	mmod32_refresh_mont(d, modulus, rng);
}

void
mmod32_refresh(mmod32 *a, uint32_t modulus, fndsa_mask_rng *rng)
{
	if (modulus < 3 || (modulus & 1) == 0) {
		a->share[0] = 0;
		a->share[1] = 0;
		return;
	}
	mmod32_refresh_mont(a, modulus, rng);
}

uint32_t
mmod32_encode_public(uint32_t x, uint32_t modulus)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		return 0;
	}
	return mont32_encode(x, &mp);
}


void
mmod32_set_residue_shares(mmod32 *d, uint32_t s0, uint32_t s1,
	uint32_t modulus, fndsa_mask_rng *rng)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	d->share[0] = mont32_encode(s0, &mp);
	d->share[1] = mont32_encode(s1, &mp);
	mmod32_refresh_mont(d, modulus, rng);
}

void
mmod32_decode_residue_shares(uint32_t d[FNDSA_MASK_SHARES],
	const mmod32 *a, uint32_t modulus)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d[0] = 0;
		d[1] = 0;
		return;
	}
	d[0] = mont32_decode(a->share[0], &mp);
	d[1] = mont32_decode(a->share[1], &mp);
}

void
mmod32_mul_public_mont(mmod32 *d, const mmod32 *a,
	uint32_t b_mont, uint32_t modulus)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	b_mont = fndsa_mask_reduce_u64(b_mont, modulus);
	d->share[0] = mont32_mul(a->share[0], b_mont, &mp);
	d->share[1] = mont32_mul(a->share[1], b_mont, &mp);
}

void
mi32_share(mi32 *d, int32_t x, fndsa_mask_rng *rng)
{
	uint32_t r = fndsa_mask_rng_next_u32(rng);
	d->share[1] = r;
	d->share[0] = (uint32_t)x - r;
}

int32_t
mi32_open_compat(const mi32 *a)
{
	return (int32_t)(a->share[0] + a->share[1]);
}

/*
 * Sign-extend a compact shared integer into integral Q32.64.
 * The 32-bit arithmetic shares are first converted with the same independently
 * re-shared word-sliced A2B gadget as the 160/192-bit paths; neither operand is
 * ever represented as (a_i, 0).
 */
MASK_NOINLINE void
mi32_to_mrx160(mrx160 *d, const mi32 *a, fndsa_mask_rng *rng)
{
	uint32_t z0[1], z1[1];
	mb160 b;

	arith_to_bool_words(z0, z1,
		&a->share[0], &a->share[1], 1, rng);
	memset(&b, 0, sizeof b);
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		uint32_t z = j == 0 ? z0[0] : z1[0];
		uint32_t sw = (uint32_t)(0 - (z >> 31));
		b.share[j][2] = z;
		b.share[j][3] = sw;
		b.share[j][4] = sw;
	}
	mrx160_b2a(d, &b, rng);
	mrx160_refresh(d, rng);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
	mask_clear(&b, sizeof b);
}

/* Extract the integral Q32.64 word into a compact arithmetic 32-bit sharing. */
MASK_NOINLINE void
mrx160_to_mi32(mi32 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mb160 b;
	uint32_t z0[1], z1[1];
	mrx160_a2b(&b, a, rng);
	psi_b2a_words(z0, z1, &b.share[0][2], &b.share[1][2], 1, rng);
	d->share[0] = z0[0];
	d->share[1] = z1[0];
	mask_clear(&b, sizeof b);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
}

/* Convert one Boolean-shared bit into arithmetic Montgomery shares. */
static MASK_NOINLINE void
mbit_to_mmod32(mmod32 *d, const mbit *a, uint32_t modulus,
	fndsa_mask_rng *rng)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}

	/* Correlated Boolean/arithmetic sharings of an independent random bit. */
	uint32_t r = fndsa_mask_rng_next_bit(rng);
	uint32_t rb0 = fndsa_mask_rng_next_bit(rng);
	uint32_t rb1 = rb0 ^ r;
	uint32_t ar1 = fndsa_mask_rng_mod(rng, modulus);
	uint32_t ar0 = mont32_sub(r, ar1, modulus);
	uint32_t c = ((a->share[0] & 1) ^ rb0)
		^ ((a->share[1] & 1) ^ rb1);

	/* c is a one-time-padded public bit: c ? (1-r) : r. */
	uint32_t nr0 = mont32_sub(1, ar0, modulus);
	uint32_t nr1 = mont32_sub(0, ar1, modulus);
	uint32_t ctl = (uint32_t)0 - (c & 1);
	ar0 ^= ctl & (ar0 ^ nr0);
	ar1 ^= ctl & (ar1 ^ nr1);
	d->share[0] = mont32_encode(ar0, &mp);
	d->share[1] = mont32_encode(ar1, &mp);
	mmod32_refresh_mont(d, modulus, rng);

	value_barrier_u32(&d->share[0]);
	value_barrier_u32(&d->share[1]);
	mask_clear(&r, sizeof r);
	mask_clear(&rb0, sizeof rb0);
	mask_clear(&rb1, sizeof rb1);
	mask_clear(&ar0, sizeof ar0);
	mask_clear(&ar1, sizeof ar1);
	mask_clear(&nr0, sizeof nr0);
	mask_clear(&nr1, sizeof nr1);
	mask_clear(&ctl, sizeof ctl);
}

MASK_NOINLINE void
mi32_to_mmod32(mmod32 *d, const mi32 *a, uint32_t modulus,
	fndsa_mask_rng *rng)
{
	mont32_params mp;
	if (!mont32_get_params(&mp, modulus)) {
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}

	/* Obtain Boolean sharings of S=a0+a1, its carry, and signed low word. */
	uint32_t a0[2] = { a->share[0], 0 };
	uint32_t a1[2] = { a->share[1], 0 };
	uint32_t z0[2], z1[2];
	arith_to_bool_words(z0, z1, a0, a1, 2, rng);
	mbit carry = { { z0[1] & 1, z1[1] & 1 } };
	mbit sign = { { (z0[0] >> 31) & 1, (z1[0] >> 31) & 1 } };

	/* a0+a1 - (carry+sign)*2^32 is the represented signed integer. */
	d->share[0] = mont32_encode(a->share[0], &mp);
	d->share[1] = mont32_encode(a->share[1], &mp);
	uint32_t k = fndsa_mask_reduce_u64(UINT64_C(0x100000000), modulus);
	uint32_t km = mont32_encode(k, &mp);
	mmod32 bc, bs, tc, ts;
	mbit_to_mmod32(&bc, &carry, modulus, rng);
	mbit_to_mmod32(&bs, &sign, modulus, rng);
	mmod32_mul_public_mont(&tc, &bc, km, modulus);
	mmod32_mul_public_mont(&ts, &bs, km, modulus);
	mmod32_sub(d, d, &tc, modulus);
	mmod32_sub(d, d, &ts, modulus);
	mmod32_refresh_mont(d, modulus, rng);

	mask_clear(a0, sizeof a0);
	mask_clear(a1, sizeof a1);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
	mask_clear(&carry, sizeof carry);
	mask_clear(&sign, sizeof sign);
	mask_clear(&bc, sizeof bc);
	mask_clear(&bs, sizeof bs);
	mask_clear(&tc, sizeof tc);
	mask_clear(&ts, sizeof ts);
	mask_clear(&k, sizeof k);
	mask_clear(&km, sizeof km);
}

/* Add a public 32-bit constant to a Boolean sharing. */
static void
mb32_add_public(uint32_t d[2], const uint32_t a[2], uint32_t c,
	fndsa_mask_rng *rng)
{
	uint32_t c1 = 0;
	bool_add_words(&d[0], &d[1], &a[0], &a[1], &c, &c1, 1, rng);
}

/* d <- ctl ? b : a for Boolean-shared 32-bit words. */
static void
mb32_select(uint32_t d[2], const uint32_t a[2], const uint32_t b[2],
	const mbit *ctl, fndsa_mask_rng *rng)
{
	uint32_t x0 = a[0] ^ b[0];
	uint32_t x1 = a[1] ^ b[1];
	uint32_t m0 = (uint32_t)0 - (ctl->share[0] & 1);
	uint32_t m1 = (uint32_t)0 - (ctl->share[1] & 1);
	uint32_t t0, t1;
	fndsa_mask_bool_and_words(&t0, &t1, &x0, &x1, &m0, &m1, 1, rng);
	d[0] = a[0] ^ t0;
	d[1] = a[1] ^ t1;
	mask_clear(&x0, sizeof x0);
	mask_clear(&x1, sizeof x1);
	mask_clear(&m0, sizeof m0);
	mask_clear(&m1, sizeof m1);
	mask_clear(&t0, sizeof t0);
	mask_clear(&t1, sizeof t1);
}

MASK_NOINLINE void
mi32_center_mod_public(mi32 *d, const mi32 *a, uint32_t modulus,
	fndsa_mask_rng *rng)
{
	uint32_t s[2], sm[2], x[2], xh[2], xm[2], y[2];
	uint32_t z0, z1;
	mbit ctl;
	if (modulus < 3 || (modulus & 1) == 0
		|| modulus >= UINT32_C(0x40000000))
	{
		d->share[0] = 0;
		d->share[1] = 0;
		return;
	}
	uint32_t h1 = (modulus + 1) >> 1;

	/* The caller guarantees 0 <= share0+share1 < 2*modulus. */
	arith_to_bool_words(&s[0], &s[1],
		&a->share[0], &a->share[1], 1, rng);

	/* First conditional subtraction: x = s - [s >= modulus]*modulus. */
	mb32_add_public(sm, s, (uint32_t)0 - modulus, rng);
	ctl.share[0] = ((sm[0] >> 31) ^ 1) & 1;
	ctl.share[1] = (sm[1] >> 31) & 1;
	mb32_select(x, s, sm, &ctl, rng);

	/* x > (modulus-1)/2 iff x-h1 is non-negative. */
	mb32_add_public(xh, x, (uint32_t)0 - h1, rng);
	ctl.share[0] = ((xh[0] >> 31) ^ 1) & 1;
	ctl.share[1] = (xh[1] >> 31) & 1;
	mb32_add_public(xm, x, (uint32_t)0 - modulus, rng);
	mb32_select(y, x, xm, &ctl, rng);

	psi_b2a_words(&z0, &z1, &y[0], &y[1], 1, rng);
	uint32_t r = fndsa_mask_rng_next_u32(rng);
	d->share[0] = z0 + r;
	d->share[1] = z1 - r;
	value_barrier_u32(&d->share[0]);
	value_barrier_u32(&d->share[1]);

	mask_clear(s, sizeof s);
	mask_clear(sm, sizeof sm);
	mask_clear(x, sizeof x);
	mask_clear(xh, sizeof xh);
	mask_clear(xm, sizeof xm);
	mask_clear(y, sizeof y);
	mask_clear(&z0, sizeof z0);
	mask_clear(&z1, sizeof z1);
	mask_clear(&r, sizeof r);
	mask_clear(&h1, sizeof h1);
	mask_clear(&ctl, sizeof ctl);
}

void
mi32_center_mod_m(mi32 *d, const mi32 *a, fndsa_mask_rng *rng)
{
	mi32_center_mod_public(d, a, UINT32_C(226523137), rng);
}

static void
mbit_and_internal(mbit *d, const mbit *a, const mbit *b,
	fndsa_mask_rng *rng)
{
	uint32_t z0, z1;
	fndsa_mask_bool_and_words(&z0, &z1,
		&a->share[0], &a->share[1],
		&b->share[0], &b->share[1], 1, rng);
	d->share[0] = z0 & 1;
	d->share[1] = z1 & 1;
	mask_clear(&z0, sizeof z0);
	mask_clear(&z1, sizeof z1);
}

/* Boolean-shared absolute value and the predicate |a| <= 840. */
static void
mi32_abs_b(uint32_t w[2], mbit *in_range, const mi32 *a,
	fndsa_mask_rng *rng)
{
	uint32_t z0, z1, t0, t1, s0, s1;
	arith_to_bool_words(&z0, &z1,
		&a->share[0], &a->share[1], 1, rng);
	s0 = (z0 >> 31) & 1;
	s1 = (z1 >> 31) & 1;
	t0 = z0 ^ ((uint32_t)0 - s0);
	t1 = z1 ^ ((uint32_t)0 - s1);
	bool_add_words(&w[0], &w[1], &t0, &t1, &s0, &s1, 1, rng);

	/* Zero-extend and subtract 841; a borrow means |a| <= 840. */
	uint32_t x0[2] = { w[0], 0 };
	uint32_t x1[2] = { w[1], 0 };
	uint32_t c0[2] = { (uint32_t)0 - UINT32_C(841), UINT32_MAX };
	uint32_t c1[2] = { 0, 0 };
	uint32_t y0[2], y1[2];
	bool_add_words(y0, y1, x0, x1, c0, c1, 2, rng);
	in_range->share[0] = (y0[1] >> 31) & 1;
	in_range->share[1] = (y1[1] >> 31) & 1;

	mask_clear(&z0, sizeof z0);
	mask_clear(&z1, sizeof z1);
	mask_clear(&t0, sizeof t0);
	mask_clear(&t1, sizeof t1);
	mask_clear(&s0, sizeof s0);
	mask_clear(&s1, sizeof s1);
	mask_clear(x0, sizeof x0);
	mask_clear(x1, sizeof x1);
	mask_clear(c0, sizeof c0);
	mask_clear(c1, sizeof c1);
	mask_clear(y0, sizeof y0);
	mask_clear(y1, sizeof y1);
}

static void
mi32_square_add(mrx160 *sum, const mi32 *a, fndsa_mask_rng *rng)
{
	mrx160 x, y;
	mi32_to_mrx160(&x, a, rng);
	mrx160_mul_wide(&y, &x, &x, rng);
	mrx160_add(sum, sum, &y);
	mask_clear(&x, sizeof x);
	mask_clear(&y, sizeof y);
}

static void
mrx160_nonnegative_le_public(mbit *d, const mrx160 *a, uint32_t bound,
	fndsa_mask_rng *rng)
{
	mb160 b;
	uint32_t c0[5] = {
		0, 0, (uint32_t)0 - (bound + 1), UINT32_MAX, UINT32_MAX
	};
	uint32_t c1[5] = { 0, 0, 0, 0, 0 };
	uint32_t z0[5], z1[5];
	mrx160_a2b(&b, a, rng);
	bool_add_words(z0, z1, b.share[0], b.share[1], c0, c1, 5, rng);
	d->share[0] = (z0[4] >> 31) & 1;
	d->share[1] = (z1[4] >> 31) & 1;
	mask_clear(&b, sizeof b);
	mask_clear(c0, sizeof c0);
	mask_clear(c1, sizeof c1);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
}

MASK_NOINLINE void
mi32_poly_sqnorm_le(mbit *accept, unsigned logn,
	const mi32 *a, const mi32 *b, uint32_t bound,
	fndsa_mask_rng *rng)
{
	if (accept == NULL || a == NULL || rng == NULL || logn > 10) {
		if (accept != NULL) {
			accept->share[0] = 0;
			accept->share[1] = 0;
		}
		return;
	}
	mrx160 sum;
	memset(&sum, 0, sizeof sum);
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mi32_square_add(&sum, &a[i], rng);
		if (b != NULL) {
			mi32_square_add(&sum, &b[i], rng);
		}
	}
	mrx160_nonnegative_le_public(accept, &sum, bound, rng);
	mask_clear(&sum, sizeof sum);
}

MASK_NOINLINE void
mi32_signature_is_acceptable(mbit *accept, unsigned logn,
	const mi32 *s0, const mi32 *s1, size_t encoded_len,
	fndsa_mask_rng *rng)
{
	static const uint32_t SQBETA[11] = {
		0, 101498, 208714, 428865, 892039, 1852696,
		3842630, 7959734, 16468416, 34034726, 70265242
	};
	if (accept == NULL || s0 == NULL || s1 == NULL || rng == NULL
		|| (logn != 9 && logn != 10)
		|| encoded_len > (size_t)UINT32_MAX / 8)
	{
		if (accept != NULL) {
			accept->share[0] = 0;
			accept->share[1] = 0;
		}
		return;
	}

	size_t n = (size_t)1 << logn;
	mrx160 sum;
	memset(&sum, 0, sizeof sum);
	mbit range = { { 1, 0 } };
	uint32_t len0[2] = { 0, 0 };
	uint32_t len1[2] = { 0, 0 };
	for (size_t i = 0; i < n; i ++) {
		uint32_t aw[2];
		mbit ok, next;

		mi32_abs_b(aw, &ok, &s0[i], rng);
		mbit_and_internal(&next, &range, &ok, rng);
		range = next;
		mi32_square_add(&sum, &s0[i], rng);
		mask_clear(aw, sizeof aw);
		mask_clear(&ok, sizeof ok);
		mask_clear(&next, sizeof next);

		mi32_abs_b(aw, &ok, &s1[i], rng);
		mbit_and_internal(&next, &range, &ok, rng);
		range = next;
		mi32_square_add(&sum, &s1[i], rng);

		/* comp_encode uses 9 + floor(|s1[i]|/128) bits. */
		uint32_t q0[2] = { aw[0] >> 7, 0 };
		uint32_t q1[2] = { aw[1] >> 7, 0 };
		uint32_t t0[2], t1[2];
		bool_add_words(t0, t1, len0, len1, q0, q1, 2, rng);
		uint32_t c0[2] = { 9, 0 };
		uint32_t c1[2] = { 0, 0 };
		bool_add_words(len0, len1, t0, t1, c0, c1, 2, rng);
		mask_clear(aw, sizeof aw);
		mask_clear(&ok, sizeof ok);
		mask_clear(&next, sizeof next);
		mask_clear(q0, sizeof q0);
		mask_clear(q1, sizeof q1);
		mask_clear(t0, sizeof t0);
		mask_clear(t1, sizeof t1);
		mask_clear(c0, sizeof c0);
		mask_clear(c1, sizeof c1);
	}

	mbit norm_ok, len_ok, both;
	mrx160_nonnegative_le_public(&norm_ok, &sum, SQBETA[logn], rng);
	uint32_t cap = (uint32_t)encoded_len * 8;
	uint32_t c0[2] = { (uint32_t)0 - (cap + 1), UINT32_MAX };
	uint32_t c1[2] = { 0, 0 };
	uint32_t z0[2], z1[2];
	bool_add_words(z0, z1, len0, len1, c0, c1, 2, rng);
	len_ok.share[0] = (z0[1] >> 31) & 1;
	len_ok.share[1] = (z1[1] >> 31) & 1;
	mbit_and_internal(&both, &range, &norm_ok, rng);
	mbit_and_internal(accept, &both, &len_ok, rng);

	mask_clear(&sum, sizeof sum);
	mask_clear(&range, sizeof range);
	mask_clear(&norm_ok, sizeof norm_ok);
	mask_clear(&len_ok, sizeof len_ok);
	mask_clear(&both, sizeof both);
	mask_clear(len0, sizeof len0);
	mask_clear(len1, sizeof len1);
	mask_clear(c0, sizeof c0);
	mask_clear(c1, sizeof c1);
	mask_clear(z0, sizeof z0);
	mask_clear(z1, sizeof z1);
	mask_clear(&cap, sizeof cap);
}

/* Exact arithmetic right shift by one; conversion prevents per-share rounding. */
void
mrx160_half(mrx160 *d, const mrx160 *a, fndsa_mask_rng *rng)
{
	mrx160_half_internal(d, a, rng);
}
