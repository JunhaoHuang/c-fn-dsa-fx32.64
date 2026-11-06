/*
 * Masked split-FFT and LDL operations for the Q32.64 research backend.
 *
 * The loop structure and arithmetic order deliberately mirror sign_xpoly.c.
 * Public constants use sharewise constant multiplication followed by one
 * masked rounding step. General secret products keep the exact 192-bit
 * fallback. Recursive LDL and iXadj use separate 160-bit helpers only at
 * Lean-proved call sites.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#if defined FNDSA_MASK_RANGE_CHECKS
#include <stdio.h>
#include <stdlib.h>
#endif

#include "sign_inner.h"
#include "mask_xpoly.h"

extern const rx fndsa_RNJ[];
extern const uint64_t fndsa_IRNJ[];
extern const rx fndsa_X_PLUS_1_SQUARED_OVER_X[];
extern const rx fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[];

static void
mpublic_rx(mrx160 *d, const uint32_t x[3])
{
	uint32_t s = (uint32_t)(0 - (x[2] >> 31));
	memset(d, 0, sizeof *d);
	d->share[0][0] = x[0];
	d->share[0][1] = x[1];
	d->share[0][2] = x[2];
	d->share[0][3] = s;
	d->share[0][4] = s;
}

static void
mpublic_i32(mrx160 *d, int32_t x)
{
	uint32_t r[3];
	r[0] = 0;
	r[1] = 0;
	r[2] = (uint32_t)x;
	mpublic_rx(d, r);
}

static void
mmul_public_rx(mrx160 *d, const mrx160 *a, const uint32_t c[3],
	fndsa_mask_rng *rng)
{
	mrx160_mul_public_wide(d, a, c, rng);
}

#if defined FNDSA_MASK_RANGE_CHECKS
static fndsa_mask_range_stats range_stats;

static void
raw96_abs(uint32_t d[3], const uint32_t a[3])
{
	uint32_t neg = a[2] >> 31;
	uint64_t cc = neg;
	uint32_t m = (uint32_t)0 - neg;
	for (size_t i = 0; i < 3; i ++) {
		uint64_t z = (uint64_t)(a[i] ^ m) + cc;
		d[i] = (uint32_t)z;
		cc = z >> 32;
	}
}

/* Test-only exact 96x96->192 oracle. Production builds contain neither
 * reconstruction nor this branch. */
static void
range_check_product(fndsa_mask_range_site site,
	const mrx160 *a, const mrx160 *b)
{
	uint32_t ar[3], br[3], am[3], bm[3], p[6] = { 0 };
	mrx160_open_rx(ar, a);
	mrx160_open_rx(br, b);
	raw96_abs(am, ar);
	raw96_abs(bm, br);
	for (size_t i = 0; i < 3; i ++) {
		uint64_t cc = 0;
		for (size_t j = 0; j < 3; j ++) {
			uint64_t z = (uint64_t)am[i] * bm[j]
				+ p[i + j] + cc;
			p[i + j] = (uint32_t)z;
			cc = z >> 32;
		}
		p[i + 3] = (uint32_t)cc;
	}
	unsigned bits = 0;
	for (size_t i = 6; i > 0; i --) {
		uint32_t w = p[i - 1];
		if (w != 0) {
#if defined(__GNUC__) || defined(__clang__)
			bits = (unsigned)((i - 1) * 32 + 32
				- (size_t)__builtin_clz(w));
#else
			bits = (unsigned)((i - 1) * 32);
			while (w != 0) {
				bits ++;
				w >>= 1;
			}
#endif
			break;
		}
	}
	range_stats.calls[site] ++;
	if (bits > range_stats.max_product_bits[site]) {
		range_stats.max_product_bits[site] = bits;
	}
	if (bits > 159) {
		fprintf(stderr,
			"masked ffSampling range violation: site=%u bits=%u\n",
			(unsigned)site, bits);
		abort();
	}
	memset(ar, 0, sizeof ar); memset(br, 0, sizeof br);
	memset(am, 0, sizeof am); memset(bm, 0, sizeof bm);
	memset(p, 0, sizeof p);
}

void
fndsa_mask_range_stats_reset(void)
{
	memset(&range_stats, 0, sizeof range_stats);
}

void
fndsa_mask_range_stats_get(fndsa_mask_range_stats *stats)
{
	*stats = range_stats;
}
#else
#define range_check_product(site, a, b)   ((void)0)
#endif

/*
 * mrx160_mul() consumes 290 fresh 32-bit words; mrx160_mul_wide()
 * consumes 909. Discard the fixed difference so this optimization preserves
 * the mask-RNG state expected by the previous wide-path call graph. The
 * discarded words are fresh and are never reused as masks.
 */
#define TREE_FAST_RNG_PAD_WORDS   619u
static void
pad_tree_fast_rng(fndsa_mask_rng *rng)
{
	for (unsigned i = 0; i < TREE_FAST_RNG_PAD_WORDS; i ++) {
		(void)fndsa_mask_rng_next_u32(rng);
	}
}

static void
mmul_tree_fast(mrx160 *d, const mrx160 *a, const mrx160 *b,
	unsigned site, fndsa_mask_rng *rng)
{
	(void)site;
	range_check_product((fndsa_mask_range_site)site, a, b);
#if defined FNDSA_MASK_TREE_FORCE_WIDE
	mrx160_mul_wide(d, a, b, rng);
#else
	mrx160_mul(d, a, b, rng);
	pad_tree_fast_rng(rng);
#endif
}

static void
mmul_public_rx_tree(mrx160 *d, const mrx160 *a, const uint32_t c[3],
	unsigned site, fndsa_mask_rng *rng)
{
	mrx160 mc;
	mpublic_rx(&mc, c);
	(void)site;
	range_check_product((fndsa_mask_range_site)site, a, &mc);
#if defined FNDSA_MASK_TREE_FORCE_WIDE
	mrx160_mul_public_wide(d, a, c, rng);
#else
	mrx160_mul_public_fast(d, a, c, rng);
	pad_tree_fast_rng(rng);
#endif
	memset(&mc, 0, sizeof mc);
}
static void
mmul_public_ufrac(mrx160 *d, const mrx160 *a, uint64_t c,
	fndsa_mask_rng *rng)
{
	mrx160_mul_public_u64_fast(d, a, c, rng);
}
static void
mcomplex_mul(mrx160 *dre, mrx160 *dim,
	const mrx160 *are, const mrx160 *aim,
	const mrx160 *bre, const mrx160 *bim, fndsa_mask_rng *rng)
{
	mrx160 m0, m1, m2, m3;
	mrx160_add(&m0, are, aim);
	mrx160_add(&m1, bre, bim);
	mrx160_mul_wide(&m0, &m0, &m1, rng);
	mrx160_mul_wide(&m2, are, bre, rng);
	mrx160_mul_wide(&m3, aim, bim, rng);
	mrx160_sub(&m0, &m0, &m2);
	mrx160_sub(dre, &m2, &m3);
	mrx160_sub(dim, &m0, &m3);
	memset(&m0, 0, sizeof m0);
	memset(&m1, 0, sizeof m1);
	memset(&m2, 0, sizeof m2);
	memset(&m3, 0, sizeof m3);
}

static void
mi32_add_local(mi32 *d, const mi32 *a, const mi32 *b)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j] = a->share[j] + b->share[j];
	}
}

static void
mi32_sub_local(mi32 *d, const mi32 *a, const mi32 *b)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j] = a->share[j] - b->share[j];
	}
}

static void
mi32_double_local(mi32 *d)
{
	for (size_t j = 0; j < FNDSA_MASK_SHARES; j ++) {
		d->share[j] <<= 1;
	}
}

void
mxpoly_fft_i32_selfadj(unsigned logn, mrx160 *d, mi32 *a,
	fndsa_mask_rng *rng)
{
	if (logn == 1) {
		mi32_to_mrx160(&d[0], &a[0], rng);
		return;
	}
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t qn = hn >> 1;

	for (unsigned lk = 0; lk < (logn - 2); lk ++) {
		size_t k = (size_t)1 << lk;
		size_t dk = k << 1;
		for (size_t i = hn - dk; i > 0; i -= dk) {
			for (size_t j = 0; j < k; j ++) {
				mi32_add_local(&a[i + j + k], &a[i + j + k],
					&a[i + j - k]);
			}
		}
		for (size_t j = 0; j < k; j ++) {
			mi32_double_local(&a[j + k]);
		}
	}
	for (size_t j = 0; j < hn; j ++) {
		mi32_to_mrx160(&d[j], &a[j], rng);
	}

	for (size_t j = 0; j < qn; j ++) {
		static const uint32_t SQRT2[3] = {
			UINT32_C(4089235721), UINT32_C(1779033703), 1
		};
		mrx160 z, x;
		mmul_public_rx(&z, &d[j + qn], SQRT2, rng);
		mrx160_copy(&x, &d[j]);
		mrx160_sub(&d[j + qn], &x, &z);
		mrx160_add(&d[j], &x, &z);
		memset(&z, 0, sizeof z);
		memset(&x, 0, sizeof x);
	}

	for (unsigned lm = 1; lm < (logn - 1); lm ++) {
		unsigned lk = logn - 2 - lm;
		size_t k = (size_t)1 << lk;
		size_t m = (size_t)1 << lm;
		for (size_t j = 0; j < m; j ++) {
			size_t mj = j << (lk + 1);
			const uint32_t *prz = fndsa_RNJ[j + m];
			for (size_t i = 0; i < k; i ++) {
				mrx160 z, x;
				mmul_public_rx(&z, &d[mj + i + k], prz, rng);
				mrx160_copy(&x, &d[mj + i]);
				mrx160_sub(&d[mj + i + k], &x, &z);
				mrx160_add(&d[mj + i], &x, &z);
				memset(&z, 0, sizeof z);
				memset(&x, 0, sizeof x);
			}
		}
	}
}

void
mxpoly_fft_i32(unsigned logn, mrx160 *d, const mi32 *a,
	mi32 *scratch, fndsa_mask_rng *rng)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	mi32 *cc = scratch;
	mi32 *aa = scratch + n;

	memcpy(cc, a, n * sizeof *cc);
	aa[0] = cc[0];
	mi32_double_local(&aa[0]);
	for (size_t j = 1; j < hn; j ++) {
		mi32_sub_local(&aa[j], &cc[j], &cc[n - j]);
	}
	mxpoly_fft_i32_selfadj(logn, d, aa, rng);

	aa[0] = cc[hn];
	mi32_double_local(&aa[0]);
	for (size_t j = 1; j < hn; j ++) {
		mi32_add_local(&aa[j], &cc[hn + j], &cc[hn - j]);
	}
	mxpoly_fft_i32_selfadj(logn, d + hn, aa, rng);

	for (size_t j = 0; j < n; j ++) {
		mrx160_half(&d[j], &d[j], rng);
	}
	memset(scratch, 0, (n + hn) * sizeof *scratch);
}

void
mxpoly_add(unsigned logn, mrx160 *a, const mrx160 *b)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mrx160_add(&a[i], &a[i], &b[i]);
	}
}

void
mxpoly_sub(unsigned logn, mrx160 *a, const mrx160 *b)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mrx160_sub(&a[i], &a[i], &b[i]);
	}
}

void
mxpoly_neg(unsigned logn, mrx160 *a)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		mrx160_neg(&a[i], &a[i]);
	}
}

void
mxpoly_mul_fft(unsigned logn, mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	for (size_t i = 0; i < hn; i ++) {
		mrx160 re, im;
		mcomplex_mul(&re, &im, &a[i], &a[i + hn],
			&b[i], &b[i + hn], rng);
		mrx160_copy(&a[i], &re);
		mrx160_copy(&a[i + hn], &im);
		memset(&re, 0, sizeof re);
		memset(&im, 0, sizeof im);
	}
}

void
mxpoly_div_selfadj_spec_fft(unsigned logn, mrx160 *a,
	const mrx160 *b, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	for (size_t i = 0; i < hn; i ++) {
		mrx160 inv, re, im;
		/* Mirror xpoly_div_selfadj_spec_fft exactly: inverse once, then
		 * apply two independently rounded Q32.64 multiplications. */
		mrx160_inv_spec(&inv, &b[i], rng);
		mrx160_mul_wide(&re, &a[i], &inv, rng);
		mrx160_mul_wide(&im, &a[i + hn], &inv, rng);
		mrx160_copy(&a[i], &re);
		mrx160_copy(&a[i + hn], &im);
		memset(&inv, 0, sizeof inv);
		memset(&re, 0, sizeof re);
		memset(&im, 0, sizeof im);
	}
}

void
mxpoly_mul_iXadj_fft(unsigned logn, mrx160 *d,
	const mrx160 *a, const mrx160 *b, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	for (size_t i = 0; i < hn; i ++) {
		mrx160 bt, w, x;
		mmul_public_rx(&bt, &b[i],
			fndsa_X_PLUS_1_SQUARED_OVER_X[hn + i], rng);
		mrx160_half(&bt, &bt, rng);
		mmul_public_rx(&w, &a[i + hn],
			fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[hn + i], rng);
		mrx160_add(&w, &w, &a[i]);
		mrx160_mul_wide(&x, &bt, &w, rng);
		mrx160_copy(&d[i], &x);
		mmul_public_rx(&w, &a[i],
			fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[hn + i], rng);
		mrx160_sub(&w, &a[i + hn], &w);
		mrx160_mul_wide(&x, &bt, &w, rng);
		mrx160_copy(&d[i + hn], &x);
		memset(&bt, 0, sizeof bt);
		memset(&w, 0, sizeof w);
		memset(&x, 0, sizeof x);
	}
}

void
mxpoly_mul_iXadj_tree_fft(unsigned logn, mrx160 *d,
	const mrx160 *a, const mrx160 *b, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	for (size_t i = 0; i < hn; i ++) {
		mrx160 bt, w, x;
		mmul_public_rx_tree(&bt, &b[i],
			fndsa_X_PLUS_1_SQUARED_OVER_X[hn + i],
			FNDSA_MASK_RANGE_IXADJ_THETA_B, rng);
		mrx160_half(&bt, &bt, rng);
		mmul_public_rx_tree(&w, &a[i + hn],
			fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[hn + i],
			FNDSA_MASK_RANGE_IXADJ_MU_A, rng);
		mrx160_add(&w, &w, &a[i]);
		mmul_tree_fast(&x, &bt, &w,
			FNDSA_MASK_RANGE_IXADJ_FINAL, rng);
		mrx160_copy(&d[i], &x);
		mmul_public_rx_tree(&w, &a[i],
			fndsa_X_MINUS_1_OVER_X_PLUS_1_TIMES_I[hn + i],
			FNDSA_MASK_RANGE_IXADJ_MU_A, rng);
		mrx160_sub(&w, &a[i + hn], &w);
		mmul_tree_fast(&x, &bt, &w,
			FNDSA_MASK_RANGE_IXADJ_FINAL, rng);
		mrx160_copy(&d[i + hn], &x);
		memset(&bt, 0, sizeof bt);
		memset(&w, 0, sizeof w);
		memset(&x, 0, sizeof x);
	}
}

void
mxpoly_ldl_fft(unsigned logn, mrx160 *d, mrx160 *e,
	const mrx160 *a0, const mrx160 *b1, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	if (logn == 1) {
		mrx160 t, u;
		mrx160_div_spec(&t, &b1[0], &a0[0], rng);
		mrx160_half(&e[0], &t, rng);
		mrx160_mul_wide(&u, &e[0], &b1[0], rng);
		mrx160_sub(&d[0], &a0[0], &u);
		memset(&t, 0, sizeof t);
		memset(&u, 0, sizeof u);
		return;
	}
	for (size_t i = 0; i < hn; i ++) {
		mrx160 t, u;
		mmul_public_rx(&t, &a0[i],
			fndsa_X_PLUS_1_SQUARED_OVER_X[hn + i], rng);
		mrx160_div_spec(&e[i], &b1[i], &t, rng);
		mrx160_mul_wide(&u, &e[i], &b1[i], rng);
		mrx160_sub(&d[i], &a0[i], &u);
		memset(&t, 0, sizeof t);
		memset(&u, 0, sizeof u);
	}
}

void
mxpoly_ldl_tree_fft(unsigned logn, mrx160 *d, mrx160 *e,
	const mrx160 *a0, const mrx160 *b1, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	if (logn == 1) {
		mrx160 t, u;
		mrx160_div_spec(&t, &b1[0], &a0[0], rng);
		mrx160_half(&e[0], &t, rng);
		mmul_tree_fast(&u, &e[0], &b1[0],
			FNDSA_MASK_RANGE_LDL_SCHUR, rng);
		mrx160_sub(&d[0], &a0[0], &u);
		memset(&t, 0, sizeof t);
		memset(&u, 0, sizeof u);
		return;
	}
	for (size_t i = 0; i < hn; i ++) {
		mrx160 t, u;
		mmul_public_rx_tree(&t, &a0[i],
			fndsa_X_PLUS_1_SQUARED_OVER_X[hn + i],
			FNDSA_MASK_RANGE_LDL_THETA_A0, rng);
		mrx160_div_spec(&e[i], &b1[i], &t, rng);
		mmul_tree_fast(&u, &e[i], &b1[i],
			FNDSA_MASK_RANGE_LDL_SCHUR, rng);
		mrx160_sub(&d[i], &a0[i], &u);
		memset(&t, 0, sizeof t);
		memset(&u, 0, sizeof u);
	}
}

void
mxpoly_split_fft(unsigned logn, mrx160 *f, mrx160 *f1,
	fndsa_mask_rng *rng)
{
	if (logn == 1) {
		mrx160_copy(&f1[0], &f[1]);
		return;
	}
	if (logn == 2) {
		mrx160 tr, ti, ur, ui, t0, t1;
		mrx160_add(&tr, &f[0], &f[1]);
		mrx160_sub(&ur, &f[0], &f[1]);
		mrx160_add(&ti, &f[2], &f[3]);
		mrx160_sub(&ui, &f[2], &f[3]);
		mrx160_half(&f[0], &tr, rng);
		mrx160_half(&f[1], &ti, rng);
		mmul_public_ufrac(&t0, &ur, fndsa_x_SIN[1], rng);
		mmul_public_ufrac(&t1, &ui, fndsa_x_SIN[1], rng);
		mrx160_add(&ur, &t0, &t1);
		mrx160_sub(&ui, &t1, &t0);
		mrx160_half(&f1[0], &ur, rng);
		mrx160_half(&f1[1], &ui, rng);
		memset(&tr, 0, sizeof tr); memset(&ti, 0, sizeof ti);
		memset(&ur, 0, sizeof ur); memset(&ui, 0, sizeof ui);
		memset(&t0, 0, sizeof t0); memset(&t1, 0, sizeof t1);
		return;
	}
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;
	for (size_t i = 0; i < qn; i += 2) {
		size_t j = i << 1;
		mrx160 tr, ti, ur, ui, t0, t1, t2, t3;
		uint64_t sr = fndsa_x_SIN[hn - 1 - (i >> 1)];
		uint64_t si = fndsa_x_SIN[qn + (i >> 1)];

		mrx160_add(&tr, &f[j], &f[j + 1]);
		mrx160_sub(&ur, &f[j], &f[j + 1]);
		mrx160_add(&ti, &f[j + hn], &f[j + 1 + hn]);
		mrx160_sub(&ui, &f[j + hn], &f[j + 1 + hn]);
		mrx160_half(&f[j], &tr, rng);
		mrx160_half(&f[j + hn], &ti, rng);
		mmul_public_ufrac(&t0, &ur, sr, rng);
		mmul_public_ufrac(&t1, &ui, si, rng);
		mmul_public_ufrac(&t2, &ur, si, rng);
		mmul_public_ufrac(&t3, &ui, sr, rng);
		mrx160_add(&ur, &t0, &t1);
		mrx160_sub(&ui, &t3, &t2);
		mrx160_half(&f1[i], &ur, rng);
		mrx160_half(&f1[i + qn], &ui, rng);

		mrx160_add(&tr, &f[j + 2], &f[j + 3]);
		mrx160_sub(&ur, &f[j + 2], &f[j + 3]);
		mrx160_add(&ti, &f[j + 2 + hn], &f[j + 3 + hn]);
		mrx160_sub(&ui, &f[j + 2 + hn], &f[j + 3 + hn]);
		mrx160_half(&f[j + 2], &tr, rng);
		mrx160_half(&f[j + 2 + hn], &ti, rng);
		mmul_public_ufrac(&t0, &ur, si, rng);
		mmul_public_ufrac(&t1, &ui, sr, rng);
		mmul_public_ufrac(&t2, &ur, sr, rng);
		mmul_public_ufrac(&t3, &ui, si, rng);
		mrx160_sub(&ur, &t1, &t0);
		mrx160_add(&ui, &t2, &t3);
		mrx160_neg(&ui, &ui);
		mrx160_half(&f1[i + 1], &ur, rng);
		mrx160_half(&f1[i + 1 + qn], &ui, rng);
		memset(&tr, 0, sizeof tr); memset(&ti, 0, sizeof ti);
		memset(&ur, 0, sizeof ur); memset(&ui, 0, sizeof ui);
		memset(&t0, 0, sizeof t0); memset(&t1, 0, sizeof t1);
		memset(&t2, 0, sizeof t2); memset(&t3, 0, sizeof t3);
	}
	for (size_t i = 1; i < hn; i ++) {
		mrx160_copy(&f[i], &f[i << 1]);
	}
}

void
mxpoly_split_selfadj_fft(unsigned logn, mrx160 *a0, mrx160 *b1,
	const mrx160 *a, fndsa_mask_rng *rng)
{
	if (logn == 1) {
		mrx160_copy(&a0[0], &a[0]);
		mpublic_i32(&b1[0], 0);
		return;
	}
	if (logn == 2) {
		mrx160 t0, t1;
		mrx160_add(&t0, &a[0], &a[1]);
		mrx160_sub(&t1, &a[0], &a[1]);
		mrx160_half(&a0[0], &t0, rng);
		mmul_public_ufrac(&b1[0], &t1, fndsa_IRNJ[1], rng);
		memset(&t0, 0, sizeof t0); memset(&t1, 0, sizeof t1);
		return;
	}
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;
	for (size_t i = 0; i < qn; i += 2) {
		mrx160 t0, t1;
		mrx160_add(&t0, &a[(i << 1)], &a[(i << 1) + 1]);
		mrx160_sub(&t1, &a[(i << 1)], &a[(i << 1) + 1]);
		mrx160_half(&a0[i], &t0, rng);
		mmul_public_ufrac(&b1[i], &t1, fndsa_IRNJ[qn + i], rng);
		mrx160_add(&t0, &a[(i << 1) + 2], &a[(i << 1) + 3]);
		mrx160_sub(&t1, &a[(i << 1) + 3], &a[(i << 1) + 2]);
		mrx160_half(&a0[i + 1], &t0, rng);
		mmul_public_ufrac(&b1[i + 1], &t1,
			fndsa_IRNJ[qn + i + 1], rng);
		memset(&t0, 0, sizeof t0); memset(&t1, 0, sizeof t1);
	}
}

void
mxpoly_merge_fft(unsigned logn, mrx160 *f,
	const mrx160 *f0, const mrx160 *f1, fndsa_mask_rng *rng)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;
	if (logn == 1) {
		mrx160_copy(&f[0], &f0[0]);
		mrx160_copy(&f[hn], &f1[0]);
		return;
	}
	if (logn == 2) {
		mrx160 br, bi;
		mrx160_sub(&br, &f1[0], &f1[1]);
		mrx160_add(&bi, &f1[0], &f1[1]);
		mmul_public_ufrac(&br, &br, fndsa_x_SIN[1], rng);
		mmul_public_ufrac(&bi, &bi, fndsa_x_SIN[1], rng);
		mrx160_add(&f[0], &f0[0], &br);
		mrx160_sub(&f[1], &f0[0], &br);
		mrx160_add(&f[2], &f0[1], &bi);
		mrx160_sub(&f[3], &f0[1], &bi);
		memset(&br, 0, sizeof br); memset(&bi, 0, sizeof bi);
		return;
	}
	for (size_t i = 0; i < qn; i += 2) {
		uint64_t sr = fndsa_x_SIN[hn - 1 - (i >> 1)];
		uint64_t si = fndsa_x_SIN[qn + (i >> 1)];
		mrx160 t0, t1, t2, t3, br, bi;
		mmul_public_ufrac(&t0, &f1[i], sr, rng);
		mmul_public_ufrac(&t1, &f1[i + qn], si, rng);
		mmul_public_ufrac(&t2, &f1[i], si, rng);
		mmul_public_ufrac(&t3, &f1[i + qn], sr, rng);
		mrx160_sub(&br, &t0, &t1);
		mrx160_add(&bi, &t2, &t3);
		mrx160_add(&f[(i << 1)], &f0[i], &br);
		mrx160_sub(&f[(i << 1) + 1], &f0[i], &br);
		mrx160_add(&f[(i << 1) + hn], &f0[i + qn], &bi);
		mrx160_sub(&f[(i << 1) + 1 + hn], &f0[i + qn], &bi);

		mmul_public_ufrac(&t0, &f1[i + 1], sr, rng);
		mmul_public_ufrac(&t1, &f1[i + 1 + qn], si, rng);
		mmul_public_ufrac(&t2, &f1[i + 1], si, rng);
		mmul_public_ufrac(&t3, &f1[i + 1 + qn], sr, rng);
		mrx160_add(&br, &t2, &t3);
		mrx160_sub(&bi, &t0, &t1);
		mrx160_sub(&f[(i << 1) + 2], &f0[i + 1], &br);
		mrx160_add(&f[(i << 1) + 3], &f0[i + 1], &br);
		mrx160_add(&f[(i << 1) + 2 + hn], &f0[i + 1 + qn], &bi);
		mrx160_sub(&f[(i << 1) + 3 + hn], &f0[i + 1 + qn], &bi);
		memset(&t0, 0, sizeof t0); memset(&t1, 0, sizeof t1);
		memset(&t2, 0, sizeof t2); memset(&t3, 0, sizeof t3);
		memset(&br, 0, sizeof br); memset(&bi, 0, sizeof bi);
	}
}
