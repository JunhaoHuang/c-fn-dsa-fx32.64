/*
 * Two-share Q32.64 dynamic ffSampling backend.
 *
 * This mirrors sign_sampler.c's left/right recursion and public traversal.
 * Secret leaves, Gram values, centers and sampled integers remain shared.
 * Only the already documented SamplerZ/BerExp compatibility decisions open.
 */

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#include "mask_tree.h"

#if defined(__GNUC__) || defined(__clang__)
#define MASK_TREE_NOINLINE __attribute__((noinline))
#else
#define MASK_TREE_NOINLINE
#endif

static const uint32_t RX_IQ_M[3] = {
	UINT32_C(3827689952), UINT32_C(349496), 0
};
static const uint32_t RX_QQ_M[3] = { 0, 0, UINT32_C(151019521) };
static const uint32_t RX_INV_SIGMA_M[11][3] = {
	{ 0, 0, 0 },
	{ UINT32_C(3795937681), UINT32_C(29658807), 0 },
	{ UINT32_C(1221258620), UINT32_C(29249701), 0 },
	{ UINT32_C(166159018), UINT32_C(28857070), 0 },
	{ UINT32_C(1018012490), UINT32_C(28296685), 0 },
	{ UINT32_C(2540388031), UINT32_C(27767728), 0 },
	{ UINT32_C(749324652), UINT32_C(27267368), 0 },
	{ UINT32_C(3857871696), UINT32_C(26793116), 0 },
	{ UINT32_C(1458413264), UINT32_C(26342780), 0 },
	{ UINT32_C(3863772309), UINT32_C(25914413), 0 },
	{ UINT32_C(1662103065), UINT32_C(25506287), 0 }
};

static void
tree_clear(void *buf, size_t len)
{
	volatile uint8_t *p = (volatile uint8_t *)buf;
	while (len -- > 0) {
		*p ++ = 0;
	}
}

static void
tree_public_rx(mrx160 *d, const uint32_t x[3])
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
tree_mul_public(mrx160 *d, const mrx160 *a, const uint32_t x[3],
	fndsa_mask_rng *rng)
{
	mrx160_mul_public_fast(d, a, x, rng);
}

size_t
fndsa_masked_dyn_tmp_size(unsigned logn)
{
	if (logn != 9 && logn != 10) {
		return 0;
	}
	size_t n = (size_t)1 << logn;
	return (2 * n) * sizeof(mrx160) + n * sizeof(mi32) + 32;
}

int
fndsa_mask_tree_bind(fndsa_mask_tree_state *st, unsigned logn,
	const fndsa_mask_sampler_source *source, fndsa_mask_rng *mask_rng,
	void *tmp, size_t tmp_len)
{
	size_t need = fndsa_masked_dyn_tmp_size(logn);
	if (st == NULL || source == NULL || source->next_u8 == NULL
		|| source->logn != logn || mask_rng == NULL || tmp == NULL
		|| need == 0 || tmp_len < need)
	{
		if (st != NULL) {
			memset(st, 0, sizeof *st);
		}
		return 0;
	}
	uintptr_t p = ((uintptr_t)tmp + 31u) & ~(uintptr_t)31u;
	size_t n = (size_t)1 << logn;
	size_t used = (size_t)(p - (uintptr_t)tmp)
		+ (2 * n) * sizeof(mrx160) + n * sizeof(mi32);
	if (used > tmp_len) {
		memset(st, 0, sizeof *st);
		return 0;
	}
	memset(st, 0, sizeof *st);
	st->source = *source;
	st->mask_rng = mask_rng;
	st->work = (mrx160 *)(void *)p;
	st->dz = (mi32 *)(void *)(st->work + 2 * n);
	st->tmp_base = tmp;
	st->tmp_len = tmp_len;
	st->logn = logn;
	memset(st->work, 0, (2 * n) * sizeof(mrx160) + n * sizeof(mi32));
	return 1;
}

mrx160 *
fndsa_mask_tree_work(fndsa_mask_tree_state *st)
{
	return st == NULL ? NULL : st->work;
}

mi32 *
fndsa_mask_tree_samples(fndsa_mask_tree_state *st)
{
	return st == NULL ? NULL : st->dz;
}

static void
mask_ffsampling_deepest(fndsa_mask_tree_state *st,
	unsigned flags, mrx160 *tmp)
{
	mrx160 leaf, t, z;
	mrx160_sqrt_spec(&leaf, &tmp[2], st->mask_rng);
	tree_mul_public(&leaf, &leaf, RX_INV_SIGMA_M[st->logn],
		st->mask_rng);

	unsigned j = -- st->off;
	j = ((j & 0x1C0u) >> 6) | (j & 0x038u) | ((j & 0x007u) << 6);
	j = ((j & 0x124u) >> 2) | (j & 0x092u) | ((j & 0x049u) << 2);
	j >>= (10 - st->logn);
	unsigned j2 = j + (1u << (st->logn - 1));

	mi32_to_mrx160(&t, &st->dz[j], st->mask_rng);
	tree_mul_public(&t, &t, RX_IQ_M, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		mrx160_add(&t, &t, &tmp[0]);
	}
	mrx160_sampler_z(&z, &st->source, &t, &leaf, st->mask_rng);
	mrx160_to_mi32(&st->dz[j], &z, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_NOOUT) == 0) {
		mrx160_sub(&tmp[0], &t, &z);
	}
	tree_clear(&t, sizeof t);
	tree_clear(&z, sizeof z);

	mi32_to_mrx160(&t, &st->dz[j2], st->mask_rng);
	tree_mul_public(&t, &t, RX_IQ_M, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		mrx160_add(&t, &t, &tmp[1]);
	}
	mrx160_sampler_z(&z, &st->source, &t, &leaf, st->mask_rng);
	mrx160_to_mi32(&st->dz[j2], &z, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_NOOUT) == 0) {
		mrx160_sub(&tmp[1], &t, &z);
	}
	tree_clear(&leaf, sizeof leaf);
	tree_clear(&t, sizeof t);
	tree_clear(&z, sizeof z);
}

static MASK_TREE_NOINLINE void
mask_ffsampling_left_inner(fndsa_mask_tree_state *st,
	unsigned logn, unsigned flags, mrx160 *tmp)
{
	if (logn == 1) {
		mask_ffsampling_deepest(st, flags, tmp);
		return;
	}
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t qn = n >> 2;
	mrx160 *tmp0 = tmp;
	mrx160 *tmp6 = tmp + (n >> 1);
	mrx160 *tmp9 = tmp + (3 * n >> 2);
	mrx160 *tmp12 = tmp + n;
	mrx160 *tmp18 = tmp + (3 * n >> 1);
	mrx160 *tmp21 = tmp + (7 * n >> 2);

	mxpoly_split_selfadj_fft(logn, tmp18, tmp21, tmp12, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		mxpoly_split_fft(logn, tmp0, tmp12, st->mask_rng);
	}
	memcpy(tmp6, tmp18, qn * sizeof *tmp6);
	mxpoly_ldl_tree_fft(logn - 1, tmp18, tmp9, tmp6, tmp21,
		st->mask_rng);

	mask_ffsampling_left_inner(st, logn - 1,
		flags & ~FNDSA_MASK_TREE_NOOUT, tmp12);

	memcpy(tmp18, tmp6, hn * sizeof *tmp18);
	memcpy(tmp6, tmp12, hn * sizeof *tmp6);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		memcpy(tmp12, tmp0, hn * sizeof *tmp12);
	}
	mxpoly_mul_iXadj_tree_fft(logn - 1, tmp0, tmp6, tmp21,
		st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) != 0) {
		memcpy(tmp12, tmp0, hn * sizeof *tmp12);
	} else {
		mxpoly_add(logn - 1, tmp12, tmp0);
	}

	mask_ffsampling_left_inner(st, logn - 1,
		flags & ~FNDSA_MASK_TREE_INIT0, tmp12);

	if ((flags & FNDSA_MASK_TREE_NOOUT) == 0) {
		mxpoly_sub(logn - 1, tmp12, tmp0);
		memcpy(tmp18, tmp6, hn * sizeof *tmp18);
		mxpoly_merge_fft(logn, tmp0, tmp12, tmp18, st->mask_rng);
	}
}

static MASK_TREE_NOINLINE void
mask_ffsampling_right_inner(fndsa_mask_tree_state *st,
	unsigned logn, unsigned flags, mrx160 *tmp)
{
	if (logn == 1) {
		mrx160 q;
		tree_public_rx(&q, RX_QQ_M);
		mrx160_div_spec(&tmp[2], &q, &tmp[2], st->mask_rng);
		tree_clear(&q, sizeof q);
		mask_ffsampling_deepest(st, flags, tmp);
		return;
	}
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	mrx160 *tmp0 = tmp;
	mrx160 *tmp6 = tmp + (n >> 1);
	mrx160 *tmp9 = tmp + (3 * n >> 2);
	mrx160 *tmp12 = tmp + n;
	mrx160 *tmp18 = tmp + (3 * n >> 1);
	mrx160 *tmp21 = tmp + (7 * n >> 2);

	mxpoly_split_selfadj_fft(logn, tmp18, tmp21, tmp12, st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		mxpoly_split_fft(logn, tmp0, tmp12, st->mask_rng);
	}
	mxpoly_ldl_tree_fft(logn - 1, tmp6, tmp9, tmp18, tmp21,
		st->mask_rng);
	mxpoly_neg(logn - 2, tmp9);

	mask_ffsampling_right_inner(st, logn - 1,
		flags & ~FNDSA_MASK_TREE_NOOUT, tmp12);

	memcpy(tmp18, tmp6, hn * sizeof *tmp18);
	memcpy(tmp6, tmp12, hn * sizeof *tmp6);
	if ((flags & FNDSA_MASK_TREE_INIT0) == 0) {
		memcpy(tmp12, tmp0, hn * sizeof *tmp12);
	}
	mxpoly_mul_iXadj_tree_fft(logn - 1, tmp0, tmp6, tmp21,
		st->mask_rng);
	if ((flags & FNDSA_MASK_TREE_INIT0) != 0) {
		memcpy(tmp12, tmp0, hn * sizeof *tmp12);
	} else {
		mxpoly_add(logn - 1, tmp12, tmp0);
	}

	mask_ffsampling_right_inner(st, logn - 1,
		flags & ~FNDSA_MASK_TREE_INIT0, tmp12);

	if ((flags & FNDSA_MASK_TREE_NOOUT) == 0) {
		mxpoly_sub(logn - 1, tmp12, tmp0);
		memcpy(tmp18, tmp6, hn * sizeof *tmp18);
		mxpoly_merge_fft(logn, tmp0, tmp12, tmp18, st->mask_rng);
	}
}

void
fndsa_mask_ffsampling_left(fndsa_mask_tree_state *st, unsigned flags)
{
	if (st == NULL || st->work == NULL || st->mask_rng == NULL) {
		return;
	}
	st->off = 1u << (st->logn - 1);
	mask_ffsampling_left_inner(st, st->logn, flags, st->work);
}

void
fndsa_mask_ffsampling_right(fndsa_mask_tree_state *st, unsigned flags)
{
	if (st == NULL || st->work == NULL || st->mask_rng == NULL) {
		return;
	}
	st->off = 1u << (st->logn - 1);
	mask_ffsampling_right_inner(st, st->logn, flags, st->work);
}

void
fndsa_mask_tree_clear(fndsa_mask_tree_state *st)
{
	if (st == NULL) {
		return;
	}
	if (st->tmp_base != NULL && st->tmp_len != 0) {
		tree_clear(st->tmp_base, st->tmp_len);
	}
	tree_clear(st, sizeof *st);
}
