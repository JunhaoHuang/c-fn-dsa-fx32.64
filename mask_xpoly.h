#ifndef FNDSA_MASK_XPOLY_H__
#define FNDSA_MASK_XPOLY_H__

/*
 * Internal masked polynomial/FFT backend.
 *
 * Arrays use the same split FFT representation and operation order as
 * sign_xpoly.c.  Values remain two-share mrx160 objects throughout.
 * This header is research-only and is not part of the stable FN-DSA API.
 */

#include "mask160.h"

void mxpoly_fft_i32_selfadj(unsigned logn, mrx160 *d, mi32 *a,
	fndsa_mask_rng *rng);
/* scratch must contain at least n + n/2 mi32 elements. */
void mxpoly_fft_i32(unsigned logn, mrx160 *d, const mi32 *a,
	mi32 *scratch, fndsa_mask_rng *rng);

void mxpoly_add(unsigned logn, mrx160 *a, const mrx160 *b);
void mxpoly_sub(unsigned logn, mrx160 *a, const mrx160 *b);
void mxpoly_neg(unsigned logn, mrx160 *a);

void mxpoly_mul_fft(unsigned logn, mrx160 *a, const mrx160 *b,
	fndsa_mask_rng *rng);
void mxpoly_div_selfadj_spec_fft(unsigned logn, mrx160 *a,
	const mrx160 *b, fndsa_mask_rng *rng);
void mxpoly_mul_iXadj_fft(unsigned logn, mrx160 *d,
	const mrx160 *a, const mrx160 *b, fndsa_mask_rng *rng);
void mxpoly_ldl_fft(unsigned logn, mrx160 *d, mrx160 *e,
	const mrx160 *a0, const mrx160 *b1, fndsa_mask_rng *rng);

/*
 * Dynamic-tree-only variants. Their five multiplication classes are
 * covered by the call-site theorems in FalconRange/FFSamplingFastPath.lean;
 * the general helpers above deliberately keep the 192-bit fallback.
 */
void mxpoly_mul_iXadj_tree_fft(unsigned logn, mrx160 *d,
	const mrx160 *a, const mrx160 *b, fndsa_mask_rng *rng);
void mxpoly_ldl_tree_fft(unsigned logn, mrx160 *d, mrx160 *e,
	const mrx160 *a0, const mrx160 *b1, fndsa_mask_rng *rng);

void mxpoly_split_fft(unsigned logn, mrx160 *f, mrx160 *f1,
	fndsa_mask_rng *rng);
void mxpoly_split_selfadj_fft(unsigned logn, mrx160 *a0, mrx160 *b1,
	const mrx160 *a, fndsa_mask_rng *rng);
void mxpoly_merge_fft(unsigned logn, mrx160 *f,
	const mrx160 *f0, const mrx160 *f1, fndsa_mask_rng *rng);


typedef enum {
	FNDSA_MASK_RANGE_LDL_THETA_A0,
	FNDSA_MASK_RANGE_LDL_SCHUR,
	FNDSA_MASK_RANGE_IXADJ_THETA_B,
	FNDSA_MASK_RANGE_IXADJ_MU_A,
	FNDSA_MASK_RANGE_IXADJ_FINAL,
	FNDSA_MASK_RANGE_SITE_COUNT
} fndsa_mask_range_site;

#if defined FNDSA_MASK_RANGE_CHECKS
typedef struct {
	uint64_t calls[FNDSA_MASK_RANGE_SITE_COUNT];
	unsigned max_product_bits[FNDSA_MASK_RANGE_SITE_COUNT];
} fndsa_mask_range_stats;

void fndsa_mask_range_stats_reset(void);
void fndsa_mask_range_stats_get(fndsa_mask_range_stats *stats);
#endif
#endif
