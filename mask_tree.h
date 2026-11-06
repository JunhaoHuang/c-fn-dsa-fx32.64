#ifndef FNDSA_MASK_TREE_H__
#define FNDSA_MASK_TREE_H__

/*
 * Internal dynamic ffSampling backend.  This is not a public signing API.
 * The caller supplies one contiguous temporary area; it is split into
 * 2*n mrx160 work values and n compact mi32 samples after 32-byte alignment.
 */

#include <stddef.h>

#include "mask_sampler.h"
#include "mask_xpoly.h"

#define FNDSA_MASK_TREE_INIT0   1u
#define FNDSA_MASK_TREE_NOOUT   2u

typedef struct {
	fndsa_mask_sampler_source source;
	fndsa_mask_rng *mask_rng;
	mrx160 *work;
	mi32 *dz;
	void *tmp_base;
	size_t tmp_len;
	unsigned off;
	unsigned logn;
} fndsa_mask_tree_state;

size_t fndsa_masked_dyn_tmp_size(unsigned logn);

int fndsa_mask_tree_bind(fndsa_mask_tree_state *st, unsigned logn,
	const fndsa_mask_sampler_source *source, fndsa_mask_rng *mask_rng,
	void *tmp, size_t tmp_len);

mrx160 *fndsa_mask_tree_work(fndsa_mask_tree_state *st);
mi32 *fndsa_mask_tree_samples(fndsa_mask_tree_state *st);

void fndsa_mask_ffsampling_left(fndsa_mask_tree_state *st, unsigned flags);
void fndsa_mask_ffsampling_right(fndsa_mask_tree_state *st, unsigned flags);

void fndsa_mask_tree_clear(fndsa_mask_tree_state *st);

#endif
