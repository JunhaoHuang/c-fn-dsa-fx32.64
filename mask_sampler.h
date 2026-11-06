#ifndef FNDSA_MASK_SAMPLER_H__
#define FNDSA_MASK_SAMPLER_H__

/*
 * Internal research interfaces for the standalone masked fixed-point
 * sampler gadgets.  They are not part of fndsa.h and are DEMO ONLY.
 * The byte callback is the signature PRNG; fndsa_mask_rng supplies an
 * independent mask stream and must never reuse consumed mask bits.
 */
#include <stdint.h>

#include "mask160.h"

typedef uint8_t (*fndsa_mask_next_u8)(void *ctx);

/*
 * Signature-random byte source. next_u8() must advance exactly once per
 * requested byte; it is deliberately separate from fndsa_mask_rng. logn
 * selects the public FN-DSA-512/1024 sampler constant.
 */
typedef struct {
	void *ctx;
	fndsa_mask_next_u8 next_u8;
	unsigned logn;
} fndsa_mask_sampler_source;

/*
 * Evaluate the legacy FACCT exp(-x)*ccs approximation with the same 13
 * coefficients and truncation order. x and ccs are non-negative Q0.64
 * values in mrx160 containers; d is the Boolean-shared 64-bit threshold.
 */
void mrx160_expm_p63(mb64 *d, const mrx160 *x, const mrx160 *ccs,
	fndsa_mask_rng *rng);

/*
 * Bernoulli exponential test. The threshold and final decision stay
 * Boolean-shared. To preserve the legacy byte stream, equality of each
 * compared byte is opened only to decide whether another byte is needed.
 */
void mrx160_ber_exp(mbit *accept, fndsa_mask_sampler_source *src,
	const mrx160 *x, const mrx160 *ccs, fndsa_mask_rng *rng);

/*
 * Sample an integer from the FN-DSA discrete Gaussian centered at mu.
 * isigma is the masked reciprocal standard deviation, not sigma itself.
 * Candidates remain shared; only the accept bit is opened at the rejection-
 * loop compatibility boundary. The result is integral Q32.64 sharing.
 */
void mrx160_sampler_z(mrx160 *z, fndsa_mask_sampler_source *src,
	const mrx160 *mu, const mrx160 *isigma, fndsa_mask_rng *rng);

#endif
