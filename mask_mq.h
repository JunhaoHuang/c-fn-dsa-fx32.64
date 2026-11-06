#ifndef FNDSA_MASK_MQ_H__
#define FNDSA_MASK_MQ_H__

/*
 * Research-only masked polynomial arithmetic modulo the two Falcon primes.
 *
 * Values are two canonical Montgomery shares (R = 2^32). Every linear NTT,
 * iNTT, add/sub/neg, public-pointwise and CRT map is implemented as a full
 * share-0 pass followed by a full share-1 pass. Boundary refresh is separate.
 * Secret Gram/cross products use the four-term masked multiplication gadget;
 * overall centering uses a joint masked compare/select. The demo RNG is not
 * suitable for a production side-channel claim.
 */

#include <stddef.h>
#include <stdint.h>

#include "mask160.h"

#define FNDSA_MASK_Q       UINT32_C(12289)
#define FNDSA_MASK_Q2      UINT32_C(18433)
#define FNDSA_MASK_QQ2     UINT32_C(226523137)
#define FNDSA_MASK_QQ2_H   UINT32_C(113261568)

typedef enum {
	FNDSA_MASK_MOD_Q = 12289,
	FNDSA_MASK_MOD_Q2 = 18433
} fndsa_mask_modulus;

int mmpoly_ntt(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
int mmpoly_intt(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
void mmpoly_refresh(unsigned logn, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);

/* Linear 1-NI operations. They intentionally do not refresh the result. */
void mmpoly_add(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b, fndsa_mask_modulus mod);
void mmpoly_sub(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b, fndsa_mask_modulus mod);
void mmpoly_neg(unsigned logn, mmod32 *d,
	const mmod32 *a, fndsa_mask_modulus mod);

void mmpoly_mul_ntt(unsigned logn, mmod32 *d,
	const mmod32 *a, const mmod32 *b,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
void mmpoly_mul_public_ntt(unsigned logn, mmod32 *d,
	const mmod32 *a, const uint32_t *b_mont,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);


/* Import compact signed shares and perform a public-input NTT. */
void mmpoly_from_mi32(unsigned logn, mmod32 *d, const mi32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
int mmpoly_public_ntt(unsigned logn, uint32_t *d, const uint16_t *a,
	fndsa_mask_modulus mod);

/* Public constant/product helpers preserve both shares independently. */
void mmpoly_mul_public_const(unsigned logn, mmod32 *d,
	const mmod32 *a, uint32_t c, fndsa_mask_modulus mod,
	fndsa_mask_rng *rng);
void mmpoly_public_sub_ntt(unsigned logn, mmod32 *d,
	const uint32_t *a_public_mont, const mmod32 *b,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);

/* iNTT, Montgomery decode, and one joint masked centering operation. */
int mmpoly_intt_center(unsigned logn, mi32 *d, mmod32 *a,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
/*
 * Secret Gram/cross products cannot be evaluated as only share0*share0 and
 * share1*share1: mmod32_mul() supplies both cross terms under fresh masks.
 * d must not overlap any input for these adjoint products.
 */
void mmpoly_gram_ntt(unsigned logn, mmod32 *d,
	const mmod32 *f, const mmod32 *g,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);
void mmpoly_cross_ntt(unsigned logn, mmod32 *d,
	const mmod32 *F, const mmod32 *f,
	const mmod32 *G, const mmod32 *g,
	fndsa_mask_modulus mod, fndsa_mask_rng *rng);

/*
 * CRT mapping runs independently on the two shares. Centering then operates
 * jointly on their masked sum; centering each share separately is incorrect.
 */
void mmpoly_crt_center(unsigned logn, mi32 *d,
	const mmod32 *a_q, const mmod32 *a_q2,
	fndsa_mask_rng *rng);

/* Caller-owned arrays used by the masked-key expansion hook. */
typedef struct {
	mmod32 *f_q;
	mmod32 *g_q;
	mmod32 *F_q;
	mmod32 *G_q;
	mmod32 *f_q2;
	mmod32 *g_q2;
	mmod32 *F_q2;
	mmod32 *G_q2;
} fndsa_masked_basis_mq;

int fndsa_masked_key_expand_mq_internal(
	const void *masked_key, size_t masked_key_len,
	fndsa_masked_basis_mq *basis,
	const uint8_t **vrfy_key, size_t *vrfy_key_len,
	fndsa_mask_rng *rng);

#endif
