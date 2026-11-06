#ifndef FNDSA_MASKED_H__
#define FNDSA_MASKED_H__

/*
 * Experimental first-order masked FN-DSA API.
 *
 * WARNING: DEMO ONLY.
 *
 * This interface uses the deterministic LFSR-127 masking source described
 * in the accompanying implementation.  It is intended for reproducible
 * functional experiments, not for production use.  In particular, the
 * fixed LFSR seed is not a cryptographically secure source of masks.
 * The signing entry points execute the two-share q/q-prime, FFT/LDL,
 * dynamic-sampling, final-vector, and aggregate-acceptance data path directly;
 * they do not reconstruct an ordinary private key.  Compatibility openings
 * (BerExp prefix continuation, SamplerZ acceptance, whole-attempt retry and
 * the accepted public s1 vector) remain observable by design, so this demo
 * is not a trace-oblivious or production SNI implementation.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef FNDSA_MASK_TEST_HOOKS
/* Test-only controls; absent from production objects and headers. */
void fndsa_masked_test_configure(
	unsigned mask_rng_skip_words, unsigned forced_reject_attempts);
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define FNDSA_MASK_SHARES   2

/*
 * Size of a seed-compressed masked key, or 0 for an unsupported logn.
 * The v2 object contains a 16-byte header, the encoded public key, a 16-byte
 * share-expansion seed, a 40-byte SHAKE256(f || g) hedging token, and four
 * correction polynomials modulo 226523137.  V1 containers are rejected.
 */
size_t fndsa_masked_key_size(unsigned logn);

/*
 * Temporary area required by the masked signing core.  It holds separated
 * q/q-prime shares, polynomial scratch, the masked dynamic tree, final
 * vectors and a staged signature.  The logn=10 arena is about 217 KiB and
 * is supplied by the caller; no large secret arrays are placed on the stack.
 */
size_t fndsa_masked_tmp_size(unsigned logn);

/*
 * Import an ordinary encoded signing key into the seed-compressed demo
 * format.  Only standard degrees (logn 9 and 10) are accepted.  Import is
 * a trusted provisioning operation and is outside the masking boundary.
 * Returns 1 on success and 0 on validation, size, or overlap failure.
 */
int fndsa_masked_key_import_demo(
	const void *sign_key, size_t sign_key_len,
	void *masked_key, size_t masked_key_len);

/* Explicitly overwrite a masked-key buffer through volatile stores. */
void fndsa_masked_key_clear(void *masked_key, size_t masked_key_len);

/*
 * Experimental masked standard-signature entry points. Their argument order
 * and return convention mirror fndsa_sign_temp() and
 * fndsa_sign_seeded_temp().  Passing sig == NULL returns the maximum
 * signature size after validating the masked-key container. On failure
 * they return 0 and clear the portion of the output that could be observed.
 */
size_t fndsa_sign_masked_temp(
	const void *masked_key, size_t masked_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len);

size_t fndsa_sign_masked_seeded_temp(
	const void *masked_key, size_t masked_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len);

#ifdef __cplusplus
}
#endif

#endif
