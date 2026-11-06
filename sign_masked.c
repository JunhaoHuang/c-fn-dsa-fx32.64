/*
 * Experimental seed-compressed masked-key container and signing API.
 *
 * WARNING: DEMO ONLY.
 *
 * The v2 container follows the Raccoon-style seed-plus-correction pattern.
 * Signing expands two separate shares in q and q-prime and carries them
 * through NTT/CRT, FFT/LDL, dynamic sampling, final-vector generation and
 * aggregate acceptance without reconstructing an ordinary private key.
 *
 * The fixed LFSR is predictable and repeats across top-level calls.  The
 * stored SHAKE256(f || g) hedging token and the documented BerExp/SamplerZ/
 * retry openings are explicit compatibility boundaries.  Consequently this
 * remains a reproducible research demo, not a production SNI claim.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "fndsa.h"
#include "fndsa_masked.h"
#include "inner.h"
#include "sign_inner.h"
#include "mask_tree.h"
#include "mask_xpoly.h"
#include "mask160.h"
#include "mask_mq.h"

#define MKEY_VERSION        2
#define MKEY_FLAGS_DEMO     1
#define MKEY_HEADER_LEN    16
#define MKEY_SEED_LEN      16
#define MKEY_HEDGE_LEN     40
#define MKEY_POLYS          4
#define MKEY_MOD    UINT32_C(226523137)
#define MKEY_LIMIT  UINT64_C(0xF3087012)

static const uint8_t MKEY_MAGIC[4] = { 'F', 'M', 'K', '2' };
#ifdef FNDSA_MASK_TEST_HOOKS
static unsigned mask_test_rng_skip;
static unsigned mask_test_forced_rejects;

void
fndsa_masked_test_configure(unsigned mask_rng_skip_words,
	unsigned forced_reject_attempts)
{
	mask_test_rng_skip = mask_rng_skip_words;
	mask_test_forced_rejects = forced_reject_attempts > 27
		? 27 : forced_reject_attempts;
}
#endif

/* Best-effort erasure that the compiler may not optimize into a no-op. */
static void
secure_clear(void *buf, size_t len)
{
	volatile uint8_t *d = (volatile uint8_t *)buf;
	while (len -- > 0) {
		*d ++ = 0;
	}
}

static uint32_t
dec32le(const uint8_t *src)
{
	return (uint32_t)src[0]
		| ((uint32_t)src[1] << 8)
		| ((uint32_t)src[2] << 16)
		| ((uint32_t)src[3] << 24);
}

static void
enc32le(uint8_t *dst, uint32_t x)
{
	dst[0] = (uint8_t)x;
	dst[1] = (uint8_t)(x >> 8);
	dst[2] = (uint8_t)(x >> 16);
	dst[3] = (uint8_t)(x >> 24);
}

static void
enc64le(uint8_t *dst, uint64_t x)
{
	enc32le(dst, (uint32_t)x);
	enc32le(dst + 4, (uint32_t)(x >> 32));
}

size_t
fndsa_masked_key_size(unsigned logn)
{
	if (logn != 9 && logn != 10) {
		return 0;
	}
	return (size_t)72 + FNDSA_VRFY_KEY_SIZE(logn)
		+ ((size_t)16 << logn);
}

size_t
fndsa_masked_tmp_size(unsigned logn)
{
	if (logn != 9 && logn != 10) {
		return 0;
	}
	return (size_t)216 * ((size_t)1 << logn) + 63
		+ FNDSA_SIGNATURE_SIZE(logn);
}

static int
ct_equal(const void *a, const void *b, size_t len)
{
	const uint8_t *x = (const uint8_t *)a;
	const uint8_t *y = (const uint8_t *)b;
	uint32_t r = 0;
	for (size_t i = 0; i < len; i ++) {
		r |= x[i] ^ y[i];
	}
	return r == 0;
}

static uint32_t
shake_next_u32le(shake_context *sc)
{
	uint8_t t[4];
	shake_extract(sc, t, sizeof t);
	return dec32le(t);
}

static void
/*
 * Expand the seeded share of one polynomial. The domain binds the stream
 * to the format version, degree and polynomial number. Rejection sampling
 * makes coefficients uniform modulo MKEY_MOD without using '%'.
 */
expand_share1(unsigned logn, unsigned poly, const uint8_t seed[16],
	uint32_t *d)
{
	uint8_t domain[8] = {
		'F', 'M', 'K', '2', (uint8_t)logn, (uint8_t)poly, 0, 0
	};
	shake_context sc;
	shake_init(&sc, 256);
	shake_inject(&sc, domain, sizeof domain);
	shake_inject(&sc, seed, MKEY_SEED_LEN);
	shake_flip(&sc);
	const uint64_t lim = MKEY_LIMIT;
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		uint32_t x;
		do {
			x = shake_next_u32le(&sc);
		} while ((uint64_t)x >= lim);
		d[i] = fndsa_mask_reduce_u64(x, MKEY_MOD);
	}
	secure_clear(&sc, sizeof sc);
}

static int
decode_private(unsigned logn, const uint8_t *sk,
	int8_t *f, int8_t *g, int8_t *F)
{
	size_t n = (size_t)1 << logn;
	unsigned nbits = fg_nbits(logn);
	size_t off = 1;
	off += trim_i8_decode(logn, sk + off, f, nbits);
	off += trim_i8_decode(logn, sk + off, g, nbits);
	off += trim_i8_decode(logn, sk + off, F, 8);
	if (off + 64 != FNDSA_SIGN_KEY_SIZE(logn)) {
		return 0;
	}
	int fmin = -(1 << (nbits - 1));
	uint32_t bad = 0;
	for (size_t i = 0; i < n; i ++) {
		bad |= (uint32_t)(f[i] == fmin);
		bad |= (uint32_t)(g[i] == fmin);
		bad |= (uint32_t)(F[i] == -128);
	}
	return bad == 0;
}

static int
/*
 * Derive both public h=g/f mod q and the missing fourth basis polynomial
 * G=gF/f mod q. G is centered and range-checked before it is admitted to
 * the masked-key container.
 */
derive_public_and_G(unsigned logn,
	const int8_t *f, const int8_t *g, const int8_t *F,
	uint8_t *vk, int16_t *G)
{
	size_t n = (size_t)1 << logn;
	uint16_t ft[1024], gt[1024], Ft[1024];
	uint16_t h[1024], t[1024], work[1024];
	int ok = 0;

	mqpoly_small_to_int(logn, f, ft);
	mqpoly_small_to_int(logn, g, gt);
	mqpoly_small_to_int(logn, F, Ft);
	mqpoly_int_to_ntt(logn, ft);
	mqpoly_int_to_ntt(logn, gt);
	mqpoly_int_to_ntt(logn, Ft);

	/* h = g/f mod q, in the same form used by key generation. */
	memcpy(h, gt, n * sizeof *h);
	if (!mqpoly_div_ntt(logn, h, ft, work)) {
		goto cleanup;
	}
	mqpoly_int_to_ext(logn, h);
	vk[0] = (uint8_t)logn;
	if (1 + mqpoly_encode(logn, h, vk + 1)
		!= FNDSA_VRFY_KEY_SIZE(logn))
	{
		goto cleanup;
	}

	/* G = g*F/f mod q; valid FN-DSA keys have small centered G. */
	memcpy(t, gt, n * sizeof *t);
	mqpoly_mul_ntt(logn, t, Ft);
	if (!mqpoly_div_ntt(logn, t, ft, work)) {
		goto cleanup;
	}
	mqpoly_ntt_to_int(logn, t);
	mqpoly_int_to_signed(logn, t);
	uint32_t bad = 0;
	for (size_t i = 0; i < n; i ++) {
		int32_t x = (int16_t)t[i];
		bad |= (uint32_t)(x < -127);
		bad |= (uint32_t)(x > 127);
		G[i] = (int16_t)x;
	}
	if (bad != 0) {
		goto cleanup;
	}
	ok = 1;

cleanup:
	secure_clear(ft, sizeof ft);
	secure_clear(gt, sizeof gt);
	secure_clear(Ft, sizeof Ft);
	secure_clear(h, sizeof h);
	secure_clear(t, sizeof t);
	secure_clear(work, sizeof work);
	return ok;
}

int
fndsa_masked_key_import_demo(
	const void *sign_key, size_t sign_key_len,
	void *masked_key, size_t masked_key_len)
{
	if (sign_key == NULL || masked_key == NULL || sign_key_len == 0) {
		return 0;
	}
	const uint8_t *sk = (const uint8_t *)sign_key;
	unsigned logn = sk[0] & 0x0F;
	if ((sk[0] & 0xF0) != 0x50 || (logn != 9 && logn != 10)
		|| sign_key_len != FNDSA_SIGN_KEY_SIZE(logn))
	{
		return 0;
	}
	size_t need = fndsa_masked_key_size(logn);
	if (masked_key_len < need) {
		return 0;
	}

	uintptr_t sa = (uintptr_t)sign_key;
	uintptr_t sb = sa + sign_key_len;
	uintptr_t da = (uintptr_t)masked_key;
	uintptr_t db = da + need;
	if (sa < db && da < sb) {
		return 0;
	}

	size_t n = (size_t)1 << logn;
	int8_t f[1024], g[1024], F[1024];
	int16_t G[1024];
	uint8_t vk[FNDSA_VRFY_KEY_SIZE(10)];
	uint8_t hvk[64];
	uint8_t hedge[MKEY_HEDGE_LEN];
	uint32_t share1[1024];
	int ok = 0;

	if (!decode_private(logn, sk, f, g, F)
		|| !derive_public_and_G(logn, f, g, F, vk, G))
	{
		goto cleanup;
	}
	fndsa_hashed_vrfykey_from_vrfykey(
		hvk, vk, FNDSA_VRFY_KEY_SIZE(logn));
	if (!ct_equal(hvk, sk + sign_key_len - 64, 64)) {
		goto cleanup;
	}
	/* Compatibility token used by make_seedbuf(NULL, ...). */
	shake_context hsc;
	shake_init(&hsc, 256);
	shake_inject(&hsc, sk + 1, fg_nbits(logn) << (logn - 2));
	shake_flip(&hsc);
	shake_extract(&hsc, hedge, sizeof hedge);
	secure_clear(&hsc, sizeof hsc);

	/* Serialize header || public key || expansion seed || corrections. */
	uint8_t *mk = (uint8_t *)masked_key;
	memset(mk, 0, need);
	memcpy(mk, MKEY_MAGIC, sizeof MKEY_MAGIC);
	mk[4] = MKEY_VERSION;
	mk[5] = (uint8_t)logn;
	mk[6] = FNDSA_MASK_SHARES;
	mk[7] = MKEY_FLAGS_DEMO;
	enc32le(mk + 8, (uint32_t)need);
	memcpy(mk + MKEY_HEADER_LEN, vk, FNDSA_VRFY_KEY_SIZE(logn));

	uint8_t *mseed = mk + MKEY_HEADER_LEN + FNDSA_VRFY_KEY_SIZE(logn);
	fndsa_mask_rng rng;
	fndsa_mask_rng_init_demo(&rng);
	enc64le(mseed, fndsa_mask_rng_next_u64(&rng));
	enc64le(mseed + 8, fndsa_mask_rng_next_u64(&rng));

	/*
	 * correction + expanded_share1 == secret modulo MKEY_MOD for every
	 * coefficient; only the correction share is stored explicitly.
	 */
	uint8_t *mhedge = mseed + MKEY_SEED_LEN;
	memcpy(mhedge, hedge, MKEY_HEDGE_LEN);
	uint8_t *corr = mhedge + MKEY_HEDGE_LEN;
	for (unsigned j = 0; j < MKEY_POLYS; j ++) {
		expand_share1(logn, j, mseed, share1);
		for (size_t i = 0; i < n; i ++) {
			int32_t sx;
			switch (j) {
			case 0: sx = f[i]; break;
			case 1: sx = g[i]; break;
			case 2: sx = F[i]; break;
			default: sx = G[i]; break;
			}
			uint32_t neg = (uint32_t)0 - (uint32_t)(sx < 0);
			uint32_t x = (uint32_t)sx + (MKEY_MOD & neg);
			uint32_t c = fndsa_mask_reduce_u64(
				(uint64_t)x + MKEY_MOD - share1[i], MKEY_MOD);
			enc32le(corr + ((size_t)j * n + i) * 4, c);
		}
	}
	secure_clear(&rng, sizeof rng);
	ok = 1;

cleanup:
	if (!ok) {
		secure_clear(masked_key, masked_key_len);
	}
	secure_clear(f, sizeof f);
	secure_clear(g, sizeof g);
	secure_clear(F, sizeof F);
	secure_clear(G, sizeof G);
	secure_clear(vk, sizeof vk);
	secure_clear(hvk, sizeof hvk);
	secure_clear(hedge, sizeof hedge);
	secure_clear(share1, sizeof share1);
	return ok;
}

void
fndsa_masked_key_clear(void *masked_key, size_t masked_key_len)
{
	if (masked_key != NULL) {
		secure_clear(masked_key, masked_key_len);
	}
}

static int
/* Validate container metadata and canonical correction coefficients. */
validate_mkey(const uint8_t *mk, size_t mk_len,
	unsigned *logn_out, const uint8_t **vk_out,
	const uint8_t **seed_out, const uint8_t **hedge_out,
	const uint8_t **corr_out)
{
	if (mk == NULL || mk_len < MKEY_HEADER_LEN
		|| !ct_equal(mk, MKEY_MAGIC, sizeof MKEY_MAGIC)
		|| mk[4] != MKEY_VERSION || mk[6] != FNDSA_MASK_SHARES
		|| mk[7] != MKEY_FLAGS_DEMO
		|| mk[12] != 0 || mk[13] != 0 || mk[14] != 0 || mk[15] != 0)
	{
		return 0;
	}
	unsigned logn = mk[5];
	size_t need = fndsa_masked_key_size(logn);
	if (need == 0 || mk_len != need || dec32le(mk + 8) != need) {
		return 0;
	}
	const uint8_t *vk = mk + MKEY_HEADER_LEN;
	if (vk[0] != logn) {
		return 0;
	}
	const uint8_t *seed = vk + FNDSA_VRFY_KEY_SIZE(logn);
	const uint8_t *hedge = seed + MKEY_SEED_LEN;
	const uint8_t *corr = hedge + MKEY_HEDGE_LEN;
	size_t n = (size_t)1 << logn;
	uint32_t bad = 0;
	for (size_t i = 0; i < MKEY_POLYS * n; i ++) {
		bad |= (uint32_t)(dec32le(corr + 4 * i) >= MKEY_MOD);
	}
	if (bad != 0) {
		return 0;
	}
	*logn_out = logn;
	*vk_out = vk;
	*seed_out = seed;
	*hedge_out = hedge;
	*corr_out = corr;
	return 1;
}

int
fndsa_masked_key_expand_mq_internal(
	const void *masked_key, size_t masked_key_len,
	fndsa_masked_basis_mq *basis,
	const uint8_t **vrfy_key, size_t *vrfy_key_len,
	fndsa_mask_rng *rng)
{
	unsigned logn;
	const uint8_t *vk, *seed, *hedge, *corr;
	uint32_t share1[1024];
	if (vrfy_key != NULL) {
		*vrfy_key = NULL;
	}
	if (vrfy_key_len != NULL) {
		*vrfy_key_len = 0;
	}
	if (basis == NULL || rng == NULL
		|| !validate_mkey((const uint8_t *)masked_key, masked_key_len,
			&logn, &vk, &seed, &hedge, &corr))
	{
		return 0;
	}
	mmod32 *dq[MKEY_POLYS] = {
		basis->f_q, basis->g_q, basis->F_q, basis->G_q
	};
	mmod32 *dq2[MKEY_POLYS] = {
		basis->f_q2, basis->g_q2, basis->F_q2, basis->G_q2
	};
	for (unsigned j = 0; j < MKEY_POLYS; j ++) {
		if (dq[j] == NULL || dq2[j] == NULL) {
			secure_clear(share1, sizeof share1);
			return 0;
		}
	}

	size_t n = (size_t)1 << logn;
	for (unsigned j = 0; j < MKEY_POLYS; j ++) {
		expand_share1(logn, j, seed, share1);

		/*
		 * Import each arithmetic share in a separate complete pass. M is a
		 * multiple of q and q', so independent reduction preserves the
		 * shared secret in both prime fields.
		 */
		for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
			for (size_t i = 0; i < n; i ++) {
				uint32_t v = sh == 0
					? dec32le(corr + ((size_t)j * n + i) * 4)
					: share1[i];
				dq[j][i].share[sh] = mmod32_encode_public(
					v, FNDSA_MASK_Q);
				dq2[j][i].share[sh] = mmod32_encode_public(
					v, FNDSA_MASK_Q2);
			}
		}

		/* Refresh only after both independent import passes are complete. */
		for (size_t i = 0; i < n; i ++) {
			mmod32_refresh(&dq[j][i], FNDSA_MASK_Q, rng);
			mmod32_refresh(&dq2[j][i], FNDSA_MASK_Q2, rng);
		}
	}
	secure_clear(share1, sizeof share1);
	if (vrfy_key != NULL) {
		*vrfy_key = vk;
	}
	if (vrfy_key_len != NULL) {
		*vrfy_key_len = FNDSA_VRFY_KEY_SIZE(logn);
	}
	return 1;
}

/* One-pass arena used by the explicit temporary-area API. */
typedef struct {
	uint8_t *ptr;
	uint8_t *end;
} masked_arena;

static void *
arena_take(masked_arena *a, size_t len, size_t align)
{
	uintptr_t p = (uintptr_t)a->ptr;
	uintptr_t q = (p + align - 1) & ~(uintptr_t)(align - 1);
	if (q < p || q > (uintptr_t)a->end
		|| len > (size_t)((uintptr_t)a->end - q))
	{
		return NULL;
	}
	a->ptr = (uint8_t *)(q + len);
	return (void *)q;
}

static uint8_t
masked_shake_next_u8(void *ctx)
{
	uint8_t x;
	shake_extract((shake_context *)ctx, &x, 1);
	return x;
}

static int
masked_crt_from_ntt(unsigned logn, mi32 *d, mmod32 *aq, mmod32 *aq2,
	fndsa_mask_rng *rng)
{
	if (!mmpoly_intt(logn, aq, FNDSA_MASK_MOD_Q, rng)
		|| !mmpoly_intt(logn, aq2, FNDSA_MASK_MOD_Q2, rng))
	{
		return 0;
	}
	mmpoly_crt_center(logn, d, aq, aq2, rng);
	return 1;
}

static int
masked_make_c_ntt(unsigned logn, const uint8_t nonce[40],
	const uint8_t mu[64], uint32_t *pub)
{
	uint16_t *c = (uint16_t *)(void *)pub;
	hash_to_point(logn, nonce, mu, c);
	return mmpoly_public_ntt(logn, pub, c, FNDSA_MASK_MOD_Q);
}

static void
masked_neg_i32(unsigned logn, mi32 *a)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		for (size_t sh = 0; sh < FNDSA_MASK_SHARES; sh ++) {
			a[i].share[sh] = (uint32_t)0 - a[i].share[sh];
		}
	}
}

static void
masked_attempt_clear(unsigned logn, mmod32 *mods, uint32_t *pub,
	void *tree_tmp, size_t tree_len, mi32 *iv, uint8_t *stage,
	size_t sig_len)
{
	size_t n = (size_t)1 << logn;
	secure_clear(mods, 4 * n * sizeof *mods);
	secure_clear(pub, n * sizeof *pub);
	secure_clear(tree_tmp, tree_len);
	secure_clear(iv, (7 * n / 2) * sizeof *iv);
	secure_clear(stage, sig_len);
}


/*
 * End-to-end masked standard signing core. The only opens are the aggregate
 * key/retry bits and, after acceptance, the signature component s1.
 */
static size_t
sign_masked_core(unsigned logn,
	const void *masked_key, size_t masked_key_len,
	const uint8_t *vk, size_t vk_len, const uint8_t mu[64],
	const uint8_t seedbuf[SEEDBUF_LEN],
	uint8_t *sig, void *tmp, size_t tmp_len)
{
	size_t n = (size_t)1 << logn;
	size_t sig_len = FNDSA_SIGNATURE_SIZE(logn);
	size_t tree_len = fndsa_masked_dyn_tmp_size(logn);
	masked_arena ar = { (uint8_t *)tmp, (uint8_t *)tmp + tmp_len };

	mmod32 *basis_mem = arena_take(&ar, 8 * n * sizeof(mmod32), 32);
	mmod32 *mods = arena_take(&ar, 4 * n * sizeof(mmod32), 32);
	uint32_t *pub = arena_take(&ar, n * sizeof(uint32_t), 32);
	void *tree_tmp = arena_take(&ar, tree_len, 32);
	mi32 *iv = arena_take(&ar, (7 * n / 2) * sizeof(mi32), 32);
	uint8_t *stage = arena_take(&ar, sig_len, 1);
	if (basis_mem == NULL || mods == NULL || pub == NULL
		|| tree_tmp == NULL || iv == NULL || stage == NULL)
	{
		secure_clear(tmp, tmp_len);
		return 0;
	}

	fndsa_masked_basis_mq basis = {
		basis_mem + 0 * n, basis_mem + 1 * n,
		basis_mem + 2 * n, basis_mem + 3 * n,
		basis_mem + 4 * n, basis_mem + 5 * n,
		basis_mem + 6 * n, basis_mem + 7 * n
	};
	mmod32 *m0 = mods + 0 * n;
	mmod32 *m1 = mods + 1 * n;
	mmod32 *m2 = mods + 2 * n;
	mmod32 *m3 = mods + 3 * n;
	mi32 *A = iv;
	mi32 *B = iv + n;                 /* n + n/2 elements */
	mi32 *z1 = iv + (5 * n / 2);      /* n elements */

	fndsa_mask_rng rng;
	shake_context pc;
	fndsa_mask_tree_state tree;
	memset(&pc, 0, sizeof pc);
	memset(&tree, 0, sizeof tree);
	fndsa_mask_rng_init_demo(&rng);
#ifdef FNDSA_MASK_TEST_HOOKS
	for (unsigned i = 0; i < mask_test_rng_skip; i ++) {
		(void)fndsa_mask_rng_next_u64(&rng);
	}
#endif
	const uint8_t *got_vk;
	size_t got_vk_len;
	if (!fndsa_masked_key_expand_mq_internal(masked_key, masked_key_len,
		&basis, &got_vk, &got_vk_len, &rng)
		|| got_vk != vk || got_vk_len != vk_len)
	{
		goto fail;
	}

	/* The clear signer rejects a basis whose (g,-f) norm exceeds 16822. */
	mmpoly_crt_center(logn, A, basis.f_q, basis.f_q2, &rng);
	mmpoly_crt_center(logn, B, basis.g_q, basis.g_q2, &rng);
	mbit key_ok;
	mi32_poly_sqnorm_le(&key_ok, logn, A, B, UINT32_C(16822), &rng);
	if (!mbit_open_compat(&key_ok)) {
		secure_clear(&key_ok, sizeof key_ok);
		goto fail;
	}
	secure_clear(&key_ok, sizeof key_ok);

	/* Expand and transform the eight immutable basis polynomials once. */
	if (!mmpoly_ntt(logn, basis.f_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.g_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.F_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.G_q, FNDSA_MASK_MOD_Q, &rng)
		|| !mmpoly_ntt(logn, basis.f_q2, FNDSA_MASK_MOD_Q2, &rng)
		|| !mmpoly_ntt(logn, basis.g_q2, FNDSA_MASK_MOD_Q2, &rng)
		|| !mmpoly_ntt(logn, basis.F_q2, FNDSA_MASK_MOD_Q2, &rng)
		|| !mmpoly_ntt(logn, basis.G_q2, FNDSA_MASK_MOD_Q2, &rng))
	{
		goto fail;
	}

	fndsa_mask_sampler_source source = {
		&pc, masked_shake_next_u8, logn
	};
	if (!fndsa_mask_tree_bind(&tree, logn, &source, &rng,
		tree_tmp, tree_len))
	{
		goto fail;
	}

	for (uint8_t counter = 0; counter < 27; counter ++) {
		shake_init(&pc, 256);
		shake_inject(&pc, seedbuf, SEEDBUF_LEN);
		shake_inject(&pc, &counter, 1);
		shake_flip(&pc);
		shake_extract(&pc, stage + 1, 40);
		if (!masked_make_c_ntt(logn, stage + 1, mu, pub)) {
			goto fail;
		}

		/* A = f*adj(f)+g*adj(g), recovered exactly with both primes. */
		mmpoly_gram_ntt(logn, m0, basis.f_q, basis.g_q,
			FNDSA_MASK_MOD_Q, &rng);
		mmpoly_gram_ntt(logn, m1, basis.f_q2, basis.g_q2,
			FNDSA_MASK_MOD_Q2, &rng);
		if (!masked_crt_from_ntt(logn, A, m0, m1, &rng)) {
			goto fail;
		}
		mxpoly_fft_i32_selfadj(logn, tree.work + n, A, &rng);

		/* Right target: center(c*f mod q). */
		mmpoly_mul_public_ntt(logn, m0, basis.f_q, pub,
			FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_intt_center(logn, tree.dz, m0,
			FNDSA_MASK_MOD_Q, &rng))
		{
			goto fail;
		}
		fndsa_mask_ffsampling_right(&tree, FNDSA_MASK_TREE_INIT0);
		memcpy(z1, tree.dz, n * sizeof *z1);

		/* B = F*adj(f)+G*adj(g), then (t1-z1)*B/A. */
		mmpoly_cross_ntt(logn, m0, basis.F_q, basis.f_q,
			basis.G_q, basis.g_q, FNDSA_MASK_MOD_Q, &rng);
		mmpoly_cross_ntt(logn, m1, basis.F_q2, basis.f_q2,
			basis.G_q2, basis.g_q2, FNDSA_MASK_MOD_Q2, &rng);
		if (!masked_crt_from_ntt(logn, B, m0, m1, &rng)) {
			goto fail;
		}
		/* The right tree consumes its Gram workspace; rebuild A before
		 * the division, exactly as the clear signing core does. */
		mmpoly_gram_ntt(logn, m0, basis.f_q, basis.g_q,
			FNDSA_MASK_MOD_Q, &rng);
		mmpoly_gram_ntt(logn, m1, basis.f_q2, basis.g_q2,
			FNDSA_MASK_MOD_Q2, &rng);
		if (!masked_crt_from_ntt(logn, A, m0, m1, &rng)) {
			goto fail;
		}
		mxpoly_fft_i32_selfadj(logn, tree.work + n, A, &rng);
		mxpoly_div_selfadj_spec_fft(logn, tree.work, tree.work + n, &rng);
		mxpoly_fft_i32(logn, tree.work + n, B, B, &rng);
		mxpoly_mul_fft(logn, tree.work, tree.work + n, &rng);

		/* Restore A for the left dynamic tree. */
		mmpoly_gram_ntt(logn, m0, basis.f_q, basis.g_q,
			FNDSA_MASK_MOD_Q, &rng);
		mmpoly_gram_ntt(logn, m1, basis.f_q2, basis.g_q2,
			FNDSA_MASK_MOD_Q2, &rng);
		if (!masked_crt_from_ntt(logn, A, m0, m1, &rng)) {
			goto fail;
		}
		mxpoly_fft_i32_selfadj(logn, tree.work + n, A, &rng);

		/* Left target: -center(c*F mod q). */
		mmpoly_mul_public_ntt(logn, m0, basis.F_q, pub,
			FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_intt_center(logn, tree.dz, m0,
			FNDSA_MASK_MOD_Q, &rng))
		{
			goto fail;
		}
		masked_neg_i32(logn, tree.dz);
		fndsa_mask_ffsampling_left(&tree, FNDSA_MASK_TREE_NOOUT);

		/* u0 = z0 + 3*center(cF mod q), in q'. */
		mmpoly_from_mi32(logn, m0, tree.dz, FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_mul_public_ntt(logn, m2, basis.F_q, pub,
			FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_intt_center(logn, A, m2,
			FNDSA_MASK_MOD_Q, &rng))
		{
			goto fail;
		}
		mmpoly_from_mi32(logn, m2, A, FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_mul_public_const(logn, m2, m2, UINT32_C(3),
			FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_add(logn, m0, m0, m2, FNDSA_MASK_MOD_Q2);

		/* u1 = z1 - 3*center(cf mod q), in q'. */
		mmpoly_from_mi32(logn, m1, z1, FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_mul_public_ntt(logn, m2, basis.f_q, pub,
			FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_intt_center(logn, B, m2,
			FNDSA_MASK_MOD_Q, &rng))
		{
			goto fail;
		}
		mmpoly_from_mi32(logn, m2, B, FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_mul_public_const(logn, m2, m2, UINT32_C(3),
			FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_sub(logn, m1, m1, m2, FNDSA_MASK_MOD_Q2);

		/* s1 = f*u0 + F*u1 modulo q', then one joint centering. */
		if (!mmpoly_ntt(logn, m0, FNDSA_MASK_MOD_Q2, &rng)
			|| !mmpoly_ntt(logn, m1, FNDSA_MASK_MOD_Q2, &rng))
		{
			goto fail;
		}
		mmpoly_mul_ntt(logn, m2, basis.f_q2, m0,
			FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_mul_ntt(logn, m3, basis.F_q2, m1,
			FNDSA_MASK_MOD_Q2, &rng);
		mmpoly_add(logn, m2, m2, m3, FNDSA_MASK_MOD_Q2);
		mmpoly_refresh(logn, m2, FNDSA_MASK_MOD_Q2, &rng);
		if (!mmpoly_intt_center(logn, A, m2,
			FNDSA_MASK_MOD_Q2, &rng))
		{
			goto fail;
		}

		/* s0 = center(c-s1*h mod q), with h read from the public key. */
		mmpoly_from_mi32(logn, m0, A, FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_ntt(logn, m0, FNDSA_MASK_MOD_Q, &rng)) {
			goto fail;
		}
		uint16_t *hp = (uint16_t *)(void *)pub;
		if (mqpoly_decode(logn, vk + 1, hp) != vk_len - 1) {
			goto fail;
		}
		/* FN-DSA verification keys already encode h in NTT order. */
		for (size_t i = n; i -- > 0;) {
			pub[i] = mmod32_encode_public(hp[i], FNDSA_MASK_Q);
		}
		mmpoly_mul_public_ntt(logn, m1, m0, pub,
			FNDSA_MASK_MOD_Q, &rng);
		if (!masked_make_c_ntt(logn, stage + 1, mu, pub)) {
			goto fail;
		}
		mmpoly_public_sub_ntt(logn, m2, pub, m1,
			FNDSA_MASK_MOD_Q, &rng);
		if (!mmpoly_intt_center(logn, B, m2,
			FNDSA_MASK_MOD_Q, &rng))
		{
			goto fail;
		}

		mbit accept;
		mi32_signature_is_acceptable(&accept, logn, B, A,
			sig_len - 41, &rng);
		uint32_t accepted = mbit_open_compat(&accept);
#ifdef FNDSA_MASK_TEST_HOOKS
		if (counter < mask_test_forced_rejects) {
			accepted = 0;
		}
#endif
		secure_clear(&accept, sizeof accept);
		if (accepted) {
			int16_t *opened = (int16_t *)(void *)pub;
			for (size_t i = 0; i < n; i ++) {
				opened[i] = (int16_t)mi32_open_compat(&A[i]);
			}
			if (!comp_encode(logn, opened, stage + 41, sig_len - 41)) {
				goto fail;
			}
			stage[0] = (uint8_t)((SIG_TYPE_STD << 4) + logn);
			memcpy(sig, stage, sig_len);
			secure_clear(&pc, sizeof pc);
			secure_clear(&tree, sizeof tree);
			secure_clear(&rng, sizeof rng);
			secure_clear(tmp, tmp_len);
			return sig_len;
		}

		masked_attempt_clear(logn, mods, pub, tree_tmp, tree_len,
			iv, stage, sig_len);
	}

fail:
	secure_clear(&pc, sizeof pc);
	secure_clear(&tree, sizeof tree);
	secure_clear(&rng, sizeof rng);
	secure_clear(tmp, tmp_len);
	return 0;
}

static size_t
sign_masked_real(int seeded,
	const void *masked_key, size_t masked_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	unsigned logn;
	const uint8_t *vk, *mseed, *hedge, *corr;
	if (!validate_mkey((const uint8_t *)masked_key, masked_key_len,
		&logn, &vk, &mseed, &hedge, &corr))
	{
		if (sig != NULL && max_sig_len != 0) {
			secure_clear(sig, max_sig_len);
		}
		return 0;
	}
	(void)mseed;
	(void)corr;
	size_t sig_len = FNDSA_SIGNATURE_SIZE(logn);
	if (sig == NULL) {
		return sig_len;
	}
	if (tmp == NULL || tmp_len < fndsa_masked_tmp_size(logn)
		|| max_sig_len < sig_len)
	{
		secure_clear(sig, max_sig_len < sig_len ? max_sig_len : sig_len);
		return 0;
	}
	secure_clear(sig, sig_len);

	uint8_t mu[64], hvk[64], seedbuf[SEEDBUF_LEN];
	if (id != NULL && *(const uint8_t *)id == 0xFE) {
		if (hv_len != sizeof mu) {
			goto fail;
		}
		memcpy(mu, hv, sizeof mu);
	} else {
		fndsa_hashed_vrfykey_from_vrfykey(
			hvk, vk, FNDSA_VRFY_KEY_SIZE(logn));
		if (!fndsa_compute_mu(mu, hvk, ctx, ctx_len, id, hv, hv_len)) {
			goto fail;
		}
	}
	memcpy(seedbuf, hedge, sizeof seedbuf);
	if (!make_seedbuf(seedbuf, NULL, mu,
		seeded ? seed : NULL, seeded ? seed_len : 0, "20261006"))
	{
		goto fail;
	}

	size_t r = sign_masked_core(logn, masked_key, masked_key_len,
		vk, FNDSA_VRFY_KEY_SIZE(logn), mu, seedbuf,
		(uint8_t *)sig, tmp, tmp_len);
	secure_clear(mu, sizeof mu);
	secure_clear(hvk, sizeof hvk);
	secure_clear(seedbuf, sizeof seedbuf);
	if (r == 0) {
		secure_clear(sig, sig_len);
	}
	return r;

fail:
	secure_clear(mu, sizeof mu);
	secure_clear(hvk, sizeof hvk);
	secure_clear(seedbuf, sizeof seedbuf);
	secure_clear(tmp, tmp_len);
	secure_clear(sig, sig_len);
	return 0;
}

size_t
fndsa_sign_masked_temp(
	const void *masked_key, size_t masked_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_masked_real(0, masked_key, masked_key_len,
		ctx, ctx_len, id, hv, hv_len, NULL, 0,
		sig, max_sig_len, tmp, tmp_len);
}

size_t
fndsa_sign_masked_seeded_temp(
	const void *masked_key, size_t masked_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_masked_real(1, masked_key, masked_key_len,
		ctx, ctx_len, id, hv, hv_len, seed, seed_len,
		sig, max_sig_len, tmp, tmp_len);
}
