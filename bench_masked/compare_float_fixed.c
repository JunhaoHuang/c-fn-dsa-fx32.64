/*
 * Same-process microbenchmark for burak231708/masked_FALCON's first-order
 * masked IEEE-754 gadgets and this project's first-order masked Q32.64
 * gadgets.
 *
 * This file intentionally contains no source copied from masked_FALCON.  It
 * calls that repository's public gadget API and must be linked with a local
 * checkout of that repository.  The results are research measurements only:
 * both implementations use deterministic/non-cryptographic demo randomness,
 * and their representations and currently implemented operation sets differ.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined __x86_64__ || defined __i386__
#include <x86intrin.h>
#else
#error This comparison benchmark currently requires an x86 TSC.
#endif

#include "mask160.h"
#include "mask_sampler.h"

/* Public headers from the separately checked-out masked_FALCON repository. */
#include "gadgets.h"
#include "fpr_gadgets.h"
#include "secfpr.h"

#define MAX_SAMPLES   101

static volatile uint64_t bench_sink;

typedef void (*bench_proc)(void *ctx, unsigned reps);

static inline uint64_t
read_tsc(void)
{
	_mm_lfence();
	return __rdtsc();
}

static int
cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t *)a;
	uint64_t y = *(const uint64_t *)b;
	return x < y ? -1 : x > y;
}

static double
measure(bench_proc proc, void *ctx, unsigned reps,
	unsigned samples, unsigned warmup)
{
	uint64_t tt[MAX_SAMPLES];
	if (reps == 0 || samples == 0 || samples > MAX_SAMPLES
		|| (samples & 1) == 0)
	{
		fprintf(stderr, "invalid benchmark parameters\n");
		exit(EXIT_FAILURE);
	}
	for (unsigned i = 0; i < warmup + samples; i ++) {
		uint64_t begin = read_tsc();
		proc(ctx, reps);
		uint64_t end = read_tsc();
		if (i >= warmup) {
			tt[i - warmup] = end - begin;
		}
	}
	qsort(tt, samples, sizeof tt[0], cmp_u64);
	return (double)tt[samples >> 1] / (double)reps;
}

typedef struct {
	maskedb_t fa;
	maskedb_t fb;
	maskedb_t fhalf;
	maskedb_t ffloor;
	maskedb_t fsigma;
	maskedb_t falea;
	maskedb_t fout;
	maskeda_t fexp_a;
	maskedb_t fexp_b;

	fndsa_mask_rng rng;
	mrx160 xa;
	mrx160 xb;
	mrx160 xrefresh;
	mrx160 xout;
	mb160 xbool;
	mrx160 xhalf;
	mrx160 xfloor;
	mrx160 xsigma;
	mrx160 xisigma;
	mb64 xword;
	mbit xbit;
	uint64_t byte_state;
	unsigned byte_pos;
} bench_ctx;

static uint64_t
double_bits(double x)
{
	uint64_t r;
	memcpy(&r, &x, sizeof r);
	return r;
}

static double
bits_double(uint64_t x)
{
	double r;
	memcpy(&r, &x, sizeof r);
	return r;
}

static void
init_ctx(bench_ctx *c)
{
	uint32_t a[3] = { 0, 0, 5 };
	uint32_t b[3] = { 0, 0, 3 };

	memset(c, 0, sizeof *c);
	srand(1);
	MaskB(c->fa, double_bits(5.0));
	MaskB(c->fb, double_bits(3.0));
	MaskB(c->fhalf, double_bits(0.5));
	MaskB(c->ffloor, double_bits(5.5));
	MaskB(c->fsigma, double_bits(1.5));
	MaskB(c->falea, UINT64_C(0x3141592653589793));
	MaskA(c->fexp_a, UINT64_C(12345), UINT64_C(1) << 16);
	MaskB(c->fexp_b, UINT64_C(12345));

	fndsa_mask_rng_init_demo(&c->rng);
	mrx160_share_rx(&c->xa, a, &c->rng);
	mrx160_share_rx(&c->xb, b, &c->rng);
	mrx160_copy(&c->xrefresh, &c->xa);
	mrx160_a2b(&c->xbool, &c->xa, &c->rng);
	uint32_t half[3] = { 0, UINT32_C(0x80000000), 0 };
	uint32_t floorv[3] = { 0, UINT32_C(0x80000000), 5 };
	uint32_t sigma[3] = { 0, UINT32_C(0x80000000), 1 };
	uint32_t isigma[3] = {
		UINT32_C(0xAAAAAAAA), UINT32_C(0xAAAAAAAA), 0
	};
	mrx160_share_rx(&c->xhalf, half, &c->rng);
	mrx160_share_rx(&c->xfloor, floorv, &c->rng);
	mrx160_share_rx(&c->xsigma, sigma, &c->rng);
	mrx160_share_rx(&c->xisigma, isigma, &c->rng);
	c->byte_state = UINT64_C(0xD31947A6C85E20BF);
}

static void
proc_float_mask(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		MaskB(c->fout, double_bits(5.0));
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_fixed_mask(void *v, unsigned reps)
{
	bench_ctx *c = v;
	uint32_t a[3] = { 0, 0, 5 };
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_share_rx(&c->xout, a, &c->rng);
	}
	bench_sink ^= c->xout.share[0][0] ^ c->xout.share[1][4];
}

static void
proc_float_refresh(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		RefreshXOR_64(c->fout, c->fa, MASKSIZE);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_fixed_refresh(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_refresh(&c->xrefresh, &c->rng);
	}
	bench_sink ^= c->xrefresh.share[0][0] ^ c->xrefresh.share[1][4];
}

static void
proc_float_add(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprAdd(c->fout, c->fa, c->fb);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_fixed_add(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_add(&c->xout, &c->xa, &c->xb);
	}
	bench_sink ^= c->xout.share[0][0] ^ c->xout.share[1][4];
}

static void
proc_float_mul(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprMul(c->fout, c->fa, c->fb);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_fixed_mul(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul(&c->xout, &c->xa, &c->xb, &c->rng);
	}
	bench_sink ^= c->xout.share[0][0] ^ c->xout.share[1][4];
}

static void
proc_fixed_mul_wide(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul_wide(&c->xout, &c->xa, &c->xb, &c->rng);
	}
	bench_sink ^= c->xout.share[0][0] ^ c->xout.share[1][4];
}

static void
proc_float_a2b16(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		A2B(c->fout, c->fexp_a, UINT64_C(1) << 16);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_fixed_a2b160(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_a2b(&c->xbool, &c->xa, &c->rng);
	}
	bench_sink ^= c->xbool.share[0][0] ^ c->xbool.share[1][4];
}

static void
proc_float_b2a16(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		B2A(c->fexp_a, c->fexp_b, UINT64_C(1) << 16, MASKSIZE);
	}
	bench_sink ^= c->fexp_a[0] ^ c->fexp_a[1];
}

static void
proc_fixed_b2a160(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_b2a(&c->xout, &c->xbool, &c->rng);
	}
	bench_sink ^= c->xout.share[0][0] ^ c->xout.share[1][4];
}

static uint8_t
fixed_ber_next_u8(void *v)
{
	bench_ctx *c = v;
	static const uint64_t ALEA = UINT64_C(0x3141592653589793);
	unsigned n = c->byte_pos ++ & 7;
	return (uint8_t)(ALEA >> (56 - 8 * n));
}

static uint8_t
fixed_sampler_next_u8(void *v)
{
	bench_ctx *c = v;
	uint64_t x = c->byte_state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	c->byte_state = x;
	return (uint8_t)x;
}

static void
proc_fixed_floor(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_floor(&c->xout, &c->xfloor, &c->rng);
	}
	bench_sink ^= c->xout.share[0][2] ^ c->xout.share[1][2];
}

static void
proc_fixed_inv(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_inv_spec(&c->xout, &c->xa, &c->rng);
	}
	bench_sink ^= c->xout.share[0][1] ^ c->xout.share[1][1];
}

static void
proc_fixed_exp(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_expm_p63(&c->xword, &c->xhalf, &c->xhalf, &c->rng);
	}
	bench_sink ^= c->xword.share[0][0] ^ c->xword.share[1][1];
}

static void
proc_fixed_berexp(void *v, unsigned reps)
{
	bench_ctx *c = v;
	fndsa_mask_sampler_source src = { c, fixed_ber_next_u8, 9 };
	for (unsigned i = 0; i < reps; i ++) {
		c->byte_pos = 0;
		mrx160_ber_exp(&c->xbit, &src,
			&c->xhalf, &c->xhalf, &c->rng);
	}
	bench_sink ^= c->xbit.share[0] ^ c->xbit.share[1];
}

static void
proc_fixed_samplerz(void *v, unsigned reps)
{
	bench_ctx *c = v;
	fndsa_mask_sampler_source src = { c, fixed_sampler_next_u8, 9 };
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_sampler_z(&c->xout, &src,
			&c->xhalf, &c->xisigma, &c->rng);
	}
	bench_sink ^= c->xout.share[0][2] ^ c->xout.share[1][2];
}

static void
proc_fixed_inv_samplerz(void *v, unsigned reps)
{
	bench_ctx *c = v;
	fndsa_mask_sampler_source src = { c, fixed_sampler_next_u8, 9 };
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_inv_spec(&c->xout, &c->xsigma, &c->rng);
		mrx160_sampler_z(&c->xout, &src,
			&c->xhalf, &c->xout, &c->rng);
	}
	bench_sink ^= c->xout.share[0][2] ^ c->xout.share[1][2];
}

static void
proc_float_floor(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprFloor(c->fout, c->ffloor);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_float_inv(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprInv(c->fout, c->fa);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_float_exp(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecApproxExp(c->fout, c->fhalf, c->fhalf);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_float_berexp(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		bench_sink ^= SecFprBerExp(c->fout,
			c->fhalf, c->fhalf, c->falea);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
proc_float_samplerz(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SamplerZ(c->fout, c->fhalf, c->fsigma);
	}
	bench_sink ^= c->fout[0] ^ c->fout[1];
}

static void
validate(bench_ctx *c)
{
	uint64_t u;
	uint32_t x[3];

	SecFprAdd(c->fout, c->fa, c->fb);
	UnmaskB(&u, c->fout);
	if (bits_double(u) != 8.0) {
		fprintf(stderr, "masked_FALCON add validation failed\n");
		exit(EXIT_FAILURE);
	}
	SecFprMul(c->fout, c->fa, c->fb);
	UnmaskB(&u, c->fout);
	if (bits_double(u) != 15.0) {
		fprintf(stderr, "masked_FALCON multiply validation failed\n");
		exit(EXIT_FAILURE);
	}
	mrx160_add(&c->xout, &c->xa, &c->xb);
	mrx160_open_rx(x, &c->xout);
	if (x[0] != 0 || x[1] != 0 || x[2] != 8) {
		fprintf(stderr, "fixed-point add validation failed\n");
		exit(EXIT_FAILURE);
	}
	mrx160_mul(&c->xout, &c->xa, &c->xb, &c->rng);
	mrx160_open_rx(x, &c->xout);
	if (x[0] != 0 || x[1] != 0 || x[2] != 15) {
		fprintf(stderr, "fixed-point multiply validation failed\n");
		exit(EXIT_FAILURE);
	}
	SecFprFloor(c->fout, c->ffloor);
	UnmaskB(&u, c->fout);
	if (bits_double(u) != 5.0) {
		fprintf(stderr, "masked_FALCON floor validation failed\n");
		exit(EXIT_FAILURE);
	}
	mrx160_floor(&c->xout, &c->xfloor, &c->rng);
	mrx160_open_rx(x, &c->xout);
	if (x[0] != 0 || x[1] != 0 || x[2] != 5) {
		fprintf(stderr, "fixed-point floor validation failed\n");
		exit(EXIT_FAILURE);
	}
	uint64_t fb = SecFprBerExp(c->fout,
		c->fhalf, c->fhalf, c->falea);
	fndsa_mask_sampler_source src = { c, fixed_ber_next_u8, 9 };
	c->byte_pos = 0;
	mrx160_ber_exp(&c->xbit, &src,
		&c->xhalf, &c->xhalf, &c->rng);
	if (fb != mbit_open_compat(&c->xbit)) {
		fprintf(stderr, "BerExp decision validation failed\n");
		exit(EXIT_FAILURE);
	}
}

static double
show(const char *name, bench_proc proc, bench_ctx *c, unsigned reps)
{
	double v = measure(proc, c, reps, 101, 20);
	printf("%-34s %12.2f TSC-ticks/op\n", name, v);
	return v;
}

int
main(void)
{
	bench_ctx c;
	double fadd, xadd, fmul, xmul, xmulw;
	double ffloor, xfloor, finv, xinv, fexp, xexp;
	double fber, xber, fsamp, xsamp, xpipeline;

	init_ctx(&c);
	validate(&c);
	printf("masked_FALCON float vs masked Q32.64 -- DEMO ONLY\n");
	printf("Counter: invariant-TSC reference ticks (not hardware core cycles)\n");
	printf("Both sides: first order / two shares; same process and compiler.\n\n");

	printf("Representation support operations (different widths)\n");
	show("float MaskB (64-bit Boolean)", proc_float_mask, &c, 2000);
	show("fixed share (160-bit arithmetic)", proc_fixed_mask, &c, 1000);
	show("float refresh (64-bit Boolean)", proc_float_refresh, &c, 1000);
	show("fixed refresh (160-bit arithmetic)", proc_fixed_refresh, &c, 1000);
	show("float A2B (16-bit exponent)", proc_float_a2b16, &c, 100);
	show("fixed A2B (160-bit value)", proc_fixed_a2b160, &c, 100);
	show("float B2A (16-bit exponent)", proc_float_b2a16, &c, 200);
	show("fixed B2A (160-bit value)", proc_fixed_b2a160, &c, 200);

	printf("\nNumerical operations\n");
	fadd = show("float add (IEEE-754)", proc_float_add, &c, 20);
	xadd = show("fixed add (Q32.64)", proc_fixed_add, &c, 4000);
	fmul = show("float multiply (IEEE-754)", proc_float_mul, &c, 10);
	xmul = show("fixed multiply (160-bit fast)", proc_fixed_mul, &c, 100);
	xmulw = show("fixed multiply (192-bit fallback)",
		proc_fixed_mul_wide, &c, 30);
	printf("\nDirect-operation ratios (float/fixed; lower time is faster)\n");
	printf("add:                  %.2fx\n", fadd / xadd);
	printf("multiply vs fast:     %.2fx\n", fmul / xmul);
	printf("multiply vs fallback: %.2fx\n", fmul / xmulw);

	printf("\nCorresponding floor, inverse, exponential and sampling gadgets\n");
	ffloor = show("float floor", proc_float_floor, &c, 10);
	xfloor = show("fixed floor", proc_fixed_floor, &c, 100);
	finv = show("float inverse", proc_float_inv, &c, 1);
	xinv = show("fixed inverse", proc_fixed_inv, &c, 5);
	fexp = show("float ApproxExp", proc_float_exp, &c, 1);
	xexp = show("fixed expm_p63", proc_fixed_exp, &c, 5);
	fber = show("float BerExp", proc_float_berexp, &c, 1);
	xber = show("fixed BerExp", proc_fixed_berexp, &c, 5);
	fsamp = measure(proc_float_samplerz, &c, 1, 31, 5);
	xsamp = measure(proc_fixed_samplerz, &c, 1, 31, 5);
	xpipeline = measure(proc_fixed_inv_samplerz, &c, 1, 31, 5);
	printf("%-34s %12.2f TSC-ticks/op\n", "float SamplerZ(sigma)", fsamp);
	printf("%-34s %12.2f TSC-ticks/op\n", "fixed SamplerZ(isigma)", xsamp);
	printf("%-34s %12.2f TSC-ticks/op\n",
		"fixed inverse+SamplerZ(sigma)", xpipeline);

	printf("\nCorresponding-gadget ratios (float/fixed)\n");
	printf("floor:                    %.2fx\n", ffloor / xfloor);
	printf("inverse:                  %.2fx\n", finv / xinv);
	printf("ApproxExp/expm_p63:       %.2fx\n", fexp / xexp);
	printf("BerExp:                   %.2fx\n", fber / xber);
	printf("SamplerZ precomputed i-s: %.2fx\n", fsamp / xsamp);
	printf("SamplerZ full pipeline:   %.2fx\n", fsamp / xpipeline);

	printf("\nbenchmark sink: %llu\n",
		(unsigned long long)bench_sink);
	return 0;
}
