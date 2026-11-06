/*
 * Same-host microbenchmark for three first-order/two-share implementations:
 * burak231708/masked_FALCON, kengyuchen/masked-falcon-power-analysis and
 * this project's Q32.64/160-bit gadgets.
 *
 * Third-party implementations are compiled from separate pinned checkouts and
 * linked through their public headers. Reported values are invariant-TSC
 * reference ticks, not PMU core cycles. All random sources are demo/test
 * sources and their cost is part of the measured implementation.
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
#define rand64 cc_rand64
#include "mask_fpr_nshares.h"
#include "mask_rand.h"
#include "uint128.h"
#undef rand64
#include "gadgets.h"
#include "fpr_gadgets.h"

#define BENCH_SHARES  2
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
	fndsa_mask_rng rng;
	mrx160 xa;
	mrx160 xb;
	mrx160 xout;

	uint64_t ka[BENCH_SHARES];
	uint64_t kb[BENCH_SHARES];
	uint64_t kout[BENCH_SHARES];
	uint8_t ksa[BENCH_SHARES];
	uint8_t ksb[BENCH_SHARES];
	uint16_t kea[BENCH_SHARES];
	uint16_t keb[BENCH_SHARES];
	UINT128 kma[BENCH_SHARES];
	UINT128 kmb[BENCH_SHARES];

	maskedb_t bfa;
	maskedb_t bfb;
	maskedb_t bfout;
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
share_bool64(uint64_t d[BENCH_SHARES], uint64_t x)
{
	d[1] = cc_rand64();
	d[0] = x ^ d[1];
}

static void
share_kengyu_mul(uint64_t x, uint8_t s[BENCH_SHARES],
	uint16_t e[BENCH_SHARES], UINT128 m[BENCH_SHARES])
{
	UINT128 v;
	uint64_t mant = (x & ((UINT64_C(1) << 52) - 1))
		| (UINT64_C(1) << 52);
	v.value[0] = 0;
	v.value[1] = 0;
	v.value[2] = (uint32_t)(mant >> 32);
	v.value[3] = (uint32_t)mant;

	m[1] = rand128();
	m[0] = UINT128Sub(v, m[1]);
	e[1] = rand16();
	e[0] = (uint16_t)((x >> 52) & 0x7FF) - e[1];
	s[1] = rand8();
	s[0] = (uint8_t)(x >> 63) ^ s[1];
}

static void
init_ctx(bench_ctx *c)
{
	uint32_t a[3] = { 0, 0, 5 };
	uint32_t b[3] = { 0, 0, 3 };
	uint64_t ba = double_bits(5.0);
	uint64_t bb = double_bits(3.0);

	memset(c, 0, sizeof *c);
	fndsa_mask_rng_init_demo(&c->rng);
	mrx160_share_rx(&c->xa, a, &c->rng);
	mrx160_share_rx(&c->xb, b, &c->rng);

	init_prng_table(UINT64_C(1));
	share_bool64(c->ka, ba);
	share_bool64(c->kb, bb);
	share_kengyu_mul(ba, c->ksa, c->kea, c->kma);
	share_kengyu_mul(bb, c->ksb, c->keb, c->kmb);

	srand(1);
	MaskB(c->bfa, ba);
	MaskB(c->bfb, bb);
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
proc_kengyu_add(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprAdd_nshares(c->ka, c->kb, c->kout, BENCH_SHARES);
	}
	bench_sink ^= c->kout[0] ^ c->kout[1];
}

static void
proc_kengyu_mul(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprMul_nshares(c->ksa, c->kea, c->kma,
			c->ksb, c->keb, c->kmb, c->kout, BENCH_SHARES);
	}
	bench_sink ^= c->kout[0] ^ c->kout[1];
}

static void
proc_burak_add(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprAdd(c->bfout, c->bfa, c->bfb);
	}
	bench_sink ^= c->bfout[0] ^ c->bfout[1];
}

static void
proc_burak_mul(void *v, unsigned reps)
{
	bench_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		SecFprMul(c->bfout, c->bfa, c->bfb);
	}
	bench_sink ^= c->bfout[0] ^ c->bfout[1];
}

static void
validate(bench_ctx *c)
{
	uint32_t x[3];
	uint64_t u;

	mrx160_add(&c->xout, &c->xa, &c->xb);
	mrx160_open_rx(x, &c->xout);
	if (x[0] != 0 || x[1] != 0 || x[2] != 8) {
		fprintf(stderr, "fixed add validation failed\n");
		exit(EXIT_FAILURE);
	}
	mrx160_mul(&c->xout, &c->xa, &c->xb, &c->rng);
	mrx160_open_rx(x, &c->xout);
	if (x[0] != 0 || x[1] != 0 || x[2] != 15) {
		fprintf(stderr, "fixed multiply validation failed\n");
		exit(EXIT_FAILURE);
	}

	SecFprAdd_nshares(c->ka, c->kb, c->kout, BENCH_SHARES);
	u = c->kout[0] ^ c->kout[1];
	if (bits_double(u) != 8.0) {
		fprintf(stderr, "Chen-Chen add validation failed\n");
		exit(EXIT_FAILURE);
	}
	SecFprMul_nshares(c->ksa, c->kea, c->kma,
		c->ksb, c->keb, c->kmb, c->kout, BENCH_SHARES);
	u = c->kout[0] ^ c->kout[1];
	if (bits_double(u) != 15.0) {
		fprintf(stderr, "Chen-Chen multiply validation failed\n");
		exit(EXIT_FAILURE);
	}

	SecFprAdd(c->bfout, c->bfa, c->bfb);
	UnmaskB(&u, c->bfout);
	if (bits_double(u) != 8.0) {
		fprintf(stderr, "Burak add validation failed\n");
		exit(EXIT_FAILURE);
	}
	SecFprMul(c->bfout, c->bfa, c->bfb);
	UnmaskB(&u, c->bfout);
	if (bits_double(u) != 15.0) {
		fprintf(stderr, "Burak multiply validation failed\n");
		exit(EXIT_FAILURE);
	}
}

static double
show(const char *name, bench_proc proc, bench_ctx *c, unsigned reps)
{
	double v = measure(proc, c, reps, 101, 20);
	printf("%-38s %12.2f TSC-ticks/op\n", name, v);
	return v;
}

int
main(void)
{
	bench_ctx c;
	double xadd, xmul, xmulw, kadd, kmul, badd, bmul;

	init_ctx(&c);
	validate(&c);
	printf("Three masked arithmetic implementations -- DEMO ONLY\n");
	printf("Counter: invariant-TSC reference ticks (not PMU core cycles)\n");
	printf("All sides: first order / two shares; inputs 5 and 3.\n\n");

	xadd = show("fixed add (Q32.64)", proc_fixed_add, &c, 4000);
	kadd = show("Chen-Chen add (IEEE-754)", proc_kengyu_add, &c, 20);
	badd = show("Burak add (IEEE-754)", proc_burak_add, &c, 20);
	xmul = show("fixed multiply (160-bit fast)", proc_fixed_mul, &c, 100);
	xmulw = show("fixed multiply (192-bit fallback)",
		proc_fixed_mul_wide, &c, 30);
	kmul = show("Chen-Chen multiply (IEEE-754)", proc_kengyu_mul, &c, 20);
	bmul = show("Burak multiply (IEEE-754)", proc_burak_mul, &c, 10);

	printf("\nratios\n");
	printf("Chen/fixed add:            %.2fx\n", kadd / xadd);
	printf("Burak/fixed add:           %.2fx\n", badd / xadd);
	printf("Chen/fixed fast multiply:  %.2fx\n", kmul / xmul);
	printf("Burak/fixed fast multiply: %.2fx\n", bmul / xmul);
	printf("Chen/fixed wide multiply:  %.2fx\n", kmul / xmulw);
	printf("Burak/fixed wide multiply: %.2fx\n", bmul / xmulw);
	printf("Burak/Chen add:            %.2fx\n", badd / kadd);
	printf("Burak/Chen multiply:       %.2fx\n", bmul / kmul);
	printf("\nbenchmark sink: %llu\n",
		(unsigned long long)bench_sink);
	return 0;
}
