/*
 * Cycle benchmarks for the experimental two-share masking backend.
 *
 * On Linux x86, the default opens PERF_COUNT_HW_CPU_CYCLES through
 * perf_event_open(), so the kernel selects the correct PMU on AMD and Intel.
 * This requires root, CAP_PERFMON, or kernel.perf_event_paranoid=0. Other
 * targets use their architecture cycle counter. Use "-tsc" for the userspace
 * reference counter; TSC ticks are not necessarily CPU core cycles.
 *
 * WARNING: the masking backend and its deterministic LFSR are DEMO ONLY.
 */

#if defined __linux__
#define _GNU_SOURCE
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined __linux__
#include <errno.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

#include "fndsa.h"
#include "fndsa_masked.h"
#include "sign_inner.h"
#include "mask160.h"
#include "mask_mq.h"
#include "mask_sampler.h"
#include "mask_xpoly.h"
#include "mask_tree.h"

static int use_tsc;
static volatile uint64_t bench_sink;

#if defined __linux__ && (defined __x86_64__ || defined __i386__)
static int perf_cycles_fd = -1;

static void
counter_init(void)
{
	if (use_tsc) {
		return;
	}
	struct perf_event_attr attr;
	memset(&attr, 0, sizeof attr);
	attr.type = PERF_TYPE_HARDWARE;
	attr.size = sizeof attr;
	attr.config = PERF_COUNT_HW_CPU_CYCLES;
	attr.disabled = 1;
	attr.pinned = 1;
	attr.exclude_kernel = 1;
	attr.exclude_hv = 1;
	perf_cycles_fd = (int)syscall(__NR_perf_event_open,
		&attr, 0, -1, -1, 0);
	if (perf_cycles_fd < 0) {
		fprintf(stderr,
			"cannot open hardware CPU-cycle event: %s\n",
			strerror(errno));
		fprintf(stderr,
			"run as root/CAP_PERFMON, or set kernel.perf_event_paranoid=0\n");
		exit(EXIT_FAILURE);
	}
	if (ioctl(perf_cycles_fd, PERF_EVENT_IOC_RESET, 0) != 0
		|| ioctl(perf_cycles_fd, PERF_EVENT_IOC_ENABLE, 0) != 0)
	{
		fprintf(stderr, "cannot enable hardware CPU-cycle event: %s\n",
			strerror(errno));
		close(perf_cycles_fd);
		exit(EXIT_FAILURE);
	}
}

static uint64_t
read_perf_cycles(void)
{
	uint64_t x;
	if (read(perf_cycles_fd, &x, sizeof x) != (ssize_t)sizeof x) {
		fprintf(stderr, "cannot read hardware CPU-cycle event: %s\n",
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	return x;
}

static void
counter_close(void)
{
	if (perf_cycles_fd >= 0) {
		(void)ioctl(perf_cycles_fd, PERF_EVENT_IOC_DISABLE, 0);
		close(perf_cycles_fd);
		perf_cycles_fd = -1;
	}
}
#else
static void counter_init(void) { }
static void counter_close(void) { }
#endif

#if defined __x86_64__ || defined _M_X64 || defined __i386__ || defined _M_IX86

#include <immintrin.h>
#ifdef _MSC_VER
#ifndef __rdpmc
#define __rdpmc   __readpmc
#endif
#else
#include <x86intrin.h>
#endif

#if defined __GNUC__ || defined __clang__
__attribute__((target("sse2")))
#endif
static inline uint64_t
read_counter(void)
{
	if (use_tsc) {
		_mm_lfence();
		return __rdtsc();
	}
#if defined __linux__
	return read_perf_cycles();
#else
	_mm_lfence();
	return __rdpmc(0x40000001);
#endif
}

#elif defined __aarch64__ && (defined __GNUC__ || defined __clang__)

static inline uint64_t
read_counter(void)
{
	uint64_t x;
	if (use_tsc) {
		__asm__ __volatile__ (
			"dsb sy\n\tmrs %0, cntvct_el0" : "=r" (x) : : "memory");
	} else {
		__asm__ __volatile__ (
			"dsb sy\n\tmrs %0, pmccntr_el0" : "=r" (x) : : "memory");
	}
	return x;
}

#elif defined __riscv && defined __riscv_xlen && __riscv_xlen >= 64

static inline uint64_t
read_counter(void)
{
	uint64_t x;
	if (use_tsc) {
		__asm__ __volatile__ ("rdtime %0" : "=r" (x) : : "memory");
	} else {
		__asm__ __volatile__ ("rdcycle %0" : "=r" (x) : : "memory");
	}
	return x;
}

#else
#error Architecture not supported by speed_masked
#endif

static inline void
compiler_barrier(void)
{
#if defined __GNUC__ || defined __clang__
	__asm__ __volatile__ ("" : : : "memory");
#else
	volatile uint64_t x = bench_sink;
	(void)x;
#endif
}

static int
cmp_u64(const void *a, const void *b)
{
	uint64_t x = *(const uint64_t *)a;
	uint64_t y = *(const uint64_t *)b;
	return x < y ? -1 : x > y;
}

typedef void (*bench_proc)(void *ctx, unsigned reps);

static double
measure(bench_proc proc, void *ctx, unsigned reps,
	unsigned samples, unsigned warmup)
{
	uint64_t tt[101];
	if (reps == 0 || samples == 0 || samples > 101 || (samples & 1) == 0) {
		fprintf(stderr, "invalid benchmark parameters\n");
		exit(EXIT_FAILURE);
	}
	for (unsigned i = 0; i < warmup + samples; i ++) {
		compiler_barrier();
		uint64_t begin = read_counter();
		proc(ctx, reps);
		uint64_t end = read_counter();
		compiler_barrier();
		if (i >= warmup) {
			tt[i - warmup] = end - begin;
		}
	}
	qsort(tt, samples, sizeof tt[0], cmp_u64);
	return (double)tt[samples >> 1] / (double)reps;
}

static void
print_bench(const char *name, bench_proc proc, void *ctx, unsigned reps)
{
	printf("%-38s %13.2f %s/op\n",
		name, measure(proc, ctx, reps, 101, 20),
		use_tsc ? "TSC-ticks" : "cycles");
}

static void
print_slow_bench(const char *name, bench_proc proc, void *ctx)
{
	printf("%-38s %13.2f %s/op\n",
		name, measure(proc, ctx, 1, 11, 2),
		use_tsc ? "TSC-ticks" : "cycles");
}
/* Modular polynomial benchmark state. */
typedef struct {
	unsigned logn;
	size_t n;
	fndsa_mask_rng rng;
	mmod32 *q_a;
	mmod32 *q_b;
	mmod32 *q_d;
	mmod32 *q2_a;
	mmod32 *q2_b;
	mmod32 *q2_d;
	mi32 *crt;
	mi32 raw;
	uint16_t *plain_q_a;
	uint16_t *plain_q_b;
	uint16_t *plain_q2_a;
	uint16_t *plain_q2_b;
	int16_t *plain_crt_q;
	int16_t *plain_crt_q2;
	int32_t *plain_crt_d;
} mq_backend_ctx;

static void
init_mq_backend_ctx(mq_backend_ctx *c, unsigned logn)
{
	memset(c, 0, sizeof *c);
	c->logn = logn;
	c->n = (size_t)1 << logn;
	c->q_a = malloc(c->n * sizeof *c->q_a);
	c->q_b = malloc(c->n * sizeof *c->q_b);
	c->q_d = malloc(c->n * sizeof *c->q_d);
	c->q2_a = malloc(c->n * sizeof *c->q2_a);
	c->q2_b = malloc(c->n * sizeof *c->q2_b);
	c->q2_d = malloc(c->n * sizeof *c->q2_d);
	c->crt = malloc(c->n * sizeof *c->crt);
	c->plain_q_a = malloc(c->n * sizeof *c->plain_q_a);
	c->plain_q_b = malloc(c->n * sizeof *c->plain_q_b);
	c->plain_q2_a = malloc(c->n * sizeof *c->plain_q2_a);
	c->plain_q2_b = malloc(c->n * sizeof *c->plain_q2_b);
	c->plain_crt_q = malloc(c->n * sizeof *c->plain_crt_q);
	c->plain_crt_q2 = malloc(c->n * sizeof *c->plain_crt_q2);
	c->plain_crt_d = malloc(c->n * sizeof *c->plain_crt_d);
	if (c->q_a == NULL || c->q_b == NULL || c->q_d == NULL
		|| c->q2_a == NULL || c->q2_b == NULL || c->q2_d == NULL
		|| c->crt == NULL || c->plain_q_a == NULL
		|| c->plain_q_b == NULL || c->plain_q2_a == NULL
		|| c->plain_q2_b == NULL || c->plain_crt_q == NULL
		|| c->plain_crt_q2 == NULL || c->plain_crt_d == NULL)
	{
		fprintf(stderr, "masked MQ benchmark allocation failure\n");
		exit(EXIT_FAILURE);
	}
	fndsa_mask_rng_init_demo(&c->rng);
	for (size_t i = 0; i < c->n; i ++) {
		uint32_t x = (uint32_t)(17 * i + 3);
		uint32_t y = (uint32_t)(29 * i + 5);
		mmod32_share(&c->q_a[i], x, FNDSA_MASK_Q, &c->rng);
		mmod32_share(&c->q_b[i], y, FNDSA_MASK_Q, &c->rng);
		mmod32_share(&c->q2_a[i], x, FNDSA_MASK_Q2, &c->rng);
		mmod32_share(&c->q2_b[i], y, FNDSA_MASK_Q2, &c->rng);
		c->plain_q_a[i] = (uint16_t)(x % FNDSA_MASK_Q);
		c->plain_q_b[i] = (uint16_t)(y % FNDSA_MASK_Q);
		c->plain_q2_a[i] = (uint16_t)(x % FNDSA_MASK_Q2);
		c->plain_q2_b[i] = (uint16_t)(y % FNDSA_MASK_Q2);
		int32_t z = (int32_t)(i % 2001) - 1000;
		c->plain_crt_q[i] = (int16_t)z;
		c->plain_crt_q2[i] = (int16_t)z;
	}
	mqpoly_ext_to_int(logn, c->plain_q_a);
	mqpoly_ext_to_int(logn, c->plain_q_b);
	mq2poly_unsigned_to_mod(logn, c->plain_q2_a);
	mq2poly_unsigned_to_mod(logn, c->plain_q2_b);
	c->raw.share[0] = FNDSA_MASK_QQ2 - 1;
	c->raw.share[1] = FNDSA_MASK_QQ2_H + 1;
}

static void
free_mq_backend_ctx(mq_backend_ctx *c)
{
	free(c->plain_crt_d);
	free(c->plain_crt_q2);
	free(c->plain_crt_q);
	free(c->plain_q2_b);
	free(c->plain_q2_a);
	free(c->plain_q_b);
	free(c->plain_q_a);
	free(c->crt);
	free(c->q2_d);
	free(c->q2_b);
	free(c->q2_a);
	free(c->q_d);
	free(c->q_b);
	free(c->q_a);
	memset(c, 0, sizeof *c);
}

static void proc_plain_ntt_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mqpoly_int_to_ntt(c->logn, c->plain_q_a);
	}
	bench_sink ^= c->plain_q_a[0];
}

static void proc_plain_intt_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mqpoly_ntt_to_int(c->logn, c->plain_q_a);
	}
	bench_sink ^= c->plain_q_a[0];
}

static void proc_plain_ntt_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mq2poly_NTT(c->logn, c->plain_q2_a);
	}
	bench_sink ^= c->plain_q2_a[0];
}

static void proc_plain_intt_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mq2poly_iNTT(c->logn, c->plain_q2_a);
	}
	bench_sink ^= c->plain_q2_a[0];
}

static void proc_plain_mul_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mqpoly_mul_ntt(c->logn, c->plain_q_a, c->plain_q_b);
	}
	bench_sink ^= c->plain_q_a[0];
}

static void proc_plain_mul_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mq2poly_mul_ntt(c->logn, c->plain_q2_a, c->plain_q2_b);
	}
	bench_sink ^= c->plain_q2_a[0];
}

static void proc_plain_gram_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mqpoly_muladj_x2_ntt(c->logn, c->plain_q_a, c->plain_q_b);
	}
	bench_sink ^= c->plain_q_a[0];
}

static void proc_plain_gram_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mq2poly_muladj_x2_ntt(c->logn, c->plain_q2_a, c->plain_q2_b);
	}
	bench_sink ^= c->plain_q2_a[0];
}

static void proc_plain_crt(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mq2mqpoly_CRT(c->logn, c->plain_crt_d,
			c->plain_crt_q, c->plain_crt_q2);
	}
	bench_sink ^= (uint32_t)c->plain_crt_d[0];
}

static void proc_mq_ntt_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		(void)mmpoly_ntt(c->logn, c->q_a, FNDSA_MASK_MOD_Q, &c->rng);
	}
	bench_sink ^= c->q_a[0].share[0];
}

static void proc_mq_intt_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		(void)mmpoly_intt(c->logn, c->q_a, FNDSA_MASK_MOD_Q, &c->rng);
	}
	bench_sink ^= c->q_a[0].share[1];
}

static void proc_mq_ntt_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		(void)mmpoly_ntt(c->logn, c->q2_a, FNDSA_MASK_MOD_Q2, &c->rng);
	}
	bench_sink ^= c->q2_a[0].share[0];
}

static void proc_mq_intt_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		(void)mmpoly_intt(c->logn, c->q2_a, FNDSA_MASK_MOD_Q2, &c->rng);
	}
	bench_sink ^= c->q2_a[0].share[1];
}

static void proc_mq_mul_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_mul_ntt(c->logn, c->q_d, c->q_a, c->q_b,
			FNDSA_MASK_MOD_Q, &c->rng);
	}
	bench_sink ^= c->q_d[0].share[0];
}

static void proc_mq_mul_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_mul_ntt(c->logn, c->q2_d, c->q2_a, c->q2_b,
			FNDSA_MASK_MOD_Q2, &c->rng);
	}
	bench_sink ^= c->q2_d[0].share[0];
}

static void proc_mq_gram_q(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_gram_ntt(c->logn, c->q_d, c->q_a, c->q_b,
			FNDSA_MASK_MOD_Q, &c->rng);
	}
	bench_sink ^= c->q_d[0].share[1];
}

static void proc_mq_gram_q2(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_gram_ntt(c->logn, c->q2_d, c->q2_a, c->q2_b,
			FNDSA_MASK_MOD_Q2, &c->rng);
	}
	bench_sink ^= c->q2_d[0].share[1];
}

static void proc_mq_crt(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_crt_center(c->logn, c->crt, c->q_a, c->q2_a, &c->rng);
	}
	bench_sink ^= c->crt[0].share[0];
}

static void proc_mq_center(void *v, unsigned reps)
{
	mq_backend_ctx *c = v;
	mi32 d;
	for (unsigned i = 0; i < reps; i ++) {
		mi32_center_mod_m(&d, &c->raw, &c->rng);
	}
	bench_sink ^= d.share[1];
}

static void
print_mq_backend_group(unsigned logn)
{
	mq_backend_ctx c;
	char name[64];
	init_mq_backend_ctx(&c, logn);
	printf("\nMasked modular polynomial backend (FN-DSA-%u)\n", 1u << logn);
	printf("Unmasked Falcon baseline on the same host\n");
	snprintf(name, sizeof name, "plain NTT q logn=%u", logn);
	print_slow_bench(name, proc_plain_ntt_q, &c);
	snprintf(name, sizeof name, "plain iNTT q logn=%u", logn);
	print_slow_bench(name, proc_plain_intt_q, &c);
	snprintf(name, sizeof name, "plain NTT q-prime logn=%u", logn);
	print_slow_bench(name, proc_plain_ntt_q2, &c);
	snprintf(name, sizeof name, "plain iNTT q-prime logn=%u", logn);
	print_slow_bench(name, proc_plain_intt_q2, &c);
	snprintf(name, sizeof name, "plain pointwise q logn=%u", logn);
	print_slow_bench(name, proc_plain_mul_q, &c);
	snprintf(name, sizeof name, "plain pointwise q-prime logn=%u", logn);
	print_slow_bench(name, proc_plain_mul_q2, &c);
	snprintf(name, sizeof name, "plain Gram q logn=%u", logn);
	print_slow_bench(name, proc_plain_gram_q, &c);
	snprintf(name, sizeof name, "plain Gram q-prime logn=%u", logn);
	print_slow_bench(name, proc_plain_gram_q2, &c);
	snprintf(name, sizeof name, "plain CRT logn=%u", logn);
	print_slow_bench(name, proc_plain_crt, &c);
	printf("Masked Raccoon-style sharewise backend\n");
	snprintf(name, sizeof name, "masked NTT q logn=%u", logn);
	print_slow_bench(name, proc_mq_ntt_q, &c);
	snprintf(name, sizeof name, "masked iNTT q logn=%u", logn);
	print_slow_bench(name, proc_mq_intt_q, &c);
	snprintf(name, sizeof name, "masked NTT q-prime logn=%u", logn);
	print_slow_bench(name, proc_mq_ntt_q2, &c);
	snprintf(name, sizeof name, "masked iNTT q-prime logn=%u", logn);
	print_slow_bench(name, proc_mq_intt_q2, &c);
	snprintf(name, sizeof name, "masked pointwise q logn=%u", logn);
	print_slow_bench(name, proc_mq_mul_q, &c);
	snprintf(name, sizeof name, "masked pointwise q-prime logn=%u", logn);
	print_slow_bench(name, proc_mq_mul_q2, &c);
	snprintf(name, sizeof name, "masked Gram q logn=%u", logn);
	print_slow_bench(name, proc_mq_gram_q, &c);
	snprintf(name, sizeof name, "masked Gram q-prime logn=%u", logn);
	print_slow_bench(name, proc_mq_gram_q2, &c);
	snprintf(name, sizeof name, "masked CRT+center logn=%u", logn);
	print_slow_bench(name, proc_mq_crt, &c);
	print_bench("masked center one coefficient", proc_mq_center, &c, 20);
	free_mq_backend_ctx(&c);
}


typedef struct {
	fndsa_mask_rng rng;
	mrx160 a;
	mrx160 b;
	mrx160 d;
	mb160 ba;
	mrx160 fx;
	mrx160 fccs;
	mrx160 fsqrt;
	mb64 d64;
	mbit bit;
	uint8_t stream_byte;
	uint64_t stream_state;
	unsigned stream_count;
	unsigned sampler_logn;
	uint32_t x[3];
	uint32_t y[3];
	uint32_t opened[3];
	mmod32 ma;
	mmod32 mb;
	mmod32 md;
	uint32_t modulus;
} gadget_ctx;

static void
init_gadget_ctx(gadget_ctx *c)
{
	memset(c, 0, sizeof *c);
	c->x[0] = UINT32_C(0x12345678);
	c->x[1] = UINT32_C(0x00030000);
	c->x[2] = 2;
	c->y[0] = UINT32_C(0x9ABCDEF0);
	c->y[1] = UINT32_C(0x00010000);
	c->y[2] = 1;
	fndsa_mask_rng_init_demo(&c->rng);
	uint32_t fs[3] = { 0, 0, 10000 };
	mrx160_share_rx(&c->fsqrt, fs, &c->rng);
	mrx160_share_rx(&c->a, c->x, &c->rng);
	mrx160_share_rx(&c->b, c->y, &c->rng);
	mrx160_a2b(&c->ba, &c->a, &c->rng);
	uint32_t fx[3] = { 0, UINT32_C(0x40000000), 0 };
	uint32_t fc[3] = { 0, UINT32_C(0xC0000000), 0 };
	mrx160_share_rx(&c->fx, fx, &c->rng);
	mrx160_share_rx(&c->fccs, fc, &c->rng);
	c->stream_byte = 0;
	c->stream_state = UINT64_C(0xD31947A6C85E20BF);
	c->sampler_logn = 9;
}

static void
proc_rng_init(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		fndsa_mask_rng_init_demo(&c->rng);
	}
	bench_sink ^= c->rng.s[0];
}

static void
proc_rng_u64(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint64_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_mask_rng_next_u64(&c->rng);
	}
	bench_sink ^= x;
}

static void
proc_rng_u32(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint32_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_mask_rng_next_u32(&c->rng);
	}
	bench_sink ^= x;
}

static void
proc_rng_bit(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint32_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_mask_rng_next_bit(&c->rng);
	}
	bench_sink ^= x;
}

static void
proc_rng_mod(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint32_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_mask_rng_mod(&c->rng, c->modulus);
	}
	bench_sink ^= x;
}

static void
proc_mrx_share(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_share_rx(&c->d, c->x, &c->rng);
	}
	bench_sink ^= c->d.share[0][0];
}

static void
proc_mrx_sqr(void *v, unsigned reps)
{
        gadget_ctx *c = v;
        for (unsigned i = 0; i < reps; i ++) {
                mrx160_sqr(&c->d, &c->a, &c->rng);
        }
        bench_sink ^= c->d.share[1][0];
}


static void
proc_mrx_open(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint32_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_open_rx(c->opened, &c->a);
		x ^= c->opened[i % 3];
	}
	bench_sink ^= x;
}

static void
proc_mrx_copy(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_copy(&c->d, &c->a);
	}
	bench_sink ^= c->d.share[1][1];
}

static void
proc_mrx_add(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_add(&c->d, &c->a, &c->b);
	}
	bench_sink ^= c->d.share[0][2];
}

static void
proc_mrx_sub(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_sub(&c->d, &c->a, &c->b);
	}
	bench_sink ^= c->d.share[1][2];
}

static void
proc_mrx_neg(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_neg(&c->d, &c->a);
	}
	bench_sink ^= c->d.share[0][3];
}

static void
proc_mrx_refresh(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_refresh(&c->a, &c->rng);
	}
	bench_sink ^= c->a.share[1][3];
}

static void
proc_mrx_a2b(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_a2b(&c->ba, &c->a, &c->rng);
	}
	bench_sink ^= c->ba.share[0][4];
}

static void
proc_mrx_b2a(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_b2a(&c->d, &c->ba, &c->rng);
	}
	bench_sink ^= c->d.share[1][4];
}

static void
proc_mrx_mul(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul(&c->d, &c->a, &c->b, &c->rng);
	}
	bench_sink ^= c->d.share[0][0];
}

static void
proc_mrx_floor(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_floor(&c->d, &c->a, &c->rng);
	}
	bench_sink ^= c->d.share[0][2];
}

static void
proc_mrx_division(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_div_spec(&c->d, &c->a, &c->b, &c->rng);
	}
	bench_sink ^= c->d.share[0][1];
}

static void
proc_mrx_sqrt(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_sqrt_spec(&c->d, &c->fsqrt, &c->rng);
	}
	bench_sink ^= c->d.share[1][1];
}

static void
proc_mrx_inverse(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_inv_spec(&c->d, &c->a, &c->rng);
	}
	bench_sink ^= c->d.share[1][1];
}

static void
proc_mrx_mul_u64_hi(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul_u64_hi(&c->d64, &c->fx, &c->fccs, &c->rng);
	}
	bench_sink ^= c->d64.share[0][0];
}

static void
proc_mrx_expm(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_expm_p63(&c->d64, &c->fx, &c->fccs, &c->rng);
	}
	bench_sink ^= c->d64.share[1][1];
}

static uint8_t
bench_constant_u8(void *v)
{
	gadget_ctx *c = v;
	c->stream_count ++;
	return c->stream_byte;
}

static uint8_t
bench_random_u8(void *v)
{
	gadget_ctx *c = v;
	uint64_t x = c->stream_state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	c->stream_state = x;
	c->stream_count ++;
	return (uint8_t)x;
}

static void
proc_mrx_ber_accept(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	fndsa_mask_sampler_source src = { c, bench_constant_u8, 9 };
	c->stream_byte = 0;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_ber_exp(&c->bit, &src, &c->fx, &c->fccs, &c->rng);
	}
	bench_sink ^= c->bit.share[0] ^ c->stream_count;
}

static void
proc_mrx_ber_reject(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	fndsa_mask_sampler_source src = { c, bench_constant_u8, 9 };
	c->stream_byte = UINT8_C(0xFF);
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_ber_exp(&c->bit, &src, &c->fx, &c->fccs, &c->rng);
	}
	bench_sink ^= c->bit.share[1] ^ c->stream_count;
}

static void
proc_mrx_sampler_z(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	fndsa_mask_sampler_source src = { c, bench_random_u8, c->sampler_logn };
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_sampler_z(&c->d, &src, &c->a, &c->fccs, &c->rng);
	}
	bench_sink ^= c->d.share[0][2] ^ c->stream_count;
}

static void
proc_mrx_mul_wide(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul_wide(&c->d, &c->a, &c->b, &c->rng);
	}
	bench_sink ^= c->d.share[1][0];
}

static void
proc_mrx_mul_public_fast(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul_public_fast(&c->d, &c->a, c->y, &c->rng);
	}
	bench_sink ^= c->d.share[0][1];
}

static void
proc_mrx_mul_public_wide(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mrx160_mul_public_wide(&c->d, &c->a, c->y, &c->rng);
	}
	bench_sink ^= c->d.share[1][1];
}

static void
init_mod_ctx(gadget_ctx *c, uint32_t modulus)
{
	c->modulus = modulus;
	mmod32_share(&c->ma, 1234 % modulus, modulus, &c->rng);
	mmod32_share(&c->mb, 5678 % modulus, modulus, &c->rng);
}

static void
proc_mod_share(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmod32_share(&c->md, (1234u + i) % c->modulus,
			c->modulus, &c->rng);
	}
	bench_sink ^= c->md.share[0];
}

static void
proc_mod_open(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	uint32_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= mmod32_open(&c->ma, c->modulus);
	}
	bench_sink ^= x;
}

static void
proc_mod_add(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmod32_add(&c->md, &c->ma, &c->mb, c->modulus);
	}
	bench_sink ^= c->md.share[0];
}

static void
proc_mod_sub(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmod32_sub(&c->md, &c->ma, &c->mb, c->modulus);
	}
	bench_sink ^= c->md.share[1];
}

static void
proc_mod_mul(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmod32_mul(&c->md, &c->ma, &c->mb, c->modulus, &c->rng);
	}
	bench_sink ^= c->md.share[0];
}

static void
proc_mod_refresh(void *v, unsigned reps)
{
	gadget_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mmod32_refresh(&c->ma, c->modulus, &c->rng);
	}
	bench_sink ^= c->ma.share[1];
}

typedef struct {
	unsigned logn;
	size_t n;
	size_t hn;
	fndsa_mask_rng rng;
	mi32 *coeff;
	mi32 *scratch;
	mrx160 *fft;
	mrx160 *a0;
	mrx160 *b1;
	mrx160 *ld;
	mrx160 *le;
	void *tree_tmp;
	size_t tree_tmp_len;
	fndsa_mask_tree_state tree;
	uint64_t stream_state;
	uint64_t stream_count;
} backend_ctx;

static uint8_t
backend_next_u8(void *v)
{
	backend_ctx *c = v;
	uint64_t x = c->stream_state;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	c->stream_state = x;
	c->stream_count ++;
	return (uint8_t)x;
}

static void
init_backend_ctx(backend_ctx *c, unsigned logn)
{
	memset(c, 0, sizeof *c);
	c->logn = logn;
	c->n = (size_t)1 << logn;
	c->hn = c->n >> 1;
	c->coeff = malloc(c->n * sizeof *c->coeff);
	c->scratch = malloc((c->n + c->hn) * sizeof *c->scratch);
	c->fft = malloc(c->n * sizeof *c->fft);
	c->a0 = malloc(c->hn * sizeof *c->a0);
	c->b1 = malloc(c->hn * sizeof *c->b1);
	c->ld = malloc(c->hn * sizeof *c->ld);
	c->le = malloc(c->hn * sizeof *c->le);
	c->tree_tmp_len = fndsa_masked_dyn_tmp_size(logn);
	c->tree_tmp = malloc(c->tree_tmp_len);
	if (c->coeff == NULL || c->scratch == NULL || c->fft == NULL
		|| c->a0 == NULL || c->b1 == NULL || c->ld == NULL
		|| c->le == NULL || c->tree_tmp == NULL)
	{
		fprintf(stderr, "backend benchmark allocation failure\n");
		exit(EXIT_FAILURE);
	}
	fndsa_mask_rng_init_demo(&c->rng);
	for (size_t i = 0; i < c->n; i ++) {
		mi32_share(&c->coeff[i], (int32_t)(i % 31) - 15, &c->rng);
	}
	rx ra, rb;
	for (size_t i = 0; i < c->hn; i ++) {
		rx_of32(ra, 100 + (int32_t)(i & 15));
		rx_of32(rb, (int32_t)(i % 7) - 3);
		mrx160_share_rx(&c->a0[i], ra, &c->rng);
		mrx160_share_rx(&c->b1[i], rb, &c->rng);
	}
	c->stream_state = UINT64_C(0xA31E7C5D942B680F) ^ logn;
	fndsa_mask_sampler_source source = {
		c, backend_next_u8, logn
	};
	if (!fndsa_mask_tree_bind(&c->tree, logn, &source, &c->rng,
		c->tree_tmp, c->tree_tmp_len))
	{
		fprintf(stderr, "backend benchmark tree bind failure\n");
		exit(EXIT_FAILURE);
	}
}

static void
prepare_backend_tree(backend_ctx *c)
{
	mrx160 *w = fndsa_mask_tree_work(&c->tree);
	mi32 *dz = fndsa_mask_tree_samples(&c->tree);
	memset(w, 0, 2 * c->n * sizeof *w);
	memset(dz, 0, c->n * sizeof *dz);
	for (size_t i = 0; i < c->hn; i ++) {
		w[c->n + i].share[0][2] = 10000;
	}
}

static void
free_backend_ctx(backend_ctx *c)
{
	fndsa_mask_tree_clear(&c->tree);
	free(c->tree_tmp);
	free(c->le);
	free(c->ld);
	free(c->b1);
	free(c->a0);
	free(c->fft);
	free(c->scratch);
	free(c->coeff);
	memset(c, 0, sizeof *c);
}

static void
proc_backend_fft(void *v, unsigned reps)
{
	backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mxpoly_fft_i32(c->logn, c->fft, c->coeff,
			c->scratch, &c->rng);
	}
	bench_sink ^= c->fft[0].share[0][0];
}

static void
proc_backend_ldl(void *v, unsigned reps)
{
	backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		mxpoly_ldl_fft(c->logn, c->ld, c->le,
			c->a0, c->b1, &c->rng);
	}
	bench_sink ^= c->ld[0].share[1][0] ^ c->le[0].share[0][0];
}

static void
proc_backend_tree_left(void *v, unsigned reps)
{
	backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		prepare_backend_tree(c);
		fndsa_mask_ffsampling_left(&c->tree, FNDSA_MASK_TREE_NOOUT);
	}
	bench_sink ^= c->tree.work[0].share[0][0] ^ c->stream_count;
}

static void
proc_backend_tree_right(void *v, unsigned reps)
{
	backend_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		prepare_backend_tree(c);
		fndsa_mask_ffsampling_right(&c->tree, FNDSA_MASK_TREE_INIT0);
	}
	bench_sink ^= c->tree.work[0].share[1][0] ^ c->stream_count;
}

static void
print_backend_group(unsigned logn)
{
	backend_ctx c;
	char name[64];
	init_backend_ctx(&c, logn);
	printf("\nMasked dynamic backend (FN-DSA-%u)\n", 1u << logn);
	snprintf(name, sizeof name, "masked FFT logn=%u", logn);
	print_slow_bench(name, proc_backend_fft, &c);
	snprintf(name, sizeof name, "masked LDL logn=%u", logn);
	print_slow_bench(name, proc_backend_ldl, &c);
	snprintf(name, sizeof name, "masked ffSampling left NOOUT logn=%u", logn);
	print_slow_bench(name, proc_backend_tree_left, &c);
	snprintf(name, sizeof name, "masked ffSampling right INIT0 logn=%u", logn);
	print_slow_bench(name, proc_backend_tree_right, &c);
	free_backend_ctx(&c);
}

typedef struct {
	unsigned logn;
	size_t sk_len;
	size_t vk_len;
	size_t mk_len;
	size_t sig_len;
	size_t tmp_plain_len;
	size_t tmp_masked_len;
	uint8_t *sk;
	uint8_t *vk;
	uint8_t *mk;
	uint8_t *sig;
	uint8_t *tmp_plain;
	uint8_t *tmp_masked;
	size_t n;
	fndsa_mask_rng stage_rng;
	mmod32 *stage_basis_mem;
	mmod32 *stage_work;
	uint32_t *stage_h;
	mi32 *stage_vec;
	fndsa_masked_basis_mq stage_basis;
	uint8_t sign_seed[40];
} key_ctx;

static void
init_key_ctx(key_ctx *c, unsigned logn)
{
	memset(c, 0, sizeof *c);
	c->logn = logn;
	c->n = (size_t)1 << logn;
	c->sk_len = FNDSA_SIGN_KEY_SIZE(logn);
	c->vk_len = FNDSA_VRFY_KEY_SIZE(logn);
	c->mk_len = fndsa_masked_key_size(logn);
	c->sig_len = FNDSA_SIGNATURE_SIZE(logn);
	c->tmp_plain_len = FNDSA_SIGN_TMP_SIZE(logn);
	c->tmp_masked_len = fndsa_masked_tmp_size(logn);
	c->sk = malloc(c->sk_len);
	c->vk = malloc(c->vk_len);
	c->mk = malloc(c->mk_len);
	c->sig = malloc(c->sig_len);
	c->tmp_plain = malloc(c->tmp_plain_len);
	c->tmp_masked = malloc(c->tmp_masked_len);
	c->stage_basis_mem = malloc(8 * c->n * sizeof *c->stage_basis_mem);
	c->stage_work = malloc(4 * c->n * sizeof *c->stage_work);
	c->stage_h = malloc(c->n * sizeof *c->stage_h);
	c->stage_vec = malloc(4 * c->n * sizeof *c->stage_vec);
	if (c->sk == NULL || c->vk == NULL || c->mk == NULL
		|| c->sig == NULL || c->tmp_plain == NULL
		|| c->tmp_masked == NULL || c->stage_basis_mem == NULL
		|| c->stage_work == NULL || c->stage_h == NULL
		|| c->stage_vec == NULL)
	{
		fprintf(stderr, "allocation failure\n");
		exit(EXIT_FAILURE);
	}

	c->stage_basis = (fndsa_masked_basis_mq) {
		c->stage_basis_mem + 0 * c->n, c->stage_basis_mem + 1 * c->n,
		c->stage_basis_mem + 2 * c->n, c->stage_basis_mem + 3 * c->n,
		c->stage_basis_mem + 4 * c->n, c->stage_basis_mem + 5 * c->n,
		c->stage_basis_mem + 6 * c->n, c->stage_basis_mem + 7 * c->n
	};
	uint8_t key_seed[48];
	for (size_t i = 0; i < sizeof key_seed; i ++) {
		key_seed[i] = (uint8_t)(17 * logn + 29 * i + 3);
	}
	for (size_t i = 0; i < sizeof c->sign_seed; i ++) {
		c->sign_seed[i] = (uint8_t)(0xA5 ^ logn ^ i);
	}
	fndsa_keygen_seeded(logn, key_seed, sizeof key_seed, c->sk, c->vk);
	if (!fndsa_masked_key_import_demo(
		c->sk, c->sk_len, c->mk, c->mk_len))
	{
		fprintf(stderr, "masked-key import failure\n");
		exit(EXIT_FAILURE);
	}
}

static void
free_key_ctx(key_ctx *c)
{
	free(c->stage_vec);
	free(c->stage_h);
	free(c->stage_work);
	free(c->stage_basis_mem);
	fndsa_masked_key_clear(c->mk, c->mk_len);
	free(c->sk);
	free(c->vk);
	free(c->mk);
	free(c->sig);
	free(c->tmp_plain);
	free(c->tmp_masked);
	memset(c, 0, sizeof *c);
}

static void
proc_key_expand(void *v, unsigned reps)
{
	key_ctx *c = v;
	unsigned x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		const uint8_t *vk;
		size_t vk_len;
		x ^= (unsigned)fndsa_masked_key_expand_mq_internal(
			c->mk, c->mk_len, &c->stage_basis,
			&vk, &vk_len, &c->stage_rng);
		x ^= (unsigned)vk_len;
	}
	bench_sink ^= x ^ c->stage_basis_mem[0].share[0];
}

static int
prepare_signing_stages(key_ctx *c)
{
	const uint8_t *vk;
	size_t vk_len;
	fndsa_mask_rng_init_demo(&c->stage_rng);
	if (!fndsa_masked_key_expand_mq_internal(
		c->mk, c->mk_len, &c->stage_basis, &vk, &vk_len,
		&c->stage_rng))
	{
		return 0;
	}
	mmod32 *bp[8] = {
		c->stage_basis.f_q, c->stage_basis.g_q,
		c->stage_basis.F_q, c->stage_basis.G_q,
		c->stage_basis.f_q2, c->stage_basis.g_q2,
		c->stage_basis.F_q2, c->stage_basis.G_q2
	};
	for (unsigned j = 0; j < 8; j ++) {
		fndsa_mask_modulus mod = j < 4
			? FNDSA_MASK_MOD_Q : FNDSA_MASK_MOD_Q2;
		if (!mmpoly_ntt(c->logn, bp[j], mod, &c->stage_rng)) {
			return 0;
		}
	}
	uint16_t *hp = (uint16_t *)(void *)c->stage_h;
	if (mqpoly_decode(c->logn, vk + 1, hp) != vk_len - 1) {
		return 0;
	}
	for (size_t i = c->n; i -- > 0;) {
		c->stage_h[i] = mmod32_encode_public(hp[i], FNDSA_MASK_Q);
	}
	for (size_t i = 0; i < 4 * c->n; i ++) {
		mi32_share(&c->stage_vec[i], 0, &c->stage_rng);
	}
	return 1;
}

static void
proc_final_vectors(void *v, unsigned reps)
{
	key_ctx *c = v;
	mmod32 *w0 = c->stage_work + 0 * c->n;
	mmod32 *w1 = c->stage_work + 1 * c->n;
	mmod32 *w2 = c->stage_work + 2 * c->n;
	mmod32 *w3 = c->stage_work + 3 * c->n;
	mi32 *u0 = c->stage_vec + 0 * c->n;
	mi32 *u1 = c->stage_vec + 1 * c->n;
	mi32 *s1 = c->stage_vec + 2 * c->n;
	mi32 *s0 = c->stage_vec + 3 * c->n;
	for (unsigned i = 0; i < reps; i ++) {
		mmpoly_from_mi32(c->logn, w0, u0, FNDSA_MASK_MOD_Q2,
			&c->stage_rng);
		mmpoly_from_mi32(c->logn, w1, u1, FNDSA_MASK_MOD_Q2,
			&c->stage_rng);
		(void)mmpoly_ntt(c->logn, w0, FNDSA_MASK_MOD_Q2,
			&c->stage_rng);
		(void)mmpoly_ntt(c->logn, w1, FNDSA_MASK_MOD_Q2,
			&c->stage_rng);
		mmpoly_mul_ntt(c->logn, w2, c->stage_basis.f_q2, w0,
			FNDSA_MASK_MOD_Q2, &c->stage_rng);
		mmpoly_mul_ntt(c->logn, w3, c->stage_basis.F_q2, w1,
			FNDSA_MASK_MOD_Q2, &c->stage_rng);
		mmpoly_add(c->logn, w2, w2, w3, FNDSA_MASK_MOD_Q2);
		mmpoly_refresh(c->logn, w2, FNDSA_MASK_MOD_Q2, &c->stage_rng);
		(void)mmpoly_intt_center(c->logn, s1, w2,
			FNDSA_MASK_MOD_Q2, &c->stage_rng);

		mmpoly_from_mi32(c->logn, w0, s1, FNDSA_MASK_MOD_Q,
			&c->stage_rng);
		(void)mmpoly_ntt(c->logn, w0, FNDSA_MASK_MOD_Q,
			&c->stage_rng);
		mmpoly_mul_public_ntt(c->logn, w1, w0, c->stage_h,
			FNDSA_MASK_MOD_Q, &c->stage_rng);
		mmpoly_neg(c->logn, w1, w1, FNDSA_MASK_MOD_Q);
		(void)mmpoly_intt_center(c->logn, s0, w1,
			FNDSA_MASK_MOD_Q, &c->stage_rng);
	}
	bench_sink ^= s0[0].share[0] ^ s1[0].share[1];
}

static void
proc_signature_acceptance(void *v, unsigned reps)
{
	key_ctx *c = v;
	mi32 *s1 = c->stage_vec + 2 * c->n;
	mi32 *s0 = c->stage_vec + 3 * c->n;
	mbit accept = { { 0, 0 } };
	for (unsigned i = 0; i < reps; i ++) {
		mi32_signature_is_acceptable(&accept, c->logn, s0, s1,
			c->sig_len - 41, &c->stage_rng);
	}
	bench_sink ^= accept.share[0] ^ accept.share[1];
}

static void
proc_key_size(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_masked_key_size(c->logn);
	}
	bench_sink ^= x;
}

static void
proc_tmp_size(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_masked_tmp_size(c->logn);
	}
	bench_sink ^= x;
}

static void
proc_key_import(void *v, unsigned reps)
{
	key_ctx *c = v;
	unsigned x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= (unsigned)fndsa_masked_key_import_demo(
			c->sk, c->sk_len, c->mk, c->mk_len);
	}
	bench_sink ^= x ^ c->mk[16];
}

static void
proc_key_validate(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_sign_masked_seeded_temp(
			c->mk, c->mk_len, NULL, 0,
			FNDSA_HASH_ID_RAW, "test", 4,
			c->sign_seed, sizeof c->sign_seed,
			NULL, 0, NULL, 0);
	}
	bench_sink ^= x;
}

static void
proc_sign_plain(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_sign_seeded_temp(
			c->sk, c->sk_len, NULL, 0,
			FNDSA_HASH_ID_RAW, "test", 4,
			c->sign_seed, sizeof c->sign_seed,
			c->sig, c->sig_len, c->tmp_plain, c->tmp_plain_len);
	}
	bench_sink ^= x ^ c->sig[1];
}

static void
proc_sign_masked_seeded(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_sign_masked_seeded_temp(
			c->mk, c->mk_len, NULL, 0,
			FNDSA_HASH_ID_RAW, "test", 4,
			c->sign_seed, sizeof c->sign_seed,
			c->sig, c->sig_len, c->tmp_masked, c->tmp_masked_len);
	}
	bench_sink ^= x ^ c->sig[2];
}

static void
proc_sign_masked_random(void *v, unsigned reps)
{
	key_ctx *c = v;
	size_t x = 0;
	for (unsigned i = 0; i < reps; i ++) {
		x ^= fndsa_sign_masked_temp(
			c->mk, c->mk_len, NULL, 0,
			FNDSA_HASH_ID_RAW, "test", 4,
			c->sig, c->sig_len, c->tmp_masked, c->tmp_masked_len);
	}
	bench_sink ^= x ^ c->sig[3];
}

static void
proc_key_clear(void *v, unsigned reps)
{
	key_ctx *c = v;
	for (unsigned i = 0; i < reps; i ++) {
		fndsa_masked_key_clear(c->mk, c->mk_len);
	}
	bench_sink ^= c->mk[c->mk_len - 1];
}

static double
measure_key(bench_proc proc, key_ctx *ctx, unsigned reps)
{
	return measure(proc, ctx, reps, 51, 10);
}

static void
print_key_bench(const char *name, bench_proc proc, key_ctx *ctx,
	unsigned reps)
{
	printf("%-38s %13.2f %s/op\n",
		name, measure_key(proc, ctx, reps),
		use_tsc ? "TSC-ticks" : "cycles");
}

static void
print_mod_group(gadget_ctx *c, uint32_t modulus)
{
	char name[64];
	init_mod_ctx(c, modulus);
	printf("\nModulo %u gadgets\n", modulus);
	snprintf(name, sizeof name, "mmod32_share mod %u", modulus);
	print_bench(name, proc_mod_share, c, 1000);
	snprintf(name, sizeof name, "mmod32_open mod %u", modulus);
	print_bench(name, proc_mod_open, c, 2000);
	snprintf(name, sizeof name, "mmod32_add mod %u", modulus);
	print_bench(name, proc_mod_add, c, 2000);
	snprintf(name, sizeof name, "mmod32_sub mod %u", modulus);
	print_bench(name, proc_mod_sub, c, 2000);
	snprintf(name, sizeof name, "mmod32_mul mod %u", modulus);
	print_bench(name, proc_mod_mul, c, 500);
	snprintf(name, sizeof name, "mmod32_refresh mod %u", modulus);
	print_bench(name, proc_mod_refresh, c, 1000);
}

static void
print_key_group(unsigned logn)
{
	key_ctx c;
	char name[64];
	init_key_ctx(&c, logn);
	printf("\nMasked-key and signing API (FN-DSA-%u)\n", 1u << logn);

	snprintf(name, sizeof name, "fndsa_masked_key_size logn=%u", logn);
	print_key_bench(name, proc_key_size, &c, 2000);
	snprintf(name, sizeof name, "fndsa_masked_tmp_size logn=%u", logn);
	print_key_bench(name, proc_tmp_size, &c, 2000);
	snprintf(name, sizeof name, "masked key import logn=%u", logn);
	print_key_bench(name, proc_key_import, &c, 1);
	snprintf(name, sizeof name, "masked basis expand logn=%u", logn);
	fndsa_mask_rng_init_demo(&c.stage_rng);
	print_key_bench(name, proc_key_expand, &c, 1);
	if (!prepare_signing_stages(&c)) {
		fprintf(stderr, "masked signing-stage setup failure\n");
		exit(EXIT_FAILURE);
	}
	snprintf(name, sizeof name, "masked final vectors logn=%u", logn);
	print_key_bench(name, proc_final_vectors, &c, 1);
	snprintf(name, sizeof name, "masked aggregate accept logn=%u", logn);
	print_key_bench(name, proc_signature_acceptance, &c, 1);
	snprintf(name, sizeof name, "masked key validate logn=%u", logn);
	print_key_bench(name, proc_key_validate, &c, 1);
	snprintf(name, sizeof name, "plain seeded sign logn=%u", logn);
	print_key_bench(name, proc_sign_plain, &c, 1);
	snprintf(name, sizeof name, "masked seeded sign logn=%u", logn);
	print_key_bench(name, proc_sign_masked_seeded, &c, 1);
	snprintf(name, sizeof name, "masked random sign logn=%u", logn);
	print_key_bench(name, proc_sign_masked_random, &c, 1);
	snprintf(name, sizeof name, "masked key clear logn=%u", logn);
	print_key_bench(name, proc_key_clear, &c, 1);

	free_key_ctx(&c);
}

static void
usage(void)
{
	fprintf(stderr, "usage: speed_masked [ -tsc ]\n");
	fprintf(stderr,
		"default: hardware core cycles (usually requires OS permission)\n");
	fprintf(stderr,
		"-tsc: userspace reference-counter ticks, not core cycles\n");
	exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
	if (argc == 2 && strcmp(argv[1], "-tsc") == 0) {
		use_tsc = 1;
	} else if (argc != 1) {
		usage();
	}

	counter_init();

	gadget_ctx c;
	init_gadget_ctx(&c);

	printf("Experimental masking gadget benchmark -- DEMO ONLY\n");
	printf("Counter unit: %s\n",
		use_tsc ? "TSC/reference ticks" : "hardware core cycles");

	printf("\nMask RNG\n");
	print_bench("mask RNG init", proc_rng_init, &c, 2000);
	print_bench("mask RNG next_u64", proc_rng_u64, &c, 2000);
	print_bench("mask RNG next_u32", proc_rng_u32, &c, 4000);
	print_bench("mask RNG next_bit", proc_rng_bit, &c, 8000);
	c.modulus = 12289;
	print_bench("mask RNG mod 12289", proc_rng_mod, &c, 1000);
	c.modulus = 18433;
	print_bench("mask RNG mod 18433", proc_rng_mod, &c, 1000);
	c.modulus = UINT32_C(226523137);
	print_bench("mask RNG mod 226523137", proc_rng_mod, &c, 1000);

	printf("\n160-bit arithmetic and conversion gadgets\n");
	print_bench("mrx160_share_rx", proc_mrx_share, &c, 1000);
	print_bench("mrx160_open_rx", proc_mrx_open, &c, 2000);
	print_bench("mrx160_copy", proc_mrx_copy, &c, 4000);
	print_bench("mrx160_add", proc_mrx_add, &c, 2000);
	print_bench("mrx160_sub", proc_mrx_sub, &c, 2000);
	print_bench("mrx160_neg", proc_mrx_neg, &c, 2000);
	print_bench("mrx160_refresh", proc_mrx_refresh, &c, 1000);
	print_bench("mrx160_a2b", proc_mrx_a2b, &c, 200);
	print_bench("mrx160_b2a", proc_mrx_b2a, &c, 200);
	print_bench("mrx160_mul fast", proc_mrx_mul, &c, 100);
	print_bench("mrx160_mul_public_fast", proc_mrx_mul_public_fast, &c, 100);
        print_bench("mrx160_sqr", proc_mrx_sqr, &c, 100);
	print_bench("mrx160_mul_wide", proc_mrx_mul_wide, &c, 30);
	print_bench("mrx160_mul_public_wide", proc_mrx_mul_public_wide, &c, 30);
	print_bench("mrx160_floor", proc_mrx_floor, &c, 100);
	print_bench("mrx160_inv_spec", proc_mrx_inverse, &c, 5);
	print_bench("mrx160_div_spec", proc_mrx_division, &c, 5);
	print_bench("mrx160_sqrt_spec", proc_mrx_sqrt, &c, 3);
	print_bench("mrx160_mul_u64_hi", proc_mrx_mul_u64_hi, &c, 100);
	print_bench("mrx160_expm_p63", proc_mrx_expm, &c, 5);
	print_bench("mrx160_BerExp accept path", proc_mrx_ber_accept, &c, 5);
	print_bench("mrx160_BerExp reject path", proc_mrx_ber_reject, &c, 5);
	c.sampler_logn = 9;
	print_bench("mrx160_SamplerZ logn=9", proc_mrx_sampler_z, &c, 1);
	c.sampler_logn = 10;
	print_bench("mrx160_SamplerZ logn=10", proc_mrx_sampler_z, &c, 1);

	print_mod_group(&c, 12289);
	print_mod_group(&c, 18433);
	print_mod_group(&c, UINT32_C(226523137));

	print_backend_group(9);
	print_backend_group(10);
	print_mq_backend_group(9);
	print_mq_backend_group(10);

	print_key_group(9);
	print_key_group(10);

	printf("\nbenchmark sink: %llu\n",
		(unsigned long long)bench_sink);
	counter_close();
	return 0;
}
