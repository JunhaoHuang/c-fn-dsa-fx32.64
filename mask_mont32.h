#ifndef FNDSA_MASK_MONT32_H__
#define FNDSA_MASK_MONT32_H__

#include <stdint.h>

static inline int32_t
mask_mont32_signed_word(uint32_t word)
{
	return (int32_t)((int64_t)(word ^ UINT32_C(0x80000000))
		- INT64_C(2147483648));
}

static inline int32_t
mask_mont32_redc(int64_t product, uint32_t modulus, uint32_t q0i)
{
	int32_t correction = mask_mont32_signed_word((uint32_t)product * q0i);
	int64_t sum = product + (int64_t)correction * modulus;
	return mask_mont32_signed_word((uint32_t)((uint64_t)sum >> 32));
}

static inline uint32_t
mask_mont32_cadd(int32_t value, uint32_t modulus)
{
	uint32_t word = (uint32_t)value;
	return word + (modulus & ((uint32_t)0 - (word >> 31)));
}

static inline uint32_t
mask_mont32_add(uint32_t left, uint32_t right, uint32_t modulus)
{
	uint32_t word = left + right - modulus;
	return word + (modulus & ((uint32_t)0 - (word >> 31)));
}

static inline uint32_t
mask_mont32_sub(uint32_t left, uint32_t right, uint32_t modulus)
{
	uint32_t word = left - right;
	return word + (modulus & ((uint32_t)0 - (word >> 31)));
}

#endif
