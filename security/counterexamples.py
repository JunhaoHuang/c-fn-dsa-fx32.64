#!/usr/bin/env python3
"""Reduced exhaustive checks for the repaired first-order masking gadgets.

The script first reproduces the three distributions that broke 1-NI in the
old code. It then checks the invariant used by the repaired proofs: every
multi-share partial sum is masked before the first term and remains uniform.
This is an algebraic regression checker, not a compiler-level verifier.
"""

from collections import Counter
from fractions import Fraction


def total_variation(left, right):
    nl = sum(left.values())
    nr = sum(right.values())
    assert nl == nr
    keys = set(left) | set(right)
    return Fraction(sum(abs(left[k] - right[k]) for k in keys), 2 * nl)


def legacy_ring_cross_probe(secret_a, secret_b):
    observed = Counter()
    for a0 in range(2):
        a1 = secret_a ^ a0
        for b0 in range(2):
            b1 = secret_b ^ b0
            observed[(a0 & b1) ^ (a1 & b0)] += 1
    return observed


def legacy_mi32_probe(secret):
    observed = Counter()
    for a0 in range(2):
        a1 = secret ^ a0
        observed[a0 ^ a1] += 1
    return observed


def legacy_mmod_probe(secret_a, secret_b, q=5):
    observed = Counter()
    for a1 in range(q):
        a0 = (secret_a - a1) % q
        for b1 in range(q):
            b0 = (secret_b - b1) % q
            for r in range(q):
                z1 = a0 * b1 + a1 * b0 + a1 * b1 + q - r
                observed[z1] += 1
    return observed


def check_repaired_ring():
    uniform_bit = Counter({0: 1, 1: 1})
    for secret_a in range(2):
        for secret_b in range(2):
            for a0 in range(2):
                a1 = secret_a ^ a0
                for b0 in range(2):
                    b1 = secret_b ^ b0
                    probes = [Counter() for _ in range(7)]
                    for r in range(2):
                        p00 = a0 & b0
                        p01 = a0 & b1
                        p10 = a1 & b0
                        p11 = a1 & b1
                        probes[0][r] += 1
                        probes[1][r] += 1
                        probes[2][r ^ p00] += 1
                        probes[3][r ^ p01] += 1
                        probes[4][r ^ p01 ^ p10] += 1
                        probes[5][r ^ p01 ^ p10 ^ p11] += 1
                        probes[6][r ^ p00] += 1
                    assert all(p == uniform_bit for p in probes)


def check_repaired_mi32():
    for secret in range(2):
        for a0 in range(2):
            a1 = secret ^ a0
            out0 = Counter()
            out1 = Counter()
            for rx in range(2):
                for ry in range(2):
                    for refresh in range(2):
                        t0 = (a0 ^ rx) ^ (a1 ^ ry)
                        t1 = rx ^ ry
                        out0[t0 ^ refresh] += 1
                        out1[t1 ^ refresh] += 1
            assert out0 == Counter({0: 4, 1: 4})
            assert out1 == Counter({0: 4, 1: 4})


def check_repaired_mmod(q=5):
    uniform_q = Counter({x: 1 for x in range(q)})
    for secret_a in range(q):
        for secret_b in range(q):
            for a0 in range(q):
                a1 = (secret_a - a0) % q
                for b0 in range(q):
                    b1 = (secret_b - b0) % q
                    probes = [Counter() for _ in range(6)]
                    for r in range(q):
                        z0 = r
                        z1 = (-r) % q
                        probes[0][z0] += 1
                        probes[1][z1] += 1
                        z0 = (z0 + a0 * b0) % q
                        probes[2][z0] += 1
                        z1 = (z1 + a0 * b1) % q
                        probes[3][z1] += 1
                        z1 = (z1 + a1 * b0) % q
                        probes[4][z1] += 1
                        z1 = (z1 + a1 * b1) % q
                        probes[5][z1] += 1
                    assert all(p == uniform_q for p in probes)


def main():
    ring_tv = total_variation(
        legacy_ring_cross_probe(0, 0), legacy_ring_cross_probe(1, 0))
    mi32_tv = total_variation(
        legacy_mi32_probe(0), legacy_mi32_probe(1))
    mod_tv = total_variation(
        legacy_mmod_probe(0, 0), legacy_mmod_probe(1, 0))
    assert ring_tv != 0 and mi32_tv != 0 and mod_tv != 0

    check_repaired_ring()
    check_repaired_mi32()
    check_repaired_mmod()

    print("legacy ring_mul distinguisher: TV =", ring_tv)
    print("legacy mi32_to_mrx160 distinguisher: TV =", mi32_tv)
    print("legacy mmod32_mul distinguisher: TV =", mod_tv)
    print("repaired ring, mi32 A2B, and mod-q masked accumulators: uniform")


if __name__ == "__main__":
    main()
