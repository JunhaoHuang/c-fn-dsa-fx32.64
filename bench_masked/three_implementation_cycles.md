# Three-way masked arithmetic cycle comparison

> Superseded for the repaired implementation by
> `three_implementation_reassessment_2026-10-31.md`. The older multiplier
> values below predate the final random-first first-order repair and must not
> be used as current performance numbers.

## Scope

This report compares three first-order, two-share implementations in one x86-64
process:

1. the local Q32.64 implementation with 160-bit arithmetic shares;
2. Chen--Chen masked IEEE-754 from
   `kengyuchen/masked-falcon-power-analysis`, commit
   `24be736ad5544f02f2c11f217de708740d988c8e`;
3. Burak's masked IEEE-754 from `burak231708/masked_FALCON`, commit
   `4c103fb774adfb41edc406d8d52283b9b74418c9`.

The benchmark checks that every implementation reconstructs `5 + 3 = 8` and
`5 * 3 = 15` before measuring it. All sources use Clang 18.1.3, `-O2`, and
`-fno-lto`. Runs are pinned to logical CPU 2 on an AMD Ryzen 9 9950X3D.

The counter is the invariant x86 TSC. Therefore the values below are
**TSC reference ticks, not PMU core cycles**. They are suitable for a same-host
relative comparison, but must not be mixed with Cortex-M4 cycle counts.

## Median of three complete runs

| Operation | Local fixed Q32.64 | Chen--Chen float | Burak float |
|---|---:|---:|---:|
| add | 7.80 | 6,579.00 | 25,174.35 |
| multiply, local 160-bit fast path | 1,227.65 | 3,311.00 | 15,690.70 |
| multiply, local 192-bit fallback | 3,812.67 | 3,311.00 | 15,690.70 |

All entries are TSC ticks per operation.

Derived from the medians:

- Chen--Chen add / local add: 843.46x.
- Burak add / local add: 3,227.48x.
- Burak add / Chen--Chen add: 3.83x.
- Chen--Chen multiply / local 160-bit multiply: 2.70x.
- Burak multiply / local 160-bit multiply: 12.78x.
- Burak multiply / Chen--Chen multiply: 4.74x.
- Chen--Chen multiply / local 192-bit fallback: 0.87x; equivalently, the local
  fallback costs 1.15x Chen--Chen's multiply.
- Burak multiply / local 192-bit fallback: 4.12x.

## Raw runs

| Run | fixed add | Chen add | Burak add | fixed mul 160 | fixed mul 192 | Chen mul | Burak mul |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7.80 | 6,604.80 | 25,167.90 | 1,237.97 | 3,812.67 | 3,311.00 | 15,656.30 |
| 2 | 7.82 | 6,579.00 | 25,195.85 | 1,227.65 | 3,885.77 | 3,334.65 | 15,690.70 |
| 3 | 7.80 | 6,566.10 | 25,174.35 | 1,212.60 | 3,743.87 | 3,308.85 | 15,733.70 |

## Interpretation and limitations

The fixed-point addition is a five-limb share-wise linear operation. Floating
addition must also align exponents, shift a secret mantissa, normalize, round,
and repack IEEE-754, so the gap is expected and does not represent equivalent
work.

For multiplication, the local proven-range 160-bit path is the fastest. The
192-bit exact fallback is slightly slower than Chen--Chen's floating multiply,
so reducing fallback frequency or optimizing its A2B/B2A and wide multiply is
the local priority. Chen--Chen is substantially faster than Burak for both
operations in this build.

Random generation is included in each measured gadget, but none of the three
random sources is a production-equivalent CSPRNG: the local code uses its demo
LFSR, Chen--Chen uses a repeating test table, and Burak uses `rand()`. Thus this
comparison measures the submitted/demo gadget implementations, not a complete
production masking stack with an entropy backend.

The [Chen--Chen TCHES 2024 paper](https://doi.org/10.46586/tches.v2024.i2.276-303)
separately reports actual Cortex-M4 results for its two-share implementation:
17,154 cycles for masked addition and 7,134 cycles for masked multiplication
(unmasked baselines: 487 and 308 cycles).
Those published M4 values confirm the relative ordering inside that work, but
are deliberately not combined numerically with the x86 TSC table above.

## Reproduction

The source is `bench_masked/compare_kengyu_fixed.c`. It links the two pinned
third-party checkouts into one process, renaming Chen--Chen's `rand64` symbol to
avoid colliding with Burak's symbol. The comparison uses 101 timing samples and
20 warm-up samples per operation; the per-run value is the sample median.
