# masked_FALCON A2B/B2A comparison

Reference repository: <https://github.com/burak231708/masked_FALCON>, commit
`4c103fb774adfb41edc406d8d52283b9b74418c9` (2025-01-02). The repository
contains no paper citation, release tag, or license. Consequently, this project
does not copy its source code; it independently re-expresses only the
mathematical two-share Psi construction.

## Algorithm assessment

The reference A2B recursively combines arithmetic shares with `SecAdd`. Its
addition is the same broad prefix-adder approach already used by this project's
word-sliced Kogge-Stone A2B, so it offers no new asymptotic advantage. Moreover,
UBSan stops the reference test at `src/gadgets.c:161`: its 64-bit `SecAdd`
shifts a `uint64_t` by 64 bits. The 128-bit variant also reaches invalid shift
counts. This A2B was therefore not imported.
The reference `test_gadgets.c` also calls `B2A(m5, ...)` but then unmasks
`m4`, so it does not validate the produced B2A output. The local adaptation is
therefore accepted only on this project's independent differential tests.

The useful reference idea is the first-order, two-share B2A construction based
on `Psi(x,y) = (x XOR y) - y`. With two fresh masks r1 and r2, it converts a
Boolean sharing without a masked prefix addition. The implementation in
`mask160.c` performs exact multi-limb subtraction modulo 2^160 or 2^192, uses
fresh random words, compiler value barriers, noinline boundaries, and explicit
zeroization of all temporaries.

For 160 bits, the previous B2A consumed 90 fresh 32-bit words: five for the
arithmetic output mask plus 85 for the Kogge-Stone secure-AND network. The Psi
version consumes ten fresh 32-bit words. For 192 bits the corresponding change
is from 108 to 12 words.

## Performance

Host: AMD Ryzen 9 9950X3D, Ubuntu clang 18.1.3, `-O2`, mask module
`-fno-lto`. Values are medians after 20 warm-up samples and are invariant-TSC
ticks, not hardware core cycles.

The primary comparison alternated the exact previous production B2A and the
new production B2A in one process, with 101 samples and 2,000 conversions per
sample. This controls for the substantial CPU-frequency drift seen between
independent runs.

| Implementation | TSC ticks/op | Relative |
|---|---:|---:|
| Previous prefix-adder B2A | 713.13 | 1.000x |
| Production Psi-B2A, including barriers and zeroization | 138.46 | 5.150x faster |

A2B is unchanged. Complete independent-run logs are preserved as
[`mask_gadget_cycles.txt`](mask_gadget_cycles.txt) and
[`mask_gadget_tsc_after_psi_b2a.txt`](mask_gadget_tsc_after_psi_b2a.txt).
Because their absolute TSC rates differ, downstream results are normalized to
the unchanged A2B in each run: fast multiplication improves from 2.123 to
1.385 A2B-equivalents (1.533x), and wide multiplication from 6.922 to 4.380
A2B-equivalents (1.580x).

A stripped exploratory Psi-B2A without production zeroization and barriers
measured 48.98 ticks versus 756.95 ticks for the previous function (15.46x),
but that number is not used as the production claim.

## Validation

- 10,000 randomized 96-bit sign-extended round trips in the isolated prototype.
- Existing 2,000-case full 160-bit A2B/B2A carry-chain differential test.
- Existing Q32.64 edge, fast multiply, 192-bit fallback, signing compatibility,
  and FN-DSA-512/1024 verification tests.
- ASan and UBSan clean on the complete `test_masked` suite.
- Generated x86-64 assembly check: the multi-limb borrow uses `seta`; branches
  depend only on the public limb count, loop counters, and RNG-reservoir state.

This remains a DEMO ONLY implementation. The deterministic LFSR and observable
signing control flow still prevent a production side-channel security claim;
generated assembly and physical leakage require separate review.
