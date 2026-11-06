# Sampler Boolean gadgets: first-order MIMO-SNI composition map

## Scope and assumptions

This note covers the Boolean nonlinear path shared by `mask160.c` and
`mask_sampler.c`. Claims marked `1-SNI*` or `MIMO-1-SNI*` are conditional on
the source-level 32-bit word-probing model from [MODEL.md](MODEL.md), fresh
independent uniform masks, and the compiler preserving the `noinline` and value
barriers. The bundled deterministic LFSR does not meet the randomness premise.
This is not a glitch, transition-leakage, fault, or machine-code proof.

## Single nonlinear boundary

`fndsa_mask_bool_and_words()` is the only Boolean AND implementation used by
the fixed-point and sampler gadgets. For every 32-bit slice it consumes three
independent random words `rx`, `ry`, and `r`:

```text
ax0 = x0 XOR rx       ax1 = x1 XOR rx       [Rin-x]
ay0 = y0 XOR ry       ay1 = y1 XOR ry       [Rin-y]
z0  = r                                                   [Rout]
z1  = r                                                   [Rout]
z0 ^= ax0 AND ay0
z1 ^= ax0 AND ay1
z1 ^= ax1 AND ay0
z1 ^= ax1 AND ay1
```

`Rin-x` and `Rin-y` are independent even when the caller supplies the same
sharing twice, a shifted view of one sharing, or two correlated logical values.
They therefore cut alias and fan-out edges before a cross-product is formed.
`Rout` exists before the first product is accumulated. Each product is placed
in a temporary, value-barriered, and added separately; no C RHS combines two
cross-products. Reconstruction gives `(x0 XOR x1) AND (y0 XOR y1)`.

`fndsa_mask_bool_refresh_words()` is the common Boolean output refresh. Each
word and each invocation consume an independent random word.

## Composition notation

The graphs below use:

- `F`: share-level fan-out;
- `L`: sharewise public linear operation;
- `A[Rin,Rout]`: the common AND, including two input refreshes and its fresh
  output mask;
- `R`: an explicit independent output refresh;
- `O`: an intentional opening boundary.

### OR, bit select, and nonzero reduction

```text
OR:       a --F--> L(a XOR b) ---------------------+
          b --F--> A[Rin,Rout](a,b) ---------------+--> L --> R --> out

SELECT:   a --F--> L(a XOR b) --+
          b ---------------------+--> A[Rin,Rout](a XOR b, ctl) --+
          a --F----------------------------------------------------+--> L --> R

NONZERO32:
          x --F--> x, L(x >> 16) --> A[Rin,Rout] --> L(OR)
                    ... >> 8,4,2,1; every stage repeats A ... --> R --> bit
```

The input refreshes inside every `A` cut the repeated/shifted-sharing fan-out.
`mbit_or()`, `mbit_select()`, and `mbit_nonzero32()` explicitly refresh their
returned bit, so their output share can be simulated independently of internal
probes under the first-order budget.

### 64-bit select and secret right shift

```text
MB64_SELECT, independently for each 32-bit word:
  a --F--> L(a XOR b) --+
  b ---------------------+--> A[Rin,Rout](a XOR b, ctl-mask) --+
  a --F---------------------------------------------------------+--> L --> Rword

SECRET_RSHIFT:
  x0 --> [L shift 1, MB64_SELECT, R] --> x1
     --> [L shift 2, MB64_SELECT, R] --> x2
     --> [L shift 4, MB64_SELECT, R] --> x3
     --> [L shift 8, MB64_SELECT, R] --> x4
     --> [L shift 16, MB64_SELECT, R] --> x5
     --> [L shift 32, MB64_SELECT, R] --> out
```

The two output words of every select receive independent AND masks and an
independent explicit word refresh. Thus a word probe on one output does not
fix the mask of the other output. Every stage cuts the fan-out between the
unshifted state and its publicly shifted view before the next stage.

### Byte comparator with two outputs

For each public bit `w_i`, the comparator maintains shared states `e` (all
higher bits equal) and `l` (threshold is lower):

```text
                    +--> A[Rin,Rout](e, diff) --> take --+
e --F---------------+                                    +--> OR --> R --> l'
                    +--> A[Rin,Rout](e, NOT diff) -------------> e'

diff --F------------+-------------------------------------------+

final e --> R_eq --> eq
final l --> R_lt --> lt
```

The fan-out of `e` and `diff` enters distinct AND invocations; each invocation
independently refreshes both inputs. `R_eq` and `R_lt` use distinct random
words. Consequently the pair `(eq,lt)` has an explicit first-order MIMO output
isolation boundary instead of inheriting a claim from one AND. The branch on
`w_i` is public.

## Conditional result and remaining boundary

Under the stated ideal-randomness and source-level assumptions:

- the common AND and refresh are `1-SNI*`;
- bit OR/select/nonzero and 64-bit select/secret shift are conditionally
  `1-SNI*` after the marked output refreshes;
- `mb64_compare_byte()` is conditionally `MIMO-1-SNI*` for its two shared
  outputs because `R_eq` and `R_lt` are independent.

This result ends at the comparator outputs. `mrx160_ber_exp()` deliberately
opens `eq` after each byte (`O_eq`) to preserve lazy comparison, and
`mrx160_sampler_z()` opens the acceptance bit (`O_accept`). Prefix length,
accept/reject, and retry count therefore remain explicit trace leakage; neither
caller is promoted to end-to-end SNI or trace-oblivious security.

## Current compiler audit

With Ubuntu Clang 18.1.3, `-O2 -fno-lto`, `mask_sampler.o` retains undefined
references to both common functions (`nm -u mask_sampler.o`); it does not carry
a private inlined AND implementation. `objdump -d
--disassemble=fndsa_mask_bool_and_words mask160.o` shows three RNG draws before
the input-dependent `and` sequence and four separate `and`/accumulate steps.
Branches in that symbol depend on the public word count or RNG-reservoir
position, not on `x` or `y`. This observation is compiler/version specific.

## Verification obligations

`test_mask_bool_nonlinear()` checks independent inputs, the same sharing used
as both operands, shifted correlated fan-out, in-place aliasing, refresh
correctness, and actual share changes. The full `test_masked` suite checks
BerExp byte consumption, SamplerZ decisions, dynamic sampling, and FN-DSA-512
and FN-DSA-1024 signing compatibility. These are functional regressions, not
leakage proofs. Production review still requires assembly inspection and a
cryptographic randomness backend.
