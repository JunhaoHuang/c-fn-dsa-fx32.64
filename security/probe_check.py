#!/usr/bin/env python3
"""Exhaustive single-probe NI / SNI / PINI / MIMO classifier for the two-share
gadgets of this repository, evaluated on reduced word widths.

Model (see MODEL.md / DEFINITIONS.md): two shares, t = 1, one probe on one
intermediate value ("wire"), ideal fresh randomness.  Every wire that the C
code materialises is recorded by the gadget model below.  For each wire and
each possible secret assignment, the script computes the exact distribution
of the wire conditioned on a chosen set of input shares, and then decides:

  0-share   : the wire's distribution is identical for every secret value
              (simulatable without reading any input share);
  PINI      : there is one share index i such that conditioning on share i of
              *every* input makes the distribution secret-independent;
  NI        : there is one share index per input (possibly different per
              input) with the same property;
  FAIL      : no such choice exists -> first-order flaw.

A gadget is 1-NI if every wire is at least NI; 1-SNI if, in addition, every
output share wire is 0-share; 1-PINI if every wire is PINI and output share j
is simulatable from index-j input shares only.  For a two-output gadget the
MIMO-1-SNI condition additionally requires that the pair (out1 share i, out2
share j) is jointly 0-share for every i, j.

Run:  python3 -I security/probe_check.py
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product


class Trace:
    """Records every intermediate value of one gadget execution."""

    def __init__(self):
        self.wires = {}

    def w(self, name, value):
        self.wires[name] = value
        return value


def distributions(gadget, secrets_space, share_space, rand_space, n_inputs):
    """Return {wire: {(secrets, cond_key): Counter}} for each conditioning
    choice.  cond_key = ('none',) | ('pini', i) | ('ni', (i_1..i_n))."""
    conds = [("none",)]
    conds += [("pini", i) for i in range(2)]
    conds += [("ni", idx) for idx in product(range(2), repeat=n_inputs)]
    table = defaultdict(lambda: defaultdict(Counter))
    for sec in secrets_space:
        for sp in share_space:
            shares = share_space_to_shares(sec, sp)
            for rnd in rand_space:
                tr = Trace()
                gadget(tr, shares, rnd)
                for name, val in tr.wires.items():
                    for c in conds:
                        if c[0] == "none":
                            key = ()
                        elif c[0] == "pini":
                            key = tuple(shares[k][c[1]] for k in range(n_inputs))
                        else:
                            key = tuple(shares[k][c[1][k]] for k in range(n_inputs))
                        table[name][c][(sec, key, val)] += 1
    return table, conds


_share_builder = None


def share_space_to_shares(sec, sp):
    return _share_builder(sec, sp)


def secret_independent(counter_for_cond, n_inputs_unused):
    """Check that for every conditioning key the normalised distribution of
    the wire is the same across all secret values."""
    by_key = defaultdict(dict)
    totals = defaultdict(Counter)
    for k, cnt in counter_for_cond.items():
        if len(k) == 3:
            sec, key, val = k
            by_key[key].setdefault(sec, Counter())[val] += cnt
    for key, per_sec in by_key.items():
        ref = None
        for sec, c in per_sec.items():
            tot = sum(c.values())
            norm = {v: Fraction(n, tot) for v, n in c.items()}
            if ref is None:
                ref = norm
            elif norm != ref:
                return False
    return True


def classify(gadget, name, secrets_space, share_space, rand_space, n_inputs,
             outputs, pairs=()):
    table, conds = distributions(gadget, secrets_space, share_space,
                                 rand_space, n_inputs)
    worst = "0-share"
    order = {"0-share": 0, "PINI": 1, "NI": 2, "FAIL": 3}
    report = {}
    for wire in table:
        cls = "FAIL"
        if secret_independent(table[wire][("none",)], n_inputs):
            cls = "0-share"
        elif any(secret_independent(table[wire][("pini", i)], n_inputs)
                 for i in range(2)):
            cls = "PINI"
        elif any(secret_independent(table[wire][c], n_inputs)
                 for c in conds if c[0] == "ni"):
            cls = "NI"
        report[wire] = cls
        if order[cls] > order[worst]:
            worst = cls
    out_ok = all(report[o] == "0-share" for o in outputs)
    sni = worst != "FAIL" and out_ok
    pini = all(report[w] in ("0-share", "PINI") for w in report)
    mimo = True
    for a, b in pairs:
        mimo = mimo and report.get(a) == "0-share" and report.get(b) == "0-share"
    print(f"{name:34s} worst-wire={worst:7s} outputs-0-share={out_ok!s:5s} "
          f"1-NI={worst != 'FAIL'!s:5s} 1-SNI={sni!s:5s} 1-PINI={pini!s:5s}")
    bad = [w for w, c in report.items() if c == "FAIL"]
    if bad:
        print("    FAIL wires:", ", ".join(bad))
    ni_only = [w for w, c in report.items() if c == "NI"]
    if ni_only:
        print("    NI-only wires (not PINI):", ", ".join(ni_only))
    return report


# ---------------------------------------------------------------- helpers
def bool_shares(sec, sp):
    return [(s ^ r, r) for s, r in zip(sec, sp)]


def arith_shares_mod(m):
    def f(sec, sp):
        return [((s - r) % m, r) for s, r in zip(sec, sp)]
    return f


def run(name, gadget, n_inputs, secret_vals, share_vals, rand_vals, builder,
        outputs, pairs=()):
    global _share_builder
    _share_builder = builder
    secrets_space = list(product(secret_vals, repeat=n_inputs))
    share_space = list(product(share_vals, repeat=n_inputs))
    return classify(gadget, name, secrets_space, share_space, rand_vals,
                    n_inputs, outputs, pairs)


# ---------------------------------------------------------------- gadgets
# 1. fndsa_mask_bool_and_words, one 1-bit lane.
def g_and(tr, sh, rnd):
    (x0, x1), (y0, y1) = sh
    rx, ry, r = rnd
    ax0 = tr.w("ax0", x0 ^ rx); ax1 = tr.w("ax1", x1 ^ rx)
    ay0 = tr.w("ay0", y0 ^ ry); ay1 = tr.w("ay1", y1 ^ ry)
    a = tr.w("a=r", r); b = tr.w("b=r", r)
    t = tr.w("t00", ax0 & ay0); a = tr.w("a^t00", a ^ t)
    t = tr.w("t01", ax0 & ay1); b = tr.w("b^t01", b ^ t)
    t = tr.w("t10", ax1 & ay0); b = tr.w("b^t01^t10", b ^ t)
    t = tr.w("t11", ax1 & ay1); b = tr.w("z1", b ^ t)
    tr.w("z0", a)


# 1b. A plain ISW-style AND without input refresh, for comparison.
def g_and_norefresh(tr, sh, rnd):
    (x0, x1), (y0, y1) = sh
    r = rnd[0]
    tr.w("t00", x0 & y0); tr.w("t01", x0 & y1)
    tr.w("t10", x1 & y0); tr.w("t11", x1 & y1)
    tr.w("z0", r ^ (x0 & y0))
    tr.w("z1", r ^ (x0 & y1) ^ (x1 & y0) ^ (x1 & y1))


# 1c. AND applied to the same sharing twice (fan-out/alias), with and
# without the input refresh of fndsa_mask_bool_and_words.
def g_and_alias_norefresh(tr, sh, rnd):
    (x0, x1), = sh
    r = rnd[0]
    tr.w("t01=x0&x1", x0 & x1)
    tr.w("z0", r ^ (x0 & x0))
    tr.w("z1", r ^ (x0 & x1) ^ (x1 & x0) ^ (x1 & x1))


def g_and_alias_refreshed(tr, sh, rnd):
    (x0, x1), = sh
    rx, ry, r = rnd
    ax0 = tr.w("ax0", x0 ^ rx); ax1 = tr.w("ax1", x1 ^ rx)
    ay0 = tr.w("ay0", x0 ^ ry); ay1 = tr.w("ay1", x1 ^ ry)
    tr.w("t01", ax0 & ay1); tr.w("t10", ax1 & ay0)
    tr.w("z0", r ^ (ax0 & ay0))
    tr.w("z1", r ^ (ax0 & ay1) ^ (ax1 & ay0) ^ (ax1 & ay1))


# 2. Boolean refresh.
def g_refresh(tr, sh, rnd):
    (x0, x1), = sh
    r = rnd[0]
    tr.w("z0", x0 ^ r); tr.w("z1", x1 ^ r)


# 3. psi_b2a_words (BCZ18 C_2), k-bit words modulo 2^k.
def make_b2a(k):
    m = 1 << k

    def g(tr, sh, rnd):
        (x0, x1), = sh
        r1, r2 = rnd
        y0 = tr.w("y0", x0 ^ r1 ^ r2); y1 = tr.w("y1", x1 ^ r1); y2 = tr.w("y2", r2)
        xi = tr.w("psi1.xor", y0 ^ y1)
        p = tr.w("p", (xi - y1) % m)
        z0 = tr.w("z0", y0 ^ p)
        xi2 = tr.w("psi2.xor", y0 ^ y2)
        z1 = tr.w("z1", (xi2 - y2) % m)
        tr.w("d0", y1 ^ y2); tr.w("d1", (z0 ^ z1) % m)
    return g


# 4. ring_mul random-first accumulator, modulo 2^k (same structure as mmod32).
def make_ring_mul(m):
    def g(tr, sh, rnd):
        (a0, a1), (b0, b1) = sh
        r = rnd[0]
        t0 = tr.w("t0=r", r); t1 = tr.w("t1=-r", (-r) % m)
        p = tr.w("p00", a0 * b0 % m); t0 = tr.w("z0", (t0 + p) % m)
        p = tr.w("p01", a0 * b1 % m); t1 = tr.w("t1+p01", (t1 + p) % m)
        p = tr.w("p10", a1 * b0 % m); t1 = tr.w("t1+p01+p10", (t1 + p) % m)
        p = tr.w("p11", a1 * b1 % m); tr.w("z1", (t1 + p) % m)
    return g


# 4b. ring_mul(a, a): the aliasing hazard (single logical input).
def make_ring_sqr_alias(m):
    def g(tr, sh, rnd):
        (a0, a1), = sh
        r = rnd[0]
        tr.w("p01=a0*a1", a0 * a1 % m)
        tr.w("z0", (r + a0 * a0) % m)
        tr.w("z1", (-r + 2 * a0 * a1 + a1 * a1) % m)
    return g


# 4c. mrx160_sqr: two independent refreshes, then ring_mul.
def make_sqr_refreshed(m):
    def g(tr, sh, rnd):
        (a0, a1), = sh
        r, s, u = rnd
        x0 = tr.w("x0", (a0 + r) % m); x1 = tr.w("x1", (a1 - r) % m)
        y0 = tr.w("y0", (a0 + s) % m); y1 = tr.w("y1", (a1 - s) % m)
        tr.w("p01=x0*y1", x0 * y1 % m); tr.w("p10=x1*y0", x1 * y0 % m)
        tr.w("z0", (u + x0 * y0) % m)
        tr.w("z1", (-u + x0 * y1 + x1 * y0 + x1 * y1) % m)
    return g


# 5. mbit_to_mmod32 with q = 5 (Boolean bit -> arithmetic mod q).
def make_bit_to_mod(q):
    def g(tr, sh, rnd):
        (a0, a1), = sh
        r, rb0, ar1, rr = rnd
        rb1 = tr.w("rb1", rb0 ^ r)
        ar0 = tr.w("ar0", (r - ar1) % q)
        c0 = tr.w("a0^rb0", (a0 & 1) ^ rb0)
        c1 = tr.w("a1^rb1", (a1 & 1) ^ rb1)
        tr.w("COMPILER:a0^a1", (a0 ^ a1) & 1)   # what clang -O2 materialises
        c = tr.w("c", c0 ^ c1)
        nr0 = tr.w("nr0", (1 - ar0) % q); nr1 = tr.w("nr1", (-ar1) % q)
        s0 = tr.w("sel0", nr0 if c else ar0); s1 = tr.w("sel1", nr1 if c else ar1)
        tr.w("d0", (s0 + rr) % q); tr.w("d1", (s1 - rr) % q)
    return g


# 6. mb64_compare_byte on a 1-bit "byte", public w, two outputs (eq, lt).
# AND gadgets are modelled with fresh randomness per invocation; only their
# outputs are recorded as wires here (their internals are checked in #1).
def make_compare(w_pub, nbits):
    def AND(tr, tag, x, y, rnd):
        rx, ry, r = rnd
        ax0, ax1 = x[0] ^ rx, x[1] ^ rx
        ay0, ay1 = y[0] ^ ry, y[1] ^ ry
        z0 = r ^ (ax0 & ay0)
        z1 = r ^ (ax0 & ay1) ^ (ax1 & ay0) ^ (ax1 & ay1)
        tr.w(tag + ".z0", z0); tr.w(tag + ".z1", z1)
        return (z0, z1)

    def g(tr, sh, rnd):
        (t0, t1), = sh
        it = iter(rnd)
        e = (next(it), next(it) ^ 1)
        l = (next(it), next(it))
        for i in range(nbits - 1, -1, -1):
            b = ((t0 >> i) & 1, (t1 >> i) & 1)
            wi = (w_pub >> i) & 1
            diff = (b[0] ^ wi, b[1])
            ndiff = (diff[0] ^ 1, diff[1])
            take = AND(tr, f"take{i}", e, diff, (next(it), next(it), next(it)))
            if wi == 0:
                tt = AND(tr, f"or{i}.and", l, take, (next(it), next(it), next(it)))
                o = (l[0] ^ take[0] ^ tt[0], l[1] ^ take[1] ^ tt[1])
                rr = next(it)
                l = (tr.w(f"l{i}.0", o[0] ^ rr), tr.w(f"l{i}.1", o[1] ^ rr))
            ne = AND(tr, f"ne{i}", e, ndiff, (next(it), next(it), next(it)))
            e = ne
        r_eq = next(it); r_lt = next(it)
        eq = (tr.w("eq0", e[0] ^ r_eq), tr.w("eq1", e[1] ^ r_eq))
        lt = (tr.w("lt0", l[0] ^ r_lt), tr.w("lt1", l[1] ^ r_lt))
        # MIMO pairs: one share of each output observed jointly.
        for i in range(2):
            for j in range(2):
                tr.w(f"pair(eq{i},lt{j})", (eq[i], lt[j]))
    return g


# 7. mb160_select (one lane): d = a ^ AND(a ^ b, ctl).
def g_select(tr, sh, rnd):
    (a0, a1), (b0, b1), (c0, c1) = sh
    rx, ry, r = rnd
    x0 = tr.w("x0", a0 ^ b0); x1 = tr.w("x1", a1 ^ b1)
    ax0, ax1 = x0 ^ rx, x1 ^ rx
    ay0, ay1 = c0 ^ ry, c1 ^ ry
    t0 = tr.w("t0", r ^ (ax0 & ay0))
    t1 = tr.w("t1", r ^ (ax0 & ay1) ^ (ax1 & ay0) ^ (ax1 & ay1))
    tr.w("d0", a0 ^ t0); tr.w("d1", a1 ^ t1)


def main():
    bits = [0, 1]
    w2 = range(4)
    print("== Boolean gadgets (1-bit lanes) ==")
    run("bool_and_words", g_and, 2, bits, bits, list(product(bits, repeat=3)),
        bool_shares, ["z0", "z1"])
    run("ISW-AND no refresh, indep. inputs", g_and_norefresh, 2, bits, bits,
        [(0,), (1,)], bool_shares, ["z0", "z1"])
    run("ISW-AND no refresh, AND(x,x) alias", g_and_alias_norefresh, 1, bits,
        bits, [(0,), (1,)], bool_shares, ["z0", "z1"])
    run("bool_and_words AND(x,x) alias", g_and_alias_refreshed, 1, bits,
        bits, list(product(bits, repeat=3)), bool_shares, ["z0", "z1"])
    run("bool_refresh_words", g_refresh, 1, bits, bits, [(0,), (1,)],
        bool_shares, ["z0", "z1"])
    run("mb160_select lane", g_select, 3, bits, bits,
        list(product(bits, repeat=3)), bool_shares, ["d0", "d1"])
    print("== Conversions ==")
    run("psi_b2a_words (2-bit)", make_b2a(2), 1, w2, w2,
        list(product(w2, repeat=2)), bool_shares, ["d0", "d1"])
    run("mbit_to_mmod32 (q=5)", make_bit_to_mod(5), 1, bits, bits,
        list(product(bits, bits, range(5), range(5))), bool_shares,
        ["d0", "d1"])
    print("== Multiplications (Z/4 stands for Z/2^160 and Z/q) ==")
    run("ring_mul / mmod32_mul (mod 4)", make_ring_mul(4), 2, w2, w2,
        [(r,) for r in w2], arith_shares_mod(4), ["z0", "z1"])
    run("ring_mul(a,a) alias (mod 4)", make_ring_sqr_alias(4), 1, w2, w2,
        [(r,) for r in w2], arith_shares_mod(4), ["z0", "z1"])
    run("mrx160_sqr refreshed (mod 4)", make_sqr_refreshed(4), 1, w2, w2,
        list(product(w2, repeat=3)), arith_shares_mod(4), ["z0", "z1"])
    print("== Two-output comparator (1-bit byte, MIMO pairs) ==")
    for w in range(2):
        rnd_len = 4 + (3 + 3 + 1 + 3) + 2
        rnd = list(product(bits, repeat=rnd_len))
        rep = run(f"mb64_compare_byte w={w}", make_compare(w, 1), 1, bits, bits,
                  rnd, bool_shares, ["eq0", "eq1", "lt0", "lt1"],
                  pairs=[("eq0", "lt0"), ("eq0", "lt1"),
                         ("eq1", "lt0"), ("eq1", "lt1")])
        mimo = all(rep[f"pair(eq{i},lt{j})"] == "0-share"
                   for i in range(2) for j in range(2))
        print(f"    MIMO-1-SNI (joint eq/lt share pairs 0-share): {mimo}")


if __name__ == "__main__":
    main()
