# 低层 gadget 的 NI/SNI 证明与反例

除非明确标为失败，本文件的正面结论都以 [MODEL.md](MODEL.md) 的
IDEAL-RNG 和 32 位 word-probing 模型为前提。

## 1. Boolean refresh：`fndsa_mask_bool_refresh_words`

对每个字，代码计算

```text
z0 = x0 XOR r
z1 = x1 XOR r
```

其中 `r` 是新鲜均匀 32 位字。

**正确性。** `z0 XOR z1 = x0 XOR x1`。

**1-SNI 模拟。**

- 一个输出探针：`z0` 或 `z1` 都是均匀字，模拟器不读取输入份额；
- 一个内部探针：`r` 可直接均匀采样；`x0 XOR r` 或 `x1 XOR r` 也可直接
  均匀采样。即使把 XOR 的一个操作数单独作为机器临时量探测，它至多涉及
  一个输入份额。

所以该 gadget 是条件 1-SNI。

## 2. Arithmetic refresh：`mrx160_refresh` 与 `mmod32_refresh`

`mrx160_refresh` 计算

```text
z0 = x0 + r mod 2^160
z1 = x1 - r mod 2^160
```

`mmod32_refresh` 在 `Z/qZ` 中计算相同结构。

**正确性。** 两个输出份额之和不变。

**1-SNI 模拟。** 任一输出份额在相应环上均匀；输出探针不需要输入。
多字加减中的单个进位最多依赖一个输入份额和随机数，可由一个输入份额模拟。
模 `q` 版本要求 `fndsa_mask_rng_mod()` 真正均匀采样 `[0,q)`。

## 3. 两份额 Boolean AND：`fndsa_mask_bool_and_words`

每个 32 位字首先用独立随机字 `rx,ry` 重分享两个输入，再用独立 `r` 初始化
输出：

```text
ax0 = x0 XOR rx       ax1 = x1 XOR rx
ay0 = y0 XOR ry       ay1 = y1 XOR ry
z0 = r                 z1 = r
z0 ^= ax0 AND ay0
z1 ^= ax0 AND ay1
z1 ^= ax1 AND ay0
z1 ^= ax1 AND ay1
```

**正确性。** 输入 refresh 不改变 `x0 XOR x1` 或 `y0 XOR y1`；展开四项可得
`z0 XOR z1 = (x0 XOR x1) AND (y0 XOR y1)`。

**1-SNI 模拟。** `rx,ry,r` 均可直接采样；任一 `ax_i` 或 `ay_i` 是均匀字。
即使两个逻辑输入是同一 sharing、移位视图或相关秘密，`rx` 与 `ry` 仍使送入
交叉项的两组 sharing 具有独立掩码。单个原始乘积是两个独立均匀字的公开函数，
其分布与秘密无关；`z0,z1` 及每个部分和从首项之前就含 `r`，单独均匀。因此
一个内部或输出 word probe 可在不读取秘密，或至多读取允许的一个输入份额下
模拟。

实现为 `noinline`；输入 refresh、输出初始化、每个乘积和每次累加均设值屏障，
且没有一个 RHS 同时包含两个交叉项。结论仍限理想随机数和源码 word-probing
模型，不是任意编译器或转换泄漏模型的机器码定理。
## 4. 前缀加法：`bool_add_words`

Kogge-Stone 网络的 XOR、公开移位和复制是逐份额线性操作；每个 generate/
propagate 非线性节点使用上一节的 1-SNI AND。网络拓扑只依赖公开宽度。
按 SNI gadget 的顺序/并行组合规则，内部单探针至多需要每个输入的一个份额，
因此该函数是 1-NI。

该函数本身不声明 SNI。最低结果位没有进位：

```text
z0[bit 0] = x0[bit 0] XOR y0[bit 0]
```

当输入份额被视为 gadget 的固定输入时，单独观察这个输出需要输入份额，违反
SNI 对“仅输出探针使用零个输入份额”的要求。调用者若需要 SNI，必须像
`arith_to_bool_words()` 一样在输出后执行新鲜 refresh。

## 5. A2B：`arith_to_bool_words`、`mrx160_a2b`、`a2b192`

对算术输入 `(a0,a1)`，实现先生成

```text
X = (a0 XOR rx, rx)
Y = (a1 XOR ry, ry)
```

再用 Boolean masked adder 求 `X+Y`，最后 Boolean refresh。

**正确性。** `X`、`Y` 分别重构为 `a0`、`a1`，因此输出重构为
`a0+a1 mod 2^k`。

**1-SNI 模拟。**

1. 重分享阶段的任一值要么均匀，要么只涉及一个算术输入份额；
2. 前缀加法的 1-NI 模拟器至多索取 `X`、`Y` 的一个份额，而每个这样的份额
   本身由独立 `rx`/`ry` 一次一密，因此可不读取原始算术输入模拟；
3. 最终 refresh 使任一输出份额均匀。
因此，任何 A2B 调用点都不得把算术份额直接编码为 `(a_i,0)`。当前
`mi32_to_mrx160()` 已改为调用同一安全重分享路径，见第 10 节。

## 6. B2A：`psi_b2a_words`、`mrx160_b2a`、`b2a192`

实现是 BCZ18 Algorithm 1 的 `n=2` 递归 `C2` 特化：先用两个独立随机字
把 Boolean 输入扩为三个分享，再在两个 `C1` 叶上计算
`Psi(x,y) = (x XOR y) - y`，最后组合两个算术输出份额。

BCZ18 给出 ISW probing 安全证明；FALCON floor 论文的 gadget 表把该 B2A
列为 `t-SNI`。特化到两个份额得到条件 1-SNI：一个内部探针至多由一个输入
份额模拟；单个输出份额由 `r1,r2` 随机化，可由零输入份额模拟。

多字 `Psi` 的进位链必须作为 32 位 word-probe 程序审查。当前代码的每个输出
字都有值屏障，但这仍不是编译器语义的机械化证明。

## 7. 线性算术：`copy/add/sub/neg`

这些函数只把 share `j` 与同编号 share `j` 组合。任一内部或输出探针可由
每个输入的至多一个份额直接模拟，所以是 1-NI。

它们通常不是独立 1-SNI gadget，因为输出 share `j` 是输入 share `j` 的确定
函数，没有新鲜输出随机化。在线性层之后增加正确的 arithmetic refresh 可获得
SNI 输出隔离，但 refresh 不能修复此前不满足 NI 的非线性层。

## 8. `floor`、`half` 与 `mrx160_to_mi32`

`mrx160_floor()` 的结构是：

```text
1-SNI A2B -> 公开固定投影/符号扩展 -> 1-SNI B2A -> 1-SNI refresh
```

投影对每个 Boolean 份额独立执行。A2B 的单份额是均匀随机字；保留的字、
固定零字以及由一个随机符号位复制出的全 0/全 1 字，其分布都与秘密无关。
因此内部单探针可模拟，最终输出由 B2A/refresh 隔离。`mrx160_floor()` 是条件
1-SNI。

`mrx160_half()` 使用相同结构，只把投影替换为公开的一位算术右移，也是条件
1-SNI。`mrx160_to_mi32()` 先做 SNI A2B，再对公开整数 word 做 SNI B2A，
所以也是条件 1-SNI。

## 9. 修复后的 `ring_mul` 与 `ring_mul192`

旧版本在随机化前生成 `p01+p10`，其最低位存在一阶区分器。当前实现已经改为：

```text
t0 = r;  t1 = -r;
t0 += p00;
t1 += p01; t1 += p10; t1 += p11;
```

随机 accumulator 在首个交叉项之前建立，每个 `pij` 单独只依赖每个输入的一个
份额，每个 accumulator 部分和都是均匀环元素。因而两个环乘在 IDEAL-RNG
word-probing 模型下为 1-SNI；完整模拟器、进位链处理和输出证明见
[REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)。

`mrx160_mul()`、`mrx160_mul_wide()` 与 `mrx160_mul_u64_hi()` 由该 SNI 环乘、
安全 A2B/B2A 和公开截断组合，也获得条件 1-SNI。该双输入结论不允许把同一
sharing 直接作为快速乘法的两个逻辑输入；平方必须调用 `mrx160_sqr()`。后者
先产生两个独立 SNI refresh，再调用双输入乘法，因此获得条件 1-SNI。完整的
别名反例和组合证明见 [REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)。

### 9.1 公开常数乘法：`mrx160_mul_public_fast/wide`

公开常数乘法不调用 `ring_mul`/`ring_mul192`。对公开常数 `c`，两个算术份额
分别计算 `a0*c` 与 `a1*c`；该层是 1-NI，因为任一中间字只依赖一个输入份额
和公开值。之后只在 share 0 加半 ULP，并依次执行已证明的 SNI A2B、公开固定
右移/符号扩展、SNI B2A 和最终 arithmetic refresh。

因此在 IDEAL-RNG、源码级 32 位 word-probing 模型以及 fast 路径精确乘积落在
signed-160 范围的前提下，`mrx160_mul_public_fast()` 为条件 1-SNI。
`mrx160_mul_public_wide()` 先对唯一秘密输入做一次 SNI A2B、Boolean 符号扩展和
SNI B2A192，再逐份额乘公开常数，最后采用同一安全截断链，亦为条件 1-SNI。
两者支持 `d == a`，因为写输出前已完成所有输入读取。

固定丢弃的 5/288 个随机字只保持与旧通用 290/909 字预算相同；这些字不得复用，
也不参与安全模拟。Lean 的 `reconstruct_mulPublic`、
`roundShift64_add_pow2_160` 与逐调用点范围定理给出功能/宽度证明。该局部结论
不能自动把 FFT、split/merge 或 LDL 的多输出图提升为 MIMO-SNI；这些组件仍为
`REVIEW`。

## 10. 修复后的 `mi32_to_mrx160`

旧 `(a0,0)+(a1,0)` 路径已删除。当前函数直接调用
`arith_to_bool_words(...,&a.share[0],&a.share[1],1,...)`，先独立重分享两个
算术份额，再执行 masked add 和输出 refresh；随后逐份额 Boolean 符号扩展、
BCZ18 B2A 和 arithmetic refresh。该组合为条件 1-SNI，证明见
[REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)。

## 11. 修复后的 Montgomery `mmod32_mul`

旧版本把三个交叉积同时放入一个原始 64 位 `z1`。当前 `mmod32` 使用
`R=2^32` 的 canonical Montgomery 份额，并从 `z0=rR`、`z1=-rR mod q`
开始。每个 `pij` 单独执行 REDC，再单独加入已掩码 accumulator。REDC 内部
只依赖当前两个输入份额；accumulator 是均匀 Montgomery residue 的平移。
末尾独立 refresh，因此条件 1-SNI 结论保持成立。完整范围和模拟证明见
[MONTGOMERY_PROOFS.md](MONTGOMERY_PROOFS.md)。
旧结构的三个区分器仍保留在 [counterexamples.py](counterexamples.py) 中作为
负面回归；同一脚本同时穷举验证修复后的一位环、1-bit A2B 和 `q=5` masked
accumulator 均匀性。

## 12. sampler Boolean 非线性层与 MIMO 输出

`mask_and_word()` 和 `mbit_and()` 仅是共享
`fndsa_mask_bool_and_words()` 的窄包装；OR、select、nonzero、64 位秘密移位和
byte comparator 不再保留第二套交叉项公式。OR/select/nonzero 和每个 64 位
select 返回前显式刷新；comparator 的 `eq`、`lt` 两个输出分别使用独立随机字
刷新。

同一状态进入多个非线性节点的扇出，由每次共享 AND 内部的两组独立输入 refresh
切断。byte comparator 的 `e` 与 `diff` 扇出、六级秘密移位的未移位/移位状态
以及所有输出 refresh 的组合图和条件 MIMO-1-SNI 模拟说明见
[SAMPLER_BOOLEAN_MIMO.md](SAMPLER_BOOLEAN_MIMO.md)。

该结论到共享 comparator 输出为止。BerExp 随后显式打开逐字节 `eq`，SamplerZ
打开接受位，因此二者仍分别是 `REVIEW + TRACE-OPEN` 与 `TRACE-OPEN`，不能因
低层修复而提升为端到端 SNI。
## 13. 参考资料

- Barthe et al., *Strong Non-Interference and Type-Directed Higher-Order Masking*,
  [IACR ePrint 2015/506](https://eprint.iacr.org/2015/506)；
- Ishai, Sahai, Wagner, *Private Circuits: Securing Hardware against Probing Attacks*,
  [CRYPTO 2003](https://iacr.org/archive/crypto2003/27290462/27290462.pdf)；
- Bettale, Coron, Zeitoun, *Improved High-Order Conversion From Boolean to Arithmetic
  Masking*, [TCHES 2018(2)](https://tches.iacr.org/index.php/TCHES/article/view/873)；
- Schneider et al., *Efficiently Masking Binomial Sampling at Arbitrary Orders*,
  [PKC 2019 preprint](https://www.h2020prometheus.eu/sites/default/files/2022-06/Schneider2019_Chapter_EfficientlyMaskingBinomialSamp.pdf)；
- Berthet et al., *Masked Computation of the Floor Function and Its Application to the
  FALCON Signature*, [IACR CiC](https://cic.iacr.org/p/1/4/9/pdf)。该文 Table 1
  明确列出 SecAnd、SecAdd、A2B、B2A 与 refresh 的 NI/SNI 前提，并用
  MIMO-SNI 组合图证明其浮点 floor；这些结论只有在本实现确实实现相同安全
  子 gadget 与组合条件时才能移植。
