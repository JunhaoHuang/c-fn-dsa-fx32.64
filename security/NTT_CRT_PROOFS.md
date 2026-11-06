# 掩码 NTT、iNTT、CRT 与中心化的 NI/SNI 论证

本文覆盖 `mask_mq.c`、`mask160.c` 中由它调用的 Montgomery 域转换、
`mi32_center_mod_m()` 和真实密钥展开钩子。安全模型和状态标记沿用
[MODEL.md](MODEL.md)：两个一阶份额、一个 word probe、源级顺序执行，以及每次
请求均返回新鲜独立均匀随机量的 **IDEAL-RNG**。

固定 LFSR 不满足最后一个前提，因此本文所有 `1-NI*`/`1-SNI*` 结论都是条件
结论；可执行程序继续标记为 `DEMO ONLY`。

## 1. 表示不变量

对公开素数 `p`（12289 或 18433），每个 `mmod32` 在接口边界满足：

```text
0 <= share[i] < p
(share[0] + share[1]) mod p = x*R mod p,  R=2^32.
```

`mq_add_share()`、`mq_sub_share()`、`mq_mul_public_share()` 和
`mq_decode_share()` 都只接收一个份额。它们从类型和调用图上都不可能读取另
一个份额。`mmod32_set_residue_shares()` 仍可供独立 gadget 测试使用；多项式和密
钥路径则显式执行完整的 share-0 遍历，再执行完整的 share-1 遍历，最后在单独的
边界遍历中加入 `rR,-rR`。编码和解码是模 p 的线性双射，因此保持重构值。单份额
线性映射满足 1-NI；带末尾 refresh 的导入在 IDEAL-RNG 下满足条件 1-SNI。

## 2. NTT 与 iNTT

每个前向蝶形对固定份额执行：

```text
y1 = x1
y2 = MontMul(x2, public_twiddle)
out0 = y1 + y2 mod p
out1 = y1 - y2 mod p.
```

逆蝶形执行同一个公开线性映射，其中第一输出乘公开 `1/2`，第二输出乘 Falcon
`iGM` 表中已经包含相应缩放的公开逆根。循环层数、地址和 twiddle 索引仅依赖
公开 `logn`。

### Signed lazy 范围与规范化

当前蝶形只对单个份额使用 signed REDC 和不带约简的加减；
内存槽仍是 `uint32_t`，负的中间值按模 `2^32` 编码，读取时显式还原为
`int32_t`。这不是跨份额中心化，也不改变公开接口的 canonical 约定。

- 前向变换开始时 `0 <= x < p`。第 k 层之后
  `|x| < (k+1)*p`，因为每个 twiddle 乘积的 signed REDC 在 `(-p,p)`。
  对全部支持的 `logn <= 10`，有 `|x| < 11*p <= 202763`；
  蝶形加减不会发生 signed-32 溢出。
- twiddle 至多为 p，所有前向乘法以及出口规范化的乘积满足
  `|t| < 11*p^2 <= 3737530379 < p*2^31`，满足 signed REDC 的前置条件。
  该乘积可能超过 `INT32_MAX`，因此必须在乘法之前提升到 `int64_t`。
- 前向出口用 `REDC(x*(R mod p))` 保持原 Montgomery 值并缩到 `(-p,p)`，
  再执行一次 cadd，恢复 canonical；这是每个系数一次额外 REDC，而非
  每个蝶形执行乘法结果及两次加减的 canonical 修正。
- 逆向变换在每层对和、差分别执行 signed Montgomery 乘法。
  初始输入 canonical，后续每个输出均在 `(-p,p)`；乘法前
  `|x1 +/- x2| < 2*p`、`|t| < 2*p^2 < p*2^31`。
  因此仅在最后对每个输出执行 cadd，无需额外 REDC。

每个份额的出口规范化完成后，才执行独立的边界 refresh。
公开 NTT 使用同一 lazy 范围及出口规则。公开常数点乘在两个 share-major
遍历中复用一次编码的常数，避免逐系数重新构建参数及重复 Barrett 约简，
最后在独立遍历中 refresh。

源码的最外层循环是 `for (sh = 0; sh < 2; sh++)`：share 0 的全部层和系数先完成，
然后才处理 share 1；刷新位于这两个遍历之后。因此整次变换是对每个输入份额独立
应用同一公开线性映射。任一内部探针可由一个输入份额和公开表模拟，故变换内部满
足 1-NI。变换返回前，每个系数使用独立模 p 随机量刷新；任一输出份额均匀且不需
要输入份额即可模拟，故输出边界满足条件 1-SNI。

这是多输出结论所必需的限制：每个系数必须消耗独立随机量。自伴位置不能复制同
一份额对；`mmpoly_gram_ntt()` 和 `mmpoly_cross_ntt()` 对每个输出位置独立计算并
刷新。

## 3. 点乘、Gram 与 cross

秘密点乘逐系数调用已有条件 1-SNI 的 `mmod32_mul()`。公开点乘是两个显式、独立
的 share-major Montgomery 乘法遍历，随后逐系数刷新。

Gram 和 cross 分别计算：

```text
f[i]*f[n-1-i] + g[i]*g[n-1-i]
F[i]*f[n-1-i] + G[i]*g[n-1-i].
```

把 Gram/cross “分开两个 share”不能解释为仅计算 `a0*b0` 与 `a1*b1`，因为
`(a0+a1)(b0+b1)` 还包含 `a0*b1+a1*b0`。实现中的每个 `mmod32_mul()` 明确计算
00、01、10、11 四项并在累加前注入新鲜随机量；只有乘法后的线性加法才逐份额。

共轭索引是公开的。两个秘密乘积各自使用独立随机量；线性求和之后再用第三个独立
随机量刷新最终输出。因此单输出组合满足条件 1-SNI。接口规定 `d` 不得与输入重
叠，防止写回改变后续公开索引读取，也避免给组合证明引入未建模的扇出。

## 4. 逐份额 CRT

令 `q=12289`、`q'=18433`、`M=q*q'=226523137`。对份额编号 i，先分别解码
`a_i mod q` 和 `b_i mod q'`，再计算：

```text
u_i = (b_i-a_i)*3 mod q'
c_i = a_i + q*u_i,       0 <= c_i < M.
```

实现以 share 为最外层循环：先生成所有 `c0`，再生成所有 `c1`，之后才进入中心
化遍历；CRT 映射阶段不存在跨份额读取。这里 3 是公开的 `q^(-1) mod q'`。由中国
剩余定理同构的加法保持性：

```text
CRT(a0,b0) + CRT(a1,b1)
  = CRT(a0+a1, b0+b1) mod M.
```

所以 `c0,c1` 是目标系数的模 M 算术分享。实现从未形成 `c0+c1`；Lean 定理
`crt_component_range`、`crt_share_sum_range` 和
`twice_crt_modulus_lt_pow2_29` 证明概念上的整数和满足：

```text
0 <= c0+c1 < 2*M < 2^29.
```

这也是后续 32 位符号检测有效的关键条件。

## 5. 安全中心化

`mi32_center_mod_m()` 首先用已证明的 word-sliced A2B 将 `(c0,c1)` 转成布尔
分享，然后在固定 Kogge-Stone 网络中执行两次公开常数减法：

```text
x = s - [s >= M]*M
y = x - [x > (M-1)/2]*M.
```

因为 `s<2^29`，`s-M` 和 `x-(H+1)` 都在有符号 32 位范围，最高位恰好等于比较
所需的借位/符号。选择器用 ISW AND 处理掩码控制位；控制位不打开，也不进入 C
分支。最后通过 Psi-B2A 返回算术分享并再次刷新。

组合使用的 A2B、前缀加法、ISW 选择和 Psi-B2A 均沿用
[GADGET_PROOFS.md](GADGET_PROOFS.md) 的条件证明。因此该中心化 gadget 在
IDEAL-RNG 和源级求值模型下为条件 1-SNI。Lean 定理 `centerM_range` 证明输出
严格处于：

```text
[-113261568, +113261568].
```

逐份额中心化不是环同态并会得到错误结果，生产路径明确禁止这种做法。

## 6. 真实密钥展开钩子

`fndsa_masked_key_expand_mq_internal()` 先展开种子份额，然后对每个多项式执行两
次完整遍历：share 0 只读取修正份额，share 1 只读取种子展开份额。每个遍历分别把
该份额编码到 q 和 q' 的 Montgomery 域，不重构系数。两个遍历全部完成后，才对每
个 q/q' 系数独立刷新。若种子展开在 IDEAL-RNG 模型中均匀，则修正份额也均匀；
因此内部编码为 1-NI，刷新后的输出边界为条件 1-SNI。固定 demo 种子不满足生产随
机性前提。

## 7. 打开边界与实现约束

该模块不调用 `mmod32_open()`、`mrx160_open_rx()` 或兼容调整位接口。测试代码
可以开值进行差分比较，生产模块不能。必须继续保持：

- `mask_mq.c` 和非线性 gadget 关闭 LTO；
- 模数、`logn`、循环边界和表索引为公开值；
- 不使用秘密 `/`、`%`、数据相关循环或秘密地址；
- q/q' 两套份额在 CRT 前都保持 canonical Montgomery 表示；
- 每次变换边界、每个 Gram/cross 输出和每次中心化使用新鲜随机量；
- 真实生产声明必须用密码学随机源替换固定 LFSR，并重新审查目标汇编。
