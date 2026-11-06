# 掩码安全性定义与组合规则汇编（NI / SNI / MIMO-SNI / PINI）

本文件汇总本仓库审计所依赖的全部形式化定义、组合定理和审计方法，
作为 [REVERIFY_2026-11-05.md](REVERIFY_2026-11-05.md) 的参考基准。
所有定义都以本仓库的参数实例化：两个份额、`t = 1`、一个探针读取一个
32 位 C 中间字；随机数前提沿用 [MODEL.md](MODEL.md) 第 5 节。

## 1. 探测模型

- **ISW `t`-probing 模型**（Ishai–Sahai–Wagner 2003）：攻击者选定至多 `t`
  个内部线路（wire）并得到其精确值；gadget 安全当且仅当这 `t` 个值的联合
  分布与秘密无关。
- **软件 word-probing 实例**：线路即 C 源码中每个被物化的 32 位中间量，
  包括多字运算的进位/借位，以及编译器为同一表达式生成的任何临时寄存器。
  后者意味着证明必须在目标汇编上复核（见第 6 节）。
- **模拟器表述**：下面所有性质都用“能否由一个只读取部分输入份额的模拟器
  完美生成探针分布”表述。`t = 1` 时只需对每一条线路单独验证其边缘分布。

## 2. `t`-NI（Non-Interference，Barthe et al. 2016）

gadget `G` 是 `t`-NI，若对任意 `t_int` 个内部探针与 `t_out` 个输出探针，
`t_int + t_out <= t`，存在模拟器仅凭每个输入 sharing 的至多 `t_int + t_out`
个份额即可完美模拟全部观察。

`t = 1` 实例：每条线路（含输出份额）至多依赖每个输入的一个份额，
或其分布与秘密无关。

## 3. `t`-SNI（Strong Non-Interference）

`G` 是 `t`-SNI，若在同样的探针预算下，模拟器只需每个输入的 **`t_int`**
个份额（输出探针不消耗输入份额）。

`t = 1` 实例等价于两条：
1. 任一内部线路至多依赖每个输入的一个份额；
2. 任一 **单个输出份额** 的分布与所有输入份额无关（通常由新鲜随机量保证）。

逐份额线性 gadget（`z_j = x_j + y_j`）满足 1 不满足 2，所以是 1-NI 而非
1-SNI。

## 4. MIMO-`t`-SNI（Cassiers–Standaert 2019）

对有 `m` 组输出的 gadget：对任意 `t_int` 个内部探针和每组输出上的探针集合
`O_1..O_m`，若 `t_int + max_i |O_i| <= t`，则可由 `t_int` 个输入份额模拟。

与“每组输出各自 SNI”的差别：MIMO 允许攻击者 **同时** 在每组输出各放
`t_out` 个探针。`t = 1` 实例：`t_int = 0` 时攻击者可同时观察 `eq` 的一个
份额和 `lt` 的一个份额，二者的联合分布必须与秘密无关；这要求两组输出使用
**相互独立** 的输出随机量。

## 5. `t`-PINI（Probe-Isolating Non-Interference，Cassiers–Standaert 2020）

`G` 是 `t`-PINI，若对任意 `t_int` 个内部探针和索引集合为 `O` 的输出份额
（`|O| + t_int <= t`），存在索引集合 `A`，`|A| <= t_int`，使全部观察可由
**所有输入** 的索引落在 `O ∪ A` 内的份额模拟。

要点：
- 份额索引对所有输入统一。`a_0 * b_1` 这种跨索引依赖违反 PINI，但不违反
  NI/SNI；
- 逐份额线性操作天然 PINI；
- PINI gadget 的任意（顺序、并行、带扇出）组合仍是 PINI，从而 `t`-probing
  安全。这是本仓库 FFT/LDL/split/merge 这类大扇出线性层的正确组合框架，
  而不是 MIMO-SNI；
- 单输入 `t`-SNI gadget 必然 `t`-PINI（A 取 SNI 模拟器用到的份额）；
- 两输入 SNI gadget 若每个输入前各有一个独立 SNI refresh，则组合为 PINI。

## 6. 本仓库采用的 `t = 1` 组合定理

把 gadget 分三类：

| 类 | 例子 | 单探针请求 |
|---|---|---|
| `L` 逐份额仿射 | `mrx160_add/sub/neg`、位移、NTT 蝶形、CRT 映射、FFT 蝶形的加减 | 输出份额 `j` → 各输入份额 `j` |
| `S` 输出刷新且 PINI | AND(含输入刷新)、refresh、A2B、B2A、`mul_wide`、`mul_public_*`、`sqr`、select/nonzero | 内部探针 → 至多一个统一索引；输出探针 → 0 份额 |
| `M` 双输入非 PINI 的 SNI | `ring_mul`、`ring_mul192`、`mrx160_mul`、`mrx160_mul_u64_hi`、`mmod32_mul` | 内部探针可能要 `a_0` 与 `b_1`；输出探针 → 0 份额 |

**定理（单探针组合）。** 设电路由上述三类 gadget 组成，随机数新鲜独立。
若对每个 `M` 类调用 `M(a, b)`，`a` 与 `b` 在调用图中 **没有仅经 `L` 类
gadget 可达的公共祖先 sharing**（等价：`a` 或 `b` 至少有一个是“份额新鲜”
的，即它的每个份额都含有一个上游 `S`/`M` 输出的新鲜随机量），则整个电路
是 1-NI（因而 1-probing 安全）；若电路输出来自 `S`/`M` 类 gadget，则还是
1-SNI。

证明梗概：单个探针落在唯一 gadget 中。请求集合沿 `L` 类向后传播时份额
索引不变，遇到 `S`/`M` 输出即终止（0 份额）。唯一可能同时请求同一
sharing 两个份额的情形是 `M` 类的 `a_0, b_1` 经 `L` 路径回溯到同一祖先；
这正是前提排除的情形（`ring_mul(a, a)` 反例）。因此每个 primary sharing
至多被请求一个份额，模拟器可完美生成。

**推论（PINI 视角）。** 在相同前提下，每个 `M` 调用连同其份额新鲜的上游
`S` gadget 构成一个 1-PINI 子 gadget，整个电路是 1-PINI 的组合，无需为
多输出线性层（FFT、LDL、split/merge、NTT 蝶形）单独证明 MIMO-SNI。

## 7. 审计方法（可复用 checklist）

1. **逐 gadget 单线路分类**：对每条 C 中间量判定 `0-share`/`PINI`/`NI`/
   `FAIL`。可执行版本见 [probe_check.py](probe_check.py)，对缩小位宽实例
   穷举全部秘密、份额分解与随机数。
2. **别名/扇出规则**：对每个 `M` 类调用检查两操作数的 `L` 祖先；同一
   sharing、其移位视图、或仅经加减得到的两个值不能同时作为两个逻辑输入。
3. **输出隔离**：声明 SNI 的 gadget，其每个输出份额必须含有本 gadget 内部
   新鲜随机量；多输出 gadget 每组输出用不同随机字（MIMO 条件）。
4. **打开边界**：所有 `*_open*` 调用列入泄漏契约，并单独论证所打开的位与
   秘密的统计关系（例如 BerExp 的 `eq` 在 PRNG 字节保密时是与秘密独立的
   Bernoulli(1/256)，而 SamplerZ 的 accept 概率依赖秘密 `x`）。
5. **明文边界**：所有未分享的秘密相关局部变量（候选样本、符号位、PRNG
   字节、hedge token）都是单探针即可读取的值，必须单列。
6. **机器码复核**：对任何在同一表达式/相邻语句中同时读取 `share[0]` 与
   `share[1]` 并与随机量异或/相加的代码，用 `objdump -d` 确认编译器没有先
   合并两个份额。值屏障必须放在每个份额与随机量结合之后、两者再结合之前。
7. **随机数账目**：每个 gadget 的随机字消费固定且只依赖公开参数；不复用。
8. **模型外泄漏**：transition/Hamming-distance（两份额先后经过同一寄存器
   或总线）、glitch、cache、时序。两份额软件实现中 transition 泄漏是最常见
   的一阶失效来源，源码级 NI/SNI 结论不覆盖它。

## 8. 参考文献

- Ishai, Sahai, Wagner. *Private Circuits.* CRYPTO 2003.
- Barthe, Belaïd, Dupressoir, Fouque, Grégoire, Strub, Zucchini.
  *Strong Non-Interference and Type-Directed Higher-Order Masking.* CCS 2016,
  ePrint 2015/506.
- Cassiers, Standaert. *Towards Globally Optimized Masking: From Low
  Randomness to Low Noise Rate (MIMO-SNI).* TCHES 2019(2).
- Cassiers, Standaert. *Trivially and Efficiently Composing Masked Gadgets
  with Probe Isolating Non-Interference.* IEEE TIFS 2020.
- Bettale, Coron, Zeitoun. *Improved High-Order Conversion From Boolean to
  Arithmetic Masking.* TCHES 2018(2).
- Schneider, Paglialonga, Oder, Güneysu. *Efficiently Masking Binomial
  Sampling at Arbitrary Orders.* PKC 2019.
- Berthet, et al. *Masked Computation of the Floor Function and Its
  Application to the FALCON Signature.* IACR CiC 2024.
- Howe, Prest, Ricosset, Rossi. *Isochronous Gaussian Sampling.* PQCrypto
  2020（SamplerZ 接受/重试次数泄漏的论证）。
- Guerreau, Martinelli, Nicolas, Wallet. *The Hidden Parallelepiped Is Back
  Again: Power Analysis Attacks on Falcon.* TCHES 2022(3)（符号位泄漏攻击）。
- Zhang, Lin, Liu, Shen, Zhang. *Improved Power Analysis Attacks on Falcon.*
  EUROCRYPT 2023（BaseSampler/`z0` 泄漏攻击）。
- Balasch, Gierlichs, Grosso, Reparaz, Standaert. *On the Cost of Lazy
  Engineering for Masked Software Implementations.* CARDIS 2014（transition
  泄漏使阶数减半）。
