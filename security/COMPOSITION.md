# 复合组件与端到端传播分析

本文件把局部证明沿当前调用图传播。三个原始 1-NI 反例已修复，但局部 1-SNI
不会自动给多输出 FFT/LDL、可观察采样循环或完整签名调用图带来端到端安全性。

## 1. 已闭合的底层路径

以下组件在 IDEAL-RNG、源码级一阶 word-probing 模型下已有条件 1-SNI 证明：

- `ring_mul()`、`ring_mul192()`；
- `mrx160_mul()`、`mrx160_mul_wide()`、`mrx160_mul_u64_hi()`、`mrx160_sqr()`；
- `mrx160_mul_public_fast()`、`mrx160_mul_public_wide()`、`mrx160_mul_public_u64_fast()`；
- `mi32_to_mrx160()` 与 `mrx160_to_mi32()`；
- `mmod32_mul()`；
- A2B、B2A、refresh、`floor`、`half`。

证明见 [REPAIRED_PROOFS.md](REPAIRED_PROOFS.md) 和
[GADGET_PROOFS.md](GADGET_PROOFS.md)。固定 demo LFSR 仍使这些结论不能作为
当前可执行实例的无条件信息论声明。

公开常数乘法的局部 1-SNI 来自逐份额公开线性乘法、SNI A2B/B2A 和最终
refresh；fast 站点另有 Lean signed-160 范围证明。该结论不跨越 FFT/LDL 的
多输出扇出，因此第 5 节状态不变。

## 2. inverse、division 与 sqrt

这三个函数原先因 `ring_mul` 反例必然失败；该阻断已消除。它们的秘密规格化、
桶形移位、masked select 和所有 Newton 乘法现在都可由 NI/SNI 基础 gadget
构造。

当前状态改为 `REVIEW`，而不是 `FAIL-1-NI`：还需要为完整五轮/八轮迭代画出
组合图，检查同一 sharing 的扇出、每个线性 NI 节点后的使用次数，以及最终输出
是否确实来自 SNI 节点。在此证明完成前，不把局部乘法 SNI 直接升级为整个迭代
函数的正式 SNI 声明。

## 3. 模 `q/q'` 层

| 函数 | 当前结论 | 说明 |
|---|---:|---|
| `mmod32_add`, `mmod32_sub` | `1-NI*` | 逐份额公开模数运算 |
| `mmod32_refresh` | `1-SNI*` | 新鲜均匀模 `q` refresh |
| `mmod32_mul` | `1-SNI*` | 随机优先 Montgomery accumulator、单项 REDC、最终 refresh |
| `mmod32_share` | `BOUNDARY` | 明文模值分享 |
| `mmod32_open` | `OPEN` | 显式重构 |

当前 masked NTT/iNTT、Gram/cross、CRT 与中心化已经接入签名核心。模乘的
已知一阶反例不再阻断该路径，但 butterfly 扇出、CRT 多输出和跨阶段复用仍需
完整的 MIMO-SNI 组合证明。

## 4. sampler

`mrx160_expm_p63()` 的底层高位乘法已修复，状态从必然失败改为 `REVIEW`。
Horner 链是单输出顺序组合，预计可以由 SNI 乘法与 NI 线性层闭合，但仍需逐节点
证明。

`mrx160_ber_exp()` 和 `mrx160_sampler_z()` 仍不能给出严格 trace-oblivious
SNI 声明：

- sampler Boolean 非线性节点已统一到输入隔离、随机优先的共享 AND；OR/select/nonzero 与 64 位 select 显式刷新，byte comparator 的 `eq/lt` 具有独立 MIMO 输出 refresh；
- BerExp 打开逐字节相等位控制懒比较；
- SamplerZ 打开接受位并暴露重试次数；
- BaseSampler 候选和签名 PRNG 字节仍以明文局部变量存在。

上述 Boolean 子图的条件 MIMO-1-SNI 组合见
[SAMPLER_BOOLEAN_MIMO.md](SAMPLER_BOOLEAN_MIMO.md)；该结论不跨越随后显式打开的
`eq` 或接受位。

这些是独立于本次三个算子修复的剩余边界。

## 5. FFT、split/merge 与 LDL

旧 `mi32_to_mrx160` 和 wide 乘法反例已经消除。因此以下函数不再标为
`FAIL-1-NI`，而统一标为 `REVIEW`：

- `mxpoly_fft_i32_selfadj()`、`mxpoly_fft_i32()`；
- `mxpoly_mul_fft()`、`mxpoly_div_selfadj_spec_fft()`；
- `mxpoly_mul_iXadj_fft()`、`mxpoly_ldl_fft()`；
- split/merge 组合。

逐系数 `mxpoly_add/sub/neg` 仍为 1-NI。剩余证明义务是 MIMO-SNI：FFT
butterfly、split/merge 和 LDL 都会让一个输入 sharing 流向多个输出。单个乘法
1-SNI 不能替代扇出图证明；必须按 FALCON floor 论文的组合图条件标记 copy、
Split、NI 节点和 SNI refresh。

## 6. 动态 ffSampling

动态树不再继承本次三个底层反例，但仍依赖尚未闭合的 FFT/LDL MIMO 证明和
SamplerZ 的显式接受/重试边界，所以状态为 `REVIEW + TRACE-OPEN`。

递归深度、位反转索引、左右遍历和 flags 是公开控制流。`tree_clear()` 只负责
返回后的擦除，不提供执行期 probing 安全性。

## 7. 掩码密钥与顶层签名

V2 签名入口不再包含 `reconstruct_legacy_key()` 或明文签名器调用。它分别展开
q/q-prime 份额，依次执行 NTT、Gram/cross、整体 CRT 中心化、FFT/LDL、两棵
动态树、最终 `s0/s1` 和掩码验收。`s0` 不打开；接受后只打开签名本来公开的
`s1`。密钥范数、BerExp/SamplerZ 控制和整体尝试结果在规定边界打开。

这消除了旧兼容桥的直接明文反例，但不自动给公开 API 一个端到端 SNI 定理：
固定 LFSR 不满足理想随机前提；40 字节 hedging token 是未掩码可信边界；
FFT/LDL 尚缺完整 MIMO-SNI 图证明；BerExp、SamplerZ 和总重试轨迹可观察。
因此入口状态是 `REVIEW + TRACE-OPEN`，不是 `OPEN`，也不是已证明的 `1-SNI`。

## 8. 当前最大声明

在理想新鲜随机数和源码级一阶 word-probing 模型下，可以声明已单独证明的
基础算子及其 fixed-point 包装器为条件 1-SNI。完整签名已经是值路径掩码实现，
但当前最大端到端结论仍是“研究原型，组合证明未闭合且轨迹显式开放”，不能
声明完整 sampler、FFT/LDL 动态树或公开 masked 签名 API 为 SNI。剩余义务是
trace 泄漏契约、MIMO-SNI 组合证明、机器码审计和真实随机掩码后端。
