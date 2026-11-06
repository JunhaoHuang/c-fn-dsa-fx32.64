# Q32.64/160 位掩码实现的 NI/SNI 安全审计

本目录审计当前工作树中的 `mask160.c`、`mask_sampler.c`、
`mask_xpoly.c`、`mask_tree.c`、`mask_mq.c` 和 `sign_masked.c`。结论针对两个份额、
一阶探测（`t=1`）。

## 2026-11-05 复核

[REVERIFY_2026-11-05.md](REVERIFY_2026-11-05.md) 用单探针组合定理和
[probe_check.py](probe_check.py) 穷举重新校验了全部组件（定义见
[DEFINITIONS.md](DEFINITIONS.md)）。要点：原 `REVIEW` 项全部闭合为条件
`1-SNI*`/`1-PINI*`，FFT/LDL 不需要 MIMO-SNI；但发现 `mbit_to_mmod32()` 在
Clang 18 `-O2` 下生成未掩码的 `share[0]^share[1]`（`MACHINE-FAIL`，经
`mi32_to_mmod32`/`mmpoly_from_mi32` 影响签名核心），以及 SamplerZ 明文
BaseSampler 候选/符号位是实际攻击面。下文矩阵以该文件第 4 节为准。

## 当前结论

此前发现的三个源码级 1-NI 反例已经修复：

1. `ring_mul()` 和 `ring_mul192()` 现在先建立 `r/-r` accumulator，再逐个
   加入交叉乘积；不再生成未掩码的 `p01+p10`。
2. `mi32_to_mrx160()` 现在复用独立重分享的 word-sliced A2B；不再把输入
   份额编码为 `(ai,0)`。
3. `mmod32` 现在长期保持 canonical Montgomery 份额；每个交叉乘积单独
   REDC 并逐项加入已带模 `q` 随机量的 accumulator，不再生成同时包含三个
   交叉项的原始 `z1`。

在 IDEAL-RNG 和源码级 32 位 word-probing 模型下，这三个 gadget 及
`mrx160_mul()`、`mrx160_mul_wide()`、`mrx160_mul_u64_hi()` 可条件性声明
1-SNI。严格模拟证明见 [REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)。

当前 demo 仍不能声明无条件端到端 NI/SNI，原因是：

- 固定 LFSR 按要求继续保留为 `DEMO ONLY`，不满足新鲜、独立、未知随机数前提；
- sampler Boolean AND 已统一并给出条件 MIMO-SNI 图，但编译后机器码与转换泄漏仍待审计；
- BerExp/SamplerZ 仍打开相等/接受判定位并暴露重试轨迹；
- FFT/LDL/split/merge 的多输入多输出组合尚未完成 MIMO-SNI 图证明；
- V2 已移除明文私钥重构并接通真实签名值路径，但 FFT/LDL MIMO-SNI、机器码
  证明和顶层可观察控制流仍未闭合。

## 状态说明

| 标记 | 含义 |
|---|---|
| `1-SNI*` | 在理想随机数、源级 word-probing 模型下有一阶 SNI 模拟证明 |
| `1-NI*` | 在同一模型下有一阶 NI 模拟证明，但没有完整 SNI 输出隔离证明 |
| `REVIEW` | 已知底层反例已消除，但完整复合/MIMO/机器码证明尚未闭合 |
| `TRACE-OPEN` | 主动打开判定位或存在秘密相关可观察控制流 |
| `OPEN` | 显式重构秘密 |
| `BOUNDARY` | 接收明文秘密的分享/导入边界 |
| `PUBLIC` | 只处理公开元数据、常数或擦除 |

星号表示当前固定 LFSR 使理想随机数前提在 demo 可执行实例中不成立。

## 组件结论矩阵

| 组件 | 当前结论 | 依据 |
|---|---:|---|
| `fndsa_mask_bool_and_words` | `1-SNI*` | 输入独立刷新、随机优先逐项累加的两份额 AND |
| Boolean/arithmetic refresh | `1-SNI*` | 标准两份额 refresh |
| `bool_add_words` | `1-NI*` | 安全前缀网络；函数自身没有最终输出 refresh |
| `mrx160_a2b`, `a2b192` | `1-SNI*` | 独立重分享、masked add、最终 Boolean refresh |
| `mrx160_b2a`, `b2a192` | `1-SNI*` | BCZ18 两份额 `C2`/Psi 构造 |
| `mrx160_copy/add/sub/neg` | `1-NI*` | 逐份额线性操作 |
| `mrx160_floor`, `mrx160_half`, `mrx160_to_mi32` | `1-SNI*` | SNI 转换与公开位操作组合 |
| `mi32_to_mrx160` | `1-SNI*` | 修复为安全 A2B、Boolean 符号扩展、B2A、refresh |
| `ring_mul`, `ring_mul192` | `1-SNI*` | `r/-r` 在任何交叉项累加之前建立 |
| `mrx160_mul`, `mrx160_mul_wide`, `mrx160_mul_u64_hi` | `1-SNI*` | SNI 环乘与 SNI 转换/截断组合；快速乘法禁止同 sharing 双输入 |
| `mrx160_mul_public_fast/wide/u64_fast` | `1-SNI*` | 公开逐份额乘法 1-NI + SNI 转换/截断/refresh；fast 受 Lean 调用点范围约束 |
| `mrx160_sqr` | `1-SNI*` | 两个独立 refresh 切断自乘扇出后调用双输入 SNI 乘法 |
| `mmod32_add/sub` | `1-NI*` | 逐份额公开模数运算 |
| `mmod32_mul` | `1-SNI*` | 每个单份额乘积单独 Montgomery REDC 并加入已掩码 accumulator |
| `mmpoly_add/sub/neg` | `1-NI*` | 逐份额公开线性 NTT 域运算 |
| `mmpoly_ntt/intt` | 条件 `1-SNI*` | 线性 1-NI 变换和逐系数独立边界刷新 |
| `mmpoly_mul_ntt`, `mmpoly_mul_public_ntt` | `1-SNI*` | SNI 秘密乘法或公开线性乘法后刷新 |
| `mmpoly_gram_ntt`, `mmpoly_cross_ntt` | 条件 `1-SNI*` | 独立乘法、公开共轭索引和最终刷新；禁止输出别名 |
| `mmpoly_crt_center`, `mi32_center_mod_m` | 条件 `1-SNI*` | 逐份额 CRT 同态与安全 A2B/比较/选择/B2A |
| `mi32_to_mmod32`, `mmpoly_intt_center` | 条件 `1-SNI*` | A2B64 提取进位/符号并联合中心化；依赖转换组合条件 |
| `mi32_poly_sqnorm_le`, `mi32_signature_is_acceptable` | `REVIEW` | 值保持掩码并只打开聚合位；长链/MIMO 组合证明待闭合 |
| `mrx160_inv_spec`, `mrx160_div_spec`, `mrx160_sqrt_spec` | `REVIEW` | 旧乘法反例已消除；完整迭代组合证明待单独闭合 |
| sampler 内部 AND/OR/select/nonzero | 条件 `1-SNI*` | 统一输入隔离 AND 与显式输出 refresh；源码 word-probing 模型 |
| sampler byte comparator | 条件 `MIMO-1-SNI*` | `eq/lt` 独立 refresh；扇出图见 `SAMPLER_BOOLEAN_MIMO.md` |
| `mrx160_expm_p63` | `REVIEW` | 底层乘法已修复；完整 Horner 组合证明待闭合 |
| `mrx160_ber_exp` | `REVIEW + TRACE-OPEN` | 值路径待复核，逐字节相等位显式打开 |
| `mrx160_sampler_z` | `TRACE-OPEN` | 接受位和重试次数可观察 |
| `mxpoly_add/sub/neg` | `1-NI*` | 逐系数线性操作 |
| FFT、split/merge、LDL | `REVIEW` | 底层三个反例已修复；MIMO-SNI/扇出证明待闭合 |
| 动态 `ffSampling` | `REVIEW + TRACE-OPEN` | 依赖未闭合 MIMO 组件和 SamplerZ 泄漏边界 |
| `*_open*` | `OPEN` | 显式重构；生产签名核心仅在已列出的兼容/公开边界调用 |
| `*_share*`, `fndsa_masked_key_import_demo` | `BOUNDARY` | 明文秘密在函数内部出现 |
| `fndsa_sign_masked_*` | `REVIEW + TRACE-OPEN` | 真实掩码值路径；固定 LFSR、hedge 边界及可观察采样/重试阻止端到端 SNI 声明 |

## 文档索引

- [MODEL.md](MODEL.md)：探测模型、NI/SNI/MIMO-SNI 定义和证明边界；
- [DEFINITIONS.md](DEFINITIONS.md)：NI/SNI/MIMO-SNI/PINI 形式定义、`t=1` 组合定理与审计 checklist；
- [REVERIFY_2026-11-05.md](REVERIFY_2026-11-05.md)：全组件复核矩阵、乘法调用点别名审计与新发现；
- [probe_check.py](probe_check.py)：缩小位宽的单探针 NI/SNI/PINI/MIMO 穷举分类器；
- [SAMPLER_BOOLEAN_MIMO.md](SAMPLER_BOOLEAN_MIMO.md)：统一 Boolean AND、refresh、扇出与双输出 comparator 的 MIMO-SNI 组合图；
- [REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)：三个修复算子的严格模拟证明；
- [MONTGOMERY_PROOFS.md](MONTGOMERY_PROOFS.md)：Montgomery 范围、正确性与 1-SNI 证明；
- [NTT_CRT_PROOFS.md](NTT_CRT_PROOFS.md)：掩码 NTT/iNTT、CRT 与中心化的组合证明；
- [GADGET_PROOFS.md](GADGET_PROOFS.md)：基础 gadget 与旧反例背景；
- [COMPOSITION.md](COMPOSITION.md)：复合路径和仍待闭合的安全义务；
- [TRACEABILITY.md](TRACEABILITY.md)：逐函数结论；
- [REMEDIATION.md](REMEDIATION.md)：剩余端到端工作；
- [counterexamples.py](counterexamples.py)：旧区分器与新均匀 accumulator 的穷举检查。

## 证明强度

本目录给出源码级模拟证明和缩小模型穷举，不是 C 编译器或目标机器码的机械化
证明。它不覆盖 glitch、transition/Hamming-distance、缓存、分支预测、故障、
高阶或多迹攻击。`noinline` 与 value barrier 是必要的编译约束，但每个发布
目标仍须单独检查汇编。
