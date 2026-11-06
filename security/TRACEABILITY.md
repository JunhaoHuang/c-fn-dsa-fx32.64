# 逐函数安全结论追踪表

本表覆盖当前五个掩码实现源文件。状态定义见 [README.md](README.md)。静态 helper
若只实现公开整数运算，其安全性取决于调用者传入的是单份额值、公开值还是已经
组合的秘密相关值；这类函数标为 `LOCAL`，不单独授予共享 gadget 证明。

## `mask160.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `value_barrier_u32`, `mask_clear` | `PUBLIC` | 编译器约束/擦除；不产生 masking 安全性 |
| `lfsr_step_u64` | `PUBLIC` | 确定性状态更新，不是密码学随机源 |
| `fndsa_mask_rng_init_demo` | `FAIL assumption` | 固定、跨调用重复的公开 KAT 状态 |
| `fndsa_mask_rng_next_u64/u32/bit/mod` | `FAIL assumption` | 消费顺序不复用，但输出不满足未知独立随机前提 |
| `mul64_hi`, `barrett_reciprocal`, `barrett_reduce_u64`, `fndsa_mask_reduce_u64` | `PUBLIC/LOCAL` | 对公开模数恒定结构；秘密安全取决于调用数据流 |
| `random160`, `random192` | `LOCAL` | 理想 RNG 下生成均匀环元素；demo 中前提失败 |
| `add5/sub5/neg5/add6/sub6` | `LOCAL` | 用于单份额时保持 NI；不能安全组合两个秘密份额 |
| `mul5`, `mul6` | `LOCAL` | 单独乘两个输入份额可由各一个份额模拟；交叉积组合另行分析 |
| `public_constant_words`, `mask_rng_discard_words` | `PUBLIC/LOCAL` | 公开符号扩展；固定消费的随机字不得复用 |
| `fndsa_mask_bool_and_words` | `1-SNI*` | 输入独立刷新、随机优先逐项累加的两份额 AND |
| `fndsa_mask_bool_refresh_words` | `1-SNI*` | Boolean refresh |
| `lshift_words`, `bool_asr64`, `bool_sign_extend_96`, `bool_asr64_192_to_160` | `1-NI*` | 公开固定逐份额位变换；单独无 SNI 输出 refresh |
| `bool_add_words` | `1-NI*` | 安全前缀网络；最低输出位说明其本身不 SNI |
| `arith_to_bool_words`, `mrx160_a2b`, `a2b192` | `1-SNI*` | 安全重分享 + NI 加法 + 输出 refresh |
| `psi_words` | `LOCAL` | BCZ Psi 叶；安全结论由完整 B2A 构造承担 |
| `psi_b2a_words`, `mrx160_b2a`, `b2a192` | `1-SNI*` | BCZ18 `C2` 特化 |
| `mrx160_share_rx` | `BOUNDARY` | 接收并暂存明文 Q32.64 |
| `mrx160_open_rx` | `OPEN` | 显式重构 |
| `mrx160_copy/add/sub/neg` | `1-NI*` | 逐份额线性 |
| `mrx160_refresh` | `1-SNI*` | 算术 refresh |
| `ring_mul`, `ring_mul192` | `1-SNI*` | `r/-r` 在首个交叉项之前建立；所有部分和保持均匀 |
| `mrx160_mul`, `mrx160_mul_wide`, `mrx160_mul_u64_hi` | `1-SNI*` | SNI 环乘与安全 A2B/B2A、公开截断组合；快速乘法禁止同 sharing 双输入 |
| `mrx160_mul_public_fast`, `mrx160_mul_public_u64_fast` | `1-SNI*` | 公开逐份额乘法为 1-NI，随后 SNI A2B/B2A/refresh；需 signed-160 调用点定理 |
| `mrx160_mul_public_wide` | `1-SNI*` | 单个秘密输入安全扩至 192 位后逐份额乘公开常数，再安全截断 |
| `mrx160_sqr` | `1-SNI*` | 两个独立 arithmetic refresh 后调用双输入 SNI 乘法 |
| `mbit_nonzero_word`, `mb160_select` | `1-SNI*` | 由加固 AND 组成；结论限源级模型 |
| `mbit_not`, `mb160_lshift1/rshift1`, `mb160_lshift_public/rshift_public` | `1-NI*` | 逐份额公开线性/仿射操作 |
| `mrx160_public_raw` | `PUBLIC` | 公开常数注入 share 0 |
| `mrx160_floor` | `1-SNI*` | 安全转换与固定投影组合 |
| `mrx160_half_internal`, `mrx160_half` | `1-SNI*` | 安全转换与固定移位组合 |
| `mrx160_inv_spec`, `mrx160_div_spec`, `mrx160_sqrt_spec` | `REVIEW` | 旧乘法反例已消除；完整迭代组合证明待闭合 |
| `mbit_open_compat`, `mb64_open_compat` | `OPEN` | 显式打开 |
| `mmod32_share` | `BOUNDARY` | 明文模值分享 |
| `mmod32_open` | `OPEN` | 显式重构 |
| `mmod32_add/sub` | `1-NI*` | 逐份额公开模数运算 |
| `mmod32_mul` | `1-SNI*` | 单项 Montgomery REDC 并逐项加入已掩码模 `q` accumulator |
| `mmod32_refresh` | `1-SNI*` | 理想均匀模 `q` 随机量下 |
| `mmod32_encode_public` | `PUBLIC` | 公开常数进入 Montgomery 域 |
| `mmod32_set_residue_shares` | `1-SNI*` | 逐份额编码后执行独立模刷新 |
| `mmod32_decode_residue_shares` | `1-NI*` | 分别解码份额，不重构秘密 |
| `mmod32_mul_public_mont` | `1-NI*` | 逐份额乘公开 Montgomery 常数 |
| `mi32_share` | `BOUNDARY` | 接收明文整数 |
| `mi32_open_compat` | `OPEN` | 显式重构 |
| `mi32_to_mrx160` | `1-SNI*` | 安全 A2B、Boolean 符号扩展、B2A、refresh |
| `mrx160_to_mi32` | `1-SNI*` | A2B 后对整数 word 做 B2A |
| `mi32_center_mod_m` | `1-SNI*` | A2B、固定前缀减法、掩码选择、Psi-B2A 与最终刷新 |
| `mbit_to_mmod32`, `mi32_to_mmod32` | 源码级 `1-SNI*`；`MACHINE-FAIL` | Clang 18 `-O2` 物化 `share[1]^share[0]`（见 REVERIFY_2026-11-05.md F1）；待加值屏障 |
| `mrx160_inv_spec`, `mrx160_div_spec`, `mrx160_sqrt_spec`（更新） | `1-SNI*` | 2026-11-05 复核：PINI 规范化 + 乘法别名规则 + 末尾 refresh |
| `mi32_poly_sqnorm_le`, `mi32_signature_is_acceptable`（更新） | `1-SNI*` | 2026-11-05 复核：输出位含 AND 新鲜随机量 |

## `mask_mq.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `mod_data`, `valid_logn` | `PUBLIC` | 仅处理公开模数和维数 |
| `mmpoly_add/sub/neg` | `1-NI*` | 逐系数线性操作；接口刻意不刷新 |
| `mmpoly_refresh` | `1-SNI*` | 每个输出系数使用独立模随机量 |
| `mmpoly_ntt/intt` | 条件 `1-SNI*` | 内部线性 1-NI；变换边界逐系数刷新 |
| `mmpoly_mul_ntt` | `1-SNI*` | 每系数调用 `mmod32_mul` |
| `mmpoly_mul_public_ntt` | `1-SNI*` | 公开线性乘法后逐系数刷新 |
| `mmpoly_gram_ntt`, `mmpoly_cross_ntt` | 条件 `1-SNI*` | 独立乘法与最终刷新；输出不得与输入重叠 |
| `mmpoly_crt_center` | 条件 `1-SNI*` | 逐份额 CRT 后调用安全整体中心化；不打开比较位 |
| `mmpoly_from_mi32` | 源码级 `1-SNI*`；`MACHINE-FAIL` | 继承 `mi32_to_mmod32` 的编译器失效 |

## `mask_sampler.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `sampler_clear` | `PUBLIC` | 擦除临时量 |
| `mrx_public_raw`, `mrx_public_u64` | `PUBLIC` | 公开常数注入 |
| `mb64_to_mrx`, `mrx_to_mb64` | `1-SNI*` | SNI B2A/A2B 与公开投影组合 |
| `mask_and_word`, `mbit_and` | `1-SNI*` | 统一调用输入隔离、随机优先、逐项屏障的共享 AND |
| `mbit_set_public` | `1-SNI*` | 公开 bit 的新鲜 Boolean sharing |
| `mbit_not`, `mb64_rshift_public` | `1-NI*` | 逐份额公开操作 |
| `mbit_or/select/nonzero32`, `mb64_select/secret_rshift` | 条件 `1-SNI*` | 非线性节点统一 AND，返回边界独立 refresh；见 `SAMPLER_BOOLEAN_MIMO.md` |
| `mb64_compare_byte` | `MIMO-1-SNI*` | `probe_check.py` 穷举确认 `(eq_i, lt_j)` 联合 0-share |
| `mrx_half_positive` | `1-SNI*` | A2B、固定右移、B2A、refresh |
| `mrx160_expm_p63` | `1-SNI*` | 2026-11-05 复核：Horner 链每步 `y` 份额新鲜，输出经 A2B 刷新 |
| `mrx160_ber_exp` | `REVIEW + TRACE-OPEN` | 底层乘法已修复；逐字节相等位仍显式打开 |
| `source_next_u64/u16` | `BOUNDARY` | 明文组装签名 PRNG 输出；不分享随机候选 |
| `mrx160_sampler_z` | `TRACE-OPEN + BOUNDARY-CRITICAL` | 接受位可观察；`lo/hi/z0/b/zi` 明文（见 REVERIFY_2026-11-05.md F2） |
| `fndsa_masked_gaussian0_helper` | `BOUNDARY` | 外部明文 BaseSampler helper，候选不掩码 |

## `mask_xpoly.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `mpublic_rx`, `mpublic_i32` | `PUBLIC` | 公开常数注入 |
| `mmul_public_ufrac` | `1-SNI*` | 调用 `mrx160_mul_public_u64_fast`；split/merge 范围由 Lean 定理覆盖 |
| `mmul_public_rx`, `mmul_public_rx_tree` | `1-SNI*` | general 路径调用 public-wide；树内已证明站点调用 public-fast 并保持旧 RNG 预算 |
| `mcomplex_mul` | `REVIEW` | 三个 SNI wide 乘法；双输出 MIMO 组合待证明 |
| `mi32_add_local/sub_local/double_local` | `1-NI*` | 逐份额线性 |
| `mxpoly_fft_i32_selfadj`, `mxpoly_fft_i32` | `1-PINI*` | 2026-11-05 复核：逐份额蝶形 + PINI 转换；不需要 MIMO-SNI |
| `mxpoly_add/sub/neg` | `1-NI*` | 逐系数线性 |
| `mxpoly_mul_fft`, `mcomplex_mul` | `1-PINI*` | 2026-11-05 复核：`mul_wide` 入口双 A2B，输出为线性组合 |
| `mxpoly_div_selfadj_spec_fft` | `1-SNI*` | 2026-11-05 复核：`inv_spec` 输出份额新鲜 |
| `mxpoly_mul_iXadj_fft`, `mxpoly_mul_iXadj_tree_fft` | `1-SNI*` | 2026-11-05 复核：`bt` 份额新鲜 |
| `mxpoly_ldl_fft`, `mxpoly_ldl_tree_fft` | `1-PINI*`（`e` 为 `1-SNI*`） | 2026-11-05 复核：`d/e` 不需要 MIMO-SNI |
| `mxpoly_split_fft`, `mxpoly_split_selfadj_fft`, `mxpoly_merge_fft` | `1-PINI*` | 逐份额线性 + 单输入 SNI 公开乘/half |

## `mask_tree.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `tree_clear` | `PUBLIC` | 擦除 |
| `tree_public_rx` | `PUBLIC` | 公开常数注入 |
| `tree_mul_public` | `1-SNI*` | `RX_INV_SIGMA`/`RX_IQ` 调用 public-fast，范围由 Q0.64 定理覆盖 |
| `fndsa_masked_dyn_tmp_size` | `PUBLIC` | 公开大小计算 |
| `fndsa_mask_tree_bind` | `PUBLIC/BOUNDARY` | workspace/公开参数绑定，不证明其中秘密值安全 |
| `fndsa_mask_tree_work`, `fndsa_mask_tree_samples` | `PUBLIC` | 返回工作区指针 |
| `mask_ffsampling_deepest` | `REVIEW + TRACE-OPEN` | 底层转换已修复；仍依赖 SamplerZ 泄漏边界 |
| `mask_ffsampling_left_inner`, `mask_ffsampling_right_inner` | `REVIEW + TRACE-OPEN` | LDL/MIMO 待证明并继承叶子 trace 边界 |
| `fndsa_mask_ffsampling_left/right` | `REVIEW + TRACE-OPEN` | 完整动态递归入口 |
| `fndsa_mask_tree_clear` | `PUBLIC` | 执行后擦除，不修复执行期泄漏 |

## `sign_masked.c`

| 函数 | 状态 | 备注 |
|---|---:|---|
| `secure_clear`, `dec32le`, `enc32le`, `enc64le`, `ct_equal` | `PUBLIC/LOCAL` | 编码、比较、擦除 helper |
| `fndsa_masked_key_size`, `fndsa_masked_tmp_size` | `PUBLIC` | 公开参数大小 |
| `shake_next_u32le`, `expand_share1` | `LOCAL` | 确定性份额展开；安全依赖不可预测 seed |
| `decode_private`, `derive_public_and_G` | `BOUNDARY` | 明文私钥处理 |
| `fndsa_masked_key_import_demo` | `BOUNDARY + FAIL assumption` | 明文导入；固定 mseed 破坏隐藏分享前提 |
| `fndsa_masked_key_clear` | `PUBLIC` | 擦除 |
| `validate_mkey` | `PUBLIC/LOCAL` | 容器检查，不提供 masking 证明 |
| `fndsa_masked_key_expand_mq_internal` | `1-SNI*` | 不重构私钥；分别约简/编码 q、q' 份额并刷新 |
| `arena_take`, `masked_attempt_clear` | `PUBLIC/LOCAL` | 公开 arena 规划和显式擦除 |
| `masked_make_c_ntt` | `PUBLIC/LOCAL` | 公开 nonce/message challenge 的 NTT |
| `masked_crt_from_ntt`, `masked_neg_i32` | `REVIEW` | 已证明底层 gadget 的签名连接层；组合证明待闭合 |
| `sign_masked_core` | `REVIEW + TRACE-OPEN` | 真实 q/q-prime、FFT/LDL、采样、最终向量和聚合验收路径；无明文私钥重构 |
| `sign_masked_real` | `BOUNDARY + TRACE-OPEN` | 未掩码 hedge token、公开消息哈希和整体尝试轨迹边界 |
| `fndsa_sign_masked_temp`, `fndsa_sign_masked_seeded_temp` | `REVIEW + TRACE-OPEN` | 调用真实掩码核心；接受后只打开公开 s1，s0 不打开 |

### 分组行的完整函数名

上表斜线缩写分别覆盖以下完整符号：`fndsa_mask_rng_next_u32`、`fndsa_mask_rng_next_bit`、`fndsa_mask_rng_mod`、`mrx160_add`、`mrx160_sub`、`mrx160_neg`、`mb160_rshift1`、`mb160_rshift_public`、`mmod32_sub`、`mbit_select`、`mbit_nonzero32`、`mb64_secret_rshift`、`mb64_compare_byte`、`source_next_u16`、`mi32_sub_local`、`mi32_double_local`、`mxpoly_sub`、`mxpoly_neg` 和 `fndsa_mask_ffsampling_right`。它们继承所在分组行的结论。

## 更新规则

任何被标为 `FAIL-1-NI` 的函数修改后，都应同时更新：

1. 其局部模拟证明或反例；
2. 本表；
3. [COMPOSITION.md](COMPOSITION.md) 中所有可达调用者；
4. 可执行反例/缩小位宽穷举；
5. 目标架构汇编审查记录。
