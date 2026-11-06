# 恢复可组合 NI/SNI 的修复顺序

下面顺序按“阻断最多上层证明”的程度排列。前一项未完成时，不应为受影响的
上层组件撰写正面 SNI 声明。

## 当前状态

`ring_mul()`、`ring_mul192()`、`mi32_to_mrx160()` 和 `mmod32_mul()` 的
源码级修复已经完成；对应证明见 [REPAIRED_PROOFS.md](REPAIRED_PROOFS.md)。
下面保留其设计要求作为审计基线。固定 LFSR 替换和复合 gadget 的 MIMO 组合
证明仍未完成；明文兼容桥已由真实掩码签名核心替代。

## P0：替换 demo 随机源

1. 顶层调用从真实熵或经过明确安全假设的独立密码学 PRNG 获得 mask seed；
2. mask RNG 与签名 SHAKE 流严格域分离；
3. 不在每次调用重置为相同状态；
4. reservoir 只负责批量消费，从不重复发放 bit/word；
5. 压缩掩码密钥的展开种子不可预测，且不同对象独立；
6. 文档区分“真随机下的信息论 NI/SNI”和“CSPRNG 下的计算安全归约”。

固定 LFSR 可以保留为 KAT 专用实现，但不能进入默认或安全基准构建。

## 已完成：160/192 位环乘

不要先保存 `p01+p10`。安全的两份额顺序应让每个跨份额累加器从第一条指令
起就被全宽随机量遮蔽，例如：

```text
r  <- uniform ring element
z0 <- r
z0 <- z0 + a0*b0

z1 <- -r
z1 <- z1 + a0*b1
z1 <- z1 + a1*b0
z1 <- z1 + a1*b1
```

每一步之间设置编译器值屏障；不要把三个乘积写在同一 C 表达式中。每个单独
乘积只使用每个输入的一个份额，`z0`/`z1` 的每个部分和都由均匀 `r` 遮蔽。
对 160 位和 192 位分别实现、证明和审查汇编。

修复后需增加以下测试：

- 代数正确性与现有 fixed-point oracle 逐位一致；
- 自动扫描汇编，确认交叉项相加前已有随机累加器；
- 对缩小位宽版本穷举所有秘密、分享和单探针分布；
- 对每个输出份额验证分布与秘密无关。

## 已完成：`mi32_to_mrx160`

复用安全 A2B 模式：

```text
arith_to_bool_words(z0,z1,&a0,&a1,1,rng)
公开符号扩展到 160 位 Boolean sharing
mrx160_b2a(...)
refresh
```

不得再次构造 `(a0,0)`、`(a1,0)`。增加回归检查：任何 A2B 调用点若有一侧
Boolean share 恒为零，应被静态检查拒绝，公开常数编码除外。

## 已完成：Montgomery `mmod32_mul`

`mmod32` 的两个份额长期保存 `x*R mod q`，其中 `R=2^32`。在该域中先注入
新鲜均匀的 `rR`，再逐项 REDC 和累加：

```text
z0 = rR
z0 = z0 + REDC(a0*b0) mod q
z1 = -rR mod q
z1 = z1 + REDC(a0*b1) mod q
z1 = z1 + REDC(a1*b0) mod q
z1 = z1 + REDC(a1*b1) mod q
```

每个部分和实际存回 `[0,q)` 的 canonical Montgomery residue 并设值屏障；
不得先组合多个交叉项。`share/open` 是唯一域转换边界。正确性、范围和 SNI
证明见 [MONTGOMERY_PROOFS.md](MONTGOMERY_PROOFS.md)。
## 已完成：统一 Boolean 非线性 gadget

`mask_sampler.c` 的 word/bit AND、OR、select、nonzero、64 位秘密移位和 byte
comparator 已统一调用 `fndsa_mask_bool_and_words()`。共享原语对两个输入分别做
独立 refresh，以切断别名和扇出，再用第三个随机字初始化输出 accumulator；四个
交叉项逐项计算、逐项加入，并在每一步设置值屏障，不再存在单个 RHS 聚合交叉项。

OR、bit select、nonzero 和 64 位 select 在返回前显式 refresh；byte comparator
的 `eq` 与 `lt` 使用不同随机字独立刷新。完整扇出、refresh 和 trace-open 边界见
[SAMPLER_BOOLEAN_MIMO.md](SAMPLER_BOOLEAN_MIMO.md)。结论仍依赖理想随机数和
源级 word-probing 模型；固定 LFSR 与机器级泄漏不在正面声明内。
## P1：重新验证固定点复合 gadget

低层修复后按以下顺序恢复证明：

1. `mul`、`mul_wide`、`mul_u64_hi`；
2. inverse、division、sqrt 的每轮 Newton/迭代；
3. expm；
4. BerExp；
5. SamplerZ。

每个证明要同时覆盖：

- 值探测 NI/SNI；
- 输出份额的 SNI 隔离；
- 随机字唯一消费；
- 编译后中间寄存器；
- 已公开判定位的显式泄漏契约。

## P1：处理 BerExp 与 SamplerZ 控制流

若目标仅是论文中明确的“值 masking + 兼容泄漏”，应把以下观察写进正式安全
接口，而不是称为完整 SNI：

- BerExp 相等前缀长度；
- SamplerZ 接受/拒绝位；
- 重试次数；
- 顶层签名重试次数。

若目标是 trace-oblivious，则 BerExp 必须固定读取并比较 8 个字节，不打开中间
相等位；SamplerZ 需要固定批次/可证明的有界采样策略或其他恒定轨迹构造。仅把
接受位分享后在循环条件处打开，不是 trace 安全。

## P2：FFT、split/merge、LDL 的组合证明

1. 为公开常数快速乘和宽乘建立不同的安全类型，防止公开常数意外被重分享后
   进入通用不安全路径；
2. 为 split/merge 画 MIMO-SNI 组合图，显式表示 copy、扇出和 refresh；
3. FFT 每层证明 twiddle 乘、butterfly 和份额重用满足组合条件；
4. LDL 对 division、乘法和两个输出 `d/e` 建立 MIMO-SNI 证明；
5. 动态树按递归深度归纳，公开遍历与秘密值路径分开证明。

功能范围的 Lean 定理只能证明无溢出/舍入正确，不能替代 probing 安全证明。
建议把安全组合图用可机检的 probing 工具或 EasyCrypt/Lean 概率语义另行建模。

## 已完成（功能路径）：移除明文兼容桥

公开 masked 签名 API 的功能路径已满足下列连接要求；但只有随机源、组合证明
和机器码/物理审计也完成后，才能声称端到端 NI/SNI：

- 掩码密钥在签名边界内始终以分离份额展开；生产后端还需要不可预测种子；
- 不存在也不调用 `reconstruct_legacy_key()`；
- masked NTT/CRT、Gram、FFT/LDL、SamplerZ、最终向量和范数检查全部接通；
- 只在协议规定的最终公开边界输出签名和整体成功位；
- 所有失败路径不生成明文私钥副本。

当前 `masked` 表示已经接通的值掩码数据路径，而不是已完成的端到端安全定理；
固定 LFSR 与显式 trace 边界仍阻止生产安全声明。

## 发布门槛

建议把以下条件设为发布阻断项：

- `security/counterexamples.py` 中所有“旧实现反例”必须仍可重现，且修复后新增
  的安全核心穷举测试必须通过；
- 固定 LFSR 只能在显式 `DEMO_ONLY` 构建中出现；
- 每个非线性 gadget 有对应的 NI/SNI 定理、随机数预算和调用点清单；
- 每个目标架构保存并审查 gadget 汇编；
- 禁止 LTO 跨越 gadget 边界；
- 统计泄漏测试只作为实现检查，不替代证明。
