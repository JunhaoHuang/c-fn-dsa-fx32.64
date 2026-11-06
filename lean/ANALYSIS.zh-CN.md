# 固定点 Falcon 精度与范围证明审计

## 结论

固定点 Falcon 把两类结论分开：范围/有界性由筛选密钥不变量和算法结构推出；精度安全定理则依赖若干中间量的误差界，而端到端常数主要来自有限样本测量。测试中没有溢出、观察到 RMSE 很小，都不能直接提升为对所有合法密钥、消息、随机币和递归路径成立的定理。本目录没有把实验阈值写成公理。

## 仍属经验性的部分

1. **端到端 RMSE 门限**：作者的生产路径基准以 256 位 mpmath 为参考，对 `Gram → ffLDL → σ_i → target` 做有限密钥聚合。门限 `2^-38, 2^-43, 2^-57, 2^-40` 来自观测 RMSE 再加约 3 bit 余量；它们是回归测试门，不是全称误差上界。
2. **ffSampling 递归误差**：实验按递归层统计平均误差和观测最大误差。有限试验的“最大值”不是数学上界；严格证明需要对 split/merge、twiddle 量化、LDL 乘除和每个递归子树建立统一归纳不变量。
3. **SamplerZ 决策保真**：生产路径基准明确不覆盖采样器决策保真，而把 SamplerZ 当作理想采样器。靠近整数边界的 `floor(mu)`、BerExp 比较及拒绝循环会把小数值误差放大成离散控制流差异。有限 KAT 不能证明所有输入相同，也不能直接给出分布距离界。
4. **参数覆盖**：作者原型明确把固定点签名路径限制在 Falcon-512；Falcon-1024 仍需重新验证各递归层的精度下降。
5. **Newton 最终常数**：倒数和逆平方根的实数收敛形式可解析证明，但“5/6/8 轮在所有实际输入上达到某 bit 界”还必须计入每轮定点舍入、规格化和输入区间。

## 已由 Lean 严格证明

| 文件 | 定理 | 含义 |
|---|---|---|
| `FixedPoint.lean` | `q32_64_fits_s160` | signed Q32.64 可无损提升到 signed 160 位 |
| | `q32_64_product_fits_s192` | 两个 Q32.64 原始值的精确乘积适合 signed 192 位 |
| | `product_fits_s160_of_magnitude` | 160 位快速乘法的通用指数条件 |
| | `roundShift64_error` | `rx_mul` 的整数舍入误差至多半 ulp |
| | `qmul_real_error` | 解码后的单次 Q32.64 乘法绝对误差至多 `2^-65` |
| | `roundShift64_tie_pos/neg` | 半值严格向 `+∞`，与当前 C 代码一致 |
| | `rounded_product_fits_q32_64` | 带半 ulp 余量时舍入结果不越过 signed 96 位 |
| `MaskedRing.lean` | `reconstruct_refresh` | `s0 += r; s1 -= r` 不改变重构值 |
| | `reconstruct_mul` | 两份额 160 位环乘法重构为明文环乘积 |
| `ErrorPropagation.lean` | `add_abs_error/sub_abs_error/mul_abs_error` | 确定性绝对误差组合 |
| | `butterfly_error` | 一层 radix-2 butterfly 的全称误差界 |
| | `errorRec_two` | 深度递归误差公式 |
| `RangeBudgets.lean` | `gram00_lt_2pow17` 等 | Falcon-512 筛选密钥范围预算 |
| | `drift_lt_222000` | 在 Lemma-13 解析前提下漂移严格小于 222000 |
| | `ffsampling_*_lt_2pow*` | tweak 路径适合 18 bit，标准路径适合 21 bit |
| | `q32_64_fast_path` | `m1+m2≤31` 时原始 Q32.64 乘积适合 signed 160 位 |
| `DynamicBackend.lean` | `q32_64_half_fits` | 算术半移位保持 Q32.64 有符号范围 |
| | `q0_64_twiddle_product_fits_s160` | 公开 Q0.64 twiddle 与任意 Q32.64 输入的原始乘积适合 signed 160 位 |
| | `div_newton_*_fits_s160` | 规格化除法 Newton 局部乘法在列明幅度不变量下适合 signed 160 位 |
| | `sqrt_*_raw_lt_2pow*` / `sqrt_newton_product_fits_s160` | 平方根输入缩放与局部迭代乘法预算 |
| | `fft512/1024_public_twiddle_budget` | 公开 FFT 深度的指数预算；不替代秘密递归幅度证明 |
| | `enabled_site_exponent_sum_le_159` / `dynamic_fallback_covers_q32_products` | 可执行 fast-path 策略与未闭合调用点的 signed-192 回退 |

源码不含 `sorry` 或自定义 `axiom`；Lean/Mathlib v4.35.0 的最终 `lake build` 已成功检查 9053 个 job。C 调用点与定理/回退策略见 `CALL_SITES.md`。

## 范围无溢出的准确含义

- 96 位 Q32.64 载荷总能安全符号扩展到 160 位。
- 任意两个完整 Q32.64 载荷的未移位乘积总能放入 signed 192 位。
- 只有调用点证明原始幅度指数和满足条件时，才可使用 160 位快速环乘法。
- 无法给出该条件的调用点必须走 192 位路径，不能依赖模 `2^160` 回绕。
- 乘法舍入后的 96 位存储还需输出范围带半 ulp 余量，由 `rounded_product_fits_q32_64` 表达。

## 尚未闭合的端到端证明义务

1. 继续把 `CALL_SITES.md` 的人工映射提升为 C 实现到 Lean 模型的机器检查精化关系。
2. 为合法密钥条件下的每条秘密数据边给出精确区间，特别是所有除数下界。
3. 当前已证明除法与平方根 Newton 乘法的局部 160 位预算；仍需形式化完整 5/8 轮含舍入递推及最终误差界。
4. 当前已证明公开 Q0.64 twiddle 的原始乘积适合 signed 160 位；仍需把 twiddle 量化误差纳入 `butterfly_error`。
5. 对 ffLDL/ffSampling 做结构归纳，闭合 512/1024 的秘密操作数幅度不变量和每个 SamplerZ 输入误差界。
6. 证明 SamplerZ 离散决策稳定性，或证明实现分布到理想分布的统计/Rényi 距离。
7. 在独立密码学 mask RNG 和物理测量下完成侧信道验证；固定 LFSR 只能用于功能复现。

在这些义务完成前，可以严格声称“局部定点算术、160/192 位乘法策略及筛选密钥范围已形式化”，不能声称“论文全部经验精度分析已经被 Lean 完全替代”。
