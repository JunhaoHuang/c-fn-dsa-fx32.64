# `mmod32` Montgomery 正确性与一阶 SNI 证明

本文证明当前 `mmod32` 的 32 位 Montgomery 表示、约简界和掩码乘法。
模型与随机性前提沿用 [MODEL.md](MODEL.md)：两个算术份额、一个 32 位 word
probe，以及每次请求均返回新鲜独立均匀值的 **IDEAL-RNG**。固定 LFSR 仍只
是 `DEMO ONLY` 后端，不满足最后一个前提。

## 1. 表示和参数

设公开奇数模数满足 `3 <= q < 2^32`，并令

```text
R   = 2^32
q0i = -q^(-1) mod R
R2  = R^2 mod q.
```

`mmod32` 的两个 canonical 份额均在 `[0,q)`，重构的内部值是

```text
(s0+s1) mod q = x*R mod q.
```

由于 `q` 为奇数，`R` 在 `Z/qZ` 中可逆，所以从普通域到 Montgomery 域的
映射 `x -> x*R` 是双射。普通域中的均匀随机量经过编码后仍然均匀。
因此随机 mask 可以直接从 `[0,q)` 均匀采样并解释为 Montgomery 余数，
无需再执行一次编码；`share/mul/refresh` 现在采用这种等分布实现。
随机数请求次数不变，但固定 demo RNG 对应的具体份额会改变。

生产模数 `12289`、`18433` 和 `226523137` 使用预计算的 `q0i`、`R2`。
其他公开奇数模数使用五轮 Newton 逆元和固定 64 轮模加倍生成参数；该路径
不执行 C 除法或余数。偶数模数不受支持并返回公开零。

## 2. REDC 正确性与范围

### 2.1 通用 unsigned 后备路径

对 `0 <= t < q*R`，unsigned 后备路径计算

```text
m = (t mod R)*q0i mod R
u = (t+m*q)/R.
```

由 `q0i*q = -1 mod R`，可得 `t+m*q = 0 mod R`，所以除以 `R` 是精确的；
同时

```text
u = t*R^(-1) mod q.
```

又因为 `m < R`，

```text
0 <= t+m*q < q*R+(R-1)*q < 2*q*R,
```

故 `0 <= u < 2*q`，一次无分支条件减 `q` 足以得到 `[0,q)` 中的 canonical
结果。实现把 `t` 分成高低 32 位计算，避免直接形成可能溢出的 `t+m*q`。

所有调用点满足输入界：

- 编码：`t=x*R2`，其中 `x<R`、`R2<q`，所以 `t<q*R`；
- 乘法：两个 Montgomery 份额均小于 `q`，所以 `t<q^2<q*R`；
- 解码：`t<q<q*R`。

因此 `encode(x)=x*R mod q`，`decode(xR)=x mod q`，而
`mont_mul(xR,yR)=xyR mod q`。

### 2.2 Signed REDC 快速路径

`mask_mont32.h` 提供 Raccoon 风格的 signed REDC，前置条件为
`3 <= q <= INT32_MAX` 且 `-q*R/2 < t < q*R/2`：

```text
m = signed32((uint32(t)*q0i) mod R),  -R/2 <= m < R/2
r = (t+m*q)/R,                      -q < r < q.
```

除以 R 仍精确，且 `r = t/R mod q`。因为 `q < 2^31`，
`|t| < q*2^31 < 2^62` 且 `|m*q| <= q*2^31 < 2^62`，
有符号 64 位乘法与求和均不溢出。负值通过先转换到 unsigned 再取高字实现，
不依赖负数右移或超范围 unsigned-to-signed 强制转换。低字乘法明确使用
`uint32_t`，其回绕由 C 标准定义。

一次 `cadd(r,q)` 即得到 canonical 结果。快速路径的调用条件：

- canonical 乘法：`t=a*b < q^2 < q*R/2`；
- 解码：`0 <= t < q`；
- 任意 32 位输入编码：仅在 `R2 <= floor(q/2)` 时走 signed 路径，
  此时 `t=x*R2 < R*q/2`；否则保留 unsigned 后备路径。

生产三个模数满足编码快速路径条件。`q > INT32_MAX` 始终使用原有的
unsigned REDC；不能直接把 signed 公式扩展到任意 32 位模数。
canonical 加减在 `q <= INT32_MAX` 时使用 32 位差值的符号位修正，
更大模数保留 64 位路径。

`mmod32_set_residue_shares()` 的输入已经由支持任意 32 位输入的 encode
完整约简，因此删除编码前重复的 Barrett 约简。加、减和 refresh 不需要
Montgomery 参数，只校验公开模数，不再构建未使用的参数。

## 3. 掩码乘法正确性

设输入重构为 `A*R` 和 `B*R`。每个单项 Montgomery 乘积为

```text
pij = REDC(ai*bj) = ai*bj*R^(-1) mod q.
```

算法先直接采样一个均匀 Montgomery 余数 `rho`，再执行

```text
z0 =  rho + p00
z1 = -rho + p01 + p10 + p11              mod q.
```

于是

```text
z0+z1
 = (a0+a1)(b0+b1)R^(-1)
 = (A*R)(B*R)R^(-1)
 = A*B*R                                     mod q.
```

输出仍在 Montgomery 域；末尾 refresh 不改变重构值。

## 4. 一阶 NI/SNI 模拟

任一内部单探针分成以下类型：

1. 随机余数 `rho` 或 `-rho`：Montgomery 编码是双射，模拟器直接采样一个
   `[0,q)` 中的均匀值；
2. `ai*bj`、`m`、REDC 的高低字或 `pij`：它们只依赖每个输入的一个份额，
   模拟器读取 `ai,bj` 后运行同一公开算法；
3. `z0` 或 `z1` 的任一 canonical 部分和：它们是均匀 `rho` 的模 `q` 平移，
   因而仍然均匀；
4. 模加法修正之前的 32/64 位中间值：模拟器采样均匀 accumulator，并只读取
   当前 `pij` 所需的 `ai,bj`，即可精确重现其分布。

所以每个内部探针至多需要每个输入的一个份额，满足 1-NI。仅观察任一输出
份额时，该份额因随机优先 accumulator 和最终独立 refresh 而均匀，模拟器
不需要输入份额，满足一阶 SNI 的输出条件。因此 `mmod32_mul()` 在上述模型和
IDEAL-RNG 前提下为条件 1-SNI。

## 5. 实现约束

- Montgomery 参数和模数必须公开；
- 每个 `pij` 必须单独 REDC，随后立即加入已随机化的 accumulator；
- 不得重新形成 `p01+p10` 或三个交叉项的未掩码临时量；
- `mmod32` 不能与普通域 `uint32_t` 混用，只有 `share/open` 执行域转换；
- 每种编译器和目标仍需检查 REDC、条件减和累加顺序的生成汇编。

特别注意：lazy signed 值仅用于独立的线性 NTT 遍历；秘密乘法的四个
`pij` 以及每一步 masked accumulator 仍逐项 canonical 化，并保留原有
value barrier。不能把 NTT 的 lazy 累加照搬到随机优先的秘密乘法 gadget。
