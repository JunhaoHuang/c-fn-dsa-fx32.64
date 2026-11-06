# 三个修复算子的严格一阶 NI/SNI 证明

本文件证明当前修复后的 `ring_mul()`、`ring_mul192()`、
`mi32_to_mrx160()` 和 `mmod32_mul()`。证明采用 [MODEL.md](MODEL.md)
定义的两个份额、单个 32 位 word probe 模型。

所有定理都带有 **IDEAL-RNG** 前提：每次请求得到新鲜、独立、均匀且攻击者
未知的随机量。仓库仍按要求保留固定 LFSR 作为 `DEMO ONLY` 后端；因此下面是
gadget 公式与数据流的条件性信息论证明，不是固定 LFSR 可执行实例的无条件
安全声明。

## 1. 160 位环乘 `ring_mul`

设输入为

```text
A = a0 + a1 mod 2^160
B = b0 + b1 mod 2^160.
```

修复后的算法按以下顺序执行：

```text
r  <- uniform Z/(2^160)
t0 <- r
t1 <- -r
p  <- a0*b0;  t0 <- t0+p
p  <- a0*b1;  t1 <- t1+p
p  <- a1*b0;  t1 <- t1+p
p  <- a1*b1;  t1 <- t1+p
```

每次 `p` 与每次 `t0/t1` 更新之间都有编译器值屏障；函数为 `noinline`。

### 正确性

```text
t0 + t1
 = r + a0*b0 - r + a0*b1 + a1*b0 + a1*b1
 = (a0+a1)(b0+b1)
 = A*B                                      mod 2^160.
```

### 1-NI 内部探针模拟

任意一个内部探针属于以下类别之一：

1. `r` 或 `-r`：模拟器直接均匀采样一个环元素；
2. `p = ai*bj` 或 `mul5()` 内的乘加/进位：模拟器读取 `A` 的一个份额 `ai`
   和 `B` 的一个份额 `bj`，执行相同公开算法；
3. `t0 = r + a0*b0`：模拟器读取 `a0,b0`，均匀采样 `r` 后执行相同加法；
4. `t1 = -r + C` 的任意阶段：其中 `C` 是固定的已累加交叉项和。由于 `r`
   在整个环上均匀，`-r+C` 仍在整个环上均匀。模拟器无需知道 `C`，直接均匀
   采样当前 accumulator word；若探针位于本次 `add5()` 的内部进位链，则读取
   本次项 `ai*bj` 所需的 `ai,bj`，均匀采样进入该次加法的完整 accumulator，
   再执行相同 word 加法即可得到完全相同的联合边缘分布。

所有类别都使用每个输入至多一个份额，因此满足 1-NI。

### SNI 输出条件

仅观察输出份额时，`t0=r+a0*b0` 是均匀环元素；`t1=-r+C` 也是均匀环元素。
模拟器不读取任何输入份额即可生成任一单独输出份额。因此满足一阶 SNI 的
`t_int=0,t_out=1` 条件。结合上一节，`ring_mul()` 为条件 1-SNI。

该结论把 `A`、`B` 当作两个逻辑输入，每个内部乘积最多读取 `A` 的一个份额和
`B` 的一个份额。它不能直接用于 `ring_mul(a,a)`：此时 `a0*a1` 同时读取同一
逻辑输入的两个份额。例如只看最低位，秘密位为 0 时 `a0*a1=a0`，秘密位为 1
时该乘积恒为 0。因此，快速乘法的两个输入不得是同一个 sharing。

## 2. 192 位环乘 `ring_mul192`

该函数只把环从 `Z/(2^160)` 替换为 `Z/(2^192)`，并用六个 32 位 limb 实现
相同的随机优先累加顺序。正确性等式、内部模拟器和输出零输入模拟器完全相同，
所以 `ring_mul192()` 也是条件 1-SNI。

因此，在范围/截断正确性前提下：

- `mrx160_mul()`：SNI 环乘、公开加半 ULP、SNI A2B、公开移位、SNI B2A、
  refresh 的组合，条件 1-SNI；
- `mrx160_mul_wide()`：两次 SNI A2B/B2A、SNI 192 位环乘及安全截断的组合，
  条件 1-SNI；
- `mrx160_mul_u64_hi()`：SNI 环乘后经 SNI A2B 公开投影，条件 1-SNI。

### 专用平方 `mrx160_sqr`

`mrx160_sqr()` 不把同一 sharing 直接传给快速乘法，而是生成：

```text
u = Refresh(a; r)
v = Refresh(a; s)
out = mrx160_mul(u, v)
```

其中 `r`、`s` 新鲜、独立。两个 refresh 都是 1-SNI，故在简化组合图中切断
输入 `a` 到两个乘法输入的扇出边；`u`、`v` 是同一重构值的独立 sharing。
随后调用的双输入乘法满足上一节前提，输出再由乘法自身刷新。因此
`mrx160_sqr()` 在 IDEAL-RNG 和源码级 word-probing 模型下为条件 1-SNI，且
重构值和舍入规则严格等于 `mrx160_mul(a,a)`。

这里的 SNI 结论与 Lean 范围证明正交：范围证明负责功能上不发生意外环回绕，
本证明负责单探针分布不泄漏秘密。

## 3. `mi32_to_mrx160`

旧实现把输入构造成 `(a0,0)`、`(a1,0)`，使加法结果 share 0 的最低位等于
秘密最低位。修复后直接调用 `arith_to_bool_words()`：

```text
X = (a0 XOR rx, rx)
Y = (a1 XOR ry, ry)
Z = Refresh(BooleanAdd(X,Y)).
```

### 正确性

`X`、`Y` 分别重构为 `a0`、`a1`，所以 `Z` 重构为
`a0+a1 mod 2^32`。对 `Z` 的每个 Boolean 份额独立执行公开映射 `T`：低 64 位
补零、位 64..95 放入原 32 位字、位 96..159 复制原符号位。该映射在 Boolean
分享上是线性的：

```text
T(z0) XOR T(z1) = T(z0 XOR z1).
```

随后 BCZ18 B2A 恢复为 `Z/(2^160)` 算术分享，故结果正是符号扩展后的整数乘
`2^64`，即 integral Q32.64。

### 1-SNI 模拟

1. `a0 XOR rx`、`rx`、`a1 XOR ry`、`ry` 的单探针均为均匀字，或最多需要
   一个算术输入份额；
2. masked adder 是 1-NI，而它索取的任一 Boolean 输入份额已由 `rx/ry`
   一次一密；
3. A2B 最终 refresh 使任一 Boolean 输出份额与输入无关；
4. 公开符号扩展只对每个份额独立操作；保留字、零字和符号复制字的分布均与
   秘密无关；
5. BCZ18 B2A 和最终 arithmetic refresh 均为 1-SNI。

由顺序组合，`mi32_to_mrx160()` 为条件 1-SNI。旧的 `(ai,0)` 数据流已从
函数中完全删除。

## 4. 模公开奇数 `q` 的 Montgomery 乘法 `mmod32_mul`

`mmod32` 已从普通 canonical 域改为 canonical Montgomery 份额：若普通值为
`x`，内部重构值为 `x*R mod q`，其中 `R=2^32`。`share/open` 是唯一的域转换
边界；add、sub、refresh 和 mul 都留在 Montgomery 域。

每个交叉项现在执行独立 REDC：

```text
pij = REDC(ai*bj)
rM  <- uniform Montgomery residue
z0  =  rM + p00
z1  = -rM + p01 + p10 + p11                 mod q.
```

REDC 的每个中间量只依赖 `ai,bj` 各一个份额；每个 accumulator 部分和仍是
均匀 `rM` 的平移。因此原随机优先模拟器不变，而 Montgomery 编码是模 `q`
上的双射，普通域均匀随机量编码后仍均匀。输出再执行独立 refresh，所以
`mmod32_mul()` 继续满足条件 1-SNI。

REDC 可除性、`t<q*R`、一次条件减充分性、域转换正确性和逐类探针模拟的完整
证明见 [MONTGOMERY_PROOFS.md](MONTGOMERY_PROOFS.md)。
## 5. 缩小模型与机器码约束

[counterexamples.py](counterexamples.py) 做两组检查：

1. 复现旧结构的三个区分器，确保测试确实覆盖已知失败模式；
2. 穷举一位环、1-bit A2B 和 `q=5` 模乘，验证修复后每个带掩码 accumulator
   阶段对每组固定输入份额都保持均匀。

源码还使用 `noinline` 和每阶段 value barrier，阻止编译器在 C 优化阶段重新
形成旧的交叉项总和。发布前仍必须对每个目标和编译配置检查汇编，确认：

- 随机 accumulator 在第一个跨份额项之前建立；
- 不存在 `p01+p10` 或 `p01+p10+p11` 的未掩码寄存器；
- 模乘在每个交叉项之间完成约简/屏障；
- `mi32_to_mrx160()` 调用安全 A2B，而没有 `(ai,0)` 临时分享。

这些编译约束是把源码证明落实到具体二进制所必需的附加义务。
