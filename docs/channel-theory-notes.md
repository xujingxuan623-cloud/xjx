# 无线信道理论推导笔记

> 配套文件：`5g-channels-learning-path.md`（学习路线）、`code/ch05` `code/ch06`（可运行脚本）
> 用途：把信道建模和估计的每一步推导写清楚，供查阅和写论文时引用

---

## 第一部分：时变多径信道的数学描述

### 1.1 为什么要用两个时间变量

一个静止的、不含多径的信道只需要一个数——增益。一个含多径但不随时间变的信道需要一条曲线——每个时延上的增益。而真实的移动通信信道两样都占：**既有多径（信号走不同路径，到达时间不同），又在变化（终端在动、环境在动）**。所以要两个变量。

定义**时变冲激响应** $h(\tau, t)$：在绝对时刻 $t$ 探测信道，施加一个冲激，在延迟 $\tau$ 之后收到响应。输出与输入的关系是卷积：

$$
y(t) = \int_{-\infty}^{\infty} h(\tau, t)\, x(t - \tau)\, d\tau
$$

两个变量的地位完全不同，这是理解信道的关键：

- $\tau$：**时延**，描述"信号走了多久"，刻画多径的展开程度。它是卷积变量。
- $t$：**绝对时间**，描述"什么时候观察的"，刻画信道的变化。它是演化变量。

如果信道不随时间变，$h(\tau,t) = h(\tau)$ 退化成普通冲激响应。如果只有单径，$h(\tau,t) = a(t)\delta(\tau)$ 退化成时变增益。

**注意**：$h(\tau,t)$ 对 $\tau$ 做傅里叶变换得到**时变频率响应** $H(f,t)$：

$$
H(f,t) = \int h(\tau, t)\, e^{-j2\pi f \tau}\, d\tau
$$

这个变换对 $\tau$ 做，不动 $t$。所以 $H(f,t)$ 表示"在时刻 $t$、频率 $f$ 上信道的复增益"。后面所有频域分析都在这个函数上展开。

### 1.2 离散多径模型

实际信道可以建模成有限条路径的叠加：

$$
h(\tau, t) = \sum_{l=1}^{L} a_l(t)\, \delta\big(\tau - \tau_l(t)\big)
$$

$a_l(t)$ 是第 $l$ 条路径的复幅度，$\tau_l(t)$ 是它的时延。每条路径的复幅度可以分解成

$$
a_l(t) = \alpha_l(t)\, e^{j\theta_l(t)}
$$

其中 $\alpha_l$ 是幅度衰减（由路径长度和散射体的反射损耗决定），$\theta_l$ 是相位（由传播距离决定，$\theta_l = -2\pi f_c \tau_l + \text{const}$）。

**关键认识**：**相位对时延极其敏感**。$\theta_l = -2\pi f_c \tau_l$，载频 $f_c$ 高达 GHz 量级，而 $\tau_l$ 的微小变化（终端移动半个波长）就能让相位转 2π。这意味着 $a_l(t)$ 在移动几个波长的距离内就会剧烈变化——这就是**快衰落**的来源。相比之下 $\alpha_l$ 变化慢得多。

这个"幅度慢变、相位快变"的分离是理解衰落统计的起点。

---

## 第二部分：WSSUS 假设与二阶统计量

$h(\tau,t)$ 是随机过程，没法精确预测，只能描述统计特性。最有效的简化是 WSSUS。

### 2.1 两个假设

**广义平稳（WSS, Wide-Sense Stationary）**：信道统计特性不随时间漂移。

$$
\mathbb{E}\big[h(\tau, t+\Delta t)\, h^*(\tau, t)\big] = R_h(\tau, \Delta t)
$$

只依赖时间差 $\Delta t$，不依赖绝对时刻 $t$。物理含义：在一段足够短的时间窗内，信道"看起来是同一种随机"。

**非相关散射（US, Uncorrelated Scattering）**：不同时延的路径互不相关。

$$
\mathbb{E}\big[h(\tau_1, t)\, h^*(\tau_2, t)\big] = 0, \quad \tau_1 \neq \tau_2
$$

物理含义：产生不同时延的散射体在空间上是分离的，各自独立。

两者合起来给出**自相关函数**：

$$
R_h(\tau_1, \tau_2; \Delta t) = \mathbb{E}\big[h(\tau_1, t+\Delta t)\, h^*(\tau_2, t)\big] = R_h(\tau_1, \Delta t)\, \delta(\tau_1 - \tau_2)
$$

**适用范围**：短时间内（几毫秒到几十毫秒）WSS 成立；多数散射环境 US 成立。**失效场景**：信道统计本身在变化（如快速穿过不同环境）、或存在强相关散射（如单一强反射体主导）。做仿真时如果假设 WSSUS 而实际场景在快速变化，结果会有偏差。

### 2.2 功率延迟谱与频率相关

令 $\Delta t = 0$，得到**功率延迟谱**（PDP, Power Delay Profile）：

$$
P(\tau) = R_h(\tau, 0) = \mathbb{E}\big[|h(\tau, t)|^2\big]
$$

它描述"能量在时延轴上怎么分布"。这是最直观的信道特征，画出来就是一条随时间延展开的功率曲线。

**频率相关函数**：考察两个相距 $\Delta f$ 的频率点上的信道增益是否相关。

$$
R_H(\Delta f) = \mathbb{E}\big[H(f + \Delta f, t)\, H^*(f, t)\big]
$$

把 $H(f,t) = \int h(\tau,t) e^{-j2\pi f\tau} d\tau$ 代进去推导：

$$
\begin{aligned}
R_H(\Delta f) &= \mathbb{E}\left[\int h(\tau_1,t) e^{-j2\pi(f+\Delta f)\tau_1} d\tau_1 \int h^*(\tau_2,t) e^{j2\pi f \tau_2} d\tau_2\right] \\
&= \iint e^{-j2\pi f \tau_1} e^{-j2\pi \Delta f \tau_1} e^{j2\pi f \tau_2}\, \underbrace{\mathbb{E}[h(\tau_1,t)h^*(\tau_2,t)]}_{P(\tau_1)\delta(\tau_1-\tau_2)}\, d\tau_1 d\tau_2 \\
&= \int P(\tau)\, e^{-j2\pi \Delta f \tau}\, d\tau
\end{aligned}
$$

**结论**：频率相关函数就是功率延迟谱的傅里叶变换。

$$
\boxed{R_H(\Delta f) = \int P(\tau)\, e^{-j2\pi \Delta f \tau}\, d\tau}
$$

这个对偶关系是全部信道分析的基础。**时延轴的展开 ↔ 频率轴的相关**：时延展得越宽，频率相关性衰减越快，即频率选择性越强。

### 2.3 相干带宽

**相干带宽** $B_c$ 定义为频率相关函数降到某个阈值（通常 0.5 或 0.9）对应的 $\Delta f$。它表示"在这个频率范围内，信道看起来是平的"。

由傅里叶对偶，$P(\tau)$ 的宽度（用时延扩展度量）与 $R_H(\Delta f)$ 的宽度成反比。定义**均方根时延扩展**：

$$
\tau_{rms} = \sqrt{\frac{\sum_l P_l \tau_l^2}{\sum_l P_l} - \left(\frac{\sum_l P_l \tau_l}{\sum_l P_l}\right)^2}
$$

常用的两个经验关系（第二个对应 0.9 相关，更严格）：

$$
B_c \approx \frac{1}{5\,\tau_{rms}}, \qquad B_c \approx \frac{1}{50\,\tau_{rms}}
$$

**判据**：若信号带宽 $B \ll B_c$，信道是**平坦衰落**（频域无选择性，时域无符号间干扰）；若 $B > B_c$，信道是**频率选择性衰落**（必须用均衡或 OFDM 对付）。

**工程含义**：$\tau_{rms}$ 是选择 CP 长度的依据。CP 必须长于最大时延扩展 $\tau_{max}$，否则 OFDM 符号之间会串扰。

### 2.4 多普勒谱与时间相关

现在看时间维度。令 $\tau_1 = \tau_2$，考察信道在时间上的相关：

$$
R_h(\Delta t) = \mathbb{E}\big[h(\tau, t+\Delta t)\, h^*(\tau, t)\big]
$$

**多普勒频移的物理来源**：终端以速度 $v$ 移动，与来波方向夹角为 $\phi$，则

$$
f_d = \frac{v}{c} f_c \cos\phi = f_{d,max}\cos\phi
$$

$f_{d,max} = v f_c / c$ 是最大多普勒频移。多径从各个方向来，所以每条路径的多普勒不同，形成一个谱。

**Clarke（Jakes）模型**：假设来波在水平面内均匀分布，得到经典多普勒谱

$$
S(f) = \frac{1}{\pi f_{d,max}\sqrt{1 - (f/f_{d,max})^2}}, \qquad |f| \le f_{d,max}
$$

这个谱在 $f = \pm f_{d,max}$ 处趋于无穷——物理含义是"来波集中在水平方向时多普勒最大"，在边缘处密度堆积。

**时间相关函数**由多普勒谱的逆傅里叶变换给出：

$$
R_H(\Delta t) = \int S(f)\, e^{j2\pi f \Delta t}\, df = J_0\big(2\pi f_{d,max} \Delta t\big)
$$

其中 $J_0(\cdot)$ 是零阶第一类贝塞尔函数。这是 Clarke 模型的经典结果。

### 2.5 相干时间

**相干时间** $T_c$ 定义为时间相关函数降到某阈值对应的 $\Delta t$。这里**必须说清用哪个相关函数、哪个阈值**，因为不同组合给出的系数能差一倍。

**判据一（最常用）：复包络相关降到 0.5。** 求 $J_0(x) = 0.5$ 的首个根（数值解 $x = 1.5211$）：

$$
J_0(2\pi f_{d,max} T_c) = 0.5 \;\Rightarrow\; 2\pi f_{d,max}T_c = 1.5211 \;\Rightarrow\; T_c = \frac{0.2421}{f_{d,max}}
$$

**判据二：功率相关（即 $|J_0|^2$）降到 0.5**，等价于 $J_0 = 0.7071$，首个根 $x = 1.1266$：

$$
T_c = \frac{1.1266}{2\pi f_{d,max}} = \frac{0.1793}{f_{d,max}} = \frac{9}{16\pi f_{d,max}}
$$

**关于常见的 $0.423/f_{d,max}$：** 这个系数在教材和论文里被大量引用，但它**不是** $J_0 = 0.5$ 的解——把 $0.423$ 代回去得 $J_0(2.658) = -0.12$，相关性已经是负的了。它对应的是更宽松的包络相关或电平交叉定义，属于工程习惯取值。引用时要说清出处，不要声称它是 $J_0=0.5$ 的结果。

三种系数并列（$\Delta t = T_c$ 时的取值）：

$$
T_c = \frac{0.2421}{f_{d,max}}\ \big(J_0{=}0.5\big), \qquad
T_c = \frac{0.1793}{f_{d,max}}\ \big(|J_0|^2{=}0.5\big), \qquad
T_c \approx \frac{0.423}{f_{d,max}}\ \big(\text{包络/电平交叉}\big)
$$

另外还有一类粗略近似 $T_c \approx 1/f_{d,max}$，只在数量级估算时用，不要写进论文的定量分析。

**判据**：若符号时长 $T_{sym} \ll T_c$，信道在符号内基本不变（**慢衰落**）；若 $T_{sym} > T_c$，信道在一个符号内就变了（**快衰落**），会产生时间选择性失真。

**工程含义**：这条判据直接决定了 DM-RS 在时域上要插多密。符号太长了，就得在同一时隙内多插几个 DM-RS 来跟踪信道变化。这就是为什么高速场景下 DM-RS 的时域密度要加大。

### 2.6 两个维度的小结

| 维度 | 度量 | 谱 | 相关函数 | 相干量 | 恶化后的后果 |
|---|---|---|---|---|---|
| 时延 $\tau$ | 时延扩展 $\tau_{rms}$ | 功率延迟谱 $P(\tau)$ | $R_H(\Delta f)$ | 相干带宽 $B_c \approx 1/(5\tau_{rms})$ | 频率选择性衰落 |
| 时间 $t$ | 多普勒扩展 $f_{d,max}$ | 多普勒谱 $S(f)$ | $R_H(\Delta t) = J_0(2\pi f_d \Delta t)$ | 相干时间 $T_c = 0.2421/f_d$（$J_0{=}0.5$） | 时间选择性衰落 |

**两者是对偶的**：时延轴的展开对应频率轴的相关，时间轴的变化对应频率轴（多普勒）的展开。这个对称性是信道理论的核心结构。

---

## 第三部分：离散时间模型与 OFDM 对角化

### 3.1 采样后的信道

以采样周期 $T_s$ 对接收信号采样，$t = nT_s$，时延以采样点数计 $\tau_l = l T_s$，得到离散模型：

$$
y[n] = \sum_{l=0}^{L-1} h_l[n]\, x[n - l]
$$

$h_l[n]$ 是第 $l$ 个抽头在时刻 $n$ 的值。写成矩阵形式，$y = Hx$，其中 $H$ 是**下三角 Toeplitz 矩阵**（因为系统是因果的）。这个矩阵形式在推导均衡器时很有用。

**关键问题**：直接对 $y[n] = \sum_l h_l[n]x[n-l]$ 做检测，需要处理符号间干扰——每个接收样本都混了前面 $L-1$ 个符号。$L$ 大时（宽带系统，时延扩展对应几十上百个采样点）计算量爆炸。这就是 OFDM 要解决的问题。

### 3.2 OFDM 如何对角化信道

OFDM 的核心操作是加 CP：把每个 OFDM 符号末尾的 $N_{cp}$ 个样本复制到开头。效果是**把线性卷积变成循环卷积**。

设信道时延扩展不超过 CP 长度（$\tau_{max} \le T_{cp}$），则去掉 CP 之后的接收信号是发送符号与信道冲激响应的**循环卷积**：

$$
y[n] = \sum_{l=0}^{L-1} h_l\, x\big[(n-l) \bmod N\big], \qquad n = 0, \ldots, N-1
$$

对循环卷积做 $N$ 点 DFT：

$$
Y[k] = \text{DFT}\{y[n]\} = H[k]\, X[k], \qquad k = 0, \ldots, N-1
$$

其中信道频域响应

$$
H[k] = \sum_{l=0}^{L-1} h_l\, e^{-j2\pi kl/N}
$$

**这就是对角化**：循环卷积矩阵可以被 DFT 矩阵对角化，所以每个子载波上的关系变成一个标量乘法 $Y[k] = H[k]X[k]$，子载波之间互不干扰。检测退化成一堆独立的标量除法（或更聪明的均衡）。

**为什么是循环卷积而不是线性卷积**：CP 的作用是把"上一个符号的拖尾"补到当前符号里，使得去掉 CP 后信号的等效操作是循环的。如果不加 CP，就是线性卷积，DFT 不能对角化，会产生子载波间干扰（ICI）和符号间干扰（ISI）。

**三个条件必须同时满足**才能实现完美对角化：

1. **时延扩展在 CP 内**：$\tau_{max} \le T_{cp}$。违反则产生 ISI 和 ICI。
2. **信道在一个 OFDM 符号内不变**：符号时长 $\ll T_c$。违反则产生 ICI（时变破坏了循环卷积结构）。
3. **收发端频率同步**：载波频偏会破坏子载波正交性，产生 ICI。

**这三条正好对应了物理信道的设计约束**：第一条决定 CP 长度和 numerology（μ 越大 CP 绝对长度越短，越不适合大时延扩展场景）；第二条决定 DM-RS 的时域密度；第三条对应同步信号和相位跟踪参考信号（PT-RS）的必要性。**把这三条记住，物理信道的很多设计细节就都有解释了。**

### 3.3 加噪后的模型

考虑噪声，每个子载波上：

$$
Y[k] = H[k] X[k] + W[k]
$$

$W[k]$ 是频域噪声，若时域噪声是方差 $\sigma_w^2$ 的复高斯白噪声，则 $W[k]$ 也是独立同分布的复高斯，方差 $\sigma_w^2$。

**接收机的任务**：已知 $Y[k]$ 和 $X$ 的调制星座（也已知导频位置上的 $X$），估计出 $H[k]$，然后恢复 $X[k]$。

---

## 第四部分：信道估计

### 4.1 问题设定

在导频（DM-RS）位置，发送符号 $X_p$ 是收发双方约定好的，接收信号

$$
\mathbf{Y}_p = \mathbf{X}_p \mathbf{h} + \mathbf{w}
$$

其中 $\mathbf{X}_p = \text{diag}(X_p[0], \ldots, X_p[N_p-1])$ 是对角矩阵（每个导频符号一个），$\mathbf{h}$ 是待估的信道向量，$\mathbf{w}$ 是噪声向量（协方差 $\sigma_w^2 \mathbf{I}$）。

目标：由 $\mathbf{Y}_p$ 和 $\mathbf{X}_p$ 估计 $\mathbf{h}$。

### 4.2 LS 估计器

最小二乘不假设任何统计先验，只最小化观测残差的平方和：

$$
\hat{\mathbf{h}}_{LS} = \arg\min_{\mathbf{h}} \|\mathbf{Y}_p - \mathbf{X}_p \mathbf{h}\|^2
$$

展开目标函数并求导。记 $J(\mathbf{h}) = (\mathbf{Y}_p - \mathbf{X}_p\mathbf{h})^H(\mathbf{Y}_p - \mathbf{X}_p\mathbf{h})$，对 $\mathbf{h}^*$ 求偏导（复数域的 Wirtinger 导数）：

$$
\frac{\partial J}{\partial \mathbf{h}^*} = -\mathbf{X}_p^H (\mathbf{Y}_p - \mathbf{X}_p \mathbf{h}) = 0
$$

解得**正规方程**：

$$
\mathbf{X}_p^H \mathbf{X}_p \mathbf{h} = \mathbf{X}_p^H \mathbf{Y}_p
$$

$$
\boxed{\hat{\mathbf{h}}_{LS} = (\mathbf{X}_p^H \mathbf{X}_p)^{-1} \mathbf{X}_p^H \mathbf{Y}_p}
$$

因为 $\mathbf{X}_p$ 是对角矩阵，这退化成逐元素的除法：

$$
\hat{H}_{LS}[k] = \frac{Y_p[k]}{X_p[k]} = H[k] + \frac{W[k]}{X_p[k]}
$$

**误差分析**：估计误差是 $\mathbf{e} = \hat{\mathbf{h}}_{LS} - \mathbf{h} = \mathbf{X}_p^{-1}\mathbf{w}$，均方误差

$$
\mathbb{E}\big[\|\mathbf{e}\|^2\big] = \sigma_w^2\, \text{tr}\big[(\mathbf{X}_p^H\mathbf{X}_p)^{-1}\big] = \sigma_w^2 \sum_{k} \frac{1}{|X_p[k]|^2}
$$

**结论**：LS 的误差正比于噪声方差，反比于导频功率。**它完全不管信道结构**——即使相邻子载波的信道高度相关（因为时延扩展有限），LS 也不利用这个信息，所以每个导频上的噪声都原封不动地进了估计结果。

**优缺点**：实现极简（一次除法），不需要任何先验知识；代价是噪声放大，且每个导频独立估计，精度受限。

### 4.3 MMSE 估计器

MMSE 利用信道的统计先验：信道 $\mathbf{h}$ 是零均值随机向量，协方差 $\mathbf{R}_{HH} = \mathbb{E}[\mathbf{h}\mathbf{h}^H]$。这个矩阵由功率延迟谱决定（做傅里叶变换得到频域相关）。

MMSE 估计器最小化均方误差 $\mathbb{E}[\|\hat{\mathbf{h}} - \mathbf{h}\|^2]$，其解是条件期望（正交投影原理）：

$$
\hat{\mathbf{h}}_{MMSE} = \mathbb{E}\big[\mathbf{h} \mathbf{Y}_p^H\big]\, \mathbb{E}\big[\mathbf{Y}_p\mathbf{Y}_p^H\big]^{-1} \mathbf{Y}_p
$$

计算两个矩：

$$
\mathbb{E}\big[\mathbf{h}\mathbf{Y}_p^H\big] = \mathbb{E}\big[\mathbf{h}(\mathbf{X}_p\mathbf{h}+\mathbf{w})^H\big] = \mathbf{R}_{HH}\mathbf{X}_p^H
$$

（用到 $\mathbf{h}$ 与 $\mathbf{w}$ 独立，且 $\mathbf{w}$ 零均值）

$$
\mathbb{E}\big[\mathbf{Y}_p\mathbf{Y}_p^H\big] = \mathbf{X}_p\mathbf{R}_{HH}\mathbf{X}_p^H + \sigma_w^2 \mathbf{I}
$$

代入：

$$
\boxed{\hat{\mathbf{h}}_{MMSE} = \mathbf{R}_{HH}\mathbf{X}_p^H \big(\mathbf{X}_p\mathbf{R}_{HH}\mathbf{X}_p^H + \sigma_w^2\mathbf{I}\big)^{-1} \mathbf{Y}_p}
$$

用矩阵求逆引理可以改写成更直观的形式，把 LS 估计作为输入：

$$
\boxed{\hat{\mathbf{h}}_{MMSE} = \mathbf{R}_{HH}\big(\mathbf{R}_{HH} + \sigma_w^2(\mathbf{X}_p^H\mathbf{X}_p)^{-1}\big)^{-1} \hat{\mathbf{h}}_{LS}}
$$

**这个形式极有启发性**——它就是一个**加权滤波器**，对 LS 估计结果做平滑：

$$
\hat{\mathbf{h}}_{MMSE} = \mathbf{W}\, \hat{\mathbf{h}}_{LS}, \qquad \mathbf{W} = \mathbf{R}_{HH}\big(\mathbf{R}_{HH} + \sigma_w^2(\mathbf{X}_p^H\mathbf{X}_p)^{-1}\big)^{-1}
$$

**理解这个权重矩阵**：

- 信噪比高时（$\sigma_w^2 \to 0$），$\mathbf{W} \to \mathbf{I}$，MMSE 退化成 LS。因为此时 LS 已经很准，不需要平滑。
- 信噪比低时（$\sigma_w^2 \to \infty$），$\mathbf{W} \to \mathbf{0}$，估计值趋于先验均值（零）。因为观测太不可靠，不如不信。
- 一般情况下，权重由信道相关结构 $\mathbf{R}_{HH}$ 决定：**相关性强的子载波之间会互相"借力"**——一个子载波上的噪声会被相邻子载波的干净估计所平滑。

**标量情形验证**：单导频时，$\hat{H}_{MMSE} = \frac{\sigma_H^2}{\sigma_H^2 + \sigma_w^2/|X_p|^2}\hat{H}_{LS}$。这就是经典的维纳滤波形式，权重是"信道功率 /（信道功率 + 噪声功率）"，即信噪比越大越信任观测。

**代价**：需要知道 $\mathbf{R}_{HH}$（信道统计，从 PDP 和多普勒得到）和噪声方差 $\sigma_w^2$。实际系统中这两个量要估计，估计不准会劣化性能（失配）。**研究点**：鲁棒 MMSE、信道统计未知时的估计方法。

### 4.4 从导频到全网格：插值

导频只出现在部分 RE 上，但数据 RE 也需要信道值。所以估计完导频位置后要插值。

- **频域插值**：线性插值、DFT 插值（利用时延扩展有限 → 频域信道是低阶多项式/少数正弦的叠加）、维纳插值。
- **时域插值**：沿符号方向插值，通常线性或基于多普勒的维纳。
- **二维插值**：频域和时域同时考虑。

**关键洞察**：插值的效果取决于导频密度是否满足采样定理。频域上，信道响应 $H(f)$ 的"带宽"由时延扩展决定（是 $P(\tau)$ 的傅里叶变换），所以频域导频间隔必须小于相干带宽；时域上，信道变化的"带宽"由多普勒决定，所以时域导频间隔必须小于相干时间。**这就是 DM-RS 密度不能太稀的根本原因。**

---

## 第五部分：均衡

### 5.1 问题设定

有了信道估计 $\hat{H}[k]$，要恢复发送符号。每个子载波上：

$$
Y[k] = H[k]X[k] + W[k]
$$

目标：找均衡权重 $g[k]$，使 $\hat{X}[k] = g[k]Y[k]$ 尽量接近 $X[k]$。

MIMO 情形下这是矩阵问题：$\mathbf{y} = \mathbf{H}\mathbf{s} + \mathbf{w}$，找 $\mathbf{G}$ 使 $\hat{\mathbf{s}} = \mathbf{G}\mathbf{y}$。

### 5.2 ZF 均衡器

迫零的目标是完全消除干扰：令 $\mathbf{G}_{ZF}\mathbf{H} = \mathbf{I}$，所以

$$
\boxed{\mathbf{G}_{ZF} = (\mathbf{H}^H\mathbf{H})^{-1}\mathbf{H}^H}
$$

（Moore-Penrose 伪逆。单天线/标量情形退化成 $\hat{X}[k] = Y[k]/\hat{H}[k]$。）

代入信号模型：

$$
\hat{\mathbf{s}}_{ZF} = \mathbf{s} + (\mathbf{H}^H\mathbf{H})^{-1}\mathbf{H}^H\mathbf{w}
$$

干扰项为零，但**噪声被放大**。标量情形的噪声方差：

$$
\text{Var}\big(\hat{X}[k]\big) = \frac{\sigma_w^2}{|\hat{H}[k]|^2}
$$

**致命问题**：深衰落子载波上 $|H[k]|$ 很小，除以一个很小的数，噪声被剧烈放大。自适应调制里这些子载波本来就要被"丢弃"（分配低码率或不用），ZF 却把它们端上来，导致误码率激增。

**优点**：不需要噪声方差，实现简单。**适用**：信噪比高且信道条件好的场景。

### 5.3 MMSE 均衡器

MMSE 均衡最小化均方误差 $\mathbb{E}[\|\mathbf{G}\mathbf{y} - \mathbf{s}\|^2]$，在干扰抑制和噪声放大之间取平衡。

假设发送符号功率归一化 $\mathbb{E}[\mathbf{ss}^H] = \mathbf{I}$，噪声 $\mathbb{E}[\mathbf{ww}^H] = \sigma_w^2\mathbf{I}$。

正交投影原理：最优 $\mathbf{G}$ 满足误差与观测正交，$\mathbb{E}[(\mathbf{G}\mathbf{y}-\mathbf{s})\mathbf{y}^H] = 0$，解得

$$
\boxed{\mathbf{G}_{MMSE} = \mathbf{H}^H \big(\mathbf{H}\mathbf{H}^H + \sigma_w^2\mathbf{I}\big)^{-1}}
$$

（这个形式对应 $\hat{\mathbf{s}} = \mathbf{G}\mathbf{y}$，$N_r \times N_t$ 的 $\mathbf{H}$。若从 $\mathbf{y} = \mathbf{H}\mathbf{s}+\mathbf{w}$ 出发且 $N_t \le N_r$。）

标量/逐子载波情形的直观形式：

$$
\boxed{\hat{X}[k] = \frac{\hat{H}^*[k]}{|\hat{H}[k]|^2 + \sigma_w^2}\, Y[k]}
$$

**对比 ZF**：ZF 是 $1/\hat{H}[k]$，MMSE 是 $\hat{H}^*[k]/(|\hat{H}[k]|^2+\sigma_w^2)$。

**看两个极限**：

- **信噪比高**（$|\hat{H}|^2 \gg \sigma_w^2$）：MMSE $\approx \hat{H}^*/|\hat{H}|^2 = 1/\hat{H}$，**退化成 ZF**。因为噪声小，放大噪声无所谓。
- **深衰落**（$|\hat{H}|^2 \ll \sigma_w^2$）：MMSE $\approx \hat{H}^*/\sigma_w^2$，**增益趋于零而非无穷**。这就是关键——MMSE 不会在深衰落子载波上疯狂放大噪声，而是"明智地放弃"这些子载波。

**残差干扰**：MMSE 不追求干扰完全消除，留下的那部分是有意的，换取噪声抑制。这就是所谓的**偏差-方差权衡**。

### 5.4 估计误差的影响

前面都假设 $\hat{H}$ 准确。实际有估计误差：$\hat{H} = H + \epsilon$。

代入 ZF：$\hat{X} = Y/\hat{H} = (HX+W)/(H+\epsilon)$。误差 $\epsilon$ 会引入额外的干扰项，且当 $\hat{H}$ 接近零时（估计误差使某些子载波的估计值虚假地接近零）问题更严重。

**这解释了导频密度的取舍**：导频越稀 → 插值误差越大 → $\epsilon$ 越大；导频越密 → 开销越大 → 有效吞吐越低。**这个权衡是 5G 参考信号设计的核心，也是"物理信道设计服务于传播信道"的最直接体现。**

---

## 第六部分：3GPP 标准信道模型

### 6.1 TDL 模型结构

抽头延迟线模型把信道建模成 $L$ 个离散抽头：

$$
h(\tau, t) = \sum_{l=1}^{L} a_l(t)\, \delta(\tau - \tau_l)
$$

每个抽头的复幅度 $a_l(t)$ 是零均值复高斯过程（Rayleigh 衰落），其时间相关由多普勒谱决定。抽头的**功率和时延是固定的**（查 38.901 的表），只有幅度随时间随机变化。

**3GPP 定义的五组标准抽头**（TR 38.901 §7.7.2）：

| 模型 | 抽头数 | 场景 | 衰落分布 |
|---|---|---|---|
| TDL-A | 23 | NLOS | Rayleigh |
| TDL-B | 23 | NLOS | Rayleigh |
| TDL-C | 24 | NLOS | Rayleigh |
| TDL-D | 13 | LOS | Rician（含直射径） |
| TDL-E | 14 | LOS | Rician |

**具体抽头的时延和功率值直接查规范表格**，不要凭记忆。表格给的是相对时延（归一化）和相对功率（dB），实际使用时按 `DelaySpread` 参数缩放。缩放规则：规范表里的时延是相对某个归一化值给出的，乘上你设定的 $\tau_{rms}$ 就得到实际时延；功率保持不变。

**LOS 模型（TDL-D/E）的 Rician 特性**：第一条路径是确定性的直射径（不是随机的），其余是散射径。Rician K 因子定义为直射径功率与散射径总功率之比，决定了衰落深度——K 越大越接近 AWGN。

### 6.2 CDL 模型结构

CDL 在 TDL 基础上加了角度维度。发射天线 $s$、接收天线 $u$ 之间的信道：

$$
H_{u,s}(t, \tau) = \sum_{n=1}^{N} \sum_{m=1}^{M} \sqrt{\frac{P_n}{M}}\, G_{u,s,n,m}(t)\, \delta(\tau - \tau_n)
$$

其中 $N$ 是簇数、$M$ 是每簇的射线数、$P_n$ 是第 $n$ 簇的功率。每簇内 $M$ 条射线的角度在簇的角扩展内随机分布。

复增益 $G_{u,s,n,m}$ 包含三项：

$$
G_{u,s,n,m}(t) = \underbrace{\sqrt{G_{RX,u}(\theta_{n,m,ZOA},\phi_{n,m,AOA})}}_{\text{接收天线增益}} \cdot \underbrace{\sqrt{G_{TX,s}(\theta_{n,m,ZOD},\phi_{n,m,AOD})}}_{\text{发射天线增益}} \cdot \underbrace{e^{j2\pi f_{d,n,m}t}}_{\text{多普勒}} \cdot e^{j\Phi_{n,m}}
$$

以及阵元间的相位差（由阵列响应向量给出）。**这个相位差就是空间相关性的来源**——不同阵元看到同一簇时的到达角不同，产生不同的相位偏移。

**为什么 CDL 计算量大**：每个簇要算 $N_r \times N_t$ 个天线对的复增益，每个又涉及 $M$ 条射线的求和。$N=24$、$M=20$、$64 \times 8$ 天线的配置下，计算量是 TDL 的几万倍。

**3GPP 定义的五组标准模型**：CDL-A/B/C（NLOS）、CDL-D/E（LOS），簇数分别为 23/23/24/13/14。

### 6.3 大尺度路损模型

**自由空间路损**（Friis 公式）：

$$
P_r = P_t \frac{G_t G_r \lambda^2}{(4\pi d)^2}
$$

写成路损（dB，$d$ 用 km、$f$ 用 MHz）：

$$
PL_{FS}(dB) = 32.44 + 20\log_{10}(d_{km}) + 20\log_{10}(f_{MHz})
$$

**推导常数 32.44**：

$$
\begin{aligned}
PL &= 20\log_{10}\frac{4\pi d f}{c} \\
&= 20\log_{10}(4\pi/c) + 20\log_{10}(d_m) + 20\log_{10}(f_{Hz})
\end{aligned}
$$

代入 $c = 3\times10^8$，且 $d_m = 1000d_{km}$、$f_{Hz} = 10^6 f_{MHz}$：

$$
\begin{aligned}
PL &= 20\log_{10}(4\pi/3\times10^8) + 20\log_{10}(1000d_{km}) + 20\log_{10}(10^6 f_{MHz}) \\
&= -147.56 + 60 + 120 + 20\log_{10}d_{km} + 20\log_{10}f_{MHz} \\
&= 32.44 + 20\log_{10}d_{km} + 20\log_{10}f_{MHz}
\end{aligned}
$$

**3GPP 场景化模型**：TR 38.901 §7.4 定义了 UMa、UMi、RMa、InH 四类场景，每类分 LOS/NLOS。典型的 LOS 路损形式（UMa）：

$$
PL = \begin{cases}
PL_1 & 10\text{ m} \le d \le d_{BP} \\
PL_2 & d_{BP} \le d \le 5\text{ km}
\end{cases}
$$

**断点距离 $d_{BP}$ 的物理来源**：路损之所以在某个距离上斜率翻倍，是因为**直射径与地面反射径的干涉**。用双射线模型可以严格推导出来：

直射径长度 $d_1 = \sqrt{d^2 + (h_{BS}-h_{UT})^2}$，地面反射径长度 $d_2 = \sqrt{d^2 + (h_{BS}+h_{UT})^2}$。当 $d \gg h$ 时做一阶展开：

$$
\Delta d = d_2 - d_1 \approx \frac{2 h_{BS} h_{UT}}{d}
$$

路径差 $\Delta d$ 带来的相位差是 $2\pi\Delta d/\lambda$。令它等于 $\pi$（即半波长，第一次破坏性干涉）：

$$
\frac{2 h_{BS} h_{UT}}{d} = \frac{\lambda}{2} \;\Rightarrow\; d_{BP} = \frac{4 h_{BS} h_{UT}}{\lambda} = \frac{4 h_{BS} h_{UT} f_c}{c}
$$

**距离小于 $d_{BP}$**：路径差小于半波长，两径近似同相，相长干涉，衰减慢（LOS 段指数约 22/10 ≈ 2.2）。**距离超过 $d_{BP}$**：开始进入相消干涉区，衰减加快（指数约 4）。这就是分段公式的来源。

**注意系数是 4 而不是 $2\pi$。** 有些教材用 $2\pi h_{BS}h_{UT}f_c/c$，那是第一菲涅尔区接触地面的判据，不是 38.901 采用的写法，两者相差 $\pi/2 \approx 1.57$ 倍。**引用 38.901 时必须用 4**，否则路损曲线在断点位置会对不上规范表格。

**代入验证**：$h_{BS} = 25$ m、$h_{UT} = 1.5$ m、$f_c = 3.5$ GHz 时，

$$
d_{BP} = \frac{4 \times 25 \times 1.5 \times 3.5\times10^9}{3\times10^8} = 1750 \text{ m}
$$

### 6.4 阴影衰落

路损给的是平均损耗，实际接收功率在其周围**对数正态**波动：

$$
PL_{total} = PL_{mean} + X_\sigma, \qquad X_\sigma \sim \mathcal{N}(0, \sigma^2)
$$

$\sigma$ 是阴影衰落标准差，典型值：UMa 4 dB（LOS）/ 6 dB（NLOS），UMi 4/7.82 dB，RMa 4/8 dB，InH 3/8.03 dB。

**两点重要性质**：

**距离相关**：同一终端移动时，相邻位置的阴影衰落是相关的（同一个遮挡物会挡住一段距离），相关距离通常几十米。建模时用一阶自回归过程或指数相关函数。

**对覆盖率的影响**：链路预算里要留阴影衰落裕量。要求边缘覆盖率 95% 时，裕量约 $1.645\sigma$；要求 99% 时约 $2.33\sigma$。$\sigma$ 的方差本身也影响中断概率，且这个影响是非线性的——**阴影衰落方差的增大对中断概率的影响比均值偏移更大**，这是做覆盖分析时容易忽略的点。

---

## 第七部分：常用公式速查

**信道相关与相干量**

$$
R_H(\Delta f) = \int P(\tau)e^{-j2\pi\Delta f\tau}d\tau, \qquad B_c \approx \frac{1}{5\tau_{rms}}
$$

$$
R_H(\Delta t) = J_0(2\pi f_d \Delta t), \qquad T_c = \frac{0.2421}{f_d}, \qquad f_d = \frac{v f_c}{c}
$$

$$
\tau_{rms} = \sqrt{\frac{\sum P_l\tau_l^2}{\sum P_l} - \left(\frac{\sum P_l\tau_l}{\sum P_l}\right)^2}
$$

**OFDM 关系**

$$
Y[k] = H[k]X[k] + W[k], \qquad H[k] = \sum_{l} h_l e^{-j2\pi kl/N}
$$

对角化条件：$\tau_{max} \le T_{cp}$ 且 $T_{sym} \ll T_c$ 且频率同步。

**估计器**

$$
\hat{H}_{LS} = \frac{Y_p}{X_p}, \qquad \text{MSE}_{LS} = \frac{\sigma_w^2}{|X_p|^2}
$$

$$
\hat{H}_{MMSE} = \mathbf{R}_{HH}\big(\mathbf{R}_{HH} + \sigma_w^2(\mathbf{X}_p^H\mathbf{X}_p)^{-1}\big)^{-1}\hat{\mathbf{h}}_{LS}
$$

**均衡器**

$$
\text{ZF：}\quad \hat{X} = \frac{Y}{\hat{H}}, \qquad \text{噪声放大} \; \frac{\sigma_w^2}{|\hat{H}|^2}
$$

$$
\text{MMSE：}\quad \hat{X} = \frac{\hat{H}^*}{|\hat{H}|^2 + \sigma_w^2}Y
$$

**路损**

$$
PL_{FS} = 32.44 + 20\log_{10}d_{km} + 20\log_{10}f_{MHz}
$$

$$
d_{BP} = \frac{4 h_{BS}h_{UT}f_c}{c}
$$

---

## 附：推导中容易出错的几个地方

**$R_H(\Delta f)$ 用的是功率延迟谱的傅里叶变换，不是信道本身的傅里叶变换。** 前者是统计量（对时间取期望后的结果），后者是具体实现。仿真里要平均多个信道实现才能得到正确的 $R_H$。

**相干带宽和相干时间的系数不是绝对的，而且相干时间的常见取值里混进了一个不自洽的数。** $1/(5\tau_{rms})$ 对应 0.5 相关，$1/(50\tau_{rms})$ 对应 0.9 相关。相干时间这边要小心：$J_0=0.5$ 严格解出 $0.2421/f_d$，$|J_0|^2=0.5$ 解出 $0.1793/f_d$，而流传很广的 $0.423/f_d$ 并不是 $J_0=0.5$ 的结果（代回去得 $J_0=-0.12$）。三个数都能在文献里找到，写论文时必须说明自己用的是哪一个、对应什么判据，否则和别人结果对不上，甚至会被审稿人指出推导错误。

**MMSE 的两个等价形式别搞混。** 一个是从 $\mathbf{Y}_p$ 直接估计的闭式（含 $\mathbf{R}_{HH}\mathbf{X}_p^H(\ldots)^{-1}$），另一个是把 LS 估计作为输入的平滑形式。两者数学上等价，但后者更直观、计算上通常也更好实现（先算 LS 再做矩阵乘法）。

**多普勒谱的 Clarke 模型假设来波水平面均匀分布。** 实际场景（有主导来波方向、有仰角）不一定满足，此时用 Jakes 谱会有偏差。3GPP 的 CDL 模型不假设 Clarke 谱，而是直接规定每条射线的角度，所以更通用。

**CP 长度和相干带宽是两个不同层面的事。** CP 对付的是**符号间干扰**（时域，要求 $\tau_{max} < T_{cp}$）；均衡对付的是**符号内失真**（频域，子载波上的幅度相位起伏）。即使 CP 足够，也仍需要均衡。这两件事经常被混淆。
