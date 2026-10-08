# MATLAB 5G Toolbox 系统学习路线

> 面向对象：MATLAB 熟练、5G NR 理论零基础、目标为科研仿真
> 文档用途：按阶段推进的学习地图。每个阶段给出「理论补什么 → 工具箱怎么用 → 官方示例做什么 → 验收标准」四段式，方便自查是否真的掌握。

---

## 0. 开始之前：环境与前置准备

### 版本与工具箱

5G Toolbox 的 API 在 R2021a 有过一次大规模改版。旧写法是 `nrCarrierConfig` 配合一大堆独立函数和 `nrXxxConfig` 对象，新写法把上下行链路各自的参数收进统一的 `nrULSCHConfig` / `nrDLSCHConfig` 以及 `nrPUSCHConfig` / `nrPDSCHConfig` 这类配置对象。网上大量教程和博客还是旧版写法，遇到函数不存在时先确认版本。

建议用 R2022b 之后的版本，理由是新版示例更完整、对 Release 17 特性的覆盖更好。运行以下命令确认手上有什么：

```matlab
ver('nr')              % 5G Toolbox 版本
ver('comm')            % Communications Toolbox（很多例子依赖它）
```

需要留意的依赖关系：通信工具箱提供调制、信道编码基础模块和信道模型；DSP System Toolbox 提供滤波器和重采样（`nrOFDMModulate` 内部会用到）；Phased Array System Toolbox 在多天线和波束赋形相关示例里会用到；RF Blockset / RF Toolbox 在功率放大器和射频损伤建模时需要。做纯链路级仿真，前三个基本够。

### 规范文档：理论补课的主线

工具箱本质上是 3GPP 38 系列规范的可执行实现。读代码时对照规范，比看二手教程效率高得多。常用几本：

| 规范 | 内容 |
|---|---|
| TS 38.211 | 物理信道与调制：帧结构、numerology、参考信号序列、OFDM 基带信号生成 |
| TS 38.212 | 复用与信道编码：LDPC、Polar、CRC、速率匹配、交织 |
| TS 38.213 | 物理层过程：功率控制、链路自适应、波束管理、HARQ 时序 |
| TS 38.214 | 物理层数据传输：MCS 表、TBS 计算、资源分配、CSI 反馈 |
| TS 38.101 | 用户设备无线发射与接收指标 |
| TS 38.104 | 基站无线发射与接收指标 |

38.211 和 38.212 是必须精读的两本，前者管"信号长什么样"，后者管"比特怎么变成码字"。

### 检查官方示例库

```matlab
doc 5G Toolbox
```

打开后重点看 Examples 和 Functions 两栏。所有官方示例都可以直接在 MATLAB 里以实时脚本形式打开并逐段执行，这是最快的学习材料，比任何第三方教程都准。每个示例顶部都有 "Introduction" 段落，说清了这个例子在仿真什么场景。

---

## 阶段一：打好 NR 物理层的三块地基

**目标**：搞清 NR 和 LTE 在物理层上的关键差异，能用工具箱生成参数正确、结构正确的空载波形。

### 1.1 Numerology、帧结构与资源网格

NR 最核心的变化是取消了固定子载波间隔，改用参数集（numerology）μ 来定义。子载波间隔 Δf = 2^μ × 15 kHz，μ 取 0 到 4，对应 15/30/60/120/240 kHz。μ 越大，符号越短、时隙越短、支持的带宽越大，但覆盖越差、相位噪声越敏感。这是 FR1（Sub-6 GHz）和 FR2（毫米波）用不同 numerology 的根本原因。

配套要理解的概念：

- **帧结构**：无线帧 10 ms，子帧 1 ms 固定不变，时隙长度随 μ 变化（1 ms 里放 2^μ 个时隙）。
- **时隙与符号**：常规 CP 下一个时隙 14 个 OFDM 符号；扩展 CP 只在 μ=2 时有，一个时隙 12 个符号。
- **资源网格**：1 个资源元素（RE）是 1 个子载波 × 1 个符号；1 个资源块（RB）是 12 个子载波。频域上的分配以 RB 为粒度。
- **带宽部分（BWP）**：载波的一部分，UE 不需要处理整个载波带宽，这是 NR 降低终端复杂度的重要设计。
- **TDD 与灵活的时隙配置**：一个时隙内的符号可以是下行、上行或灵活（flexible），这让 NR 能适应各种上下行配比。

工具箱对应：

```matlab
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;   % kHz，即 mu = 1
carrier.NSizeGrid = 52;           % 频域 RB 数
carrier.NStartGrid = 0;           % 相对公共资源块的偏移
carrier.CyclicPrefix = 'normal';
carrier.NSlot = 0;                % 时隙号
carrier.NFrame = 0;               % 无线帧号
disp(carrier);                    % 查看由以上参数派生出的所有量
```

`disp(carrier)` 会打印出采样率、每时隙样本数、时隙长度等信息，把这个对象和规范里的公式逐个对上，是第一周最有价值的练习。

**验收**：能够手算任意 μ 下 1 个时隙的样本数和时长，并用 `carrier` 对象验证；能说清为什么 FR2 要用 μ=2 或 3。

### 1.2 OFDM 波形生成与资源网格映射

工具箱把 OFDM 调制拆成三层，理解这个分层是读懂后续所有示例的前提：

- `nrCarrierConfig` 描述载波，不含任何数据；
- 各个 `nrXxxConfig`（如 `nrPDSCHConfig`）描述某个信道占哪些 RE；
- `nrXxx` 函数（如 `nrPDSCH`）把比特映射到这些 RE，返回一个网格；
- `nrOFDMModulate` / `nrOFDMDemodulate` 在网格和时域波形之间转换。

```matlab
carrier = nrCarrierConfig('SubcarrierSpacing', 15, 'NSizeGrid', 52);
grid = nrResourceGrid(carrier);          % 空网格：12*NSizeGrid 行 x 每时隙符号数 列
[waveform, info] = nrOFDMModulate(carrier, grid);
```

**验收**：理解网格的行是子载波、列是 OFDM 符号；知道 `info` 里返回的采样率和 CP 长度怎么用；能解释为什么 NR 的 CP 开销在 μ 增大时反而变小（因为 CP 绝对时长随符号缩短，但比例大致固定）。

### 1.3 参考信号

参考信号是 NR 里最容易混淆的部分，因为它按用途分了好几类，而且上下行不对称。需要建立的分类：

下行有 DM-RS（解调，跟着数据走）、PT-RS（相位跟踪，抑制相位噪声，毫米波场景重要）、CSI-RS（信道状态信息测量）、SSB（同步信号块，含 PSS/SSS/PBCH）。上行有 DM-RS、PT-RS、SRS（探测参考信号，基站用来做上行信道估计）。LTE 里的 CRS（小区专用参考信号）在 NR 里被取消了，这是重要变化——NR 用 DM-RS 承担解调，用 CSI-RS 承担测量，职责分离。

工具箱里每一类都有 `nrXxx` 生成函数和对应的 `nrXxxConfig` 或属性。序列本身在 38.211 里有闭式定义（Zadoff-Chu、Gold 序列等），工具箱用 `nrPSS`、`nrSSS`、`nrPRBS` 这类函数直接实现。

**验收**：能画出 SSB 的时频位置（4 个符号、20 个 RB 宽、在频域上的候选位置），能说出 DM-RS 的前置符号（front-loaded）设计为什么有利于低时延解码。

---

## 阶段二：上下行业务信道——从比特到波形

这是工具箱使用频率最高、也最能体现科研价值的阶段。上下行的结构是对称的，学一个就能推另一个。

### 2.1 下行共享信道（PDSCH）

完整链路是：传输块 → CRC 附着 → LDPC 编码 → 速率匹配 → 加扰 → 调制 → 层映射 → 预编码 → RE 映射。

```matlab
carrier = nrCarrierConfig('SubcarrierSpacing', 30, 'NSizeGrid', 52);

pdsch = nrPDSCHConfig;
pdsch.Modulation = '16QAM';
pdsch.NumLayers = 2;
pdsch.PRBSet = 0:51;                 % 占满整个 BWP
pdsch.SymbolAllocation = [0 14];      % 起始符号 0，长度 14 个符号
pdsch.DMRS.DMRSTypeAPosition = 2;     % 第一个 DM-RS 在第 2 个符号
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.NumCDMGroupsWithoutData = 1;

% 传输块大小：R2021a 之后可以直接从配置对象算
pdsch.TargetCodeRate = 0.5;           % 目标码率，TBS 计算要用
tbs = nrTBS(pdsch);

% 生成比特并走完编码链
tb = randi([0 1], tbs, 1);
dlsch = nrDLSCH;                       % 下行共享信道编码链对象
dlsch.setTransportBlock(tb, 0);        % 第二个参数是冗余版本 rv
cw = dlsch();                          % 输出码字（已含 LDPC 编码、速率匹配）

% 加扰 + 调制 + 层映射 + 预编码 + RE 映射，一步完成，返回映射好的网格
grid = nrResourceGrid(carrier, pdsch.NumLayers);
grid = nrPDSCH(carrier, pdsch, cw) + nrPDSCHDMRS(carrier, pdsch);

% 生成波形
[waveform, info] = nrOFDMModulate(carrier, grid);
```

注：不同版本里 `nrTBS` 和 `nrDLSCH` 的调用形式略有差别，老版本 `nrTBS` 支持显式传参 `nrTBS(mod, nLayers, nPRB, nDMRS, xOverhead)`。以本机 `doc nrTBS`、`doc nrDLSCH` 为准。

几个必须搞清的点：

- **TBS 怎么算的**。`nrTBS` 背后是 38.214 里那张量化的表：先算中间量 N_info，按码率量化后查表得到近似的 TBS。这是链路自适应里最关键的一步，读一下 `nrTBS` 的源码（`edit nrTBS`）比看任何解释都清楚。
- **LDPC 基图选择**。BG1 和 BG2 的切换条件是码率和 TBS 阈值，工具箱在 `nrDLSCH` 里自动处理。
- **速率匹配**。把编码后的比特塞进实际可用的 RE 数里，涉及比特选择和交织。
- **DM-RS 的开销**。DM-RS 占掉的 RE 不能传数据，`NumCDMGroupsWithoutData` 这类参数直接影响可用 RE 数，进而影响有效码率。

**验收**：能独立跑通"随机比特 → 波形"的完整下行链路；能改动 MCS、层数、PRB 数后说清 TBS 和有效码率怎么变。

### 2.2 上行共享信道（PUSCH）

结构与下行高度对称，但有三处关键差异，正是考试和论文里常考的：

- **变换预编码（transform precoding）**。开启后 PUSCH 变成 DFT-s-OFDM，PAPR 显著降低，适合上行功率受限的场景。这个开关是 `nrPUSCHConfig` 的 `TransformPrecoding` 属性。什么时候开：小区边缘、覆盖受限。
- **PUSCH 的 DM-RS 类型与映射方式**。Type A 和 Type B 的时域位置规则不同，跟 PUSCH 的起始符号有关。
- **跳频**。上行支持时隙内和时隙间跳频，对覆盖和干扰随机化有意义。

```matlab
pusch = nrPUSCHConfig;
pusch.TransformPrecoding = true;   % 打开 DFT-s-OFDM
pusch.Modulation = 'QPSK';
pusch.NumLayers = 1;
pusch.PRBSet = (0:10)';
pusch.SymbolAllocation = [0 14];
```

**验收**：能分别生成 CP-OFDM 和 DFT-s-OFDM 的上行波形，用 `comm.CCDF` 或自己写 PAPR 统计函数比较两者的 PAPR 分布——这个对比图是理解变换预编码价值的最直观方式。

### 2.3 控制信道与广播信道

- **PDCCH**：承载 DCI（下行控制信息）。核心概念是 CCE（控制信道单元）、聚合等级（aggregation level）、搜索空间（search space）和 CORESET。工具栏函数是 `nrPDCCH`、`nrPDCCHConfig`、`nrCORESETConfig`、`nrSearchSpaceConfig`。DCI 的 Polar 编码由 `nrDCIEncode` 实现。
- **PUCCH**：上行控制信息，格式 0 到 4 各有适用场景（短 PUCCH 承载 1-2 比特，长 PUCCH 承载更多）。函数是 `nrPUCCH`、`nrPUCCH0` 到 `nrPUCCH4`。
- **PBCH 与 SSB**：`nrPBCH`、`nrSSB`。PBCH 用 Polar 编码，且携带的 MIB 信息需要理解小区搜索流程。
- **PRACH**：随机接入前导，用 Zadoff-Chu 序列。`nrPRACH`、`nrPRACHConfig`。科研做随机接入、定时提前估计时会用到。

**验收**：能生成一个完整的 SSB 并在时频图上可视化；能说清 CORESET 和搜索空间的关系。

---

## 阶段三：接收端与端到端链路仿真

发出去只是上半场，科研仿真真正的工作量在接收端。

### 3.1 完整链路闭环

官方有一个极好的起点示例，叫 `NR PDSCH Throughput`（在示例库里搜 "NR PDSCH Throughput"）。它把发送、信道、接收、误块率统计全部串起来，包括接收端的信道估计（用 DM-RS）、均衡（MMSE 或 ZF）、解调、解扰、LDPC 解码、CRC 校验。先把这一个例子吃透，再自己去改参数做实验。

接收端主要环节与函数：

| 环节 | 函数 |
|---|---|
| 信道估计（基于 DM-RS） | `nrChannelEstimate` |
| 解调 / 解扰 / 解层映射 / 均衡 | `nrPDSCHDecode`（新版把这几步合并，均衡方式可选 MMSE 或 ZF） |
| MMSE 均衡（单独调用） | `nrEqualizeMMSE` |
| LDPC 解码与 CRC 校验 | `nrDLSCHDecoder` / `nrULSCHDecoder` |
| CRC 编码与校验 | `nrCRCEncode` / `nrCRCDecode` |

注：老版本里加扰调制是分开的（`nrPDSCHScramble`、`nrPDSCHModulate`、`nrPDSCHPrecode` 等），新版被 `nrPDSCH` 和 `nrPDSCHDecode` 取代。如果你本机仍是老版本，这些函数依然存在。

### 3.2 信道模型

做链路级仿真必须用标准信道模型，否则结果无法和文献对比。工具箱提供：

- **`nrCDLChannel`**：3GPP TR 38.901 定义的分簇延迟线模型，参数包括延迟扩展、角度扩展、多普勒、阵列几何。这是最常用的。
- **`nrTDLChannel`**：抽头延迟线模型，更简单，适合快速验证和理论分析。
- **`nrPathLossConfig`** + `nrPathLoss`：路径损耗计算，系统级仿真用。

CDL 模型支持的场景（A/B/C 用于 NLOS 和 LOS 不同组合，D/E 用于毫米波）需要对照 38.901 的表来选。配置示例：

```matlab
channel = nrCDLChannel;
channel.DelayProfile = 'CDL-A';
channel.DelaySpread = 30e-9;
channel.CarrierFrequency = 3.5e9;
channel.MaximumDopplerShift = 100;
channel.SampleRate = info.SampleRate;   % 必须与实际波形采样率一致
rxWaveform = channel(txWaveform);
```

**验收**：能画出所配置信道的时域冲激响应和频率响应，能解释延迟扩展和多普勒如何影响相干带宽与相干时间，进而影响子载波间隔的选择——这条因果链是 NR numerology 设计的理论核心，也是论文里最好写的动机。

### 3.3 性能评估方法论

科研仿真与"跑通例子"的区别在于：能不能给出可发表的曲线，并解释曲线为什么长这样。

- **误块率（BLER）对 SNR 曲线**：通常扫 SNR 或噪声功率，绘 log 尺度 BLER。要点是样本量要足够——BLER 低于 1e-3 时每点至少几千个传输块，否则曲线抖动会掩盖差异。
- **吞吐量曲线**：链路自适应的实际效果体现在吞吐量随 SNR 的爬升曲线上，注意峰值吞吐受调度和码率上限限制。
- **链路与系统级的桥接**：链路级得到的是 BLER 曲线，系统级需要的是链路到系统映射（L2S mapping）。如果你的课题涉及系统级，要额外学这块。
- **可复现性**：固定随机种子（`rng`），记录所有参数，把配置对象序列化保存。仿真结果无法复现是审稿意见里最常见的批评之一。

---

## 阶段四：按研究方向深入

到这里基础已经打通，接下来按你的具体课题选方向。常见的几个：

**方向 A：多天线与大规模 MIMO**。学预编码矩阵指示（PMI）、码本、波束赋形、CSI 反馈。工具箱有 `nrCSIReportConfig`、`nrCodebook` 相关函数，配合 Phased Array System Toolbox 可以做波束图。示例搜 "Massive MIMO Hybrid Beamforming"。

**方向 B：高可靠低时延（URLLC）**。核心是缩短传输时长、提高可靠性的技术：更短的调度单元、重复传输（repetition）、更低的码率、更稳健的调制。研究点常在有限块长下的编码性能。

**方向 C：毫米波与相位噪声**。相位噪声建模、PT-RS 的作用、波束管理流程。RF 损伤建模需要 RF Toolbox。这是当前热点之一。

**方向 D：非授权频谱（NR-U）、车联网（NR-V2X）、NTN 卫星**。这些是 Release 16/17 的扩展场景，工具箱新版本有对应示例，具体支持程度以本机 doc 为准。

**方向 E：与 AI/ML 结合**。Release 18 起 3GPP 在讨论 AI 用于 CSI 反馈、波束管理、定位。科研上常用工具箱生成数据集，再用深度学习工具箱训练。这条路的好处是工具链都在 MathWorks 生态内。

---

## 学习方法建议

**读源码是本工具箱最高效的学习方式**。5G Toolbox 的很多函数是开放的 `.m` 文件，直接 `edit nrTBS`、`edit nrPDSCHIndices` 就能看到实现。规范条文读着抽象，对着实现看立刻明白。这也是工具箱相比"自己从零写"的最大价值——它给了你一份正确的参考实现。

**从示例出发而不是从函数列表出发**。工具箱有几百个函数，按文档目录线性学完不现实。更好的路径是：找一个接近你课题的官方示例 → 跑通 → 逐行读懂 → 改参数做实验 → 替换掉其中的模块做自己的研究。这个循环走几轮，比系统学完所有函数有效得多。

**建立自己的对照笔记**。每学一个模块，记三样东西：规范位置（哪本规范哪一节）、工具箱函数、一个可运行的验证例子。三样对上了才算真懂。

**善用可视化**。NR 的很多东西画出来就懂了：资源网格热力图、SSB 时频位置、CDL 信道的功率延迟谱、BLER 曲线。科研论文里绝大多数图都是这几种，早点把绘图模板做出来。

**注意版本陷阱**。R2021a 前后 API 差异大，网上教程大量是旧的。遇到函数报错先查 `ver('nr')` 和本机 `doc`，不要盲目搜解决方案。

---

## 附录 A：核心函数速查

**载波与网格**

```matlab
nrCarrierConfig      % 载波配置对象
nrResourceGrid       % 生成空资源网格
nrOFDMModulate       % 网格 → 时域波形
nrOFDMDemodulate     % 时域波形 → 网格
```

**下行**

```matlab
nrPDSCHConfig / nrPDSCH / nrPDSCHIndices / nrPDSCHDMRS / nrPDSCHDecode
nrDLSCH / nrDLSCHDecoder        % 下行共享信道编码/解码链
nrPDCCH / nrPDCCHConfig / nrPDCCHIndices / nrPDCCHDecode
nrCORESETConfig / nrSearchSpaceConfig
nrPBCH / nrSSB / nrPSS / nrSSS
```

**上行**

```matlab
nrPUSCHConfig / nrPUSCH / nrPUSCHIndices / nrPUSCHDMRS / nrPUSCHDecode
nrPUSCHPTRS / nrPUSCHPTRSIndices
nrULSCH / nrULSCHDecoder
nrPUCCH / nrPUCCH0 ~ nrPUCCH4
nrPRACH / nrPRACHConfig
nrSRS / nrSRSConfig
```

**编解码与工具**

```matlab
nrTBS                % 传输块大小计算
nrLDPCEncode / nrLDPCDecode
nrPolarEncode / nrPolarDecode
nrCRCEncode / nrCRCDecode
nrPRBS               % 伪随机序列
nrChannelEstimate    % 基于 DM-RS 的信道估计
nrEqualizeMMSE       % MMSE 均衡
nrLayerMap / nrLayerDemap
```

**信道模型**

```matlab
nrCDLChannel         % 分簇延迟线（38.901）
nrTDLChannel         % 抽头延迟线
nrPathLoss / nrPathLossConfig
```

## 附录 B：推荐的官方示例学习顺序

以下示例均可在 MATLAB 帮助浏览器中搜索标题找到，也建议用 `openExample` 命令直接打开。标题记得按你本机版本的实际名称核对。

1. **NR Downlink Carrier Waveform Generation** —— 最基础的波形生成，从空载波开始
2. **NR SSB Beam Sweeping** —— SSB 与波束扫描，理解同步信号
3. **NR PDSCH Resource Allocation and DM-RS** —— 资源分配与 DM-RS 开销
4. **NR PDSCH Throughput** —— 完整下行链路，重点吃透这一个
5. **NR PUSCH Throughput** —— 完整上行链路，注意变换预编码的对比
6. **NR PDCCH and CORESET Configuration** —— 控制信道
7. **NR Channel Estimation and Equalization** —— 接收端处理
8. **NR Cell Search and MIB and SIB1 Recovery** —— 小区搜索全流程
9. **NR CDL Channel Model** 相关示例 —— 信道建模
10. **NR PDSCH Throughput with CDL Channel** —— 端到端性能评估

## 附录 C：理论与工具的对照表

| 理论概念 | 规范 | 工具箱入口 |
|---|---|---|
| 参数集与帧结构 | 38.211 §4 | `nrCarrierConfig` |
| OFDM 基带信号 | 38.211 §5.3 | `nrOFDMModulate` |
| 物理资源映射 | 38.211 §7.3 | `nrPDSCHIndices` 等 |
| DM-RS / PT-RS | 38.211 §7.4 | `nrPDSCHDMRS`、`nrPUSCHPTRS` |
| CSI-RS / SRS | 38.211 §7.4 | `nrCSIRS`、`nrSRS` |
| SSB 与 PBCH | 38.211 §7.4.3 | `nrSSB`、`nrPBCH` |
| LDPC 编码 | 38.212 §5.3 / §7.2 | `nrDLSCH`、`nrLDPCEncode` |
| Polar 编码 | 38.212 §5.4 / §7.3 | `nrPolarEncode`、`nrDCIEncode` |
| TBS 与 MCS | 38.214 §5.1.3 | `nrTBS` |
| 信道模型 | TR 38.901 | `nrCDLChannel`、`nrTDLChannel` |

---

## 下一步

路线里的阶段一到阶段二大约对应两到四周的投入（每天一两小时）。如果你希望我把某个阶段展开成可直接运行的脚本，或者针对你的具体课题（告诉我课题方向就行）定制一份更聚焦的路线，随时说。
