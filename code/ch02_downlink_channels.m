%% ch02 下行物理信道全景
%
% 目标：把 NR 下行的四个物理信道都生成一遍，看清各自的位置、开销和作用。
%   SSB    同步信号块（PSS+SSS+PBCH）—— 终端开机第一件事
%   PDCCH  物理下行控制信道，承载 DCI —— 调度的神经中枢
%   PDSCH  物理下行共享信道，承载用户数据 —— 吞吐量来源
%   CSI-RS 信道状态信息参考信号 —— 专门用于测量
%
% 版本：目标 R2021a 之后的 API。加扰/调制/预编码在新版由 nrPDSCH 一步完成；
%       若你的版本仍是分开的（nrPDSCHScramble/nrPDSCHModulate/nrPDSCHPrecode），
%       用 doc nrPDSCH 核对后再调整。
%
% 运行：直接 F5。

clear; close all; clc;

%% ---- 公共载波配置 ----
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;      % kHz
carrier.NSizeGrid = 52;              % 频域 52 个 RB
carrier.NSlot = 0;
carrier.NFrame = 0;

fprintf('%s\n', repmat('=', 1, 78));
fprintf('载波配置：SCS = %d kHz，%d RB，时隙 %d\n', ...
        carrier.SubcarrierSpacing, carrier.NSizeGrid, carrier.NSlot);
fprintf('%s\n\n', repmat('=', 1, 78));

%% ---- 1. SSB：同步信号块 ----
% SSB 占 4 个 OFDM 符号 x 20 个 RB（240 个子载波）。
% 结构：符号1=PSS，符号2=PBCH，符号3=SSS+PBCH，符号4=PBCH
% SSB 的频域位置由 AbsoluteFrequencySSB 决定；这里放在网格底部方便画图。

ssb = nrSSBConfig;
ssb.NCellID = 321;                    % 物理小区 ID，决定 PSS/SSS 序列
ssb.BlockPattern = 'Case B';          % FR1 + SCS 30kHz 对应的 pattern
ssb.NSlot = 0;
ssb.SymbolAllocation = [2 4];         % 从符号 2 开始，占 4 个符号（Case B 的规定）
ssb.NCRBSSB = 0;                      % 相对公共资源块的偏移（简化取 0）
ssb.Power = 0;

% 单独看三个组成部分
[pssSym, pssInd] = nrPSS(ssb.NCellID);
[sssSym, sssInd] = nrSSS(ssb.NCellID);

fprintf('---- SSB ----\n');
fprintf('PSS 序列长度     : %d（占 127 个子载波）\n', length(pssSym));
fprintf('SSS 序列长度     : %d\n', length(sssSym));
fprintf('小区 ID 分解     : %d = 3*%d + %d\n', ...
        ssb.NCellID, floor(ssb.NCellID/3), mod(ssb.NCellID,3));
fprintf('                  即 N_ID^(1) = %d, N_ID^(2) = %d\n', ...
        floor(ssb.NCellID/3), mod(ssb.NCellID,3));

% 生成完整的 SSB 网格
ssbGrid = nrSSB(carrier, ssb);

%% ---- 2. PDSCH：用户数据 ----
% 完整链：TB -> CRC -> LDPC -> 速率匹配 -> 加扰 -> 调制 -> 层映射 -> 预编码 -> RE 映射
% 新版 nrPDSCH 把后半段（加扰到映射）一步完成

pdsch = nrPDSCHConfig;
pdsch.NSizeBWP = 52;
pdsch.NStartBWP = 0;
pdsch.Modulation = '16QAM';
pdsch.NumLayers = 2;
pdsch.PRBSet = (0:51).';              % 占满整个 BWP
pdsch.SymbolAllocation = [0 14];      % 起始符号 0，长度 14
pdsch.MappingType = 'A';
pdsch.TargetCodeRate = 0.4785;        % 目标码率，TBS 计算要用

% DM-RS 配置：这是 PDSCH 里最影响开销的部分
pdsch.DMRS.DMRSTypeAPosition = 2;          % 第一个 DM-RS 在第 3 个符号（0-based 的符号 2）
pdsch.DMRS.DMRSConfigurationType = 1;      % Type 1：梳状，每 RB 6 个子载波
pdsch.DMRS.DMRSLength = 1;                 % 单符号（前置）/ 双符号
pdsch.DMRS.DMRSAdditionalPosition = 0;     % 附加 DM-RS 数量，0 表示只有前置
pdsch.DMRS.NumCDMGroupsWithoutData = 1;    % 多少个 CDM 组不用于数据

% TBS：先算出这一坨数据要传多少比特
tbs = nrTBS(pdsch);
fprintf('\n---- PDSCH ----\n');
fprintf('调制方式         : %s\n', pdsch.Modulation);
fprintf('层数             : %d\n', pdsch.NumLayers);
fprintf('占用 PRB         : %d\n', numel(pdsch.PRBSet));
fprintf('目标码率         : %.4f\n', pdsch.TargetCodeRate);
fprintf('TBS              : %d 比特\n', tbs);

% 生成随机传输块并走编码链
rng(2026);                             % 固定种子，保证可复现
tb = randi([0 1], tbs, 1);
dlsch = nrDLSCH;
dlsch.setTransportBlock(tb, 0);        % 第二个参数是冗余版本 rv
cw = dlsch();                          % 输出码字

fprintf('编码后码字长度   : %d 比特\n', length(cw));
fprintf('编码码率         : %.4f\n', tbs / length(cw));

% 生成 PDSCH 网格（含加扰/调制/层映射/预编码/RE映射）
pdschGrid = nrPDSCH(carrier, pdsch, cw);

% DM-RS 单独生成
dmrsGrid = nrPDSCHDMRS(carrier, pdsch);

%% ---- 3. PDCCH：控制信道 ----
% 层次：REG（1符号 x 1RB）-> CCE（6个REG）-> 聚合等级 AL（1/2/4/8/16 个 CCE）
% 终端不知道 DCI 在哪儿，所以要在"搜索空间"里盲检；搜索空间的时频范围由 CORESET 给出。

coreset = nrCORESETConfig;
coreset.FrequencyResources = [1 1 1 1 0 0];   % 6 个 RB 组（每组 6 RB）的位图
coreset.Duration = 2;                          % 占 2 个符号
coreset.CCEREGMapping = 'interleaved';
coreset.REGCBundleSize = 6;

pdcch = nrPDCCHConfig;
pdcch.NStartBWP = 0;
pdcch.NSizeBWP = 52;
pdcch.CORESET = coreset;
pdcch.AggregationLevel = 4;            % 用 4 个 CCE
pdcch.AllocatedCandidate = 1;          % 第 1 个候选

fprintf('\n---- PDCCH ----\n');
fprintf('CORESET 时长     : %d 符号\n', coreset.Duration);
fprintf('聚合等级 AL      : %d（即 %d 个 CCE，%d 个 REG）\n', ...
        pdcch.AggregationLevel, pdcch.AggregationLevel, pdcch.AggregationLevel*6);

% DCI 载荷 -> Polar 编码 -> 速率匹配 -> 加扰 -> QPSK -> 映射
dciBits = randi([0 1], 40, 1);         % 一个典型的 DCI 大小

% 生成搜索空间的所有候选位置（这就是终端要盲检的地方）
[~, ~, pdcchInfo] = nrPDCCHSpace(carrier, pdcch);
fprintf('搜索空间候选数   : %d\n', size(pdcchInfo.Candidates, 1));

pdcchGrid = nrPDCCH(carrier, pdcch, dciBits);

%% ---- 4. CSI-RS：测量用参考信号 ----
% CSI-RS 不承载信息，专门给终端测量信道质量（CQI/PMI/RI 反馈的基础）。
% 与 DM-RS 的分工：DM-RS 管解调（跟着数据走），CSI-RS 管测量（按需配置，
% 可以稀疏、可以窄带）。这是 NR 相对 LTE 的重要设计分离——LTE 的 CRS 什么都管。

csirs = nrCSIRSConfig;
csirs.CSIRSType = {'nzp','nzp'};
csirs.RowNumber = [1 2];               % 图样行号，决定 RE 排布
csirs.Density = {'one','one'};         % 密度
csirs.NumRB = carrier.NSizeGrid;
csirs.RBOffset = 0;
csirs.SymbolLocations = {12, 13};      % 放在时隙末尾
csirs.SubcarrierLocations = {0, 0};
csirs.NID = 0;

csirsGrid = nrCSIRS(carrier, csirs);
csirsInd  = nrCSIRSIndices(carrier, csirs);

fprintf('\n---- CSI-RS ----\n');
fprintf('端口数           : %d\n', csirs.NumPorts);
fprintf('占用的符号位置   : %s\n', mat2str(cell2mat(csirs.SymbolLocations)));
fprintf('占用 RE 数       : %d\n', numel(csirsInd));

%% ---- 5. 把四个信道叠加到一张网格上 ----
% 注意：真实系统中这些信道不会在同一时隙这样整齐堆叠，
% 这里是为了可视化而人为安排的。

combined = ssbGrid + pdschGrid + dmrsGrid + pdcchGrid + csirsGrid;

figure('Color','w','Position',[100 60 1000 760]);

% 分信道可视化：用不同的透明度叠加
subplot(3,2,1);
plotGrid(ssbGrid, 'SSB (PSS/SSS/PBCH)', carrier);
subplot(3,2,2);
plotGrid(pdschGrid, 'PDSCH 数据', carrier);
subplot(3,2,3);
plotGrid(dmrsGrid, 'PDSCH DM-RS', carrier);
subplot(3,2,4);
plotGrid(pdcchGrid, 'PDCCH 控制', carrier);
subplot(3,2,5);
plotGrid(csirsGrid, 'CSI-RS 测量', carrier);
subplot(3,2,6);
plotGrid(abs(combined), '全部叠加（功率）', carrier);

%% ---- 6. 开销分析：这是理解设计权衡的关键 ----
totalRE = 12 * carrier.NSizeGrid * 14;
ssbRE   = nnz(ssbGrid);
dmrsRE  = nnz(dmrsGrid);
pdcchRE = nnz(pdcchGrid);
csirsRE = nnz(csirsGrid);

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('一个时隙内的资源开销分析\n');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('总 RE 数         : %6d\n', totalRE);
fprintf('DM-RS 占用       : %6d  (%.2f%%)\n', dmrsRE,  100*dmrsRE/totalRE);
fprintf('PDCCH 占用       : %6d  (%.2f%%)\n', pdcchRE, 100*pdcchRE/totalRE);
fprintf('CSI-RS 占用      : %6d  (%.2f%%)\n', csirsRE, 100*csirsRE/totalRE);
fprintf('SSB 占用         : %6d  (%.2f%%)\n', ssbRE,   100*ssbRE/totalRE);
fprintf('%s\n', repmat('-', 1, 78));
fprintf('参考信号+控制总开销: %.2f%%\n', ...
        100*(dmrsRE+pdcchRE+csirsRE+ssbRE)/totalRE);
fprintf('%s\n', repmat('=', 1, 78));

fprintf('\n关键权衡（记住这个，比记参数重要）：\n');
fprintf('  DM-RS 越密 -> 信道估计越准 -> 但可用数据 RE 越少 -> 吞吐越低\n');
fprintf('  这个权衡在 TBS 公式里是显式的：N_RE 已经扣掉了 DM-RS 占用的 RE\n');
fprintf('  所以"参考信号密度"本质上是在"估计精度"和"频谱效率"之间买保险\n\n');

fprintf('验收自查：\n');
fprintf('  1. 把 pdsch.DMRS.DMRSAdditionalPosition 从 0 改成 3，TBS 怎么变？为什么？\n');
fprintf('  2. 把 pdsch.Modulation 从 16QAM 改成 64QAM，TBS 怎么变？\n');
fprintf('  3. CORESET 和搜索空间的分工分别是什么？\n');
fprintf('  4. 为什么 DM-RS 要前置（front-loaded）？\n\n');

%% ---- 局部函数 ----
function plotGrid(g, ttl, carrier)
    imagesc(abs(g));
    colormap(gca, hot);
    axis xy;
    xlabel('符号'); ylabel('子载波');
    title(sprintf('%s  (占 %d RE)', ttl, nnz(g)), 'FontSize', 10);
    set(gca, 'FontSize', 8);
end
