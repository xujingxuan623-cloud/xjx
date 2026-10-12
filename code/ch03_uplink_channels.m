%% ch03 上行物理信道全景
%
% 目标：生成上行全部四个物理信道，并做本阶段最重要的对比实验
%       —— CP-OFDM vs DFT-s-OFDM 的 PAPR 分布对比。
%
%   PRACH  物理随机接入信道 —— 终端接入网络的第一个信号
%   PUCCH  物理上行控制信道 —— 承载 UCI（HARQ-ACK/SR/CSI）
%   PUSCH  物理上行共享信道 —— 承载用户数据
%   SRS    探测参考信号     —— 基站用来测上行信道
%
% 版本：目标 R2021a 之后的 API。PUCCH/PUSCH 的属性名在不同版本有过调整，
%       若报错用 doc nrPUCCHConfig / doc nrPUSCHConfig 核对。
%
% 运行：直接 F5。PAPR 统计含随机性，已固定种子。

clear; close all; clc;

rng(2026);

%% ---- 公共载波配置 ----
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid = 52;
carrier.NSlot = 0;

fprintf('%s\n', repmat('=', 1, 78));
fprintf('上行物理信道全景  (SCS = %d kHz, %d RB)\n', ...
        carrier.SubcarrierSpacing, carrier.NSizeGrid);
fprintf('%s\n\n', repmat('=', 1, 78));

%% ---- 1. PRACH：随机接入前导 ----
% PRACH 是 NR 里唯一不走标准 OFDM 调制的信道，它有自己独立的一套参数集
% 和资源映射规则。序列用 Zadoff-Chu，好处是恒幅（低 PAPR）和理想的自相关
% 特性（便于基站做定时估计）。

prach = nrPRACHConfig;
prach.FrequencyRange = 'FR1';
prach.SubcarrierSpacing = 30;
prach.SequenceIndex = 0;          % 逻辑根序列索引
prach.PreambleIndex = 0;          % 从根序列的循环移位派生出的前导号
prach.RestrictedSet = 'UnrestrictedSet';
prach.ZeroCorrelationZone = 0;    % 决定循环移位的步长，进而决定小区半径
prach.FrequencyStart = 0;

% 生成 PRACH 符号（时域）
prachWaveform = nrPRACH(carrier, prach);

fprintf('---- PRACH ----\n');
fprintf('序列类型         : Zadoff-Chu\n');
fprintf('前导格式         : %s\n', prach.Format);
fprintf('生成的样本数     : %d\n', size(prachWaveform, 1));
fprintf('循环移位配置     : ZeroCorrelationZone = %d\n', prach.ZeroCorrelationZone);

% 自相关分析：Zadoff-Chu 的自相关应该是一个尖锐的峰值
% 这正是基站能精确估计定时提前的原因
seq = prachWaveform(:,1);
xc = xcorr(seq, 'normalized');
[~, peakIdx] = max(abs(xc));

figure('Color','w','Position',[80 60 1120 800]);

subplot(3,3,1);
plot(abs(xc), 'LineWidth', 1);
xlabel('延迟（样本）'); ylabel('|自相关|（归一化）');
title('PRACH 序列自相关：尖锐峰值 = 精确定时估计的基础', 'FontSize', 10);
grid on;

subplot(3,3,2);
plot(real(seq(1:min(300,end))), 'LineWidth', 0.8);
xlabel('样本'); ylabel('幅度');
title('PRACH 时域波形（恒幅特性）', 'FontSize', 10);
grid on;

subplot(3,3,3);
histogram(abs(seq), 20);
xlabel('|样本|'); ylabel('计数');
title('幅度分布：Zadoff-Chu 恒幅', 'FontSize', 10);
grid on;

%% ---- 2. PUSCH：用户数据 + 变换预编码对比实验 ----
% 这是本脚本的核心。同一个数据，用两种方式发：
%   (a) CP-OFDM：不做变换预编码，直接 OFDM
%   (b) DFT-s-OFDM：先做 DFT 扩频，再映射到子载波，再 IFFT
% 对比它们的 PAPR 分布 —— 这解释了为什么上行边缘用户要用 DFT-s-OFDM

pusch = nrPUSCHConfig;
pusch.Modulation = '16QAM';
pusch.NumLayers = 1;                     % DFT-s-OFDM 要求单层
pusch.PRBSet = (0:23).';                 % DFT-s-OFDM 要求连续 RB
pusch.SymbolAllocation = [0 14];
pusch.MappingType = 'A';
pusch.DMRS.DMRSConfigurationType = 1;
pusch.DMRS.DMRSTypeAPosition = 2;
pusch.DMRS.DMRSLength = 1;
pusch.DMRS.DMRSAdditionalPosition = 0;
pusch.DMRS.NumCDMGroupsWithoutData = 1;
pusch.DMRS.NID = 0;

results = struct();

for tp = [false, true]
    pusch.TransformPrecoding = tp;
    if tp
        % DFT-s-OFDM 下 DM-RS 的序列生成方式不同，需要额外参数
        pusch.DMRS.DMRSConfigurationType = 1;
        key = 'DFT_s_OFDM';
        label = 'DFT-s-OFDM (Transform Precoding = on)';
    else
        key = 'CP_OFDM';
        label = 'CP-OFDM (Transform Precoding = off)';
    end

    % 生成一个完整的 PUSCH 时隙波形。
    % 为了聚焦 PAPR 对比，数据 RE 上直接用随机 QAM 符号，
    % 不跑完整的 LDPC 编码链（那部分见 ch06）。
    dataGrid  = nrPUSCH(carrier, pusch, randSymQAM(pusch, carrier));
    dmrsGrid_ = nrPUSCHDMRS(carrier, pusch);

    slotGrid = dataGrid + dmrsGrid_;
    [wf, ~] = nrOFDMModulate(carrier, slotGrid);

    % 计算 PAPR 的 CCDF
    [paprDb, ccdf] = paprCCDF(wf);

    results.(key).label = label;
    results.(key).papr  = paprDb;
    results.(key).ccdf  = ccdf;

    fprintf('\n---- PUSCH (%s) ----\n', label);
    fprintf('调制方式         : %s\n', pusch.Modulation);
    fprintf('占用 PRB         : %d (连续)\n', numel(pusch.PRBSet));
    fprintf('平均 PAPR        : %.2f dB\n', mean(paprDb));
    fprintf('PAPR @ CCDF=1e-2 : %.2f dB\n', ...
            quantile(paprDb, 0.99));
end

% ---- 画 PAPR 对比图：这是本阶段最重要的产出 ----
subplot(3,3,[4 6]);
hold on;
plot(results.CP_OFDM.ccdf,  results.CP_OFDM.papr,  'LineWidth', 2, ...
     'Color', [0.85 0.1 0.15]);
plot(results.DFT_s_OFDM.ccdf, results.DFT_s_OFDM.papr, 'LineWidth', 2, ...
     'Color', [0.1 0.3 0.85]);
hold off;
grid on;
set(gca, 'YScale', 'log');
xlabel('PAPR (dB)'); ylabel('CCDF  P(PAPR > x)');
title('CP-OFDM vs DFT-s-OFDM 的 PAPR 分布对比', 'FontSize', 11);
legend({results.CP_OFDM.label, results.DFT_s_OFDM.label}, ...
       'Location', 'southwest', 'FontSize', 9);

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('PAPR 对比结论：\n');
fprintf('  DFT-s-OFDM 的 PAPR 明显低于 CP-OFDM（典型低 3~5 dB）。\n');
fprintf('  原因：DFT 扩频把每个数据符号的能量摊到所有子载波上，\n');
fprintf('        避免了多个子载波同相叠加产生高峰值。\n');
fprintf('  代价：必须占用连续 RB，频域调度灵活性下降。\n');
fprintf('  工程结论：小区边缘/覆盖受限 -> 用 DFT-s-OFDM；\n');
fprintf('            小区中心/追求吞吐   -> 用 CP-OFDM。\n');
fprintf('%s\n\n', repmat('=', 1, 78));

%% ---- 3. PUCCH：上行控制信息 ----
% 格式 0~4 按"短/长"和"承载比特数"两个维度划分：
%   格式 0：短，1~2 比特（HARQ-ACK / SR），可放时隙末尾做快速反馈
%   格式 1：短，1~2 比特，带时域扩频
%   格式 2：长，多比特
%   格式 3：长，多比特，DFT-s-OFDM
%   格式 4：长，多比特，带正交扩频支持多用户

fprintf('---- PUCCH 格式对比 ----\n');
fprintf('（格式由 Format 属性和 SymbolAllocation 长度共同决定）\n');
pucchFormats = {'0', '1', '2', '3', '4'};

for f = pucchFormats
    try
        pc = nrPUCCHConfig;
        pc.Modulation = 'QPSK';
        pc.PRBSet = 0;
        pc.NID = 0;

        switch f{1}
            case '0'                            % 短 PUCCH，1~2 符号
                pc.Format = 0;
                pc.NumPRBs = 1;
                pc.SymbolAllocation = [13 2];
            case '1'                            % 短 PUCCH，带时域扩频
                pc.Format = 1;
                pc.NumPRBs = 1;
                pc.SymbolAllocation = [0 14];
                pc.FrequencyHopping = 'neither';
                pc.OCCLength = 4;
                pc.OCCIndex = 0;
            case '2'                            % 长 PUCCH
                pc.Format = 2;
                pc.NumPRBs = 3;
                pc.SymbolAllocation = [0 14];
                pc.FrequencyHopping = 'neither';
            case '3'                            % 长 PUCCH，DFT-s-OFDM
                pc.Format = 3;
                pc.NumPRBs = 6;
                pc.SymbolAllocation = [0 14];
                pc.FrequencyHopping = 'neither';
            case '4'                            % 长 PUCCH，正交扩频支持多用户
                pc.Format = 4;
                pc.NumPRBs = 6;
                pc.SymbolAllocation = [0 14];
                pc.FrequencyHopping = 'neither';
                pc.OCCLength = 2;
                pc.OCCIndex = 0;
        end

        [~, pucchInfo] = nrPUCCHIndices(carrier, pc);
        fprintf('  格式 %s：占用 %4d 个 RE，%2d 个符号\n', ...
                f{1}, numel(pucchInfo), pc.SymbolAllocation(2));
    catch ME
        fprintf('  格式 %s：本版本 API 不匹配（%s）—— 用 doc nrPUCCHConfig 核对\n', ...
                f{1}, ME.identifier);
    end
end
fprintf('\n');

%% ---- 4. SRS：探测参考信号 ----
% 基站通过 SRS 估计上行信道，用于上行调度和链路自适应。
% TDD 下还可利用信道互易性反推下行预编码。
% 关键概念：梳状结构（comb）、跳频图样、资源集。

srs = nrSRSConfig;
srs.NumSRSSymbols = 1;
srs.SymbolStart = 13;
srs.NumSRSPorts = 1;
srs.SRSPeriodicity = 10;               % 周期，单位时隙
srs.SRSPeriodicityOffset = 0;
srs.FrequencyStart = 0;
srs.NumRB = carrier.NSizeGrid;
srs.Bandwidth = 4;                     % 带宽配置索引
srs.BandwidthStart = 0;
srs.CombSize = 2;                      % 梳状：每 2 个子载波放一个 SRS
srs.CyclicShift = 0;
srs.NID = 0;

try
    srsGrid = nrSRS(carrier, srs);
    srsIdx  = nrSRSIndices(carrier, srs);

    subplot(3,3,7);
    plotGridSmall(srsGrid, 'SRS（梳状结构）');

    fprintf('---- SRS ----\n');
    fprintf('占用符号         : %d\n', srs.SymbolStart);
    fprintf('梳状大小 Comb    : %d（每 %d 个子载波一个 SRS）\n', ...
            srs.CombSize, srs.CombSize);
    fprintf('占用 RE 数       : %d\n', numel(srsIdx));
    fprintf('在 %d RB 带宽中的密度: %.1f%%\n', ...
            carrier.NSizeGrid, 100*numel(srsIdx)/(12*carrier.NSizeGrid));
catch ME
    fprintf('SRS 生成失败（%s）—— 属性和本版本不匹配，用 doc nrSRSConfig 核对\n', ...
            ME.identifier);
end

%% ---- 5. 上下行参考信号对照 ----
fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('上下行参考信号对照（用途分明，别记混）\n');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-10s %-8s %-8s %-30s\n', '信号', '下行', '上行', '用途');
fprintf('%-10s %-8s %-8s %-30s\n', 'DM-RS',  '有', '有', '解调（跟着数据走）');
fprintf('%-10s %-8s %-8s %-30s\n', 'PT-RS',  '有', '有', '相位跟踪（毫米波用）');
fprintf('%-10s %-8s %-8s %-30s\n', 'CSI-RS', '有', '—',  '下行信道测量');
fprintf('%-10s %-8s %-8s %-30s\n', 'SRS',    '—',  '有', '上行信道探测');
fprintf('%s\n\n', repmat('=', 1, 78));

fprintf('验收自查：\n');
fprintf('  1. PAPR 对比图里两条曲线的差距有多少 dB？为什么 DFT-s-OFDM 更低？\n');
fprintf('  2. 为什么 DFT-s-OFDM 要求连续 RB 分配？\n');
fprintf('  3. PRACH 为什么不用标准 OFDM？Zadoff-Chu 的两个好处是什么？\n');
fprintf('  4. SRS 的梳状结构有什么好处？\n\n');

%% ---- 局部函数 ----

% 生成随机的 QAM 符号（用于 PAPR 分析，替代完整的编码链）
function sym = randSymQAM(pusch, carrier)
    nRE = numel(nrPUSCHIndices(carrier, pusch));
    switch pusch.Modulation
        case 'QPSK',   M = 4;
        case '16QAM',  M = 16;
        case '64QAM',  M = 64;
        case '256QAM', M = 256;
        otherwise,     M = 4;
    end
    bits = randi([0 1], nRE*log2(M), 1);
    sym = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
end

% 计算 PAPR 的 CCDF
function [paprDb, ccdf] = paprCCDF(waveform)
    % 按 OFDM 符号分块统计瞬时功率峰值比
    x = waveform(:);
    % 用滑动窗口计算瞬时功率（窗口取一个符号长度量级）
    winLen = 64;
    nWin = floor(length(x)/winLen);
    paprDb = zeros(nWin, 1);
    avgPow = mean(abs(x).^2);
    for k = 1:nWin
        seg = x((k-1)*winLen + (1:winLen));
        peakPow = max(abs(seg).^2);
        paprDb(k) = 10*log10(peakPow / avgPow);
    end
    paprDb = sort(paprDb(:));
    ccdf = 1 - (0:length(paprDb)-1)' / length(paprDb);
end

function plotGridSmall(g, ttl)
    imagesc(abs(g)); colormap(gca, hot); axis xy;
    xlabel('符号'); ylabel('子载波');
    title(ttl, 'FontSize', 10);
    set(gca, 'FontSize', 8);
end
