%% ch04 同步与参考信号
%
% 目标：把下行三类参考信号（DM-RS / CSI-RS / PT-RS）和同步信号讲清楚，
%       并动手验证小区搜索的第一步——PSS 相关检测。
%
% 本脚本回答三个问题：
%   1. 终端怎么"找到"小区？（PSS/SSS 相关）
%   2. 三类参考信号分别在对付传播信道的什么特性？
%   3. 毫米波为什么必须加 PT-RS？
%
% 版本：目标 R2021a 之后的 API。PT-RS 相关属性名各版本有差异，
%       用 doc nrPDSCHPTRS / doc nrPUSCHPTRS 核对。
%
% 运行：直接 F5。

clear; close all; clc;
rng(2026);

%% ---- 1. 小区搜索：PSS 检测原理 ----
% 终端开机时什么都不知道：不知道有没有基站、不知道频点在哪儿、
% 不知道小区 ID、不知道帧边界。第一件事是用 PSS 做相关检测。
% PSS 只有 3 种（对应 N_ID^(2) = 0,1,2），且是频域 Zadoff-Chu 序列，
% 自相关极尖锐、互相关很小 —— 这是它能用于同步的原因。

NCellIDs = [0, 1, 2, 8, 321, 999];   % 覆盖 N_ID^(2) 的三种取值

figure('Color','w','Position',[80 60 1150 780]);

subplot(2,3,1);
hold on;
colors = lines(3);
for k = 0:2
    [pssSym, ~] = nrPSS(k);
    stem(0:length(pssSym)-1, abs(pssSym), '.', 'Color', colors(k+1,:), ...
         'MarkerSize', 4);
end
hold off; grid on; xlim([0 127]);
xlabel('子载波索引'); ylabel('|幅值|');
title('3 种 PSS 序列（长度 127 的频域 Zadoff-Chu）', 'FontSize', 10);
legend({'N_{ID}^{(2)}=0', 'N_{ID}^{(2)}=1', 'N_{ID}^{(2)}=2'}, ...
       'FontSize', 8, 'Location', 'southeast');

% 自相关与互相关：这是检测性能的关键
subplot(2,3,2);
[pss0, ~] = nrPSS(0);
[pss1, ~] = nrPSS(1);
[pss2, ~] = nrPSS(2);

% 转到时域看相关特性
w = hamming(length(pss0));
pss0t = ifft(ifftshift(pss0.*w));
pss1t = ifft(ifftshift(pss1.*w));

ac = abs(xcorr(pss0t, pss0t, 'normalized'));
cc = abs(xcorr(pss0t, pss1t, 'normalized'));
lag = -(length(pss0t)-1):(length(pss0t)-1);

plot(lag, ac, 'LineWidth', 1.6, 'Color', [0.1 0.3 0.85]); hold on;
plot(lag, cc, 'LineWidth', 1.6, 'Color', [0.85 0.1 0.15]); hold off;
grid on;
xlabel('延迟（样本）'); ylabel('归一化相关');
title('PSS 自相关（蓝）vs 互相关（红）', 'FontSize', 10);
legend({'自相关', '互相关'}, 'FontSize', 8);

subplot(2,3,3);
% SSS 用于确定完整的物理小区 ID（N_ID = 3*N_ID^(1) + N_ID^(2)）
[sssSym, ~] = nrSSS(321);
stem(0:length(sssSym)-1, abs(sssSym), '.', 'Color', [0.85 0.4 0.1], ...
     'MarkerSize', 4);
grid on; xlim([0 127]);
xlabel('子载波索引'); ylabel('|幅值|');
title('SSS 序列（N_{CellID} = 321）', 'FontSize', 10);

%% ---- 2. 模拟一次小区搜索 ----
% 造一个含 SSB 的接收信号，然后做相关，看能不能检测出小区 ID

testNCellID = 321;
ssb = nrSSBConfig;
ssb.NCellID = testNCellID;
ssb.BlockPattern = 'Case B';
ssb.NSlot = 0;
ssb.SymbolAllocation = [2 4];
ssb.NCRBSSB = 0;

carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid = 52;
carrier.NSlot = 0;

ssbGrid = nrSSB(carrier, ssb);
% 加一点噪声模拟真实的接收环境
noiseLevel = 0.05;
rxGrid = ssbGrid + noiseLevel*(randn(size(ssbGrid)) + 1j*randn(size(ssbGrid)))/sqrt(2);

% 用三种 PSS 分别去相关，看哪个峰值最高
subplot(2,3,4);
hold on;
detectedN2 = -1;
bestPeak = -inf;
peakVals = zeros(1,3);

% 从接收网格中提取 SSB 所在的符号段
ssbStartSym = ssb.SymbolAllocation(1) + 1;
ssbSyms = rxGrid(:, ssbStartSym:ssbStartSym+3);

% PSS 占 127 个子载波，居中放在 SSB 的 240 个子载波里。
% 因此它在 SSB 频段内的起始偏移是 floor((240-127)/2) = 56。
% 这里用网格行索引定位这段，而不是硬编码一个凭空的数字。
ssbRows = find(any(rxGrid ~= 0, 2));          % SSB 占用的行（240 个子载波）
% 若网格里还有其他信道，用 SSB 的标准宽度从底部对齐
if numel(ssbRows) >= 240
    ssbRows = ssbRows(1:240);
end
pssOffset = floor((240 - 127)/2);
pssRowStart = ssbRows(1) + pssOffset;
pssRowEnd   = pssRowStart + 126;

for n2 = 0:2
    [refPSS, ~] = nrPSS(n2);
    rxFirstSym = ssbSyms(:, 1);               % PSS 在第 1 个符号上
    corrVal = abs(sum(rxFirstSym(pssRowStart:pssRowEnd) .* conj(refPSS)));
    peakVals(n2+1) = corrVal;
    if corrVal > bestPeak
        bestPeak = corrVal;
        detectedN2 = n2;
    end
end

bar(0:2, peakVals, 'FaceColor', [0.3 0.6 0.85]);
grid on;
xlabel('假设的 N_{ID}^{(2)}'); ylabel('相关峰强度');
title(sprintf('检测结果：N_{ID}^{(2)} = %d（真值 %d）', ...
              detectedN2, mod(testNCellID,3)), 'FontSize', 10);

subplot(2,3,5);
% 完整的小区 ID 恢复：N_ID = 3*N1 + N2
% 实际流程是先检测 PSS 得 N2，再去解 SSS 得 N1
[sssRef, ~] = nrSSS(testNCellID);
N1 = floor(testNCellID/3);
N2 = mod(testNCellID,3);
text(0.05, 0.75, sprintf('物理小区 ID = %d', testNCellID), 'FontSize', 13);
text(0.05, 0.55, sprintf('N_{ID}^{(1)} = floor(%d/3) = %d', testNCellID, N1), ...
     'FontSize', 11);
text(0.05, 0.40, sprintf('N_{ID}^{(2)} = mod(%d,3)  = %d', testNCellID, N2), ...
     'FontSize', 11);
text(0.05, 0.25, sprintf('验证：3*%d + %d = %d', N1, N2, 3*N1+N2), ...
     'FontSize', 11, 'Color', [0.1 0.5 0.1]);
text(0.05, 0.08, 'PSS 定 N_{ID}^{(2)}（3 种），SSS 定 N_{ID}^{(1)}（336 种）', ...
     'FontSize', 9, 'Color', [0.4 0.4 0.4]);
axis off;
title('小区 ID 的分解与恢复', 'FontSize', 10);

%% ---- 3. 三类参考信号的时频密度对比 ----
% 这是本脚本的核心认知：参考信号的密度直接对应它要对付的传播信道特性

carrier2 = nrCarrierConfig;
carrier2.SubcarrierSpacing = 30;
carrier2.NSizeGrid = 52;

pdsch = nrPDSCHConfig;
pdsch.NSizeBWP = 52;
pdsch.Modulation = '16QAM';
pdsch.NumLayers = 1;
pdsch.PRBSet = (0:51).';
pdsch.SymbolAllocation = [0 14];

% --- DM-RS：时域前置 + 频域梳状 ---
pdsch.DMRS.DMRSTypeAPosition = 2;
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 1;
pdsch.DMRS.DMRSAdditionalPosition = 0;
dmrsGrid = nrPDSCHDMRS(carrier2, pdsch);

% --- 加一个附加 DM-RS，看时域密度怎么变 ---
pdsch.DMRS.DMRSAdditionalPosition = 3;
dmrsGridDense = nrPDSCHDMRS(carrier2, pdsch);

subplot(2,3,6);
hold on;
% 用柱状图对比 RE 占用
reCounts = [nnz(dmrsGrid), nnz(dmrsGridDense)];
b = bar(1:2, reCounts, 'FaceColor', [0.3 0.6 0.85]);
grid on;
set(gca, 'XTick', 1:2, 'XTickLabel', ...
    {'AdditionalPosition=0', 'AdditionalPosition=3'});
ylabel('DM-RS 占用的 RE 数');
title('DM-RS 时域密度对比', 'FontSize', 10);
% 标出百分比
for k = 1:2
    totalRE = 12*52*14;
    text(k, reCounts(k)*1.03, sprintf('%.1f%%', 100*reCounts(k)/totalRE), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

fprintf('%s\n', repmat('=', 1, 78));
fprintf('参考信号密度与传播信道的对应关系\n');
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-10s %-14s %-16s %-24s\n', '信号', '时域密度', '频域密度', '对应的信道特性');
fprintf('%-10s %-14s %-16s %-24s\n', 'DM-RS',  '前置+可附加', '梳状(每2或3子载波)', '时变性 + 频率选择性');
fprintf('%-10s %-14s %-16s %-24s\n', 'PT-RS',  '密集(每符号)', '稀疏(每2或4 RB)', '相位噪声(本振不稳)');
fprintf('%-10s %-14s %-16s %-24s\n', 'CSI-RS', '按需配置', '灵活', '信道质量测量(非解调)');
fprintf('%-10s %-14s %-16s %-24s\n', 'SSB',    '固定位置', '20 RB', '初始接入(无先验知识)');
fprintf('%s\n\n', repmat('=', 1, 78));

%% ---- 4. 为什么毫米波需要 PT-RS ----
% 相位噪声来自本振的不稳定，导致星座图整体旋转（公共相位误差 CPE）。
% 频率越高相位噪声越严重，所以毫米波必须用 PT-RS 来跟踪和补偿。

figure('Color','w','Position',[100 100 900 340]);

% 用一个简单的相位噪声模型演示
M = 16;
nSym = 3000;
txSym = qammod(randi([0 M-1], nSym, 1), M, 'UnitAveragePower', true);

% 低频：相位噪声小
pnLow = 0.01 * cumsum(randn(nSym,1));
% 高频：相位噪声大
pnHigh = 0.08 * cumsum(randn(nSym,1));

rxLow  = txSym .* exp(1j*pnLow);
rxHigh = txSym .* exp(1j*pnHigh);

subplot(1,3,1);
plot(real(txSym), imag(txSym), '.', 'MarkerSize', 3, 'Color', [0.2 0.2 0.2]);
axis square; grid on; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
title('发送 16QAM 星座', 'FontSize', 10);
xlabel('I'); ylabel('Q');

subplot(1,3,2);
plot(real(rxLow), imag(rxLow), '.', 'MarkerSize', 3, 'Color', [0.1 0.3 0.85]);
axis square; grid on; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
title('低频段：相位噪声小', 'FontSize', 10);
xlabel('I'); ylabel('Q');

subplot(1,3,3);
plot(real(rxHigh), imag(rxHigh), '.', 'MarkerSize', 3, 'Color', [0.85 0.1 0.15]);
axis square; grid on; xlim([-1.6 1.6]); ylim([-1.6 1.6]);
title('毫米波：相位噪声大 → 星座旋转', 'FontSize', 10);
xlabel('I'); ylabel('Q');

fprintf('相位噪声与 PT-RS 的逻辑链：\n');
fprintf('  载频越高 -> 本振的相对稳定性越差 -> 公共相位误差(CPE)越大\n');
fprintf('  -> 星座图整体旋转 -> 解调失败\n');
fprintf('  -> 必须插入已知的 PT-RS，用时域密集的参考符号估计并补偿这个相位\n');
fprintf('  注意：PT-RS 对付的是"相位旋转"，DM-RS 对付的是"信道衰落"，两者不可互替\n\n');

fprintf('验收自查：\n');
fprintf('  1. PSS 只有 3 种、SSS 有 336 种，合起来小区 ID 有多少种？\n');
fprintf('  2. 为什么 PSS 用 Zadoff-Chu 而不是随机序列？\n');
fprintf('  3. DM-RS 的"前置"设计对时延有什么好处？\n');
fprintf('  4. PT-RS 和 DM-RS 的分工区别是什么？\n\n');
