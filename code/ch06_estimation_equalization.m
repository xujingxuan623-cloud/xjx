%% ch06 信道估计、均衡与端到端链路
%
% 目标：把前面所有内容武器化 —— 完成一条能出性能曲线的完整链路。
%   1. LS 与 MMSE 估计器的实现与对比
%   2. ZF 与 MMSE 均衡器的实现与对比
%   3. 端到端 PDSCH 链路：比特 -> 波形 -> 信道 -> 估计 -> 均衡 -> 解码 -> BLER
%
% 这是本套脚本的重点。做完这个，你就具备了做链路级仿真的完整能力。
%
% 版本：目标 R2021a 之后的 API。
%   新版的 nrPDSCHDecode 内部集成了信道估计与均衡，可直接吃网格。
%   老版需要自己调 nrChannelEstimate + nrEqualizeMMSE。
%   本脚本采用"手动估计 + 手动均衡"的写法，因为它教学价值更高
%   （你能看清每一步），而且各版本兼容性更好。
%
% 运行：直接 F5。BLER 曲线需要几分钟。

clear; close all; clc;
rng(2026);

%% ---- 参数配置 ----
fprintf('%s\n', repmat('=', 1, 80));
fprintf('信道估计、均衡与端到端链路仿真\n');
fprintf('%s\n\n', repmat('=', 1, 80));

carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid = 52;
carrier.NCellID = 321;

pdsch = nrPDSCHConfig;
pdsch.NSizeBWP = 52;
pdsch.NStartBWP = 0;
pdsch.Modulation = '16QAM';
pdsch.NumLayers = 1;
pdsch.PRBSet = (0:51).';
pdsch.SymbolAllocation = [0 14];
pdsch.MappingType = 'A';
pdsch.TargetCodeRate = 0.4785;
pdsch.DMRS.DMRSTypeAPosition = 2;
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 1;
pdsch.DMRS.DMRSAdditionalPosition = 0;
pdsch.DMRS.NumCDMGroupsWithoutData = 1;
pdsch.DMRS.NID = 0;

tbs = nrTBS(pdsch);
fprintf('配置：16QAM，%d RB，TBS = %d 比特\n\n', numel(pdsch.PRBSet), tbs);

%% ==================== 第一部分：LS vs MMSE 估计器 ====================
% 这是理论笔记第 4.2/4.3 节的实现。

fprintf('%s\n', repmat('=', 1, 80));
fprintf('第一部分：LS 与 MMSE 信道估计器对比\n');
fprintf('%s\n', repmat('-', 1, 80));

% 造一个稀疏的导频图样（模拟 DM-RS）：每 4 个子载波一个导频
Nfft = 512;
nPilot = Nfft/4;
pilotIdx = 1:4:Nfft;
pilotIdx = pilotIdx(1:nPilot);

% 真实的信道：一个 4 抽头的频率选择性信道
trueDelays = [0 3 7 12];                  % 采样点
trueGains  = [1 0.6 0.35 0.2];            % 幅度
trueDelays = trueDelays + 1;              % MATLAB 索引从 1 开始
hTrue = zeros(Nfft, 1);
hTrue(trueDelays) = trueGains .* exp(1j*2*pi*rand(1,4));
Htrue = fft(hTrue);

% 导频符号（QPSK，单位功率）
pilotSym = exp(1j*2*pi*(randi([0 3], nPilot, 1))/4);

SNRs = -5:2.5:25;
mseLS = zeros(size(SNRs));
mseMMSE = zeros(size(SNRs));
nReal = 200;                              % 每个 SNR 点的信道实现数

for si = 1:numel(SNRs)
    snrDb = SNRs(si);
    noiseVar = 10^(-snrDb/10);
    accLS = 0; accMMSE = 0;

    for r = 1:nReal
        % 生成导频位置的接收信号
        noise = sqrt(noiseVar/2)*(randn(nPilot,1) + 1j*randn(nPilot,1));
        Yp = Htrue(pilotIdx) .* pilotSym + noise;

        % --- LS 估计：直接相除 ---
        H_LS_pilot = Yp ./ pilotSym;

        % --- MMSE 估计 ---
        % 需要信道相关矩阵 R_HH。用真实的 PDP 构造（实际中要估计）。
        % 简化做法：用一个基于时延扩展的 sinc 型相关模型
        tauGrid = (0:max(trueDelays)-1)' / Nfft;
        % 频域相关矩阵（只算导频位置，避免大矩阵）
        np = length(pilotIdx);
        R = zeros(np, np);
        for a = 1:np
            for b = 1:np
                df = (pilotIdx(a) - pilotIdx(b)) / Nfft;
                R(a,b) = sum(trueGains.^2 .* exp(1j*2*pi*df*(trueDelays-1)));
            end
        end
        R = R / mean(diag(R));            % 归一化

        % MMSE 估计：H_mmse = R (R + sigma^2 (X^H X)^-1)^-1 H_ls
        % 导频单位功率，(X^H X)^-1 = sigma^2 * I
        W = R / (R + noiseVar*eye(np));
        H_MMSE_pilot = W * H_LS_pilot;

        % 计算导频位置的估计误差（归一化 MSE）
        accLS   = accLS   + mean(abs(H_LS_pilot   - Htrue(pilotIdx)).^2);
        accMMSE = accMMSE + mean(abs(H_MMSE_pilot - Htrue(pilotIdx)).^2);
    end

    mseLS(si)   = accLS   / nReal;
    mseMMSE(si) = accMMSE / nReal;
end

figure('Color','w','Position',[80 60 1180 760]);

subplot(2,3,1);
semilogy(SNRs, mseLS, 'o-', 'LineWidth', 1.8, 'Color', [0.85 0.1 0.15]);
hold on;
semilogy(SNRs, mseMMSE, 's-', 'LineWidth', 1.8, 'Color', [0.1 0.3 0.85]);
% 理论参考：LS 的 MSE = sigma^2 / |X|^2 = 10^(-SNR/10)
semilogy(SNRs, 10.^(-SNRs/10), 'k--', 'LineWidth', 1.4);
hold off; grid on;
xlabel('SNR (dB)'); ylabel('归一化 MSE');
title('LS vs MMSE 估计误差', 'FontSize', 11);
legend({'LS (实测)', 'MMSE (实测)', 'LS 理论值 \sigma^2/|X|^2'}, ...
       'Location', 'southwest', 'FontSize', 8);

fprintf('LS 估计误差的理论值 = sigma^2 / |X_p|^2（导频单位功率时为噪声方差）\n');
fprintf('MMSE 利用信道相关性在各导频间"借力"，误差更低\n');
fprintf('在低 SNR 处 MMSE 优势最明显，高 SNR 处两者趋于接近\n\n');

%% ---- 插值：从导频到全网格 ----
% 导频稀疏，数据 RE 上的信道值要靠插值。
% 关键：插值效果取决于导频间隔是否满足采样定理。

pilotPos = pilotIdx;
allPos = (1:Nfft)';
H_LS_full_linear = interp1(pilotPos, H_LS_pilot, allPos, 'linear', 'extrap');

% DFT 插值：利用"时延扩展有限 -> 频域信道是少数正弦叠加"这一先验
% 步骤：导频 -> 时域（补零）-> 频域
H_LS_pilot_grid = zeros(Nfft, 1);
H_LS_pilot_grid(pilotPos) = H_LS_pilot;
h_est = ifft(H_LS_pilot_grid) * Nfft / nPilot;
% 只保留前 L 个抽头（利用信道稀疏性）
L_keep = 16;
h_est_trunc = zeros(size(h_est));
h_est_trunc(1:L_keep) = h_est(1:L_keep);
H_DFT_interp = fft(h_est_trunc) / Nfft * nPilot;

subplot(2,3,2);
plot(1:Nfft, 20*log10(abs(Htrue)), 'k', 'LineWidth', 1.6); hold on;
plot(1:Nfft, 20*log10(abs(H_LS_full_linear)), '--', ...
     'Color', [0.85 0.1 0.15], 'LineWidth', 1.2);
plot(1:Nfft, 20*log10(abs(H_DFT_interp)), ':', ...
     'Color', [0.1 0.3 0.85], 'LineWidth', 1.6);
hold off; grid on;
xlabel('子载波索引'); ylabel('|H| (dB)');
title('插值方法对比（低 SNR 下）', 'FontSize', 10);
legend({'真实信道', '线性插值', 'DFT 插值'}, 'Location', 'southwest', 'FontSize', 8);

errLinear = mean(abs(H_LS_full_linear - Htrue).^2);
errDFT    = mean(abs(H_DFT_interp - Htrue).^2);
fprintf('插值误差对比（同一组噪声实现）：\n');
fprintf('  线性插值 MSE : %.4f\n', errLinear);
fprintf('  DFT 插值 MSE : %.4f\n', errDFT);
fprintf('  DFT 插值利用了"信道在时域稀疏"的先验，通常优于线性插值\n\n');

%% ==================== 第二部分：ZF vs MMSE 均衡 ====================
% 理论笔记第 5.2/5.3 节的实现

fprintf('%s\n', repmat('=', 1, 80));
fprintf('第二部分：ZF 与 MMSE 均衡器对比\n');
fprintf('%s\n', repmat('-', 1, 80));

M = 16;
nSym = 20000;
txData = qammod(randi([0 M-1], nSym, 1), M, 'UnitAveragePower', true);

% 一个深衰落信道：某些子载波上 |H| 很小
Hdeep = (randn(nSym,1) + 1j*randn(nSym,1))/sqrt(2);
% 人为制造几个深零点
zeroIdx = randi(nSym, 200, 1);
Hdeep(zeroIdx) = 0.01 * Hdeep(zeroIdx);
Hdeep = Hdeep / sqrt(mean(abs(Hdeep).^2));

subplot(2,3,3);
plot(20*log10(abs(Hdeep)), '.', 'MarkerSize', 3, 'Color', [0.3 0.3 0.3]);
grid on; ylim([-45 12]);
xlabel('符号索引'); ylabel('|H| (dB)');
title('含深衰落的信道频响', 'FontSize', 10);

eqSNRs = 0:5:30;
serZF = zeros(size(eqSNRs));
serMMSE = zeros(size(eqSNRs));

for si = 1:numel(eqSNRs)
    nv = 10^(-eqSNRs(si)/10);
    errZ = 0; errM = 0;
    for trial = 1:5
        noise = sqrt(nv/2)*(randn(nSym,1) + 1j*randn(nSym,1));
        y = Hdeep .* txData + noise;

        % --- ZF：直接除以 H ---
        xZF = y ./ Hdeep;

        % --- MMSE：H*/(|H|^2 + sigma^2) ---
        xMMSE = conj(Hdeep) ./ (abs(Hdeep).^2 + nv) .* y;

        dZF = qamdemod(xZF, M, 'UnitAveragePower', true);
        dMM = qamdemod(xMMSE, M, 'UnitAveragePower', true);
        dTx = qamdemod(txData, M, 'UnitAveragePower', true);

        errZ = errZ + mean(dZF ~= dTx);
        errM = errM + mean(dMM ~= dTx);
    end
    serZF(si) = errZ/5;
    serMMSE(si) = errM/5;
end

subplot(2,3,4);
semilogy(eqSNRs, max(serZF, 1e-6), 'o-', 'LineWidth', 1.8, 'Color', [0.85 0.1 0.15]);
hold on;
semilogy(eqSNRs, max(serMMSE, 1e-6), 's-', 'LineWidth', 1.8, 'Color', [0.1 0.3 0.85]);
hold off; grid on;
xlabel('SNR (dB)'); ylabel('符号错误率 (SER)');
title('ZF vs MMSE 均衡', 'FontSize', 11);
legend({'ZF', 'MMSE'}, 'Location', 'southwest', 'FontSize', 9);

fprintf('ZF：x = y/H，噪声被放大为 sigma^2/|H|^2\n');
fprintf('    深衰落子载波(|H|很小)上噪声被剧烈放大 -> 这些符号基本全错\n');
fprintf('MMSE：x = H*/(|H|^2+sigma^2) * y\n');
fprintf('    深衰落时增益趋于 H*/sigma^2 而非无穷大 -> "明智地放弃"这些子载波\n');
fprintf('    这就是 MMSE 在低 SNR 和深衰落下明显更好的原因\n\n');

%% ==================== 第三部分：端到端 PDSCH 链路 ====================

fprintf('%s\n', repmat('=', 1, 80));
fprintf('第三部分：端到端 PDSCH 链路（这是最终目标）\n');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('链路：TB -> LDPC编码 -> 加扰调制 -> 映射 -> OFDM -> [TDL信道] ->\n');
fprintf('      OFDM解调 -> 信道估计 -> 均衡 -> 解调 -> LDPC解码 -> CRC -> 统计\n\n');

% 信道配置
fs = 30.72e6;                          % 30 kHz SCS 的典型采样率
tdl = nrTDLChannel;
tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 100e-9;
tdl.CarrierFrequency = 3.5e9;
tdl.MaximumDopplerShift = 20;
tdl.SampleRate = fs;
tdl.NumTransmitAntennas = 1;
tdl.NumReceiveAntennas = 1;

% DM-RS 的位置和取值：接收端用它做估计
dmrsIdx = nrPDSCHDMRSIndices(carrier, pdsch);
dmrsRefGrid = nrPDSCHDMRS(carrier, pdsch);
dmrsRef = dmrsRefGrid(dmrsIdx);        % 收发双方约好的参考符号

% DM-RS 落在哪些符号列上（用列号去分组，逐列做频域插值）
[subRows, symCols] = ind2sub(size(dmrsRefGrid), dmrsIdx);
dmrsCols = unique(symCols);
fprintf('DM-RS 占用 %d 个 RE，分布在 %d 个符号列：%s\n', ...
        numel(dmrsIdx), numel(dmrsCols), mat2str(dmrsCols(:).'));
fprintf('其余 %d 列的信道值靠插值得到\n\n', 14 - numel(dmrsCols));

% BLER 扫描配置
snrVec = -6:2:14;
nTxPerPoint = 60;                      % 注意：真实论文需要远大于此（几千次）
bler = nan(size(snrVec));
throughput = nan(size(snrVec));

fprintf('开始扫描（每个 SNR 点 %d 个传输块）...\n', nTxPerPoint);
fprintf('提示：真实论文的 BLER 曲线每点需几千次传输，这里为演示只跑 %d 次，\n', nTxPerPoint);
fprintf('      所以曲线的低 BLER 段会因样本不足而抖动，这本身就是个反面教材\n\n');

tic;
for si = 1:numel(snrVec)
    snrDb = snrVec(si);
    nErr = 0;
    nBitsOk = 0;

    for tb = 1:nTxPerPoint
        % ---- 发送端 ----
        tbBits = randi([0 1], tbs, 1);
        dlsch = nrDLSCH;
        dlsch.setTransportBlock(tbBits, 0);
        cw = dlsch();
        txGrid = nrPDSCH(carrier, pdsch, cw) + dmrsRefGrid;

        [txWave, ~] = nrOFDMModulate(carrier, txGrid);

        % ---- 传播信道 ----
        tdl.reset();
        [fadedWave, ~] = tdl(txWave);

        % 加噪：按目标 SNR 定噪声功率。
        % 注意：严格的做法是按每个 RE 上的信噪比来定义 SNR（因为 TDL 有
        % 功率起伏），这里按波形总功率定义，是常见简化，写论文时要说明。
        sigPow = mean(abs(fadedWave).^2);
        noisePow = sigPow * 10^(-snrDb/10);
        noise = sqrt(noisePow/2)*(randn(size(fadedWave)) + 1j*randn(size(fadedWave)));
        rxWave = fadedWave + noise;

        % ---- 接收端 ----
        rxGrid = nrOFDMDemodulate(carrier, rxWave);

        % 步骤 1：DM-RS 位置的 LS 估计 —— H = Y / X
        H_ls = rxGrid(dmrsIdx) ./ dmrsRef;

        % 步骤 2：插值到全网格。
        % 做法：按符号列分组做频域线性插值；没有 DM-RS 的符号列
        % 沿用最近的 DM-RS 列的估计（本仿真多普勒很低，一个时隙内
        % 信道基本不变，这个简化成立；高速场景必须做时域插值）。
        HestGrid = zeros(size(rxGrid));
        allRows = (1:size(rxGrid,1)).';
        for c = 1:size(rxGrid,2)
            if ismember(c, dmrsCols)
                % 本列有 DM-RS：用本列的估计量插值
                mask = (symCols == c);
                rows_c = subRows(mask);
                vals_c = H_ls(mask);
                [rows_c, sortIdx] = sort(rows_c);
                vals_c = vals_c(sortIdx);
                HestGrid(:, c) = interp1(rows_c, vals_c, allRows, ...
                                         'linear', 'extrap');
            else
                % 本列没有 DM-RS：先用最近的 DM-RS 列（后续可改为时域插值）
                [~, nearest] = min(abs(dmrsCols - c));
                rows_c = subRows(symCols == dmrsCols(nearest));
                vals_c = H_ls(symCols == dmrsCols(nearest));
                [rows_c, sortIdx] = sort(rows_c);
                vals_c = vals_c(sortIdx);
                HestGrid(:, c) = interp1(rows_c, vals_c, allRows, ...
                                         'linear', 'extrap');
            end
        end

        % 步骤 3：MMSE 均衡（逐 RE）
        % 用估计信道的平均功率作为信号功率的代理，噪声功率用设定值
        noiseVarEst = noisePow / mean(abs(txWave).^2);   % 归一化到数据 RE 上
        eqGrid = conj(HestGrid) ./ (abs(HestGrid).^2 + noiseVarEst) .* rxGrid;

        % 步骤 4：解调 + 解扰 + LDPC 解码 + CRC 校验
        % 注意：nrPDSCHDecode 的签名随版本变化。
        %   R2021a 之前：[cw, sym] = nrPDSCHDecode(carrier, pdsch, rxGrid, ...)
        %   R2021a 之后：cw = nrPDSCHDecode(carrier, pdsch, rxGrid)
        % 这里用新版；若报"输出参数过多"就去掉两边括号只留一个输出。
        cwHat = nrPDSCHDecode(carrier, pdsch, eqGrid);
        if iscell(cwHat), cwHat = cwHat{1}; end

        dlschDec = nrDLSCHDecoder;
        dlschDec.TransportBlockLength = tbs;
        dlschDec.TargetCodeRate = pdsch.TargetCodeRate;
        [~, crcErr] = dlschDec(cwHat, 0);

        if crcErr ~= 0
            nErr = nErr + 1;
        else
            nBitsOk = nBitsOk + tbs;
        end
    end

    bler(si) = nErr / nTxPerPoint;
    % 吞吐量：成功传输的比特数 / 所用时间。一个时隙 0.5 ms（30 kHz SCS）。
    throughput(si) = nBitsOk / nTxPerPoint / 0.5e-3;

    fprintf('  SNR = %+5.0f dB : BLER = %.4f\n', snrDb, bler(si));
end
elapsed = toc;

subplot(2,3,5);
semilogy(snrVec, max(bler, 1e-4), 'o-', 'LineWidth', 2, ...
         'Color', [0.1 0.3 0.85], 'MarkerFaceColor', [0.1 0.3 0.85]);
grid on; ylim([1e-4 1.2]);
xlabel('SNR (dB)'); ylabel('BLER');
title('端到端 BLER 曲线', 'FontSize', 11);

subplot(2,3,6);
plot(snrVec, throughput/1e6, 's-', 'LineWidth', 2, ...
     'Color', [0.1 0.6 0.3], 'MarkerFaceColor', [0.1 0.6 0.3]);
grid on;
xlabel('SNR (dB)'); ylabel('吞吐量 (Mbps)');
title('端到端吞吐量曲线', 'FontSize', 11);

fprintf('\n扫描完成，耗时 %.1f 秒\n', elapsed);

%% ---- 完整链路的核对点 ----
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('端到端链路的关键核对点（做仿真时逐条查）：\n');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('  1. 采样率对齐：tdl.SampleRate 必须等于 OFDM 波形采样率，\n');
fprintf('     否则信道的行为没有物理意义（这是最常见的低级错误）\n');
fprintf('  2. 采样量：BLER < 1e-3 时每点至少几千次传输，否则曲线抖动会\n');
fprintf('     掩盖真实的性能差异，得出错误结论\n');
fprintf('  3. 信道实现数：CDL/TDL 是随机的，每个 SNR 点要在多个信道实现上\n');
fprintf('     平均，否则结果不具统计意义\n');
fprintf('  4. 随机种子：固定 rng，记录全部参数，否则结果无法复现\n');
fprintf('     （无法复现是审稿意见里最常见的批评之一）\n');
fprintf('  5. 均衡器输入：确认均衡用的 H 是从 DM-RS 估计来的，不是理想值\n');
fprintf('     （本脚本为演示流程用了理想信道，真实仿真必须用估计值）\n');
fprintf('%s\n', repmat('=', 1, 80));

fprintf('\n建议的下一步：\n');
fprintf('  打开官方示例 "NR PDSCH Throughput"，对照本脚本的流程，\n');
fprintf('  看它如何用 nrChannelEstimate + nrEqualizeMMSE 做真实的估计与均衡。\n');
fprintf('  那是一个经过验证的正确实现，比任何教程都可靠。\n\n');

fprintf('验收自查：\n');
fprintf('  1. LS 和 MMSE 的 MSE 曲线在低 SNR 处差距多大？为什么？\n');
fprintf('  2. 为什么 MMSE 均衡在深衰落子载波上不放大噪声？\n');
fprintf('  3. BLER 曲线在哪个 SNR 处弯折？这个弯折点由什么决定？\n');
fprintf('  4. 把 DelaySpread 从 100ns 改到 1000ns，BLER 曲线怎么变？为什么？\n');
