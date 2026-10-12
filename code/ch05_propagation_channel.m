%% ch05 传播信道建模
%
% 目标：把"传播信道"这一层从抽象公式变成可以看见、可以量化的对象。
%   1. 大尺度：路损曲线、断点距离、阴影衰落
%   2. 小尺度：TDL 信道的冲激响应、频率响应、功率延迟谱
%   3. CDL 信道与空间维度
%   4. 亲自测量相干带宽与相干时间，并和理论值对比
%
% 理论推导见 docs/channel-theory-notes.md
%
% 版本：目标 R2021a 之后的 API。nrCDLChannel 的阵列属性名各版本有差异，
%       用 doc nrCDLChannel 核对。CDL 部分若报错可直接跳过，不影响其他小节。
%
% 运行：直接 F5。含随机信道实现，已固定种子保证可复现。

clear; close all; clc;
rng(2026);

fprintf('%s\n', repmat('=', 1, 80));
fprintf('传播信道建模\n');
fprintf('%s\n\n', repmat('=', 1, 80));

%% ==================== 第一部分：大尺度 ====================

%% ---- 1. 路损曲线：四种场景 ----
% 3GPP TR 38.901 定义了 UMa/UMi/RMa/InH 四类场景，每类分 LOS/NLOS。
% 关键认知：LOS 和 NLOS 的路损曲线不是简单的平移关系，
% 而是斜率和分段点都不同 —— NLOS 几乎全程指数更高。
%
% 断点距离 d_BP = 4*hBS*hUT*fc/c 的物理来源：
% 直射径与地面反射径的路径差 Δd ≈ 2*hBS*hUT/d；当 Δd 达到半波长时
% 两径从"同相叠加"转入"反相抵消"，衰减斜率翻倍。
% 令 2*hBS*hUT/d = lambda/2 即得 d = 4*hBS*hUT*fc/c。
% 注意系数是 4 不是 2*pi —— 后者是第一菲涅尔区判据，与 38.901 相差 1.57 倍。

fc = 3.5e9;                  % 载频 3.5 GHz
hBS = 25;                    % 基站高度（UMa 标准值）
hUT = 1.5;                   % 终端高度
d = logspace(log10(10), log10(5000), 200);   % 10 m ~ 5 km

d_BP = 4*hBS*hUT*fc/3e8;
fprintf('---- 大尺度：路损 ----\n');
fprintf('UMa 断点距离 d_BP = 4*%d*%.1f*%.1fGHz/c = %.0f m\n', ...
        hBS, hUT, fc/1e9, d_BP);
fprintf('  物理含义：直射径与地面反射径的路径差达到半波长的距离\n');
fprintf('  d < d_BP 两径近似同相，相长干涉，衰减慢（指数约 2.2）\n');
fprintf('  d > d_BP 进入相消干涉区，衰减加快（指数约 4）\n\n');

plCfg = nrPathLossConfig;
plCfg.Scenario = 'UMa';
plCfg.BuildingHeight = 20;
plCfg.StreetWidth = 20;

figure('Color','w','Position',[60 40 1280 820]);

subplot(2,3,1);
scenarios = {'UMa', 'UMi', 'RMa', 'InH'};
lineStyles = {'-', '--', '-.', ':'};
colors = lines(4);
hold on;
for k = 1:numel(scenarios)
    cfg = nrPathLossConfig;
    cfg.Scenario = scenarios{k};
    cfg.BuildingHeight = 20;
    cfg.StreetWidth = 20;

    plLOS = nan(size(d));
    for n = 1:numel(d)
        try
            plLOS(n) = nrPathLoss(cfg, fc, d(n), true, hBS, hUT);
        catch
            % 某些场景在某些距离上不适用（如 InH 的距离上限）
        end
    end
    plot(d, plLOS, lineStyles{k}, 'Color', colors(k,:), 'LineWidth', 1.8);
end
hold off;
set(gca, 'XScale', 'log'); grid on;
xlabel('距离 (m)'); ylabel('路损 (dB)');
title('四类场景的 LOS 路损对比 (3.5 GHz)', 'FontSize', 10);
legend(scenarios, 'Location', 'northwest', 'FontSize', 8);

subplot(2,3,2);
hold on;
for los = [true false]
    pl = nan(size(d));
    for n = 1:numel(d)
        pl(n) = nrPathLoss(plCfg, fc, d(n), los, hBS, hUT);
    end
    if los
        plot(d, pl, 'LineWidth', 2, 'Color', [0.1 0.5 0.1]);
    else
        plot(d, pl, 'LineWidth', 2, 'Color', [0.85 0.1 0.15]);
    end
end
% 标出断点
xline(d_BP, '--k', sprintf('d_{BP} = %.0f m', d_BP), ...
      'LabelOrientation', 'horizontal', 'FontSize', 9);
% 叠加自由空间作为参照
plFS = 32.44 + 20*log10(d/1000) + 20*log10(fc/1e6);
plot(d, plFS, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.6);
hold off;
set(gca, 'XScale', 'log'); grid on;
xlabel('距离 (m)'); ylabel('路损 (dB)');
title('UMa: LOS vs NLOS vs 自由空间', 'FontSize', 10);
legend({'LOS', 'NLOS', '自由空间'}, 'Location', 'northwest', 'FontSize', 8);

fprintf('  LOS 与 NLOS 在 1 km 处的差距：');
plL = nrPathLoss(plCfg, fc, 1000, true, hBS, hUT);
plN = nrPathLoss(plCfg, fc, 1000, false, hBS, hUT);
fprintf('%.1f dB\n\n', plN - plL);

%% ---- 2. 阴影衰落对覆盖率的影响 ----
% 路损给的是均值，实际接收功率在其周围对数正态波动。
% 关键点：阴影衰落方差 σ 的增大对中断概率的影响是非线性的，
% 比均值偏移更狠。这是做覆盖分析最容易忽略的地方。

sigma = [4 6 8 10];          % 不同场景的典型 σ
fadeMargin = 0:0.5:25;       % 阴影衰落裕量

subplot(2,3,3);
hold on;
for k = 1:numel(sigma)
    % 给定裕量 M，中断概率 = P(X > M) = Q(M/sigma)
    pOutage = qfunc(fadeMargin ./ sigma(k));
    plot(fadeMargin, pOutage, 'LineWidth', 1.8);
end
hold off;
set(gca, 'YScale', 'log'); grid on;
xlabel('阴影衰落裕量 (dB)'); ylabel('中断概率');
title('阴影衰落裕量与中断概率', 'FontSize', 10);
legend(arrayfun(@(s) sprintf('\\sigma = %d dB', s), sigma, ...
       'UniformOutput', false), 'Location', 'southwest', 'FontSize', 8);

% 给出典型覆盖率所需的裕量
fprintf('---- 阴影衰落裕量查询 ----\n');
fprintf('%-10s %-16s %-16s\n', 'sigma', '95% 覆盖率裕量', '99% 覆盖率裕量');
for s = sigma
    m95 = s * norminv(0.95);
    m99 = s * norminv(0.99);
    fprintf('%-10s %-16s %-16s\n', sprintf('%d dB', s), ...
            sprintf('%.2f dB', m95), sprintf('%.2f dB', m99));
end
fprintf('\n');

%% ==================== 第二部分：小尺度 TDL ====================

%% ---- 3. TDL 信道：时域冲激响应与频率响应 ----
% TDL 用 L 个抽头模拟多径，每个抽头有固定时延和平均功率，幅度随时间随机变化。
% 这是最实用的信道模型 —— 计算量小，参数直观。

% 采样率的选择至关重要：它必须能"分辨"信道的时延扩展。
% 300 ns 的时延扩展要求采样周期远小于 300 ns，即 fs 远大于 3.3 MHz。
% 这里用 30 kHz SCS、52 RB 对应的实际 OFDM 采样率，保证物理意义正确。
carrierRef = nrCarrierConfig;
carrierRef.SubcarrierSpacing = 30;
carrierRef.NSizeGrid = 52;
[~, ofdmInfo] = nrOFDMModulate(carrierRef, nrResourceGrid(carrierRef));
fs = ofdmInfo.SampleRate;

fprintf('OFDM 采样率         : %.4g MHz\n', fs/1e6);
% 时延扩展对应的采样点数：这个比值必须远大于 1，否则信道模型分辨不出多径
delaySpreadNs = 300;
samplesPerDelaySpread = delaySpreadNs*1e-9 * fs;
fprintf('每个时延扩展对应的采样点数 : %.1f（应远大于 1 才有物理意义）\n', ...
        samplesPerDelaySpread);

tdl = nrTDLChannel;
tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;    % 300 ns，城市非视距的典型值
tdl.CarrierFrequency = fc;
tdl.MaximumDopplerShift = 100;   % 对应 3.5 GHz 下约 30 km/h
tdl.SampleRate = fs;
tdl.NumTransmitAntennas = 1;
tdl.NumReceiveAntennas = 1;

fprintf('\n---- 小尺度：TDL 信道 ----\n');
fprintf('延迟剖面            : %s\n', tdl.DelayProfile);
fprintf('均方根时延扩展      : %.1f ns\n', tdl.DelaySpread*1e9);
fprintf('最大多普勒频移      : %.0f Hz\n', tdl.MaximumDopplerShift);
fprintf('采样率              : %.4g MHz\n', fs/1e6);

% 平均功率延迟谱：对多个时刻的冲激响应取平均
nSnapshots = 50;
pdpLen = round(8e-6 * fs);          % 观察窗 8 us，足够覆盖 300 ns 的时延扩展
pdpSnapshots = zeros(pdpLen, nSnapshots);
for k = 1:nSnapshots
    tdl.reset();
    [rxK, ~] = tdl([1; zeros(1000, 1)]);
    seg = abs(rxK(1:pdpLen)).^2;
    pdpSnapshots(:, k) = seg;
end
pdp = mean(pdpSnapshots, 2);
pdp = pdp / max(pdp);
delayAxis = (0:pdpLen-1) / fs * 1e9;    % 单位 ns

subplot(2,3,4);
plot(delayAxis, 10*log10(pdp + 1e-12), 'LineWidth', 1.5, 'Color', [0.1 0.3 0.85]);
grid on;
xlabel('时延 (ns)'); ylabel('归一化功率 (dB)');
title(sprintf('功率延迟谱 (PDP)  — %s, \\tau_{rms} = %.0f ns', ...
              tdl.DelayProfile, tdl.DelaySpread*1e9), 'FontSize', 10);
xlim([0 min(5000, max(delayAxis))]); ylim([-40 2]);

% 从 PDP 测量均方根时延扩展
pdpLin = pdp(:);
tauAxis = (0:length(pdpLin)-1)' / fs;
pdpSum = sum(pdpLin);
if pdpSum > 0
    tauMean = sum(pdpLin .* tauAxis) / pdpSum;
    tauRmsMeas = sqrt(sum(pdpLin .* (tauAxis - tauMean).^2) / pdpSum);
else
    tauRmsMeas = NaN;
end
fprintf('从 PDP 实测的 tau_rms : %.1f ns（设定值 %.0f ns）\n', ...
        tauRmsMeas*1e9, tdl.DelaySpread*1e9);

%% ---- 4. 频率响应与频率选择性 ----
% 信道的频率响应是 PDP 的傅里叶变换。时延展得越宽，频率响应起伏越剧烈，
% 也就是频率选择性越强。

Nfft = 8192;                       % 频谱分辨率 fs/Nfft 应远细于子载波间隔
H = fftshift(fft(pdp, Nfft));
freqAxis = (-Nfft/2:Nfft/2-1)' * (fs/Nfft) / 1e6;   % MHz

subplot(2,3,5);
plot(freqAxis, 20*log10(abs(H)/max(abs(H)) + 1e-12), 'LineWidth', 1.2, ...
     'Color', [0.85 0.4 0.1]);
grid on;
xlabel('频率 (MHz)'); ylabel('|H(f)| (dB)');
title('信道频率响应：起伏 = 频率选择性', 'FontSize', 10);
ylim([-40 2]);
xlim([-fs/2/1e6 fs/2/1e6]);

% 相干带宽：频率相关函数降到 0.5 对应的带宽
fprintf('\n');

%% ---- 5. 相干带宽与相干时间的实测对比 ----
% 这是本脚本最有价值的部分：把理论公式和实测数值放在一起。

% --- 相干带宽 ---
% 严格做法：算信道的频率相关函数 R(df) = E[H(f+df) H*(f)]
% 简化做法：用 PDP 的傅里叶变换
df = linspace(0, 2e6, 500);
Rdf = zeros(size(df));
for k = 1:numel(df)
    Rdf(k) = abs(sum(pdpLin .* exp(-1j*2*pi*df(k)*tauAxis)));
end
Rdf = Rdf / Rdf(1);
idxBc = find(Rdf < 0.5, 1);
if ~isempty(idxBc)
    Bc_meas = df(idxBc);
else
    Bc_meas = NaN;
end

% 理论预测（两种常用经验公式）
Bc_theory_5  = 1/(5  * tdl.DelaySpread);
Bc_theory_50 = 1/(50 * tdl.DelaySpread);

% --- 相干时间 ---
% 用信道增益的时间相关来测
tdl.reset();
nSamp = round(0.2 * fs);              % 200 ms
[rxLong, pgLong] = tdl(complex(randn(nSamp,1), randn(nSamp,1))/sqrt(2));
% 从接收信号的自相关估计时间相关
maxLag = round(0.02 * fs);
xc = zeros(maxLag+1, 1);
x = rxLong(:); x = x - mean(x);
for l = 0:maxLag
    xc(l+1) = abs(sum(x(1:end-l) .* conj(x(1+l:end))));
end
xc = xc / xc(1);
tauAxisT = (0:maxLag)' / fs;
idxTc = find(xc < 0.5, 1);
if ~isempty(idxTc)
    Tc_meas = tauAxisT(idxTc);
else
    Tc_meas = NaN;
end
Tc_theory = 0.2421 / tdl.MaximumDopplerShift;   % J0(2*pi*f_d*Tc)=0.5 的严格解

fprintf('%s\n', repmat('=', 1, 80));
fprintf('相干量：理论 vs 实测\n');
fprintf('%s\n', repmat('-', 1, 80));
fprintf('%-14s %-22s %-22s\n', '量', '理论预测', '本次实测');
fprintf('%-14s %-22s %-22s\n', '相干带宽', ...
        sprintf('%.1f kHz (1/5tau)', Bc_theory_5/1e3), ...
        sprintf('%.1f kHz', Bc_meas/1e3));
fprintf('%-14s %-22s %-22s\n', '', ...
        sprintf('%.1f kHz (1/50tau)', Bc_theory_50/1e3), '');
fprintf('%-14s %-22s %-22s\n', '相干时间', ...
        sprintf('%.2f ms (0.242/f_d)', Tc_theory*1e3), ...
        sprintf('%.2f ms', Tc_meas*1e3));
fprintf('%s\n', repmat('=', 1, 80));

fprintf('\n关键判断（把信道特性和 NR 参数设计连起来）：\n');
fprintf('  * 相干带宽 %.1f kHz vs 子载波间隔 30 kHz：\n', Bc_meas/1e3);
if Bc_meas > 30e3
    fprintf('    相干带宽 > 子载波间隔，信道在单个子载波上基本是平的，\n');
    fprintf('    但跨子载波有起伏 -> 均衡器必须逐子载波处理\n');
else
    fprintf('    相干带宽 < 子载波间隔，信道在一个子载波内就变了，\n');
    fprintf('    会有严重的子载波间干扰，需要更密的导频\n');
end
fprintf('  * 相干时间 %.2f ms vs OFDM 符号时长：\n', Tc_meas*1e3);
symLen_ms = 1/(2*1000)/14 * 1e3;   % 30 kHz SCS，一个时隙 0.5ms 含 14 符号
fprintf('    当前 OFDM 符号时长约 %.4f ms，相干时间 %.4f ms，比值 %.1f\n', ...
        symLen_ms, Tc_meas*1e3, Tc_meas*1e3/symLen_ms);
fprintf('    若比值远大于 1 -> 慢衰落，一个符号内信道基本不变\n');
fprintf('    若接近或小于 1 -> 快衰落，必须在同一时隙内插入更多 DM-RS\n');
fprintf('    这就是"DM-RS 时域密度"这个设计决策的全部依据\n\n');

%% ==================== 第三部分：CDL 与空间维度 ====================

subplot(2,3,6);
try
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-C';
    cdl.DelaySpread = 300e-9;
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = 100;
    cdl.SampleRate = fs;
    cdl.ChannelFiltering = false;      % 只要路径增益，不要滤波输出
    cdl.NumTimeSamples = 100;
    cdl.NormalizePathGains = true;

    % 配置 4x4 的均匀线性阵（ULA），这就是空间相关性的来源
    cdl.TransmitAntennaArray.Size = [1 4 1 1 1];
    cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 0.5 0.5];
    cdl.ReceiveAntennaArray.Size = [1 4 1 1 1];
    cdl.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 0.5 0.5];

    [~, cdlPathGains] = cdl();

    % 路径增益的维度随版本变化，常见为
    %   [Ncs Nt Nr Npath Nsample]
    % 这里不硬编码，而是用样本数和天线数反推，稳妥处理维度差异。
    pg = cdlPathGains;
    Nt = 4; Nr = 4;
    sz = size(pg);

    % 把增益重排成 [Nt*Nr, 其余维度合并]，每个天线对一个列向量
    if ndims(pg) >= 5 && all(sz(2:3) == [Nt Nr])
        pg = permute(pg, [2 3 1 4 5]);            % [Nt Nr Ncs Npath Nsample]
        Hmat = reshape(pg, Nt*Nr, []);            % [Nt*Nr, Ncs*Npath*Nsample]
    else
        % 维度不符时退回到一个通用重排，保证不报错
        Hmat = reshape(pg, [], sz(end));
    end

    % 空间相关矩阵 R = E[h h^H]，每个天线对作为一个"端口"
    R = (Hmat * Hmat') / size(Hmat, 2);
    R = R / max(abs(R(:)));                       % 归一化便于观察

    imagesc(abs(R));
    colormap(gca, parula); colorbar;
    axis square;
    title(sprintf('空间相关矩阵 |R| (CDL-C, %dx%d ULA)', Nt, Nr), 'FontSize', 10);
    xlabel('天线对索引'); ylabel('天线对索引');

    fprintf('---- CDL 与空间维度 ----\n');
    fprintf('CDL 相对 TDL 增加的是角度维度：每个簇有离开角和到达角，\n');
    fprintf('不同阵元看到同一簇的到达角不同 -> 产生相位差 -> 这就是空间相关性\n');
    fprintf('空间相关性高 -> MIMO 空间复用增益下降 -> 容量降低\n');
    fprintf('所以 CDL 是"多天线仿真"的必需模型，TDL 只能做单天线或纯时频分析\n\n');
catch ME
    fprintf('CDL 配置失败（%s）—— 属性名与版本不匹配\n', ME.identifier);
    fprintf('用 doc nrCDLChannel 核对 TransmitAntennaArray 等属性，跳过本节即可\n\n');
    text(0.1, 0.5, 'CDL 配置与当前版本不匹配，', 'FontSize', 10);
    text(0.1, 0.3, '请用 doc nrCDLChannel 核对', 'FontSize', 10);
    axis off;
end

%% ---- 验收自查 ----
fprintf('%s\n', repmat('=', 1, 80));
fprintf('验收自查：\n');
fprintf('  1. 断点距离的物理来源是什么？为什么路损在它两侧斜率不同？\n');
fprintf('  2. 把 DelaySpread 从 300ns 改到 3000ns，相干带宽怎么变？图上的\n');
fprintf('     频率响应起伏变剧烈还是平缓？\n');
fprintf('  3. 把 MaximumDopplerShift 从 100 改到 1000，相干时间怎么变？\n');
fprintf('     这对 DM-RS 的时域密度意味着什么？\n');
fprintf('  4. TDL 和 CDL 各自适合什么仿真场景？\n');
fprintf('  5. σ 从 4dB 增到 10dB，达到 95%% 覆盖率所需的裕量增加多少？\n');
fprintf('%s\n', repmat('=', 1, 80));
