%% ch01 参数集、帧结构与资源网格
%
% 目标：把 NR 的"时频舞台"彻底搞清楚。
%   1. 走遍 5 种参数集，看清子载波间隔、时隙长度、采样率的关系
%   2. 理解资源网格的坐标系（行=子载波，列=符号）
%   3. 搞清 BWP 偏移 NStartGrid 的作用
%
% 版本：目标 R2021a 之后的 API。若报属性不存在，用 doc nrCarrierConfig 核对。
%
% 运行：直接 F5。全程无随机性，结果可复现。

clear; close all; clc;

%% ---- 1. 五种参数集的横向对比 ----
% NR 取消了固定子载波间隔，改用 mu 定义：df = 2^mu * 15 kHz
% mu 越大 -> 符号越短 -> 对高速移动（多普勒）越宽容
%         但 CP 绝对长度也越短 -> 对大时延扩展越不宽容，覆盖越差

mu = 0:4;
scs = 15 * 2.^mu;              % 子载波间隔，单位 kHz

fprintf('%s\n', repmat('=', 1, 86));
fprintf('%-4s %-10s %-12s %-14s %-14s %-14s\n', ...
    'mu', 'SCS', '每子帧时隙', '时隙长度', '每时隙样本', '采样率');
fprintf('%s\n', repmat('-', 1, 86));

for k = 1:numel(mu)
    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = scs(k);
    carrier.NSizeGrid = 52;     % 频域占 52 个 RB
    carrier.NSlot = 0;
    carrier.NFrame = 0;

    % 调制一个空网格，从 info 里取实际采样率
    % （不硬编码采样率：它由工具箱内部的 FFT 尺寸和 SCS 一起决定）
    [~, info] = nrOFDMModulate(carrier, nrResourceGrid(carrier));

    % 时隙长度 = 1 / (每子帧时隙数 * 1000) 秒；每时隙样本数 = 采样率 * 时隙长度
    slotsPerSubframe = 2^mu(k);
    slotLen_ms = 1 / slotsPerSubframe;
    samplesPerSlot = info.SampleRate * slotLen_ms * 1e-3;

    fprintf('%-4d %-10s %-12d %-14s %-14d %-14s\n', ...
        mu(k), sprintf('%d kHz', scs(k)), slotsPerSubframe, ...
        sprintf('%.4f ms', slotLen_ms), round(samplesPerSlot), ...
        sprintf('%.4g MHz', info.SampleRate/1e6));
end
fprintf('%s\n\n', repmat('=', 1, 86));

fprintf('观察要点：\n');
fprintf('  * 时隙长度随 mu 翻倍递减，但无线帧(10ms)和子帧(1ms)的长度固定不变\n');
fprintf('  * 采样率随 mu 递增，因为要表示的带宽变宽了\n');
fprintf('  * 常规 CP 下每时隙恒为 14 个符号，但符号的绝对时长在缩短\n\n');

%% ---- 2. 一个具体载波的派生量 ----
% disp(carrier) 会打印所有由基本参数派生出来的量，逐个和 38.211 对照

carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid = 52;
carrier.NSlot = 3;
carrier.NFrame = 1;

disp('---- carrier 对象的全部属性（含派生量）----');
disp(carrier);

grid = nrResourceGrid(carrier);
fprintf('资源网格尺寸：[%d 子载波 x %d 符号]\n', size(grid,1), size(grid,2));
fprintf('  = [12 * %d RB  x  14 符号]\n', carrier.NSizeGrid);

%% ---- 3. 资源网格长什么样 ----
% 网格的坐标系是理解所有物理信道的起点：
%   行 = 子载波（频域，从低到高）
%   列 = OFDM 符号（时域，从左到右）
% 每个格子是一个 RE，1 个 RB 占连续 12 行

% 造一个可视化用的示意网格：不同 RB 涂不同颜色
visGrid = zeros(size(grid));
for rb = 0:carrier.NSizeGrid-1
    visGrid(rb*12 + (1:12), :) = rb + 1;
end

figure('Color','w','Position',[100 100 900 620]);

subplot(2,1,1);
imagesc(visGrid);
colormap(gca, parula);
xlabel('OFDM 符号索引'); ylabel('子载波索引');
title(sprintf('资源网格示意：%d 个 RB，每 RB 12 个子载波，共 14 个符号', carrier.NSizeGrid));
set(gca,'YDir','normal');
colorbar('Ticks', unique(visGrid(:))); caxis([0 carrier.NSizeGrid]);

subplot(2,1,2);
hold on;
for rb = 0:carrier.NSizeGrid-1
    y0 = rb*12;
    plot([0 14], [y0 y0], 'k-', 'LineWidth', 0.4);
    rectangle('Position',[0, y0, 14, 12], 'EdgeColor',[0.7 0.7 0.7]);
end
% 标出前两个 RB 的边界，让人看清"一个 RB 是一块 12x14 的矩形"
rectangle('Position',[0, 0, 14, 12], 'EdgeColor',[0.85 0.1 0.15], 'LineWidth', 2);
rectangle('Position',[0, 12, 14, 12], 'EdgeColor',[0.1 0.3 0.85], 'LineWidth', 2);
% 标出前两个符号的一个 RE
rectangle('Position',[0, 12, 1, 1], 'FaceColor',[0.9 0.7 0.1], 'EdgeColor','k');
text(2.2, 12.5, '一个 RE = 1 子载波 x 1 符号', 'FontSize', 10);
text(0.3, 5.5, 'RB #0', 'FontSize', 10, 'Color',[0.85 0.1 0.15]);
text(0.3, 17.5, 'RB #1', 'FontSize', 10, 'Color',[0.1 0.3 0.85]);
hold off;
axis([-0.5 18 -2 carrier.NSizeGrid*12+2]);
xlabel('OFDM 符号索引'); ylabel('子载波索引');
title('层级关系：RE -> RB -> BWP');
set(gca,'YDir','normal');

%% ---- 4. BWP 偏移：NStartGrid ----
% 载波不一定从频率最低端开始用。NStartGrid 给出这个载波相对公共参考点
% （Point A）的偏移，单位是 RB。这是 NR 引入 BWP 概念的基础——
% 终端只需要处理 BWP 内的信号，不必处理整个载波带宽。

NStart_values = [0, 10];
figure('Color','w','Position',[100 100 900 400]);

for k = 1:2
    subplot(1,2,k);
    c = carrier;
    c.NStartGrid = NStart_values(k);
    g = nrResourceGrid(c);

    % 画一张更宽的背景代表"整个载波"，把 BWP 叠上去
    bgWidth = 100;
    imagesc([1 size(g,2)], [1 bgWidth*12], zeros(bgWidth*12, size(g,2)));
    colormap(gca, [1 1 1]);
    hold on;
    y0 = NStart_values(k)*12;
    rectangle('Position', [1, y0, size(g,2), size(g,1)], ...
              'FaceColor', [0.85 0.92 1], 'EdgeColor', [0.1 0.3 0.85], ...
              'LineWidth', 2);
    text(2, y0 - 15, sprintf('NStartGrid = %d RB', NStart_values(k)), ...
         'FontSize', 10, 'Color', [0.1 0.3 0.85]);
    text(2, -60, '（上方白色区域代表载波内本 BWP 未使用的部分）', 'FontSize', 9);
    hold off;
    axis([1 size(g,2) -80 bgWidth*12]);
    set(gca, 'YDir', 'normal');
    xlabel('OFDM 符号'); ylabel('子载波（公共参考点起算）');
    title(sprintf('BWP 位置：NStartGrid = %d', NStart_values(k)));
end

%% ---- 5. OFDM 调制与解调：验证可逆 ----
% 网格 -> 时域波形 -> 网格，应该能无损还原（浮点误差范围内）

testGrid = (randn(size(grid)) + 1j*randn(size(grid))) / sqrt(2);
[waveform, info] = nrOFDMModulate(carrier, testGrid);
recovered = nrOFDMDemodulate(carrier, waveform);

err = max(abs(testGrid(:) - recovered(:)));
fprintf('\n---- OFDM 调制/解调往返验证 ----\n');
fprintf('波形长度        : %d 个样本\n', length(waveform));
fprintf('采样率          : %.6g Hz\n', info.SampleRate);
fprintf('往返最大误差    : %.3e（应接近机器精度 eps≈2.2e-16）\n', err);
% 注：info 结构体的具体字段随版本变化。如需查看全部字段，
%     在命令行执行 disp(info) 或 fieldnames(info) 即可。

if err < 1e-10
    fprintf('结论：可逆，符合预期\n\n');
else
    fprintf('结论：误差偏大，检查是否截图了 CP 或参数不一致\n\n');
end

%% ---- 6. 验收自查 ----
fprintf('%s\n', repmat('=', 1, 86));
fprintf('验收自查：\n');
fprintf('  1. 能否手算 mu=2 时一个时隙的样本数？用上面表格验证。\n');
fprintf('  2. 为什么 mu 越大越适合高速移动，却越不适合大覆盖？\n');
fprintf('     （提示：多普勒对符号时长，时延扩展对 CP 绝对长度）\n');
fprintf('  3. BWP 的作用是什么？为什么 NR 要引入它而不是像 LTE 那样用全带宽？\n');
fprintf('%s\n', repmat('=', 1, 86));
