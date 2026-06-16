%% Fig2 切片图：Sigma W / N
% 五条线固定为 0.1:0.2:0.9
% 左图：固定 h，观察 Sigma W / N 随 c 变化
% 右图：固定 c，观察 Sigma W / N 随 h 变化

clc; clear; close all;

outDir = 'example_project\1_Fig9_eta_mechanism';
cacheFile = fullfile(outDir, 'eta_phase_mechanism_data.mat');

if ~exist(cacheFile, 'file')
    error('找不到缓存文件: %s\n请先运行前面的相图脚本生成 eta_phase_mechanism_data.mat。', cacheFile);
end

S = load(cacheFile);
R = S.R;

c_range = R.c_range;
h_range = R.h_range;
avgW_mean = R.avgW_mean;

slice_vals = 0.1:0.2:0.9;
nSlice = numel(slice_vals);

h_idx = zeros(1, nSlice);
c_idx = zeros(1, nSlice);

for i = 1:nSlice
    [~, h_idx(i)] = min(abs(h_range - slice_vals(i)));
    [~, c_idx(i)] = min(abs(c_range - slice_vals(i)));
end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17.5 7.5]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

fontName = 'Arial';
fontTick = 8.5;
fontAxis = 9.5;
lineColors = lines(nSlice);

ax1 = nexttile(t, 1); hold(ax1, 'on');
for i = 1:nSlice
    plot(ax1, c_range, avgW_mean(h_idx(i), :), '-o', ...
        'LineWidth', 1.8, 'MarkerSize', 5.0, ...
        'Color', lineColors(i, :), ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', sprintf('h = %.1f', h_range(h_idx(i))));
end
xlabel(ax1, 'Core concentration c', 'FontName', fontName, 'FontSize', fontAxis, 'FontWeight', 'bold');
ylabel(ax1, '\Sigma W / N', 'FontName', fontName, 'FontSize', fontAxis, 'FontWeight', 'bold');
title(ax1, 'Fixed h: \Sigma W / N vs c', 'FontName', fontName, 'FontSize', 10.5, 'FontWeight', 'bold');
set(ax1, 'FontName', fontName, 'FontSize', fontTick, 'LineWidth', 1.0, 'Box', 'off', 'TickDir', 'out');
legend(ax1, 'Location', 'best', 'Box', 'off', 'FontName', fontName, 'FontSize', 8.5);

ax2 = nexttile(t, 2); hold(ax2, 'on');
for i = 1:nSlice
    plot(ax2, h_range, avgW_mean(:, c_idx(i)), '-o', ...
        'LineWidth', 1.8, 'MarkerSize', 5.0, ...
        'Color', lineColors(i, :), ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', sprintf('c = %.1f', c_range(c_idx(i))));
end
xlabel(ax2, 'Homogenization h', 'FontName', fontName, 'FontSize', fontAxis, 'FontWeight', 'bold');
ylabel(ax2, '\Sigma W / N', 'FontName', fontName, 'FontSize', fontAxis, 'FontWeight', 'bold');
title(ax2, 'Fixed c: \Sigma W / N vs h', 'FontName', fontName, 'FontSize', 10.5, 'FontWeight', 'bold');
set(ax2, 'FontName', fontName, 'FontSize', fontTick, 'LineWidth', 1.0, 'Box', 'off', 'TickDir', 'out');
legend(ax2, 'Location', 'best', 'Box', 'off', 'FontName', fontName, 'FontSize', 8.5);

pngFile = fullfile(outDir, 'Fig2_avgW_representative_slices_5lines.png');
pdfFile = fullfile(outDir, 'Fig2_avgW_representative_slices_5lines.pdf');
figFile = fullfile(outDir, 'Fig2_avgW_representative_slices_5lines.fig');

exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
savefig(fig, figFile);

fprintf('已保存:\n%s\n%s\n%s\n', pngFile, pdfFile, figFile);
