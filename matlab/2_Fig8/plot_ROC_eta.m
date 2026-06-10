%% =========================================================
% 计算并绘制 eta 指标在 pre-ictal 与 early 阶段的 ROC 曲线
% 数据源: ALL_CASES_window_level_eta_Comparison.csv
% =========================================================
clc; clear; close all;

%% 1. 数据加载
dataDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig8';
dataFile = fullfile(dataDir, 'ALL_CASES_window_level_eta_Comparison.csv');

if ~exist(dataFile, 'file')
    error('未找到数据文件，请检查路径: %s', dataFile);
end

T = readtable(dataFile);

targetPhases = {'pre-ictal', 'early'};
positiveClass = 'early';
modelOrder = {'Low-Order(2-nodes)', 'High-Order(3-nodes)', 'DMW-HLG'};

colors_NPG = [
     0, 160, 135;
    60,  84, 136;
   230,  75,  53
] / 255;

%% 2. 计算 ROC 与 AUC
roc_data = cell(length(modelOrder), 1);
auc_results = zeros(length(modelOrder), 1);

fprintf('====================================================\n');
fprintf('  eta ROC / AUC 分析 (pre-ictal vs early)\n');
fprintf('====================================================\n');

for m = 1:length(modelOrder)
    current_model = modelOrder{m};

    idx = strcmp(T.Model, current_model) & ismember(T.phase5, targetPhases);

    labels = T.phase5(idx);
    scores = T.eta(idx);

    valid_idx = ~isnan(scores);
    labels = labels(valid_idx);
    scores = scores(valid_idx);

    if isempty(labels)
        continue;
    end

    [X, Y, ~, AUC] = perfcurve(labels, scores, positiveClass);

    if AUC < 0.5
        scores_flipped = -scores;
        [X, Y, ~, AUC] = perfcurve(labels, scores_flipped, positiveClass);
    end

    roc_data{m}.X = X;
    roc_data{m}.Y = Y;
    auc_results(m) = AUC;

    fprintf('%-25s AUC = %.4f\n', current_model, AUC);
end
fprintf('====================================================\n');

%% 3. 绘图
fig = figure('Color', 'w', 'Position', [100, 100, 500, 500], 'Name', 'ROC_Comparison_eta_Nature');
ax = axes(fig);
hold(ax, 'on');

plot([0 1], [0 1], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);

hLines = gobjects(length(modelOrder), 1);
legend_labels = cell(length(modelOrder), 1);

for m = 1:length(modelOrder)
    if auc_results(m) == 0
        continue;
    end

    hLines(m) = plot(roc_data{m}.X, roc_data{m}.Y, '-', ...
        'Color', colors_NPG(m, :), ...
        'LineWidth', 2.5);

    legend_labels{m} = sprintf('%s (AUC = %.3f)', modelOrder{m}, auc_results(m));
end

xlim([0 1.01]);
ylim([0 1.01]);
axis square;

xlabel('False Positive Rate (1 - Specificity)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Positive Rate (Sensitivity)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
title('ROC Analysis: pre-ictal vs early', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

set(ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'TickDir', 'out', ...
    'TickLength', [0.015 0.015], ...
    'Box', 'off');

leg = legend(hLines, legend_labels, 'Location', 'southeast');
set(leg, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 11);

%% 4. 导出
export_png = fullfile(dataDir, 'Nature_Style_ROC_eta_Comparison.png');
export_pdf = fullfile(dataDir, 'Nature_Style_ROC_eta_Comparison.pdf');

exportgraphics(fig, export_png, 'Resolution', 600);
exportgraphics(fig, export_pdf, 'ContentType', 'vector');

fprintf('\n✅ eta ROC 对比图绘制完成！\n');
fprintf('已保存为 600DPI PNG: %s\n', export_png);
fprintf('已保存为 矢量图 PDF: %s\n', export_pdf);
