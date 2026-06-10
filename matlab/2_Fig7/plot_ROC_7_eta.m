%% =========================================================
% 绘制 eta 指标消融实验 ROC 曲线 (pre-ictal vs early)
% 数据源: ALL_CASES_stage5_eta_ablation.csv
% =========================================================
clc; clear; close all;

dataDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig7';
dataFile = fullfile(dataDir, 'ALL_CASES_stage5_eta_ablation.csv');

if ~exist(dataFile, 'file')
    error('未找到数据文件，请检查路径: %s', dataFile);
end

T = readtable(dataFile);

targetPhases = {'pre-ictal', 'early'};
positiveClass = 'early';

modelNames = {'Original', 'Coverage only', 'DegreeBias only'};
scoreColumns = {'eta_orig_mean', 'eta_cov_mean', 'eta_deg_mean'};

colors_NPG = [
     0, 160, 135;
    60,  84, 136;
   230,  75,  53
] / 255;

roc_data = cell(length(modelNames), 1);
auc_results = zeros(length(modelNames), 1);

fprintf('====================================================\n');
fprintf('  eta 消融实验 ROC / AUC 分析 (pre-ictal vs early)\n');
fprintf('====================================================\n');

idx = ismember(T.phase5, targetPhases);
labels = T.phase5(idx);

for m = 1:length(modelNames)
    scores = T.(scoreColumns{m})(idx);
    valid_idx = ~isnan(scores);
    cur_labels = labels(valid_idx);
    cur_scores = scores(valid_idx);

    if isempty(cur_labels)
        continue;
    end

    [X, Y, ~, AUC] = perfcurve(cur_labels, cur_scores, positiveClass);
    if AUC < 0.5
        [X, Y, ~, AUC] = perfcurve(cur_labels, -cur_scores, positiveClass);
    end

    roc_data{m}.X = X;
    roc_data{m}.Y = Y;
    auc_results(m) = AUC;

    fprintf('%-25s AUC = %.4f\n', modelNames{m}, AUC);
end
fprintf('====================================================\n');

fig = figure('Color', 'w', 'Position', [150, 150, 520, 520], 'Name', 'ROC_Ablation_eta_Nature');
ax = axes(fig);
hold(ax, 'on');

plot([0 1], [0 1], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);

hLines = gobjects(length(modelNames), 1);
legend_labels = cell(length(modelNames), 1);

for m = 1:length(modelNames)
    if auc_results(m) == 0
        continue;
    end
    hLines(m) = plot(roc_data{m}.X, roc_data{m}.Y, '-', ...
        'Color', colors_NPG(m, :), ...
        'LineWidth', 2.5);
    legend_labels{m} = sprintf('%s (AUC = %.3f)', modelNames{m}, auc_results(m));
end

xlim([0 1.01]);
ylim([0 1.01]);
axis square;

xlabel('False Positive Rate (1 - Specificity)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Positive Rate (Sensitivity)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
title('Ablation ROC: pre-ictal vs early', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

set(ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'TickDir', 'out', ...
    'TickLength', [0.015 0.015], ...
    'Box', 'off');

leg = legend(hLines, legend_labels, 'Location', 'southeast');
set(leg, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 11);

set(fig, 'PaperPositionMode', 'auto');
export_png = fullfile(dataDir, 'Nature_Style_ROC_Ablation_eta.png');
export_pdf = fullfile(dataDir, 'Nature_Style_ROC_Ablation_eta.pdf');

exportgraphics(fig, export_png, 'Resolution', 600);
exportgraphics(fig, export_pdf, 'ContentType', 'vector');

fprintf('\n✅ eta 消融实验 ROC 曲线绘制完成！\n');
fprintf('PNG 图片已保存至: %s\n', export_png);
fprintf('PDF 矢量图已保存至: %s\n', export_pdf);
