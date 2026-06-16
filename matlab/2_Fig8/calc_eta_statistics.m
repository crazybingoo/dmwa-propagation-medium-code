%% =========================================================
% 绘制 eta 指标统计结果：单面板方差图 + 右上角整体趋势小图
% =========================================================
clc; clear; close all;

%% 1. 路径设置与数据加载
dataDir = 'example_project\2_Fig8';
varFile = fullfile(dataDir, 'Stats_Variance_Comparison_eta.csv');

if ~exist(varFile, 'file')
    % 如果方差表尚未生成，则尝试从阶段总表现场生成
    stageFile = fullfile(dataDir, 'ALL_CASES_stage5_eta_Comparison.csv');
    if ~exist(stageFile, 'file')
        error('未找到方差数据文件，也未找到阶段总表，请检查路径: %s', dataDir);
    end

    T_stage = readtable(stageFile);
    modelOrder0 = {'Low-Order(2-nodes)', 'High-Order(3-nodes)', 'DMW-HLG'};
    phaseOrder0 = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};

    statsRows = [];
    for m = 1:numel(modelOrder0)
        for p = 1:numel(phaseOrder0)
            idx = strcmp(T_stage.Model, modelOrder0{m}) & strcmp(T_stage.phase5, phaseOrder0{p});
            x = T_stage.eta_mean(idx);
            x = x(~isnan(x));
            if ~isempty(x)
                statsRows = [statsRows; {modelOrder0{m}, phaseOrder0{p}, var(x, 'omitnan')}]; %#ok<AGROW>
            end
        end
    end

    T_var = cell2table(statsRows, 'VariableNames', {'Model', 'Phase', 'Variance'});
    writetable(T_var, varFile);
else
    T_var = readtable(varFile);
end

%% 2. 核心参数与配色设置
modelOrder = {'Low-Order(2-nodes)', 'High-Order(3-nodes)', 'DMW-HLG'};
phaseOrder = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
phaseShort = {'Pre', 'Early', 'Mid', 'Late', 'Post'};
insetShort = {'Low', 'High', 'DMW'};

colors_NPG = [
     0, 160, 135;
    60,  84, 136;
   230,  75,  53
] / 255;

nModel = length(modelOrder);
nPhase = length(phaseOrder);

%% 3. 数据重组与 Overall 计算
var_matrix = zeros(nPhase, nModel);
for m = 1:nModel
    for p = 1:nPhase
        idx = strcmp(T_var.Model, modelOrder{m}) & strcmp(T_var.Phase, phaseOrder{p});
        if any(idx)
            var_matrix(p, m) = T_var.Variance(idx);
        end
    end
end

overall_var = mean(var_matrix, 1);

%% 4. 初始化主图
fig = figure('Color', 'w', 'Position', [150, 150, 680, 500], 'Name', 'Variance_With_Inset_eta');
ax_main = axes(fig);
hold(ax_main, 'on');

b = bar(ax_main, 1:nPhase, var_matrix, 0.75, 'grouped', ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.8);

for m = 1:nModel
    b(m).FaceColor = colors_NPG(m, :);
    b(m).FaceAlpha = 0.9;
end

%% 5. 添加主图水平辅助线
y_ticks = yticks(ax_main);
for i = 2:length(y_ticks)
    plot(ax_main, [0.5, nPhase + 0.5], [y_ticks(i), y_ticks(i)], ...
        'Color', [0.85 0.85 0.85], 'LineStyle', '--', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end
uistack(b, 'top');

%% 6. 主图坐标轴格式化
set(ax_main, 'XTick', 1:nPhase, 'XTickLabel', phaseShort, ...
    'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.0, ...
    'TickDir', 'out', 'TickLength', [0.01 0.01], ...
    'Box', 'off', 'XColor', 'k', 'YColor', 'k');

ylabel(ax_main, 'Variance', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
ax_main.YAxis.Exponent = -4;

max_var = max(var_matrix(:));
ylim(ax_main, [0, max_var * 1.25]);

%% 7. 主图图例
leg = legend(ax_main, b, modelOrder, 'Location', 'northoutside', 'Orientation', 'horizontal');
set(leg, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 11);
ax_main.Position = [0.12, 0.12, 0.83, 0.73];

%% 8. 添加右上角小图
ax_inset = axes('Position', [0.65, 0.50, 0.25, 0.25]);
hold(ax_inset, 'on');

plot(ax_inset, 1:nModel, overall_var, '-k', 'LineWidth', 1.2, 'Color', [0.3 0.3 0.3]);

for m = 1:nModel
    scatter(ax_inset, m, overall_var(m), 50, ...
        'MarkerFaceColor', colors_NPG(m, :), ...
        'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.0);
end

set(ax_inset, 'XTick', 1:nModel, 'XTickLabel', insetShort, ...
    'FontName', 'Arial', 'FontSize', 10, 'LineWidth', 0.8, ...
    'TickDir', 'out', 'Box', 'off');
ylabel(ax_inset, 'Var', 'FontName', 'Arial', 'FontSize', 10);
title(ax_inset, 'Overall', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'normal');

y_margin = (max(overall_var) - min(overall_var)) * 0.25;
if y_margin == 0
    y_margin = max(overall_var) * 0.1;
end
ylim(ax_inset, [min(overall_var) - y_margin, max(overall_var) + y_margin]);
xlim(ax_inset, [0.5, nModel + 0.5]);
ax_inset.YAxis.Exponent = -4;

%% 9. 导出
set(fig, 'PaperPositionMode', 'auto');
export_png = fullfile(dataDir, 'Nature_Style_Variance_eta_Final.png');
export_pdf = fullfile(dataDir, 'Nature_Style_Variance_eta_Final.pdf');

exportgraphics(fig, export_png, 'Resolution', 600);
exportgraphics(fig, export_pdf, 'ContentType', 'vector');

fprintf('\n✅ eta 方差图绘制完成！\n');
fprintf('已保存至: %s\n', export_png);
