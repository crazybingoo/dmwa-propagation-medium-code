%% plot_eta_exp4_summary_noABC.m
% =========================================================
% 单独绘图脚本：实验 4 eta burden summary
%
% 改动：
%   1) 不再在标题中使用 A/B/C。
%   2) 修复 SOZ_only、SOZ_PZ_NIZ 等下划线显示问题。
%      原因：MATLAB 默认把下划线当成下标解释。
%      这里统一设置 legend/text 的 Interpreter 为 'none'。
%   3) 使用 subplot，兼容较老 MATLAB 版本。
%
% 运行前提：已经运行过 calc_eta_exp4_region_hyperedge_burden_full.m，
% 并生成以下两个 CSV：
%   ALL_CASES_exp4_stage_region_contribution.csv
%   ALL_CASES_exp4_stage_macro_role_contribution.csv
% =========================================================

clc; clear; close all;

%% 1) 路径设置
outDir = 'example_project\2_Fig5_size_adjusted';

regionCsv = fullfile(outDir, 'ALL_CASES_exp4_stage_region_contribution.csv');
macroCsv  = fullfile(outDir, 'ALL_CASES_exp4_stage_macro_role_contribution.csv');

if ~exist(regionCsv, 'file')
    error('找不到文件：%s\n请先运行 calc_eta_exp4_region_hyperedge_burden_full.m', regionCsv);
end
if ~exist(macroCsv, 'file')
    error('找不到文件：%s\n请先运行 calc_eta_exp4_region_hyperedge_burden_full.m', macroCsv);
end

Tregion = readtable(regionCsv);
TmacroRole = readtable(macroCsv);

%% 2) 顺序与配色
phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
phaseShort = {'Pre','Early','Mid','Late','Post'};

regionNames = {'SOZ','PZ','NIZ'};
macroOrder = {'SOZ_only','PZ_only','NIZ_only', ...
              'SOZ_PZ','SOZ_NIZ','PZ_NIZ','SOZ_PZ_NIZ'};
roleOrder = {'source-like','balanced','sink-like'};

regionColors = [213 94 0; 0 114 178; 0 158 115] / 255;
roleColors   = [204 121 167; 80 80 80; 86 180 233] / 255;
macroColors  = [0.10 0.10 0.10; ...
                0.35 0.35 0.35; ...
                0.60 0.60 0.60; ...
                0.80 0.35 0.20; ...
                0.20 0.45 0.75; ...
                0.35 0.60 0.30; ...
                0.50 0.30 0.70];

metricName = 'eta_abs_burden_mean';
yLabelText = 'eta burden, abs proxy';

%% 3) 绘图
fig = figure('Color','w', 'Position',[100 100 1280 400]);

% -------------------------
% Panel 1: SOZ / PZ / NIZ
% -------------------------
subplot(1,3,1); hold on;
for r = 1:numel(regionNames)
    [y,e] = group_mean_sem(Tregion, 'region', regionNames{r}, metricName, phaseOrder);
    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', regionColors(r,:), ...
        'MarkerFaceColor', regionColors(r,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.8, ...
        'CapSize', 7, ...
        'MarkerSize', 5.5);
end
format_axis(phaseShort);
ylabel(yLabelText, 'FontWeight','bold');
title('Regional eta burden', 'FontWeight','bold', 'Interpreter','none');
leg = legend(regionNames, 'Location','best');
set(leg, 'Box','off', 'Interpreter','none');

% -------------------------
% Panel 2: macro-state combination
% -------------------------
subplot(1,3,2); hold on;
for m = 1:numel(macroOrder)
    [y,e] = group_mean_sem_macro_all_roles(TmacroRole, macroOrder{m}, metricName, phaseOrder, roleOrder);
    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', macroColors(m,:), ...
        'MarkerFaceColor', macroColors(m,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.5, ...
        'CapSize', 6, ...
        'MarkerSize', 5.0);
end
format_axis(phaseShort);
ylabel(yLabelText, 'FontWeight','bold');
title('Macro-state combination', 'FontWeight','bold', 'Interpreter','none');
leg = legend(macroOrder, 'Location','eastoutside');
set(leg, 'Box','off', 'Interpreter','none', 'FontSize',8);

% -------------------------
% Panel 3: source / balanced / sink
% -------------------------
subplot(1,3,3); hold on;
for rr = 1:numel(roleOrder)
    [y,e] = group_mean_sem_role_all_macro(TmacroRole, roleOrder{rr}, metricName, phaseOrder, macroOrder);
    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', roleColors(rr,:), ...
        'MarkerFaceColor', roleColors(rr,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.8, ...
        'CapSize', 7, ...
        'MarkerSize', 5.5);
end
format_axis(phaseShort);
ylabel(yLabelText, 'FontWeight','bold');
title('Source/sink role', 'FontWeight','bold', 'Interpreter','none');
leg = legend(roleOrder, 'Location','best');
set(leg, 'Box','off', 'Interpreter','none');

% 总标题不含 ABC
% annotation(fig, 'textbox', [0 0.94 1 0.05], ...
%     'String', 'Experiment 4: regional and hyperedge-level contributors of eta', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', ...
%     'EdgeColor', 'none', ...
%     'FontName', 'Arial', ...
%     'FontSize', 13, ...
%     'FontWeight', 'bold', ...
%     'Interpreter', 'none');

%% 4) 保存
pngFile = fullfile(outDir, 'FIG_exp4_eta_burden_summary_noABC.png');
pdfFile = fullfile(outDir, 'FIG_exp4_eta_burden_summary_noABC.pdf');
figFile = fullfile(outDir, 'FIG_exp4_eta_burden_summary_noABC.fig');

saveas(fig, figFile);
print(fig, pngFile, '-dpng', '-r600');
print(fig, pdfFile, '-dpdf', '-painters');

fprintf('\n绘图完成：\n%s\n%s\n%s\n', pngFile, pdfFile, figFile);

%% =========================================================
% Local functions
%% =========================================================

function [y,e] = group_mean_sem(T, groupVar, groupValue, metricName, phaseOrder)
    y = nan(1,numel(phaseOrder));
    e = nan(1,numel(phaseOrder));
    for p = 1:numel(phaseOrder)
        idx = strcmp(cellstr(T.(groupVar)), groupValue) & strcmp(cellstr(T.phase5), phaseOrder{p});
        x = T.(metricName)(idx);
        x = x(isfinite(x));
        y(p) = mean_omitnan(x);
        e(p) = sem_omitnan(x);
    end
end

function [y,e] = group_mean_sem_macro_all_roles(T, macroState, metricName, phaseOrder, roleOrder)
    caseIDs = unique(cellstr(T.case_id), 'stable');
    Y = nan(numel(caseIDs), numel(phaseOrder));

    for i = 1:numel(caseIDs)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for r = 1:numel(roleOrder)
                idx = strcmp(cellstr(T.case_id), caseIDs{i}) & ...
                      strcmp(cellstr(T.macro_state), macroState) & ...
                      strcmp(cellstr(T.role), roleOrder{r}) & ...
                      strcmp(cellstr(T.phase5), phaseOrder{p});
                x = T.(metricName)(idx);
                if ~isempty(x) && isfinite(x(1))
                    val = val + x(1);
                    hasAny = true;
                end
            end
            if hasAny
                Y(i,p) = val;
            end
        end
    end

    y = nan(1,numel(phaseOrder));
    e = nan(1,numel(phaseOrder));
    for p = 1:numel(phaseOrder)
        x = Y(:,p);
        x = x(isfinite(x));
        y(p) = mean_omitnan(x);
        e(p) = sem_omitnan(x);
    end
end

function [y,e] = group_mean_sem_role_all_macro(T, roleName, metricName, phaseOrder, macroOrder)
    caseIDs = unique(cellstr(T.case_id), 'stable');
    Y = nan(numel(caseIDs), numel(phaseOrder));

    for i = 1:numel(caseIDs)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for m = 1:numel(macroOrder)
                idx = strcmp(cellstr(T.case_id), caseIDs{i}) & ...
                      strcmp(cellstr(T.role), roleName) & ...
                      strcmp(cellstr(T.macro_state), macroOrder{m}) & ...
                      strcmp(cellstr(T.phase5), phaseOrder{p});
                x = T.(metricName)(idx);
                if ~isempty(x) && isfinite(x(1))
                    val = val + x(1);
                    hasAny = true;
                end
            end
            if hasAny
                Y(i,p) = val;
            end
        end
    end

    y = nan(1,numel(phaseOrder));
    e = nan(1,numel(phaseOrder));
    for p = 1:numel(phaseOrder)
        x = Y(:,p);
        x = x(isfinite(x));
        y(p) = mean_omitnan(x);
        e(p) = sem_omitnan(x);
    end
end

function format_axis(phaseShort)
    xlim([0.75 numel(phaseShort)+0.25]);
    set(gca, 'XTick', 1:numel(phaseShort), 'XTickLabel', phaseShort);
    xtickangle(25);
    set(gca, 'FontName','Arial', ...
        'FontSize',10.5, ...
        'LineWidth',1.15, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top');
    ax = gca;
    ax.YAxis.Exponent = 0;
end

function y = mean_omitnan(x)
    x = x(isfinite(x));
    if isempty(x)
        y = nan;
    else
        y = mean(x);
    end
end

function y = sem_omitnan(x)
    x = x(isfinite(x));
    if numel(x) <= 1
        y = nan;
    else
        y = std(x) / sqrt(numel(x));
    end
end
