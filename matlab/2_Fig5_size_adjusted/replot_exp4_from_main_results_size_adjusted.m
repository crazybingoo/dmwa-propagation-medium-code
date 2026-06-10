%% replot_exp4_from_main_results_size_adjusted.m
% =========================================================
% 目的：
%   不封装为函数，直接调用主脚本已保存的结果 MAT，
%   使用主脚本原有绘图代码重新绘制总结图。
%
% 使用方式：
%   1) 先运行主脚本，确保结果 MAT 已生成：
%      ALL_CASES_exp4_region_hyperedge_contribution_results_size_adjusted.mat
%   2) 修改 outDir 为你的结果目录
%   3) 直接运行本脚本
% =========================================================

clc; clear; close all;

%% 1) 路径
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig5_size_adjusted';

resultMat = fullfile(outDir, 'ALL_CASES_exp4_region_hyperedge_contribution_results_size_adjusted.mat');
if ~exist(resultMat, 'file')
    error('未找到结果 MAT 文件：%s\n请先运行主脚本生成结果。', resultMat);
end

%% 2) 读取主脚本结果
S = load(resultMat);

requiredVars = {'ALL_stageRegionTable','ALL_stageMacroRoleTable'};
for i = 1:numel(requiredVars)
    if ~isfield(S, requiredVars{i})
        error('结果 MAT 中缺少变量：%s', requiredVars{i});
    end
end

ALL_stageRegionTable = S.ALL_stageRegionTable;
ALL_stageMacroRoleTable = S.ALL_stageMacroRoleTable;

%% 3) 与主脚本一致的参数
phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
phaseShort = {'Pre','Early','Mid','Late','Post'};
regionNames = {'SOZ','PZ','NIZ'};
macroOrder = {'SOZ_only','PZ_only','NIZ_only', ...
              'SOZ_PZ','SOZ_NIZ','PZ_NIZ','SOZ_PZ_NIZ'};
roleOrder = {'source-like','balanced','sink-like'};

regionColors = [213 94 0; 0 114 178; 0 158 115] / 255;
roleColors = [204 121 167; 80 80 80; 86 180 233] / 255;
macroColors = [0.10 0.10 0.10; 0.35 0.35 0.35; 0.60 0.60 0.60; ...
               0.80 0.35 0.20; 0.20 0.45 0.75; 0.35 0.60 0.30; 0.50 0.30 0.70];

%% 4) 作图：保持主脚本原有逻辑
fig = figure('Color','w', 'Position',[100 100 1280 380]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

%% A: SOZ/PZ/NIZ eta_abs_burden_per_node（扣除区域大小差异）
nexttile; hold on;
for r = 1:numel(regionNames)
    y = nan(1, numel(phaseOrder));
    e = nan(1, numel(phaseOrder));

    for p = 1:numel(phaseOrder)
        idx = strcmp(ALL_stageRegionTable.region, regionNames{r}) & ...
              strcmp(ALL_stageRegionTable.phase5, phaseOrder{p});
        x = ALL_stageRegionTable.eta_abs_burden_per_node_mean(idx);
        x = x(isfinite(x));
        y(p) = mean(x, 'omitnan');
        e(p) = std(x, 'omitnan') / sqrt(max(numel(x), 1));
    end

    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', regionColors(r,:), 'MarkerFaceColor', regionColors(r,:), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.8, 'CapSize', 7, 'MarkerSize', 5.5);
end
xlim([0.75 numel(phaseShort)+0.25]);
xticks(1:numel(phaseShort));
xticklabels(phaseShort);
xtickangle(25);
set(gca, 'FontName','Arial', 'FontSize',10.5, 'LineWidth',1.15, ...
    'Box','off', 'TickDir','out', 'Layer','top');
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('eta burden per node', 'FontWeight','bold');
title('Regional eta burden (size-adjusted)', 'FontWeight','bold');
legend(regionNames, 'Box','off', 'Location','best');

%% B: 宏观状态组合 eta_abs_burden，合并 role 后统计
nexttile; hold on;
caseIDs_macro = unique(ALL_stageMacroRoleTable.case_id, 'stable');

for m = 1:numel(macroOrder)
    Y = nan(numel(caseIDs_macro), numel(phaseOrder));

    for i = 1:numel(caseIDs_macro)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for r = 1:numel(roleOrder)
                idx = strcmp(ALL_stageMacroRoleTable.case_id, caseIDs_macro{i}) & ...
                      strcmp(ALL_stageMacroRoleTable.macro_state, macroOrder{m}) & ...
                      strcmp(ALL_stageMacroRoleTable.role, roleOrder{r}) & ...
                      strcmp(ALL_stageMacroRoleTable.phase5, phaseOrder{p});
                x = ALL_stageMacroRoleTable.eta_abs_burden_mean(idx);
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

    y = mean(Y, 1, 'omitnan');
    e = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));

    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', macroColors(m,:), 'MarkerFaceColor', macroColors(m,:), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.5, 'CapSize', 6, 'MarkerSize', 5.0);
end
xlim([0.75 numel(phaseShort)+0.25]);
xticks(1:numel(phaseShort));
xticklabels(phaseShort);
xtickangle(25);
set(gca, 'FontName','Arial', 'FontSize',10.5, 'LineWidth',1.15, ...
    'Box','off', 'TickDir','out', 'Layer','top');
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('eta burden, abs proxy', 'FontWeight','bold');
title('Macro-state combination', 'FontWeight','bold');
% legend(macroOrder, 'Box','off', 'Location','eastoutside', 'FontSize',8);
lgd = legend(macroOrder, 'Box','off', 'Location','eastoutside', 'FontSize',8);
set(lgd, 'Interpreter', 'none');

%% C: source-like / balanced / sink-like eta_abs_burden，合并 macro 后统计
nexttile; hold on;
caseIDs_role = unique(ALL_stageMacroRoleTable.case_id, 'stable');

for rr = 1:numel(roleOrder)
    Y = nan(numel(caseIDs_role), numel(phaseOrder));

    for i = 1:numel(caseIDs_role)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for m = 1:numel(macroOrder)
                idx = strcmp(ALL_stageMacroRoleTable.case_id, caseIDs_role{i}) & ...
                      strcmp(ALL_stageMacroRoleTable.role, roleOrder{rr}) & ...
                      strcmp(ALL_stageMacroRoleTable.macro_state, macroOrder{m}) & ...
                      strcmp(ALL_stageMacroRoleTable.phase5, phaseOrder{p});
                x = ALL_stageMacroRoleTable.eta_abs_burden_mean(idx);
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

    y = mean(Y, 1, 'omitnan');
    e = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));

    errorbar(1:numel(phaseOrder), y, e, '-o', ...
        'Color', roleColors(rr,:), 'MarkerFaceColor', roleColors(rr,:), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.8, 'CapSize', 7, 'MarkerSize', 5.5);
end
xlim([0.75 numel(phaseShort)+0.25]);
xticks(1:numel(phaseShort));
xticklabels(phaseShort);
xtickangle(25);
set(gca, 'FontName','Arial', 'FontSize',10.5, 'LineWidth',1.15, ...
    'Box','off', 'TickDir','out', 'Layer','top');
ax = gca;
ax.YAxis.Exponent = 0;
ylabel('eta burden, abs proxy', 'FontWeight','bold');
title('Source/sink role', 'FontWeight','bold');
legend(roleOrder, 'Box','off', 'Location','best');

%% 5) 导出
exportgraphics(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted_replot.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted_replot.pdf'), 'ContentType','vector');
savefig(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted_replot.fig'));

fprintf('\n重绘完成。\n输出目录：%s\n', outDir);
