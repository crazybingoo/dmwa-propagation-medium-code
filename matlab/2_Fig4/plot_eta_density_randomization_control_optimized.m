%%
% plot_eta_density_randomization_control_optimized.m
%
% 独立绘图脚本：读取实验 2 已保存的结果，只生成一张优化后的结果图。
%
% 依赖文件（二选一，优先读取 MAT）：
%   1) ALL_CASES_eta_density_randomization_control_results.mat
%   2) ALL_CASES_stage5_eta_control.csv
%
% 输出：
%   FIG_eta_control_one_optimized.png
%   FIG_eta_control_one_optimized.pdf
%   FIG_eta_control_one_optimized.fig
%
% 使用方法：
%   1) 先运行 calc_eta_density_randomization_control.m 得到结果文件
%   2) 修改下面 resultsDir 为实验 2 输出目录
%   3) 运行本脚本

clc; clear; close all;

%% 1) 路径设置
resultsDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig4';

% 如果当前机器没有该路径，则默认从当前工作目录读取，方便拷贝到结果文件夹后直接运行
if ~isfolder(resultsDir)
    warning('resultsDir 不存在，改用当前目录：%s', pwd);
    resultsDir = pwd;
end

outBase = fullfile(resultsDir, 'FIG_eta_control_one_optimized');

%% 2) 读取阶段级结果
ALL_stage5Table = load_stage5_results(resultsDir);

requiredVars = {'case_id','model','phase5','eta_mean'};
missingVars = setdiff(requiredVars, ALL_stage5Table.Properties.VariableNames);
if ~isempty(missingVars)
    error('ALL_stage5Table 缺少必要字段：%s', strjoin(missingVars, ', '));
end

ALL_stage5Table = normalize_table_text_columns(ALL_stage5Table, {'case_id','model','phase5'});

%% 3) 绘图参数
phaseOrder = ["pre-ictal","early","mid","late","post-ictal"];
phaseShort = ["Pre","Early","Mid","Late","Post"];

modelOrder = ["Original","DensityMatched","WeightShuffled","TopologyShuffled"];
modelLabel = ["Original","Density matched","Weight shuffled","Topology shuffled"];

% 色彩采用高对比、论文友好的配色
modelColors = [
    0.05 0.05 0.05
    0.16 0.42 0.72
    0.82 0.32 0.20
    0.24 0.56 0.28
];

lineStyles = ["-","--","-.","-"];
markerList = ["o","s","^","d"];

%% 4) 汇总 eta 与 delta eta
[etaMean, etaSem, etaN] = summarize_model_phase( ...
    ALL_stage5Table, modelOrder, phaseOrder, "eta_mean");

[deltaMean, deltaSem, deltaN] = summarize_delta_relative_to_pre( ...
    ALL_stage5Table, modelOrder, phaseOrder, "eta_mean");

fprintf('\n绘图使用病例数概览：\n');
for mm = 1:numel(modelOrder)
    fprintf('  %-18s eta n = [%s], delta n = [%s]\n', ...
        modelOrder(mm), ...
        strjoin(string(etaN(mm,:)), ', '), ...
        strjoin(string(deltaN(mm,:)), ', '));
end

%% 5) 生成一张整合图
fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[3 3 18 16], ...
    'Renderer','painters');

tl = tiledlayout(fig, 2, 1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% Panel A: eta 阶段变化
ax1 = nexttile(tl, 1);
hold(ax1, 'on');

for mm = 1:numel(modelOrder)
    h = errorbar(ax1, ...
        1:numel(phaseOrder), etaMean(mm,:), etaSem(mm,:), ...
        'LineStyle', lineStyles(mm), ...
        'Marker', markerList(mm), ...
        'Color', modelColors(mm,:), ...
        'MarkerFaceColor', modelColors(mm,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.9, ...
        'CapSize', 7, ...
        'MarkerSize', 5.8);
    h.DisplayName = modelLabel(mm);
end

format_eta_axis(ax1, phaseShort);
ylabel(ax1, '\eta', 'FontWeight','bold');
% title(ax1, 'A  Stage-wise \eta under density and randomization controls', ...
%     'FontWeight','bold', 'HorizontalAlignment','left');

lgd = legend(ax1, 'Location','northoutside', ...
    'Orientation','horizontal', ...
    'Box','off');
lgd.NumColumns = 2;

%% Panel B: delta eta relative to pre
ax2 = nexttile(tl, 2);
hold(ax2, 'on');

xDelta = 1:(numel(phaseOrder)-1);

for mm = 1:numel(modelOrder)
    h = errorbar(ax2, ...
        xDelta, deltaMean(mm,2:end), deltaSem(mm,2:end), ...
        'LineStyle', lineStyles(mm), ...
        'Marker', markerList(mm), ...
        'Color', modelColors(mm,:), ...
        'MarkerFaceColor', modelColors(mm,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.9, ...
        'CapSize', 7, ...
        'MarkerSize', 5.8);
    h.DisplayName = modelLabel(mm);
end

yline(ax2, 0, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);

format_delta_axis(ax2, phaseShort(2:end));
ylabel(ax2, '\Delta\eta relative to Pre', 'FontWeight','bold');
% title(ax2, 'B  Phase-related \Delta\eta after pre-ictal normalization', ...
%     'FontWeight','bold', 'HorizontalAlignment','left');

%% 6) 细节优化：统一 y 轴留白
pad_y_axis(ax1, etaMean, etaSem);
pad_y_axis(ax2, deltaMean(:,2:end), deltaSem(:,2:end));

% 保证标题左对齐在较老 MATLAB 中也尽量稳定
ax1.Title.Units = 'normalized';
ax1.Title.Position(1) = 0;
ax2.Title.Units = 'normalized';
ax2.Title.Position(1) = 0;

%% 7) 保存
try
    exportgraphics(fig, [outBase '.png'], 'Resolution', 600);
    exportgraphics(fig, [outBase '.pdf'], 'ContentType','vector');
catch
    warning('exportgraphics 不可用，改用 print 保存。');
    print(fig, [outBase '.png'], '-dpng', '-r600');
    print(fig, [outBase '.pdf'], '-dpdf', '-painters');
end

savefig(fig, [outBase '.fig']);

fprintf('\n优化绘图完成，已保存：\n');
fprintf('  %s.png\n', outBase);
fprintf('  %s.pdf\n', outBase);
fprintf('  %s.fig\n', outBase);

%% =========================================================
% Local functions
%% =========================================================

function T = load_stage5_results(resultsDir)

    matPath = fullfile(resultsDir, 'ALL_CASES_eta_density_randomization_control_results.mat');
    csvPath = fullfile(resultsDir, 'ALL_CASES_stage5_eta_control.csv');

    if exist(matPath, 'file') == 2
        S = load(matPath, 'ALL_stage5Table');

        if isfield(S, 'ALL_stage5Table') && ~isempty(S.ALL_stage5Table)
            T = S.ALL_stage5Table;
            fprintf('已读取 MAT：%s\n', matPath);
            return;
        else
            warning('MAT 中没有可用的 ALL_stage5Table，尝试读取 CSV。');
        end
    end

    if exist(csvPath, 'file') == 2
        T = readtable(csvPath);
        fprintf('已读取 CSV：%s\n', csvPath);
        return;
    end

    error(['未找到阶段级结果文件。请确认 resultsDir 中存在：\n' ...
           '  %s\n或\n  %s'], matPath, csvPath);
end

function T = normalize_table_text_columns(T, colNames)

    for i = 1:numel(colNames)
        col = char(colNames{i});

        if ~ismember(col, T.Properties.VariableNames)
            continue;
        end

        if iscell(T.(col))
            T.(col) = string(T.(col));
        elseif iscategorical(T.(col))
            T.(col) = string(T.(col));
        elseif ischar(T.(col))
            T.(col) = string(cellstr(T.(col)));
        end
    end
end

function [Y, E, N] = summarize_model_phase(T, modelOrder, phaseOrder, metricName)

    Y = nan(numel(modelOrder), numel(phaseOrder));
    E = nan(numel(modelOrder), numel(phaseOrder));
    N = zeros(numel(modelOrder), numel(phaseOrder));

    for mm = 1:numel(modelOrder)
        for pp = 1:numel(phaseOrder)
            idx = T.model == modelOrder(mm) & T.phase5 == phaseOrder(pp);

            x = T.(metricName)(idx);
            x = x(isfinite(x));

            [Y(mm,pp), E(mm,pp), N(mm,pp)] = mean_sem(x);
        end
    end
end

function [Dmean, Dsem, DN] = summarize_delta_relative_to_pre(T, modelOrder, phaseOrder, metricName)

    caseIDs = unique(T.case_id, 'stable');

    Dmean = nan(numel(modelOrder), numel(phaseOrder));
    Dsem  = nan(numel(modelOrder), numel(phaseOrder));
    DN    = zeros(numel(modelOrder), numel(phaseOrder));

    for mm = 1:numel(modelOrder)

        deltaMat = nan(numel(caseIDs), numel(phaseOrder));

        for cc = 1:numel(caseIDs)

            idxPre = T.case_id == caseIDs(cc) & ...
                     T.model == modelOrder(mm) & ...
                     T.phase5 == phaseOrder(1);

            xPre = T.(metricName)(idxPre);
            xPre = xPre(isfinite(xPre));

            if isempty(xPre)
                continue;
            end

            xPre = xPre(1);

            for pp = 1:numel(phaseOrder)
                idxCmp = T.case_id == caseIDs(cc) & ...
                         T.model == modelOrder(mm) & ...
                         T.phase5 == phaseOrder(pp);

                xCmp = T.(metricName)(idxCmp);
                xCmp = xCmp(isfinite(xCmp));

                if ~isempty(xCmp)
                    deltaMat(cc,pp) = xCmp(1) - xPre;
                end
            end
        end

        for pp = 1:numel(phaseOrder)
            x = deltaMat(:,pp);
            x = x(isfinite(x));
            [Dmean(mm,pp), Dsem(mm,pp), DN(mm,pp)] = mean_sem(x);
        end
    end
end

function [m, sem, n] = mean_sem(x)

    x = x(:);
    x = x(isfinite(x));
    n = numel(x);

    if n == 0
        m = nan;
        sem = nan;
    elseif n == 1
        m = x(1);
        sem = 0;
    else
        m = mean(x);
        sem = std(x, 0) / sqrt(n);
    end
end

function format_eta_axis(ax, phaseShort)

    set(ax, ...
        'FontName','Arial', ...
        'FontSize',10.5, ...
        'LineWidth',1.1, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top', ...
        'XGrid','off', ...
        'YGrid','on', ...
        'GridLineStyle',':', ...
        'GridAlpha',0.22);

    xlim(ax, [0.75 numel(phaseShort)+0.25]);
    xticks(ax, 1:numel(phaseShort));
    xticklabels(ax, phaseShort);
    ax.YAxis.Exponent = 0;
end

function format_delta_axis(ax, phaseShort)

    set(ax, ...
        'FontName','Arial', ...
        'FontSize',10.5, ...
        'LineWidth',1.1, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top', ...
        'XGrid','off', ...
        'YGrid','on', ...
        'GridLineStyle',':', ...
        'GridAlpha',0.22);

    xlim(ax, [0.75 numel(phaseShort)+0.25]);
    xticks(ax, 1:numel(phaseShort));
    xticklabels(ax, phaseShort);
    ax.YAxis.Exponent = 0;
end

function pad_y_axis(ax, Y, E)

    vals = [Y(:) - E(:); Y(:) + E(:)];
    vals = vals(isfinite(vals));

    if isempty(vals)
        return;
    end

    ymin = min(vals);
    ymax = max(vals);

    if ymin == ymax
        pad = max(abs(ymin), 1) * 0.08;
    else
        pad = (ymax - ymin) * 0.14;
    end

    ylim(ax, [ymin - pad, ymax + pad]);
end
