%% =========================================================
% 绘制 eta 指标对比实验结果（多模块横向拼接，Nature/Science 风格）
% 数据源: ALL_CASES_stage5_eta_Comparison.csv
% =========================================================
clc; clear; close all;

%% 1. 路径与数据加载
dataDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig8';
dataFile = fullfile(dataDir, 'ALL_CASES_stage5_eta_Comparison.csv');

if ~exist(dataFile, 'file')
    error('找不到文件: %s\n请确认之前的计算脚本已正确生成该文件。', dataFile);
end

T = readtable(dataFile);

phaseOrder  = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
phaseShort  = {'Pre', 'Early', 'Mid', 'Late', 'Post'};
modelOrder  = {'Low-Order(2-nodes)', 'High-Order(3-nodes)', 'DMW-HLG'};
modelTitles = {'Low-Order', 'High-Order', 'DMW-HLG'};

phaseColors = [76,  120, 168;
               89,  161, 79;
               242, 142, 43;
               225, 87,  89;
               128, 115, 172] / 255;

cases  = unique(T.case_id, 'stable');
nCase  = numel(cases);
nPhase = numel(phaseOrder);
nModel = numel(modelOrder);

%% 2. 坐标布局与数据提取
blockStarts  = [1, 6.5, 12];
blockCenters = blockStarts + 2;

allPositions = [];
groupCounter = 0;

dataVec  = [];
groupVec = [];

globalMin = inf;
globalMax = -inf;

Xmethod = cell(1, nModel);

for m = 1:nModel
    X = nan(nCase, nPhase);

    for p = 1:nPhase
        groupCounter = groupCounter + 1;
        allPositions(groupCounter) = blockStarts(m) + (p - 1); %#ok<SAGROW>

        idx = strcmp(T.Model, modelOrder{m}) & strcmp(T.phase5, phaseOrder{p});
        y   = T.eta_mean(idx);
        y   = y(~isnan(y));

        if ~isempty(y)
            globalMin = min(globalMin, min(y));
            globalMax = max(globalMax, max(y));
        end

        dataVec  = [dataVec; y(:)]; %#ok<AGROW>
        groupVec = [groupVec; repmat(groupCounter, numel(y), 1)]; %#ok<AGROW>

        for i = 1:nCase
            idxCase = idx & strcmp(T.case_id, cases{i});
            vals    = T.eta_mean(idxCase);
            vals    = vals(~isnan(vals));
            if ~isempty(vals)
                X(i, p) = vals(1);
            end
        end
    end
    Xmethod{m} = X;
end

if ~isfinite(globalMin) || ~isfinite(globalMax)
    error('没有可用数据，请检查 CSV 文件内容。');
end

%% 3. 画图初始化
fig = figure('Color', 'w', 'Position', [40, 80, 1100, 650]);
ax  = axes('Position', [0.06, 0.12, 0.92, 0.78]);
hold(ax, 'on');

yr = globalMax - globalMin;
if yr <= 0
    yr = 0.05;
end

ylim(ax, [globalMin - 0.08 * yr, globalMax + 0.75 * yr]);
xlim(ax, [0, 16.5]);
yl = ylim(ax);

%% 4. 绘制箱线图与散点
boxplot(ax, dataVec, groupVec, ...
    'Positions',   allPositions, ...
    'Symbol',      '', ...
    'Widths',      0.55, ...
    'Colors',      [0.18 0.18 0.18], ...
    'MedianStyle', 'line');

set(findobj(ax, 'Tag', 'Box'),            'Color', [0.20 0.20 0.20], 'LineWidth', 1.1);
set(findobj(ax, 'Tag', 'Median'),         'Color', [0 0 0],          'LineWidth', 1.3);
set(findobj(ax, 'Tag', 'Whisker'),        'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
set(findobj(ax, 'Tag', 'Adjacent Value'), 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);

rng(2025);
for m = 1:nModel
    for p = 1:nPhase
        x0  = blockStarts(m) + (p - 1);
        idx = strcmp(T.Model, modelOrder{m}) & strcmp(T.phase5, phaseOrder{p});
        y   = T.eta_mean(idx);
        y   = y(~isnan(y));

        if isempty(y)
            continue;
        end

        jitter = (rand(numel(y), 1) - 0.5) * 0.32;

        scatter(ax, x0 + jitter, y, 22, ...
            'MarkerFaceColor', phaseColors(p, :), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth',       0.40, ...
            'MarkerFaceAlpha', 0.92, ...
            'MarkerEdgeAlpha', 0.95);
    end
end

%% 5. 绘制显著性连线
pairIdx_adj = [1 2; 2 3; 3 4; 4 5];
pairIdx_key = [1 3; 3 5];

baseY = globalMax + 0.05 * yr;
stepY = 0.09 * yr;

for m = 1:nModel
    X    = Xmethod{m};
    xpos = blockStarts(m) + (0:4);

    for k = 1:size(pairIdx_adj, 1)
        i1 = pairIdx_adj(k, 1);
        i2 = pairIdx_adj(k, 2);

        [pval, dval] = paired_stats_pd(X(:, i1), X(:, i2));
        label = sprintf('%s\n%s', format_p_value(pval), format_d_value(dval));

        yLine = baseY + (k - 1) * stepY;
        add_sig_line_with_label(ax, xpos(i1), xpos(i2), yLine, label, 0.9, 8.5);
    end

    for k = 1:size(pairIdx_key, 1)
        i1 = pairIdx_key(k, 1);
        i2 = pairIdx_key(k, 2);

        [pval, dval] = paired_stats_pd(X(:, i1), X(:, i2));
        label = sprintf('%s\n%s', format_p_value(pval), format_d_value(dval));

        yLine = baseY + (size(pairIdx_adj, 1) + k - 1) * stepY;
        add_sig_line_with_label(ax, xpos(i1), xpos(i2), yLine, label, 1.1, 8.5);
    end
end

%% 6. 坐标轴与图例
set(ax, ...
    'FontName',   'Arial', ...
    'FontSize',   12, ...
    'LineWidth',  1.2, ...
    'Box',        'off', ...
    'TickDir',    'out', ...
    'TickLength', [0.012 0.012], ...
    'Layer',      'top');

xticks(ax, blockCenters);
xticklabels(ax, modelTitles);
ylabel(ax, '\eta', 'FontWeight', 'bold', 'FontSize', 13);

for m = 1:nModel
    plot(ax, [blockStarts(m) - 0.3, blockStarts(m) + 4.3], ...
         [yl(1) + 0.002 * yr, yl(1) + 0.002 * yr], ...
         'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
end

hLegends = gobjects(1, nPhase);
for p = 1:nPhase
    hLegends(p) = scatter(ax, nan, nan, 50, ...
        'MarkerFaceColor', phaseColors(p, :), ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
end
leg = legend(ax, hLegends, phaseShort, ...
    'Location',    'northoutside', ...
    'Orientation', 'horizontal', ...
    'Box',         'on', ...
    'FontSize',    11);
leg.Position(2) = 0.94;

%% 7. 导出
pngFile = fullfile(dataDir, 'GROUP_stage5_eta_Comparison_blocks.png');
pdfFile = fullfile(dataDir, 'GROUP_stage5_eta_Comparison_blocks.pdf');

exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

fprintf('\n✅ eta 多模块分布图绘制完成！\n已保存至:\n%s\n%s\n', pngFile, pdfFile);

%% 8. 局部函数
function [pval, dval] = paired_stats_pd(x1, x2)
    x1 = x1(:); x2 = x2(:);
    mask = ~(isnan(x1) | isnan(x2));
    x1 = x1(mask); x2 = x2(mask);

    if numel(x1) < 2
        pval = NaN; dval = NaN; return;
    end

    pval = signrank(x1, x2);

    d = x2 - x1;
    sd_d = std(d, 'omitnan');
    if sd_d == 0 || isnan(sd_d)
        dval = 0;
    else
        dval = mean(d, 'omitnan') / sd_d;
    end
end

function add_sig_line_with_label(ax, x1, x2, y, label, lw, fs)
    yr = diff(ylim(ax));
    h  = 0.016 * yr;

    plot(ax, [x1, x1, x2, x2], [y, y+h, y+h, y], ...
        'Color', [0.15 0.15 0.15], 'LineWidth', lw);

    text(ax, mean([x1, x2]), y + 1.15 * h, label, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontName',            'Arial', ...
        'FontSize',            fs);
end

function txt = format_p_value(p)
    if isnan(p)
        txt = 'p=NaN';
    elseif p < 0.001
        txt = 'p<0.001';
    elseif p < 0.01
        txt = sprintf('p=%.3f', p);
    else
        txt = sprintf('p=%.2f', p);
    end
end

function txt = format_d_value(d)
    if isnan(d)
        txt = 'd=NaN';
    else
        txt = sprintf('d=%.2f', d);
    end
end
