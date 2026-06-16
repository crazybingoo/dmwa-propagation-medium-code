clc; clear; close all;

%% =========================================================
% 单独作图脚本：基于 calc_eta.m 已经输出的结果文件进行作图
% 功能：
%   1) 绘制 5 个阶段的 Nature 风格箱线图 + 病人级散点
%   2) 显示指定 6 组比较的 p 值（配对 Wilcoxon signrank）
%   3) 显示指定 6 组比较的效应量 d 值（配对 Cohen''s d_z）
%
% 依赖输入：
%   ALL_CASES_stage5_eta.csv
%   （该文件由 calc_eta.m 运行后生成）
%
% 比较关系：
%   pre-ictal vs early
%   early vs mid
%   mid vs late
%   late vs post-ictal
%   pre-ictal vs mid
%   mid vs post-ictal
%
% 输出文件：
%   GROUP_stage5_eta_boxplot_with_p_d.png
%   GROUP_stage5_eta_boxplot_with_p_d.fig
%   GROUP_stage5_eta_pairwise_with_d.csv
%% =========================================================

%% 0) 路径设置
% 改成你的结果目录。按你截图中的目录默认填写：
resultDir = 'example_project\2_Fig1';
inFile    = fullfile(resultDir, 'ALL_CASES_stage5_eta.csv');

if ~exist(inFile, 'file')
    error('未找到结果文件: %s\n请先运行 calc_eta.m，或检查 resultDir 是否正确。', inFile);
end

%% 1) 读取 stage5 结果
T = readtable(inFile);

reqVars = {'case_id','phase5','eta_mean'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, T.Properties.VariableNames)
        error('结果表缺少必要变量: %s', reqVars{i});
    end
end

if iscell(T.case_id), T.case_id = string(T.case_id); end
if ischar(T.case_id), T.case_id = string(cellstr(T.case_id)); end
if iscell(T.phase5),  T.phase5  = string(T.phase5);  end
if ischar(T.phase5),  T.phase5  = string(cellstr(T.phase5));  end

phaseOrder  = ["pre-ictal","early","mid","late","post-ictal"];
phaseLabels = {'pre-ictal','early','mid','late','post-ictal'};
caseIDs = unique(T.case_id, 'stable');

%% 2) 组装 病人 × 阶段 矩阵（使用 eta_mean）
M = nan(numel(caseIDs), numel(phaseOrder));
for i = 1:numel(caseIDs)
    rows_i = T(T.case_id == caseIDs(i), :);
    for p = 1:numel(phaseOrder)
        idx = (rows_i.phase5 == phaseOrder(p));
        tmp = rows_i.eta_mean(idx);
        if ~isempty(tmp)
            M(i,p) = tmp(1);
        end
    end
end

%% 3) 指定比较组：四条相邻阶段 + pre-mid + mid-post
compPairs = [
    1 2;   % pre-ictal vs early
    2 3;   % early vs mid
    3 4;   % mid vs late
    4 5;   % late vs post-ictal
    1 3;   % pre-ictal vs mid
    3 5    % mid vs post-ictal
];

compNames = {
    'pre-ictal_vs_early';
    'early_vs_mid';
    'mid_vs_late';
    'late_vs_post-ictal';
    'pre-ictal_vs_mid';
    'mid_vs_post-ictal'
};

nComp = size(compPairs, 1);
pvals = nan(nComp, 1);
dvals = nan(nComp, 1);
nPair = nan(nComp, 1);
meanDiff = nan(nComp, 1);
medianDiff = nan(nComp, 1);

for k = 1:nComp
    g1 = compPairs(k,1);
    g2 = compPairs(k,2);

    x1 = M(:, g1);
    x2 = M(:, g2);
    valid = ~isnan(x1) & ~isnan(x2);

    nPair(k) = sum(valid);
    if sum(valid) >= 2
        x1v = x1(valid);
        x2v = x2(valid);
        diffVals = x2v - x1v;

        if exist('signrank', 'file') == 2
            pvals(k) = signrank(x1v, x2v);
        else
            pvals(k) = nan;
        end

        meanDiff(k)   = mean(diffVals);
        medianDiff(k) = median(diffVals);
        dvals(k)      = paired_cohens_dz(diffVals);
    end
end

pairwiseStats = table(compNames, nPair, pvals, dvals, meanDiff, medianDiff, ...
    'VariableNames', {'comparison','n_pair','p_signrank','cohens_dz','mean_diff','median_diff'});

writetable(pairwiseStats, fullfile(resultDir, 'GROUP_stage5_eta_pairwise_with_d.csv'));

%% 4) 作图数据
Y = cell(numel(phaseOrder), 1);
for p = 1:numel(phaseOrder)
    y = M(:,p);
    y = y(~isnan(y));
    Y{p} = y;
end

phaseColors = [
    76,120,168;
    89,161,79;
    242,142,43;
    225,87,89;
    128,115,172
    ] / 255;

fig = figure('Color','w', 'Position',[100 100 620 500]);
hold on;

boxplot(M, ...
    'Colors', [0.15 0.15 0.15], ...
    'Symbol', '', ...
    'Widths', 0.55, ...
    'MedianStyle', 'line');

hBox = findobj(gca, 'Tag', 'Box');
set(hBox, 'Color', [0.20 0.20 0.20], 'LineWidth', 1.2);

hMedian = findobj(gca, 'Tag', 'Median');
set(hMedian, 'Color', [0 0 0], 'LineWidth', 1.4);

hWhisker = findobj(gca, 'Tag', 'Whisker');
set(hWhisker, 'Color', [0.30 0.30 0.30], 'LineWidth', 1.0);

hAdj = findobj(gca, 'Tag', 'Adjacent Value');
set(hAdj, 'Color', [0.30 0.30 0.30], 'LineWidth', 1.0);

rng(1);
for p = 1:numel(phaseOrder)
    y = Y{p};
    if isempty(y)
        continue;
    end

    jitter = (rand(numel(y),1) - 0.5) * 0.18;
    scatter(p + jitter, y, 28, ...
        'MarkerFaceColor', phaseColors(p,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.45, ...
        'MarkerFaceAlpha', 0.90, ...
        'MarkerEdgeAlpha', 0.95);
end

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'Layer', 'top');

xticks(1:numel(phaseOrder));
xticklabels(phaseLabels);
xtickangle(25);
ylabel('\eta', 'FontWeight', 'bold');

%% 5) y 轴范围
Xall = M(:);
Xall = Xall(~isnan(Xall));

if isempty(Xall)
    yMin = 0;
    yMax = 1;
else
    yMin = min(Xall);
    yMax = max(Xall);
end

if yMin == yMax
    yRange = max(abs(yMin) * 0.10, 0.05);
else
    yRange = yMax - yMin;
end

%% 6) 显著性括号自动分层，避免重叠
levels = zeros(nComp,1);
for k = 1:nComp
    a1 = compPairs(k,1);
    b1 = compPairs(k,2);

    lvl = 1;
    while true
        conflict = false;
        for j = 1:k-1
            if levels(j) ~= lvl
                continue;
            end
            a2 = compPairs(j,1);
            b2 = compPairs(j,2);

            % 只要区间重叠或接触，就上移一层
            if ~(b1 < a2 || b2 < a1)
                conflict = true;
                break;
            end
        end

        if ~conflict
            levels(k) = lvl;
            break;
        else
            lvl = lvl + 1;
        end
    end
end

%% 7) 绘制括号 + p值 + d值
lineStep = 0.14 * yRange;
capH     = 0.03 * yRange;
textGap  = 0.015 * yRange;
yBase    = yMax + 0.08 * yRange;

for k = 1:nComp
    x1 = compPairs(k,1);
    x2 = compPairs(k,2);
    y  = yBase + (levels(k)-1) * lineStep;

    plot([x1 x1 x2 x2], [y y+capH y+capH y], '-', ...
        'Color', [0.1 0.1 0.1], 'LineWidth', 1.1);

    pStr = format_p_value(pvals(k));
    dStr = format_d_value(dvals(k));
    txt  = sprintf('%s\n%s', pStr, dStr);

    text(mean([x1 x2]), y + capH + textGap, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', ...
        'FontSize', 9.5, ...
        'Color', [0.1 0.1 0.1]);
end

yTop = yBase + max(levels) * lineStep + capH + 0.14 * yRange;
yBottom = yMin - 0.08 * yRange;
ylim([yBottom, yTop]);

ax = gca;
ax.YAxis.Exponent = 0;

%% 8) 保存输出
pngFile = fullfile(resultDir, 'GROUP_stage5_eta_boxplot_with_p_d.png');
figFile = fullfile(resultDir, 'GROUP_stage5_eta_boxplot_with_p_d.fig');

exportgraphics(fig, pngFile, 'Resolution', 600);
savefig(fig, figFile);

fprintf('\n========================================\n');
fprintf('作图完成。\n');
fprintf('读取文件: %s\n', inFile);
fprintf('输出图片: %s\n', pngFile);
fprintf('输出FIG : %s\n', figFile);
fprintf('统计表  : %s\n', fullfile(resultDir, 'GROUP_stage5_eta_pairwise_with_d.csv'));
fprintf('========================================\n');

%% ===================== 局部函数 =====================
function dz = paired_cohens_dz(diffVals)
% 配对设计 Cohen''s d_z = mean(diff) / std(diff)
    diffVals = diffVals(~isnan(diffVals) & isfinite(diffVals));

    if numel(diffVals) < 2
        dz = nan;
        return;
    end

    sd = std(diffVals, 0);
    md = mean(diffVals);

    if sd == 0
        if md == 0
            dz = 0;
        else
            dz = sign(md) * inf;
        end
    else
        dz = md / sd;
    end
end

function s = format_p_value(p)
    if isnan(p)
        s = 'p = NA';
    elseif p < 0.001
        s = 'p < 0.001';
    else
        s = sprintf('p = %.3g', p);
    end
end

function s = format_d_value(d)
    if isnan(d)
        s = 'd = NA';
    elseif isinf(d)
        if d > 0
            s = 'd = +Inf';
        else
            s = 'd = -Inf';
        end
    else
        s = sprintf('d = %.2f', d);
    end
end
