%%
clc; clear; close all;

%% =========================================================
% plot_eta_decomposition_summary.m
%
% 目的：
%   调用实验 1 生成的 ALL_CASES_stage5_eta_decomposition.csv
%   将 eta 分解指标放在一张综合图中比较。
%
% 输入：
%   1) ALL_CASES_stage5_eta_decomposition.csv
%   2) GROUP_eta_decomposition_all_stats.csv，可选
%
% 输出：
%   1) FIG_eta_decomposition_summary.png
%   2) FIG_eta_decomposition_summary.fig
%   3) FIG_eta_core_raw_eta_R_lambda1.png
%   4) FIG_eta_core_raw_eta_R_lambda1.fig
% =========================================================

%% 1) 路径设置
resultDir = 'example_project\2_Fig3';

stageFile = fullfile(resultDir, 'ALL_CASES_stage5_eta_decomposition.csv');
statFile  = fullfile(resultDir, 'GROUP_eta_decomposition_all_stats.csv');

if ~exist(stageFile, 'file')
    error('找不到文件: %s', stageFile);
end

T = readtable(stageFile);

if exist(statFile, 'file')
    Stats = readtable(statFile);
else
    Stats = table();
end

%% 2) 基本设置
phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
phaseShort = {'Pre','Early','Mid','Late','Post'};
nPhase = numel(phaseOrder);

caseIDs = unique(T.case_id, 'stable');
nCase = numel(caseIDs);

phaseColors = [
    76,120,168;
    89,161,79;
    242,142,43;
    225,87,89;
    128,115,172
    ] / 255;

metricNames = { ...
    'eta_mean', ...
    'resource_R_mean', ...
    'lambda1_mean', ...
    'spectral_gap_mean', ...
    'lambda2_ratio_mean', ...
    'spectral_PR_mean'};

metricLabels = { ...
    '\eta', ...
    'R', ...
    '\lambda_1', ...
    '\lambda_1-\lambda_2', ...
    '\lambda_2/\lambda_1', ...
    'Spectral PR'};

coreMetricNames = { ...
    'eta_mean', ...
    'resource_R_mean', ...
    'lambda1_mean'};

coreMetricLabels = { ...
    '\eta', ...
    'R = \Sigma W/N', ...
    '\lambda_1'};

%% 3) 整理成 case × phase × metric 矩阵
M = struct();

for m = 1:numel(metricNames)
    metric = metricNames{m};

    if ~ismember(metric, T.Properties.VariableNames)
        error('结果表中缺少指标: %s', metric);
    end

    X = nan(nCase, nPhase);

    for i = 1:nCase
        for p = 1:nPhase
            idx = strcmp(T.case_id, caseIDs{i}) & strcmp(T.phase5, phaseOrder{p});
            tmp = T.(metric)(idx);

            if ~isempty(tmp)
                X(i,p) = tmp(1);
            end
        end
    end

    M.(metric) = X;
end

%% 4) 病人内 z-score
% 目的：
%   不同指标量纲差异很大，直接放在一张图会看不清。
%   因此对每个病人的每个指标，在 5 个阶段内做 z-score。
%
%   X_z(i,p) = (X(i,p) - mean_i) / std_i
%
% 这样展示的是“相对 pre / ictal / post 的阶段变化模式”，不是原始绝对值。

Z = struct();

for m = 1:numel(metricNames)
    metric = metricNames{m};
    X = M.(metric);

    Xz = nan(size(X));

    for i = 1:nCase
        xi = X(i,:);
        valid = ~isnan(xi) & isfinite(xi);

        if sum(valid) >= 2
            mu = mean(xi(valid));
            sd = std(xi(valid));

            if sd > 0
                Xz(i,valid) = (xi(valid) - mu) ./ sd;
            else
                Xz(i,valid) = 0;
            end
        end
    end

    Z.(metric) = Xz;
end

%% 5) 计算均值和 SEM
meanZ = nan(numel(metricNames), nPhase);
semZ  = nan(numel(metricNames), nPhase);

for m = 1:numel(metricNames)
    metric = metricNames{m};
    Xz = Z.(metric);

    for p = 1:nPhase
        x = Xz(:,p);
        x = x(~isnan(x) & isfinite(x));

        meanZ(m,p) = mean(x);
        semZ(m,p)  = std(x) / sqrt(numel(x));
    end
end

meanRaw = struct();
semRaw = struct();

for m = 1:numel(coreMetricNames)
    metric = coreMetricNames{m};
    X = M.(metric);

    meanRaw.(metric) = nan(1,nPhase);
    semRaw.(metric) = nan(1,nPhase);

    for p = 1:nPhase
        x = X(:,p);
        x = x(~isnan(x) & isfinite(x));

        meanRaw.(metric)(p) = mean(x);
        semRaw.(metric)(p) = std(x) / sqrt(numel(x));
    end
end

%% 6) 综合图：标准化趋势 + 核心原始指标
fig = figure('Color','w', 'Position',[80 80 1350 780]);

tiledlayout(2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% 6.1 标准化趋势总览
nexttile([1 3]);
hold on;

lineColors = lines(numel(metricNames));

for m = 1:numel(metricNames)
    y = meanZ(m,:);
    e = semZ(m,:);

    errorbar( ...
        1:nPhase, y, e, ...
        '-o', ...
        'Color', lineColors(m,:), ...
        'MarkerFaceColor', lineColors(m,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.8, ...
        'CapSize', 8, ...
        'MarkerSize', 6);
end

yline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.1);

xlim([0.75 nPhase+0.25]);
xticks(1:nPhase);
xticklabels(phaseShort);

ylabel('Within-patient z-score', 'FontWeight','bold');
title('Stage-wise standardized trajectories of \eta decomposition metrics', ...
    'FontWeight','bold');

legend(metricLabels, ...
    'Location','eastoutside', ...
    'Box','off');

set(gca, ...
    'FontName','Arial', ...
    'FontSize',11, ...
    'LineWidth',1.1, ...
    'Box','off', ...
    'TickDir','out', ...
    'Layer','top');

%% 6.2 eta 原始值
for m = 1:numel(coreMetricNames)

    metric = coreMetricNames{m};

    nexttile;
    hold on;

    X = M.(metric);

    % individual thin lines
    for i = 1:nCase
        xi = X(i,:);
        if all(isnan(xi))
            continue;
        end

        plot(1:nPhase, xi, '-', ...
            'Color', [0.78 0.78 0.78], ...
            'LineWidth', 0.8);
    end

    % mean ± SEM
    y = meanRaw.(metric);
    e = semRaw.(metric);

    errorbar( ...
        1:nPhase, y, e, ...
        '-o', ...
        'Color', [0.05 0.05 0.05], ...
        'MarkerFaceColor', [0.05 0.05 0.05], ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 2.2, ...
        'CapSize', 9, ...
        'MarkerSize', 7);

    % scatter individual points
    for p = 1:nPhase
        x = X(:,p);
        x = x(~isnan(x) & isfinite(x));

        jitter = (rand(numel(x),1) - 0.5) * 0.16;

        scatter( ...
            p + jitter, x, 22, ...
            'MarkerFaceColor', phaseColors(p,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.4, ...
            'MarkerFaceAlpha', 0.85, ...
            'MarkerEdgeAlpha', 0.9);
    end

    xlim([0.75 nPhase+0.25]);
    xticks(1:nPhase);
    xticklabels(phaseShort);

    ylabel(coreMetricLabels{m}, 'FontWeight','bold');
    title(coreMetricLabels{m}, 'FontWeight','bold');

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',11, ...
        'LineWidth',1.1, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top');

    ax = gca;
    ax.YAxis.Exponent = 0;
end

exportgraphics(fig, fullfile(resultDir, 'FIG_eta_decomposition_summary.png'), 'Resolution', 600);
savefig(fig, fullfile(resultDir, 'FIG_eta_decomposition_summary.fig'));

%% 7) 单独输出核心三联图：eta / R / lambda1
fig2 = figure('Color','w', 'Position',[100 100 1150 360]);
tiledlayout(1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for m = 1:numel(coreMetricNames)

    metric = coreMetricNames{m};
    X = M.(metric);

    nexttile;
    hold on;

    for i = 1:nCase
        xi = X(i,:);
        if all(isnan(xi))
            continue;
        end

        plot(1:nPhase, xi, '-', ...
            'Color', [0.78 0.78 0.78], ...
            'LineWidth', 0.8);
    end

    y = meanRaw.(metric);
    e = semRaw.(metric);

    errorbar( ...
        1:nPhase, y, e, ...
        '-o', ...
        'Color', [0.05 0.05 0.05], ...
        'MarkerFaceColor', [0.05 0.05 0.05], ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 2.2, ...
        'CapSize', 9, ...
        'MarkerSize', 7);

    for p = 1:nPhase
        x = X(:,p);
        x = x(~isnan(x) & isfinite(x));

        jitter = (rand(numel(x),1) - 0.5) * 0.16;

        scatter( ...
            p + jitter, x, 22, ...
            'MarkerFaceColor', phaseColors(p,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.4, ...
            'MarkerFaceAlpha', 0.85, ...
            'MarkerEdgeAlpha', 0.9);
    end

    xlim([0.75 nPhase+0.25]);
    xticks(1:nPhase);
    xticklabels(phaseShort);
    xtickangle(25);

    ylabel(coreMetricLabels{m}, 'FontWeight','bold');
    title(coreMetricLabels{m}, 'FontWeight','bold');

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',11, ...
        'LineWidth',1.1, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top');

    ax = gca;
    ax.YAxis.Exponent = 0;
end

exportgraphics(fig2, fullfile(resultDir, 'FIG_eta_core_raw_eta_R_lambda1.png'), 'Resolution', 600);
savefig(fig2, fullfile(resultDir, 'FIG_eta_core_raw_eta_R_lambda1.fig'));

%% 8) 输出一个简短趋势表，方便写论文结果
summaryRows = table();

for m = 1:numel(metricNames)
    metric = metricNames{m};
    X = M.(metric);

    pre  = X(:,1);
    early = X(:,2);
    mid  = X(:,3);
    late = X(:,4);
    post = X(:,5);

    valid_mid = ~isnan(pre) & ~isnan(mid);

    if sum(valid_mid) > 0
        mid_minus_pre = mid(valid_mid) - pre(valid_mid);
        mean_mid_minus_pre = mean(mid_minus_pre);
        median_mid_minus_pre = median(mid_minus_pre);
    else
        mean_mid_minus_pre = nan;
        median_mid_minus_pre = nan;
    end

    y = nan(1,nPhase);
    for p = 1:nPhase
        x = X(:,p);
        x = x(~isnan(x) & isfinite(x));
        y(p) = mean(x);
    end

    [~, maxPhaseIdx] = max(y);
    [~, minPhaseIdx] = min(y);

    tmp = table( ...
        {metric}, ...
        y(1), y(2), y(3), y(4), y(5), ...
        {phaseOrder{maxPhaseIdx}}, ...
        {phaseOrder{minPhaseIdx}}, ...
        mean_mid_minus_pre, ...
        median_mid_minus_pre, ...
        'VariableNames', { ...
        'metric', ...
        'mean_pre', ...
        'mean_early', ...
        'mean_mid', ...
        'mean_late', ...
        'mean_post', ...
        'max_phase', ...
        'min_phase', ...
        'mean_mid_minus_pre', ...
        'median_mid_minus_pre'});

    summaryRows = [summaryRows; tmp];
end

writetable(summaryRows, fullfile(resultDir, 'SUMMARY_eta_decomposition_trend_table.csv'));

fprintf('\n绘图完成。\n');
fprintf('综合图已保存:\n%s\n', fullfile(resultDir, 'FIG_eta_decomposition_summary.png'));
fprintf('核心三联图已保存:\n%s\n', fullfile(resultDir, 'FIG_eta_core_raw_eta_R_lambda1.png'));
fprintf('趋势表已保存:\n%s\n', fullfile(resultDir, 'SUMMARY_eta_decomposition_trend_table.csv'));