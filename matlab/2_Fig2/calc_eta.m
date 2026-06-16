%% calc_eta.m
% Generate Fig. 2 effective refractive-index source data.
%
% Each seizure is assigned one PLV threshold selected from the elbow of its
% density-threshold curve. That seizure-specific threshold is then applied
% to every time window from the same seizure before constructing DMWA and
% calculating eta.
%
% Expected input for each approved local dataset:
%   data/seizure_XX_Gamma.mat, containing variable X1 (channels x samples)
%
% Outputs:
%   *_window_level_eta.csv
%   *_stage5_eta.csv
%   *_eta_results.mat
%   ALL_CASES_window_level_eta.csv
%   ALL_CASES_stage5_eta.csv
%   GROUP_stage5_eta_*.csv

clc; clear; close all;

outDir = fullfile('example_project', '2_Fig2');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fs = 1024;
winLenSec = 3;

caseInfo = build_case_info();
nCase = numel(caseInfo);

allWindowTables = cell(nCase, 1);
allStage5Tables = cell(nCase, 1);

for c = 1:nCase
    case_id = caseInfo(c).case_id;
    fileName = caseInfo(c).file_path;

    fprintf('\n==================================================\n');
    fprintf('(%d/%d) Processing %s\n', c, nCase, case_id);
    fprintf('Input file: %s\n', fileName);
    fprintf('==================================================\n');

    if ~exist(fileName, 'file')
        warning('Input file not found; skipping: %s', fileName);
        continue;
    end

    S0 = load(fileName);
    if ~isfield(S0, 'X1')
        warning('%s does not contain variable X1; skipping.', case_id);
        continue;
    end

    X1 = S0.X1;
    totalIterations = floor(size(X1, 2) / fs) - winLenSec;
    if totalIterations < 0
        warning('%s is too short to form %d-s windows; skipping.', case_id, winLenSec);
        continue;
    end

    nWindows = totalIterations + 1;
    time_axis = (1:nWindows)';
    [pre_idx, ictal_idx, post_idx] = sanitize_stage_indices( ...
        caseInfo(c).pre_idx, caseInfo(c).ictal_idx, caseInfo(c).post_idx, nWindows, case_id);

    if isempty(ictal_idx)
        warning('%s has no valid ictal windows after trimming; skipping.', case_id);
        continue;
    end

    [early_idx, mid_idx, late_idx] = split_ictal_indices(ictal_idx);

    [plvThreshold, plvThresholdCurve] = select_seizure_plv_elbow_threshold( ...
        X1, fs, winLenSec, 0:totalIterations);
    fprintf('%s: seizure-level PLV elbow threshold = %.4f\n', case_id, plvThreshold);

    all_hyperEdges_all = extract_hyperedges_for_windows( ...
        X1, fs, winLenSec, totalIterations, plvThreshold);
    W_all = build_dmwa_for_windows(all_hyperEdges_all);
    [eta_all, num_qi] = compute_eta_for_windows(W_all);

    phase5 = make_phase_labels(nWindows, pre_idx, early_idx, mid_idx, late_idx, post_idx);

    windowTable = table( ...
        repmat({case_id}, nWindows, 1), ...
        repmat({fileName}, nWindows, 1), ...
        time_axis, ...
        phase5, ...
        eta_all, ...
        num_qi, ...
        'VariableNames', {'case_id', 'file_name', 'window_idx', 'phase5', ...
        'eta', 'num_hyperedges'});

    stage5Table = summarize_case_stages(windowTable, case_id);

    allWindowTables{c} = windowTable;
    allStage5Tables{c} = stage5Table;

    writetable(windowTable, fullfile(outDir, [case_id '_window_level_eta.csv']));
    writetable(stage5Table, fullfile(outDir, [case_id '_stage5_eta.csv']));

    save(fullfile(outDir, [case_id '_eta_results.mat']), ...
        'case_id', 'fileName', 'fs', 'winLenSec', ...
        'plvThreshold', 'plvThresholdCurve', ...
        'pre_idx', 'ictal_idx', 'post_idx', ...
        'early_idx', 'mid_idx', 'late_idx', ...
        'all_hyperEdges_all', 'W_all', ...
        'eta_all', 'windowTable', 'stage5Table');

    fprintf('%s: results saved.\n', case_id);
end

validWindow = ~cellfun(@isempty, allWindowTables);
validStage5 = ~cellfun(@isempty, allStage5Tables);

if any(validWindow)
    ALL_windowTable = vertcat(allWindowTables{validWindow});
    writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_level_eta.csv'));
else
    ALL_windowTable = table();
end

if any(validStage5)
    ALL_stage5Table = vertcat(allStage5Tables{validStage5});
    writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_eta.csv'));
else
    ALL_stage5Table = table();
end

[groupDescTable, groupOmnibusTable, groupPairwiseTable] = ...
    summarize_group_stages(ALL_stage5Table, outDir);

if ~isempty(ALL_stage5Table)
    plot_group_eta_boxplot(ALL_stage5Table, outDir);
end

if ~isempty(allWindowTables{1})
    plot_example_eta_timecourse(allWindowTables{1}, caseInfo(1), outDir);
end

save(fullfile(outDir, 'ALL_CASES_eta_group_results.mat'), ...
    'caseInfo', ...
    'ALL_windowTable', ...
    'ALL_stage5Table', ...
    'groupDescTable', ...
    'groupOmnibusTable', ...
    'groupPairwiseTable');

fprintf('\nAll processing completed.\nResults saved to:\n%s\n', outDir);


function caseInfo = build_case_info()
caseIDs = arrayfun(@(ii) sprintf('seizure_%02d', ii), 1:24, 'UniformOutput', false);

preRanges = [
    1 200
    1 122
    1 195
    1 137
    1 252
    1 199
    1 239
    1 236
    1 59
    1 244
    1 150
    1 180
    1 69
    1 252
    1 154
    1 253
    1 135
    1 117
    1 255
    1 311
    1 94
    1 95
    1 162
    1 29
    ];

ictalRanges = [
    201 400
    123 233
    196 330
    138 297
    253 413
    200 285
    240 331
    237 313
    60 112
    245 470
    151 271
    181 261
    70 130
    253 346
    155 230
    254 331
    136 270
    118 300
    256 482
    312 540
    95 170
    96 195
    163 228
    30 95
    ];

postRanges = [
    401 598
    234 338
    331 477
    298 398
    414 498
    286 398
    332 498
    314 498
    113 385
    471 598
    272 404
    262 498
    131 198
    347 448
    231 324
    332 489
    271 309
    301 598
    483 598
    541 598
    171 318
    196 350
    229 310
    96 148
    ];

caseInfo = struct('case_id', {}, 'file_path', {}, 'pre_idx', {}, 'ictal_idx', {}, 'post_idx', {});
for ii = 1:numel(caseIDs)
    caseInfo(ii).case_id = caseIDs{ii};
    caseInfo(ii).file_path = fullfile('data', [caseIDs{ii} '_Gamma.mat']);
    caseInfo(ii).pre_idx = preRanges(ii, 1):preRanges(ii, 2);
    caseInfo(ii).ictal_idx = ictalRanges(ii, 1):ictalRanges(ii, 2);
    caseInfo(ii).post_idx = postRanges(ii, 1):postRanges(ii, 2);
end
end


function [pre_idx, ictal_idx, post_idx] = sanitize_stage_indices(pre_idx, ictal_idx, post_idx, nWindows, case_id)
pre_idx = unique(pre_idx(pre_idx >= 1 & pre_idx <= nWindows), 'stable');
ictal_idx = unique(ictal_idx(ictal_idx >= 1 & ictal_idx <= nWindows), 'stable');
post_idx = unique(post_idx(post_idx >= 1 & post_idx <= nWindows), 'stable');

hasOverlap = ~isempty(intersect(pre_idx, ictal_idx)) || ...
    ~isempty(intersect(pre_idx, post_idx)) || ...
    ~isempty(intersect(ictal_idx, post_idx));

if hasOverlap
    warning('%s has overlapping stage indices; resolving as pre > ictal > post.', case_id);
end

ictal_idx = setdiff(ictal_idx, pre_idx, 'stable');
post_idx = setdiff(post_idx, union(pre_idx, ictal_idx), 'stable');
end


function [early_idx, mid_idx, late_idx] = split_ictal_indices(ictal_idx)
nIctal = numel(ictal_idx);
splitPts = round(linspace(0, nIctal, 4));
early_idx = ictal_idx(1:splitPts(2));
mid_idx = ictal_idx(splitPts(2) + 1:splitPts(3));
late_idx = ictal_idx(splitPts(3) + 1:end);
end


function all_hyperEdges_all = extract_hyperedges_for_windows(X1, fs, winLenSec, totalIterations, plvThreshold)
all_hyperEdges_all = cell(totalIterations + 1, 1);

for time = 0:totalIterations
    start_point = fs * time + 1;
    end_point = fs * (time + winLenSec);
    datanew = X1(:, start_point:end_point);
    all_hyperEdges_all{time + 1} = gain_hyperEdges_23(datanew, plvThreshold);
end
end


function W_all = build_dmwa_for_windows(all_hyperEdges_all)
nWindows = numel(all_hyperEdges_all);
W_all = cell(nWindows, 1);

for tt = 1:nWindows
    W_all{tt} = build_dmwa_from_hyperedges(all_hyperEdges_all{tt});
end
end


function W = build_dmwa_from_hyperedges(hyperEdges)
nHyperedges = numel(hyperEdges);
if nHyperedges < 2
    W = zeros(nHyperedges);
    return;
end

sizes = cellfun(@numel, hyperEdges);
sizes = sizes(:);

overlapMatrix = zeros(nHyperedges, nHyperedges);
for ii = 1:nHyperedges
    for jj = ii + 1:nHyperedges
        overlapLen = numel(intersect(hyperEdges{ii}, hyperEdges{jj}));
        if overlapLen > 0
            overlapMatrix(ii, jj) = overlapLen;
            overlapMatrix(jj, ii) = overlapLen;
        end
    end
end

degree = sum(overlapMatrix > 0, 2);

sizes_i = repmat(sizes, 1, nHyperedges);
sizes_j = repmat(sizes', nHyperedges, 1);
degree_i = repmat(degree, 1, nHyperedges);
degree_j = repmat(degree', nHyperedges, 1);

coverage = overlapMatrix ./ sizes_i;
degreeSum = degree_i + degree_j;
degreeBias = degree_j ./ degreeSum;
degreeBias(degreeSum == 0) = 0.5;

W = coverage;
equalSizeMask = (sizes_i == sizes_j) & (overlapMatrix > 0);
W(equalSizeMask) = W(equalSizeMask) .* degreeBias(equalSizeMask);
W(1:nHyperedges + 1:end) = 0;
end


function [eta_all, num_qi] = compute_eta_for_windows(W_all)
nWindows = numel(W_all);
eta_all = nan(nWindows, 1);
num_qi = nan(nWindows, 1);

for tt = 1:nWindows
    W = double(W_all{tt});
    if isempty(W)
        eta_all(tt) = 0;
        num_qi(tt) = 0;
        continue;
    end

    nHyperedges = size(W, 1);
    if size(W, 2) ~= nHyperedges
        error('Window %d: DMWA matrix is not square.', tt);
    end

    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    W(1:nHyperedges + 1:end) = 0;
    num_qi(tt) = nHyperedges;

    if nHyperedges < 2
        eta_all(tt) = 0;
        continue;
    end

    sum_W = sum(W(:));
    if sum_W <= 0
        eta_all(tt) = 0;
        continue;
    end

    lambda_1 = max(abs(eig(W)));
    if ~isfinite(lambda_1) || lambda_1 <= eps
        eta_all(tt) = 0;
    else
        eta_all(tt) = (sum_W / nHyperedges) / lambda_1;
    end
end
end


function phase5 = make_phase_labels(nWindows, pre_idx, early_idx, mid_idx, late_idx, post_idx)
phase5 = repmat({''}, nWindows, 1);
phase5(pre_idx) = {'pre-ictal'};
phase5(early_idx) = {'early'};
phase5(mid_idx) = {'mid'};
phase5(late_idx) = {'late'};
phase5(post_idx) = {'post-ictal'};
phase5(cellfun(@isempty, phase5)) = {'unused'};
end


function stage5Table = summarize_case_stages(windowTable, case_id)
phaseOrder = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
nPhase = numel(phaseOrder);

phase_col = cell(nPhase, 1);
n_window_col = nan(nPhase, 1);
eta_mean = nan(nPhase, 1);
eta_std = nan(nPhase, 1);
eta_median = nan(nPhase, 1);

for pp = 1:nPhase
    idx = strcmp(windowTable.phase5, phaseOrder{pp});
    x = windowTable.eta(idx);
    x = x(~isnan(x));

    phase_col{pp} = phaseOrder{pp};
    n_window_col(pp) = sum(idx);

    if ~isempty(x)
        eta_mean(pp) = mean(x);
        eta_std(pp) = std(x);
        eta_median(pp) = median(x);
    end
end

stage5Table = table( ...
    repmat({case_id}, nPhase, 1), ...
    phase_col, ...
    n_window_col, ...
    eta_mean, ...
    eta_std, ...
    eta_median, ...
    'VariableNames', {'case_id', 'phase5', 'n_windows', ...
    'eta_mean', 'eta_std', 'eta_median'});
end


function [groupDescTable, groupOmnibusTable, groupPairwiseTable] = summarize_group_stages(ALL_stage5Table, outDir)
groupDescTable = table();
groupOmnibusTable = table();
groupPairwiseTable = table();

if isempty(ALL_stage5Table)
    return;
end

phaseOrder = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
caseIDs = unique(ALL_stage5Table.case_id, 'stable');

desc_phase = {};
desc_n = [];
desc_mean = [];
desc_std = [];
desc_sem = [];
desc_median = [];

for pp = 1:numel(phaseOrder)
    idx = strcmp(ALL_stage5Table.phase5, phaseOrder{pp});
    x = ALL_stage5Table.eta_mean(idx);
    x = x(~isnan(x));

    desc_phase{end + 1, 1} = phaseOrder{pp}; %#ok<AGROW>
    desc_n(end + 1, 1) = numel(x); %#ok<AGROW>

    if isempty(x)
        desc_mean(end + 1, 1) = nan; %#ok<AGROW>
        desc_std(end + 1, 1) = nan; %#ok<AGROW>
        desc_sem(end + 1, 1) = nan; %#ok<AGROW>
        desc_median(end + 1, 1) = nan; %#ok<AGROW>
    else
        desc_mean(end + 1, 1) = mean(x); %#ok<AGROW>
        desc_std(end + 1, 1) = std(x); %#ok<AGROW>
        desc_sem(end + 1, 1) = std(x) / sqrt(numel(x)); %#ok<AGROW>
        desc_median(end + 1, 1) = median(x); %#ok<AGROW>
    end
end

groupDescTable = table(desc_phase, desc_n, desc_mean, desc_std, desc_sem, desc_median, ...
    'VariableNames', {'phase5', 'n_case', 'mean', 'std', 'sem', 'median'});
writetable(groupDescTable, fullfile(outDir, 'GROUP_stage5_eta_descriptive_stats.csv'));

M = nan(numel(caseIDs), numel(phaseOrder));
for ii = 1:numel(caseIDs)
    rows_i = ALL_stage5Table(strcmp(ALL_stage5Table.case_id, caseIDs{ii}), :);
    for pp = 1:numel(phaseOrder)
        tmp = rows_i.eta_mean(strcmp(rows_i.phase5, phaseOrder{pp}));
        if ~isempty(tmp)
            M(ii, pp) = tmp(1);
        end
    end
end

validRows = all(~isnan(M), 2);
Mvalid = M(validRows, :);
if size(Mvalid, 1) >= 2 && exist('friedman', 'file') == 2
    p_friedman = friedman(Mvalid, 1, 'off');
else
    p_friedman = nan;
end

groupOmnibusTable = table(size(Mvalid, 1), p_friedman, ...
    'VariableNames', {'n_case_complete', 'p_friedman'});
writetable(groupOmnibusTable, fullfile(outDir, 'GROUP_stage5_eta_friedman_stats.csv'));

pair_comp = {};
pair_n = [];
pair_p = [];
pair_mean_diff = [];
pair_median_diff = [];

for p1 = 1:numel(phaseOrder) - 1
    for p2 = p1 + 1:numel(phaseOrder)
        x1 = M(:, p1);
        x2 = M(:, p2);
        valid = ~isnan(x1) & ~isnan(x2);
        d = x2(valid) - x1(valid);

        pair_comp{end + 1, 1} = [phaseOrder{p1} '_vs_' phaseOrder{p2}]; %#ok<AGROW>
        pair_n(end + 1, 1) = sum(valid); %#ok<AGROW>

        if sum(valid) >= 2 && exist('signrank', 'file') == 2
            pair_p(end + 1, 1) = signrank(x1(valid), x2(valid)); %#ok<AGROW>
        else
            pair_p(end + 1, 1) = nan; %#ok<AGROW>
        end

        if isempty(d)
            pair_mean_diff(end + 1, 1) = nan; %#ok<AGROW>
            pair_median_diff(end + 1, 1) = nan; %#ok<AGROW>
        else
            pair_mean_diff(end + 1, 1) = mean(d); %#ok<AGROW>
            pair_median_diff(end + 1, 1) = median(d); %#ok<AGROW>
        end
    end
end

groupPairwiseTable = table(pair_comp, pair_n, pair_p, pair_mean_diff, pair_median_diff, ...
    'VariableNames', {'comparison', 'n_pair', 'p_signrank', 'mean_diff', 'median_diff'});
writetable(groupPairwiseTable, fullfile(outDir, 'GROUP_stage5_eta_pairwise_signrank.csv'));
end


function plot_group_eta_boxplot(ALL_stage5Table, outDir)
phaseOrder = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
phaseColors = [
    76, 120, 168
    89, 161, 79
    242, 142, 43
    225, 87, 89
    128, 115, 172
    ] / 255;

Y = cell(numel(phaseOrder), 1);
maxLen = 0;
for pp = 1:numel(phaseOrder)
    idx = strcmp(ALL_stage5Table.phase5, phaseOrder{pp});
    y = ALL_stage5Table.eta_mean(idx);
    y = y(~isnan(y));
    Y{pp} = y;
    maxLen = max(maxLen, numel(y));
end

if maxLen == 0
    return;
end

Mplot = nan(maxLen, numel(phaseOrder));
for pp = 1:numel(phaseOrder)
    if ~isempty(Y{pp})
        Mplot(1:numel(Y{pp}), pp) = Y{pp};
    end
end

fig = figure('Color', 'w', 'Position', [100 100 520 420]);
hold on;

boxplot(Mplot, 'Colors', [0.15 0.15 0.15], 'Symbol', '', ...
    'Widths', 0.55, 'MedianStyle', 'line');

set(findobj(gca, 'Tag', 'Box'), 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
set(findobj(gca, 'Tag', 'Median'), 'Color', [0 0 0], 'LineWidth', 1.4);
set(findobj(gca, 'Tag', 'Whisker'), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
set(findobj(gca, 'Tag', 'Adjacent Value'), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

for pp = 1:numel(phaseOrder)
    y = Y{pp};
    if isempty(y)
        continue;
    end
    jitter = (rand(numel(y), 1) - 0.5) * 0.18;
    scatter(pp + jitter, y, 26, ...
        'MarkerFaceColor', phaseColors(pp, :), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.45, ...
        'MarkerFaceAlpha', 0.9, ...
        'MarkerEdgeAlpha', 0.95);
end

set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
    'Box', 'off', 'TickDir', 'out', 'Layer', 'top');
xticks(1:numel(phaseOrder));
xticklabels(phaseOrder);
xtickangle(25);
ylabel('\eta', 'FontWeight', 'bold');

allY = vertcat(Y{:});
if ~isempty(allY)
    pad = max(0.08 * (max(allY) - min(allY)), 0.01);
    ylim([min(allY) - pad, max(allY) + pad]);
    ax = gca;
    ax.YAxis.Exponent = 0;
end

exportgraphics(fig, fullfile(outDir, 'GROUP_stage5_eta_boxplot.png'), 'Resolution', 600);
savefig(fig, fullfile(outDir, 'GROUP_stage5_eta_boxplot.fig'));
close(fig);
end


function plot_example_eta_timecourse(windowTable, caseMeta, outDir)
time_axis = windowTable.window_idx;
eta_all = windowTable.eta;

fig = figure('Color', 'w', 'Position', [100, 100, 820, 300]);
ax = axes(fig);
hold(ax, 'on');
plot(ax, time_axis, eta_all, 'k-', 'LineWidth', 1.6);

ictalStart = min(caseMeta.ictal_idx);
ictalEnd = max(caseMeta.ictal_idx);
xline(ax, ictalStart, '--', 'LineWidth', 1.1, 'Color', [0.45 0.45 0.45]);
xline(ax, ictalEnd, '--', 'LineWidth', 1.1, 'Color', [0.45 0.45 0.45]);

xlabel(ax, 'Time windows', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, '\eta', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

set(ax, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.015], 'Layer', 'top');
xlim(ax, [time_axis(1), time_axis(end)]);

validEta = eta_all(~isnan(eta_all) & isfinite(eta_all));
if ~isempty(validEta)
    pad = max(0.08 * (max(validEta) - min(validEta)), 0.01);
    ylim(ax, [min(validEta) - pad, max(validEta) + pad]);
end
ax.YAxis.Exponent = 0;

yl = ylim(ax);
y_text = yl(2) - 0.04 * (yl(2) - yl(1));
text(mean(caseMeta.pre_idx), y_text, 'pre-ictal', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text(mean(caseMeta.ictal_idx), y_text, 'ictal', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text(mean(caseMeta.post_idx), y_text, 'post-ictal', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(outDir, 'seizure_01_eta_all_timecourse_nature_labeled.png'), 'Resolution', 600);
savefig(fig, fullfile(outDir, 'seizure_01_eta_all_timecourse_nature_labeled.fig'));
close(fig);
end
