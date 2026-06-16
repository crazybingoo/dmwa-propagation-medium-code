%%
clc; clear; close all;

%% =========================================================
% 实验 2：密度与随机化控制
%
% 目的：
%   检验 eta 的阶段变化是否只是由以下因素造成：
%   1) 网络密度变化
%   2) 权重分布变化
%   3) 随机拓扑结构变化
%
% 输入：
%   实验 1 输出的每个病例 MAT：
%       *_eta_decomposition_results.mat
%
%   其中需要包含：
%       W_all
%       case_id
%       pre_idx
%       early_idx
%       mid_idx
%       late_idx
%       post_idx
%
% 输出：
%   1) ALL_CASES_window_level_eta_control.csv
%   2) ALL_CASES_stage5_eta_control.csv
%   3) GROUP_eta_control_descriptive_stats.csv
%   4) GROUP_eta_control_phase_stats.csv
%   5) GROUP_eta_control_delta_original_vs_controls.csv
%   6) FIG_eta_control_by_model.png
%   7) FIG_delta_eta_control_by_model.png
%
% 模型：
%   Original:
%       原始 W
%
%   DensityMatched:
%       病例内统一二值密度，保留最强 K 条非零连接
%
%   WeightShuffled:
%       保留非零拓扑位置，随机打乱非零权重
%
%   TopologyShuffled:
%       保留非零权重分布和连接数量，随机放置到非对角位置
%
% 注意：
%   随机模型需要多次置换。nRand 越大越稳，但越慢。
%   初步建议 nRand = 30；最终论文可设为 100。
% =========================================================

%% 1) 路径设置
inputDir = 'example_project\2_Fig3';
outDir   = 'example_project\2_Fig4';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

matFiles = dir(fullfile(inputDir, '*_eta_decomposition_results.mat'));

if isempty(matFiles)
    error('在 inputDir 中没有找到 *_eta_decomposition_results.mat。请先运行实验 1。');
end

%% 2) 参数
rng(2025);

nRand = 30;               % 初步 30，论文最终建议 100
doFullSpectrum = false;   % false: 只计算 eta/R/lambda1，速度快
                           % true : 同时计算 lambda2 和 N_eff_lambda，但会慢很多

densityTargetMode = 'minPositive';
% 'minPositive'：每个病例内取最小非零二值密度，最严格，但可能较稀疏
% 'quantile'   ：每个病例内取 densityQuantile 分位数，更稳健但少数窗可能无法完全匹配

densityQuantile = 0.10;

phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
phaseShort = {'Pre','Early','Mid','Late','Post'};

modelOrder = {'Original','DensityMatched','WeightShuffled','TopologyShuffled'};

metricList = { ...
    'eta', ...
    'resource_R', ...
    'lambda1', ...
    'sum_W', ...
    'density_W', ...
    'binary_density', ...
    'num_hyperedges'};

if doFullSpectrum
    metricList = [metricList, {'lambda2','lambda2_ratio','N_eff_lambda'}];
end

%% 3) 汇总容器
allWindowTables = cell(numel(matFiles), 1);
caseSummaryTables = cell(numel(matFiles), 1);

%% =========================================================
% 主循环：逐病例处理
%% =========================================================
for f = 1:numel(matFiles)

    matPath = fullfile(matFiles(f).folder, matFiles(f).name);

    fprintf('\n==================================================\n');
    fprintf('(%d/%d) 读取: %s\n', f, numel(matFiles), matFiles(f).name);
    fprintf('==================================================\n');

    S = load(matPath);

    if ~isfield(S, 'W_all')
        warning('%s 中没有 W_all，跳过。', matFiles(f).name);
        continue;
    end

    W_all = S.W_all;
    T = numel(W_all);

    if isfield(S, 'case_id')
        case_id = S.case_id;
    else
        [~, case_id] = fileparts(matFiles(f).name);
        case_id = strrep(case_id, '_eta_decomposition_results', '');
    end

    requiredIdx = {'pre_idx','early_idx','mid_idx','late_idx','post_idx'};
    for rr = 1:numel(requiredIdx)
        if ~isfield(S, requiredIdx{rr})
            error('%s 缺少变量 %s。', matFiles(f).name, requiredIdx{rr});
        end
    end

    pre_idx   = S.pre_idx(:)';
    early_idx = S.early_idx(:)';
    mid_idx   = S.mid_idx(:)';
    late_idx  = S.late_idx(:)';
    post_idx  = S.post_idx(:)';

    %% 3.1 阶段标签
    phase5 = repmat({''}, T, 1);

    pre_idx   = pre_idx(pre_idx >= 1 & pre_idx <= T);
    early_idx = early_idx(early_idx >= 1 & early_idx <= T);
    mid_idx   = mid_idx(mid_idx >= 1 & mid_idx <= T);
    late_idx  = late_idx(late_idx >= 1 & late_idx <= T);
    post_idx  = post_idx(post_idx >= 1 & post_idx <= T);

    phase5(pre_idx)   = {'pre-ictal'};
    phase5(early_idx) = {'early'};
    phase5(mid_idx)   = {'mid'};
    phase5(late_idx)  = {'late'};
    phase5(post_idx)  = {'post-ictal'};

    emptyMask = cellfun(@isempty, phase5);
    phase5(emptyMask) = {'unused'};

    %% 3.2 计算该病例的原始二值密度，用于 DensityMatched
    origBinaryDensity = nan(T,1);
    origNumEdges = nan(T,1);
    origN = nan(T,1);

    for t = 1:T
        W = sanitize_W(W_all{t});
        n = size(W,1);

        origN(t) = n;

        if n >= 2
            K = nnz(W > 0);
            origNumEdges(t) = K;
            origBinaryDensity(t) = K / (n * (n - 1));
        else
            origNumEdges(t) = 0;
            origBinaryDensity(t) = 0;
        end
    end

    validDens = origBinaryDensity(origBinaryDensity > 0 & isfinite(origBinaryDensity));

    if isempty(validDens)
        targetDensity = 0;
    else
        switch densityTargetMode
            case 'minPositive'
                targetDensity = min(validDens);

            case 'quantile'
                targetDensity = quantile(validDens, densityQuantile);

            otherwise
                error('未知 densityTargetMode: %s', densityTargetMode);
        end
    end

    fprintf('%s: targetDensity = %.6f, nRand = %d\n', case_id, targetDensity, nRand);

    caseSummaryTables{f} = table( ...
        {case_id}, ...
        T, ...
        targetDensity, ...
        mean(origBinaryDensity, 'omitnan'), ...
        median(origBinaryDensity, 'omitnan'), ...
        min(validDens), ...
        max(validDens), ...
        'VariableNames', { ...
        'case_id', ...
        'n_windows', ...
        'target_binary_density', ...
        'mean_original_binary_density', ...
        'median_original_binary_density', ...
        'min_positive_original_binary_density', ...
        'max_original_binary_density'});

    %% 3.3 逐窗计算四类模型
    windowRows = table();

    tic_case = tic;

    for t = 1:T

        W0 = sanitize_W(W_all{t});

        %% Original
        M0 = compute_eta_metrics(W0, doFullSpectrum);

        row0 = make_window_row( ...
            case_id, t, phase5{t}, 'Original', ...
            targetDensity, nRand, 0, ...
            M0, []);

        windowRows = [windowRows; row0];

        %% DensityMatched
        W_density = density_match_topK(W0, targetDensity);
        Md = compute_eta_metrics(W_density, doFullSpectrum);

        rowD = make_window_row( ...
            case_id, t, phase5{t}, 'DensityMatched', ...
            targetDensity, nRand, 0, ...
            Md, []);

        windowRows = [windowRows; rowD];

        %% WeightShuffled
        randMetrics_weight = cell(nRand,1);

        for r = 1:nRand
            Wr = weight_shuffle_keep_topology(W0);
            randMetrics_weight{r} = compute_eta_metrics(Wr, doFullSpectrum);
        end

        [Mw_mean, Mw_sd] = average_metric_structs(randMetrics_weight, doFullSpectrum);

        rowW = make_window_row( ...
            case_id, t, phase5{t}, 'WeightShuffled', ...
            targetDensity, nRand, 1, ...
            Mw_mean, Mw_sd);

        windowRows = [windowRows; rowW];

        %% TopologyShuffled
        randMetrics_topo = cell(nRand,1);

        for r = 1:nRand
            Wr = topology_shuffle_keep_weights(W0);
            randMetrics_topo{r} = compute_eta_metrics(Wr, doFullSpectrum);
        end

        [Mt_mean, Mt_sd] = average_metric_structs(randMetrics_topo, doFullSpectrum);

        rowT = make_window_row( ...
            case_id, t, phase5{t}, 'TopologyShuffled', ...
            targetDensity, nRand, 1, ...
            Mt_mean, Mt_sd);

        windowRows = [windowRows; rowT];

        if mod(t, 50) == 0 || t == T
            fprintf('%s: window %d/%d 完成\n', case_id, t, T);
        end
    end

    fprintf('%s: 控制实验完成，用时 %.2f 秒。\n', case_id, toc(tic_case));

    allWindowTables{f} = windowRows;

    writetable(windowRows, fullfile(outDir, [case_id '_window_level_eta_control.csv']));
end

%% 4) 合并所有病例窗级结果
validWindow = ~cellfun(@isempty, allWindowTables);

if any(validWindow)
    ALL_windowTable = vertcat(allWindowTables{validWindow});
else
    ALL_windowTable = table();
end

writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_level_eta_control.csv'));

validSummary = ~cellfun(@isempty, caseSummaryTables);

if any(validSummary)
    CASE_density_summary = vertcat(caseSummaryTables{validSummary});
else
    CASE_density_summary = table();
end

writetable(CASE_density_summary, fullfile(outDir, 'CASE_density_target_summary.csv'));

%% 5) 阶段级结果
ALL_stage5Table = table();

if ~isempty(ALL_windowTable)

    caseIDs = unique(ALL_windowTable.case_id, 'stable');

    for i = 1:numel(caseIDs)
        for mm = 1:numel(modelOrder)
            for p = 1:numel(phaseOrder)

                idx = strcmp(ALL_windowTable.case_id, caseIDs{i}) & ...
                      strcmp(ALL_windowTable.model, modelOrder{mm}) & ...
                      strcmp(ALL_windowTable.phase5, phaseOrder{p});

                nWin = sum(idx);

                tmpRow = table( ...
                    {caseIDs{i}}, ...
                    {modelOrder{mm}}, ...
                    {phaseOrder{p}}, ...
                    nWin, ...
                    'VariableNames', { ...
                    'case_id', ...
                    'model', ...
                    'phase5', ...
                    'n_windows'});

                for k = 1:numel(metricList)
                    metric = metricList{k};

                    x = ALL_windowTable.(metric)(idx);
                    x = x(~isnan(x) & isfinite(x));

                    if isempty(x)
                        tmpRow.([metric '_mean']) = nan;
                        tmpRow.([metric '_std']) = nan;
                        tmpRow.([metric '_median']) = nan;
                    else
                        tmpRow.([metric '_mean']) = mean(x);
                        tmpRow.([metric '_std']) = std(x);
                        tmpRow.([metric '_median']) = median(x);
                    end
                end

                ALL_stage5Table = [ALL_stage5Table; tmpRow];
            end
        end
    end
end

writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_eta_control.csv'));

%% 6) 组水平统计
groupDescTable = table();
groupPhaseStats = table();
groupDeltaStats = table();

if ~isempty(ALL_stage5Table)

    caseIDs = unique(ALL_stage5Table.case_id, 'stable');

    stageMetricList = {};
    for k = 1:numel(metricList)
        stageMetricList{end+1} = [metricList{k} '_mean'];
    end

    %% 6.1 描述统计：model × phase × metric
    for mm = 1:numel(modelOrder)
        modelName = modelOrder{mm};

        for k = 1:numel(stageMetricList)
            metricName = stageMetricList{k};

            for p = 1:numel(phaseOrder)

                idx = strcmp(ALL_stage5Table.model, modelName) & ...
                      strcmp(ALL_stage5Table.phase5, phaseOrder{p});

                x = ALL_stage5Table.(metricName)(idx);
                x = x(~isnan(x) & isfinite(x));

                if isempty(x)
                    n_case = 0;
                    mean_x = nan;
                    std_x = nan;
                    sem_x = nan;
                    median_x = nan;
                else
                    n_case = numel(x);
                    mean_x = mean(x);
                    std_x = std(x);
                    sem_x = std(x) / sqrt(numel(x));
                    median_x = median(x);
                end

                tmp = table( ...
                    {modelName}, ...
                    {metricName}, ...
                    {phaseOrder{p}}, ...
                    n_case, ...
                    mean_x, ...
                    std_x, ...
                    sem_x, ...
                    median_x, ...
                    'VariableNames', { ...
                    'model', ...
                    'metric', ...
                    'phase5', ...
                    'n_case', ...
                    'mean', ...
                    'std', ...
                    'sem', ...
                    'median'});

                groupDescTable = [groupDescTable; tmp];
            end
        end
    end

    writetable(groupDescTable, fullfile(outDir, 'GROUP_eta_control_descriptive_stats.csv'));

    %% 6.2 每个模型内部：Friedman + pre vs 其他阶段 signrank
    for mm = 1:numel(modelOrder)
        modelName = modelOrder{mm};

        for k = 1:numel(stageMetricList)
            metricName = stageMetricList{k};

            M = nan(numel(caseIDs), numel(phaseOrder));

            for i = 1:numel(caseIDs)
                for p = 1:numel(phaseOrder)

                    idx = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                          strcmp(ALL_stage5Table.model, modelName) & ...
                          strcmp(ALL_stage5Table.phase5, phaseOrder{p});

                    tmp = ALL_stage5Table.(metricName)(idx);

                    if ~isempty(tmp)
                        M(i,p) = tmp(1);
                    end
                end
            end

            validRows = all(~isnan(M), 2);
            Mvalid = M(validRows, :);

            if size(Mvalid,1) >= 2 && exist('friedman', 'file') == 2
                p_friedman = friedman(Mvalid, 1, 'off');
            else
                p_friedman = nan;
            end

            for p = 2:numel(phaseOrder)

                x_pre = M(:,1);
                x_cmp = M(:,p);
                valid = ~isnan(x_pre) & ~isnan(x_cmp);

                if sum(valid) >= 2 && exist('signrank', 'file') == 2
                    p_signrank = signrank(x_pre(valid), x_cmp(valid));
                else
                    p_signrank = nan;
                end

                d = x_cmp(valid) - x_pre(valid);

                if isempty(d)
                    mean_diff = nan;
                    median_diff = nan;
                else
                    mean_diff = mean(d);
                    median_diff = median(d);
                end

                tmp = table( ...
                    {modelName}, ...
                    {metricName}, ...
                    {['pre-ictal_vs_' phaseOrder{p}]}, ...
                    sum(valid), ...
                    p_friedman, ...
                    p_signrank, ...
                    mean_diff, ...
                    median_diff, ...
                    'VariableNames', { ...
                    'model', ...
                    'metric', ...
                    'comparison', ...
                    'n_pair', ...
                    'p_friedman_all_phase', ...
                    'p_signrank_pre_vs_phase', ...
                    'mean_diff_phase_minus_pre', ...
                    'median_diff_phase_minus_pre'});

                groupPhaseStats = [groupPhaseStats; tmp];
            end
        end
    end

    writetable(groupPhaseStats, fullfile(outDir, 'GROUP_eta_control_phase_stats.csv'));

    %% 6.3 关键检验：Original 的 phase-pre 变化是否大于或不同于控制模型
    controlModels = {'DensityMatched','WeightShuffled','TopologyShuffled'};

    for k = 1:numel(stageMetricList)

        metricName = stageMetricList{k};

        for cm = 1:numel(controlModels)

            controlName = controlModels{cm};

            for p = 2:numel(phaseOrder)

                deltaOrig = nan(numel(caseIDs),1);
                deltaCtrl = nan(numel(caseIDs),1);

                for i = 1:numel(caseIDs)

                    idxOpre = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                              strcmp(ALL_stage5Table.model, 'Original') & ...
                              strcmp(ALL_stage5Table.phase5, 'pre-ictal');

                    idxOcmp = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                              strcmp(ALL_stage5Table.model, 'Original') & ...
                              strcmp(ALL_stage5Table.phase5, phaseOrder{p});

                    idxCpre = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                              strcmp(ALL_stage5Table.model, controlName) & ...
                              strcmp(ALL_stage5Table.phase5, 'pre-ictal');

                    idxCcmp = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                              strcmp(ALL_stage5Table.model, controlName) & ...
                              strcmp(ALL_stage5Table.phase5, phaseOrder{p});

                    xOpre = ALL_stage5Table.(metricName)(idxOpre);
                    xOcmp = ALL_stage5Table.(metricName)(idxOcmp);
                    xCpre = ALL_stage5Table.(metricName)(idxCpre);
                    xCcmp = ALL_stage5Table.(metricName)(idxCcmp);

                    if ~isempty(xOpre) && ~isempty(xOcmp)
                        deltaOrig(i) = xOcmp(1) - xOpre(1);
                    end

                    if ~isempty(xCpre) && ~isempty(xCcmp)
                        deltaCtrl(i) = xCcmp(1) - xCpre(1);
                    end
                end

                valid = ~isnan(deltaOrig) & ~isnan(deltaCtrl);

                diffDelta = deltaOrig(valid) - deltaCtrl(valid);

                if sum(valid) >= 2 && exist('signrank', 'file') == 2
                    p_delta = signrank(deltaOrig(valid), deltaCtrl(valid));
                else
                    p_delta = nan;
                end

                if isempty(diffDelta)
                    mean_diff_delta = nan;
                    median_diff_delta = nan;
                else
                    mean_diff_delta = mean(diffDelta);
                    median_diff_delta = median(diffDelta);
                end

                tmp = table( ...
                    {metricName}, ...
                    {phaseOrder{p}}, ...
                    {'Original'}, ...
                    {controlName}, ...
                    sum(valid), ...
                    p_delta, ...
                    mean(deltaOrig(valid), 'omitnan'), ...
                    mean(deltaCtrl(valid), 'omitnan'), ...
                    mean_diff_delta, ...
                    median_diff_delta, ...
                    'VariableNames', { ...
                    'metric', ...
                    'phase5', ...
                    'model_A', ...
                    'model_B', ...
                    'n_pair', ...
                    'p_signrank_delta_A_vs_B', ...
                    'mean_delta_original', ...
                    'mean_delta_control', ...
                    'mean_delta_difference_original_minus_control', ...
                    'median_delta_difference_original_minus_control'});

                groupDeltaStats = [groupDeltaStats; tmp];
            end
        end
    end

    writetable(groupDeltaStats, fullfile(outDir, 'GROUP_eta_control_delta_original_vs_controls.csv'));
end

%% 7) 绘图：eta 在四个模型下的阶段变化
if ~isempty(ALL_stage5Table)

    phaseColors = [
        76,120,168;
        89,161,79;
        242,142,43;
        225,87,89;
        128,115,172
        ] / 255;

    modelColors = [
        0.05 0.05 0.05;
        0.20 0.45 0.75;
        0.85 0.35 0.20;
        0.35 0.60 0.30
        ];

    %% 7.1 eta 原始阶段均值：四模型一张图
    fig = figure('Color','w', 'Position',[100 100 760 470]);
    hold on;

    for mm = 1:numel(modelOrder)

        modelName = modelOrder{mm};

        y = nan(1,numel(phaseOrder));
        e = nan(1,numel(phaseOrder));

        for p = 1:numel(phaseOrder)
            idx = strcmp(ALL_stage5Table.model, modelName) & ...
                  strcmp(ALL_stage5Table.phase5, phaseOrder{p});

            x = ALL_stage5Table.eta_mean(idx);
            x = x(~isnan(x) & isfinite(x));

            y(p) = mean(x);
            e(p) = std(x) / sqrt(numel(x));
        end

        errorbar( ...
            1:numel(phaseOrder), y, e, ...
            '-o', ...
            'Color', modelColors(mm,:), ...
            'MarkerFaceColor', modelColors(mm,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 2.0, ...
            'CapSize', 8, ...
            'MarkerSize', 6);
    end

    xlim([0.75 numel(phaseOrder)+0.25]);
    xticks(1:numel(phaseOrder));
    xticklabels(phaseShort);

    ylabel('\eta', 'FontWeight','bold');
    title('\eta under density and randomization controls', 'FontWeight','bold');

    legend(modelOrder, 'Location','best', 'Box','off');

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',11, ...
        'LineWidth',1.2, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top');

    ax = gca;
    ax.YAxis.Exponent = 0;

    exportgraphics(fig, fullfile(outDir, 'FIG_eta_control_by_model.png'), 'Resolution', 600);
    savefig(fig, fullfile(outDir, 'FIG_eta_control_by_model.fig'));

    %% 7.2 delta eta = phase - pre：四模型一张图
    fig2 = figure('Color','w', 'Position',[100 100 760 470]);
    hold on;

    caseIDs = unique(ALL_stage5Table.case_id, 'stable');

    for mm = 1:numel(modelOrder)

        modelName = modelOrder{mm};

        deltaMat = nan(numel(caseIDs), numel(phaseOrder)-1);

        for i = 1:numel(caseIDs)

            idxPre = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                     strcmp(ALL_stage5Table.model, modelName) & ...
                     strcmp(ALL_stage5Table.phase5, 'pre-ictal');

            xPre = ALL_stage5Table.eta_mean(idxPre);

            if isempty(xPre)
                continue;
            end

            for p = 2:numel(phaseOrder)

                idxCmp = strcmp(ALL_stage5Table.case_id, caseIDs{i}) & ...
                         strcmp(ALL_stage5Table.model, modelName) & ...
                         strcmp(ALL_stage5Table.phase5, phaseOrder{p});

                xCmp = ALL_stage5Table.eta_mean(idxCmp);

                if ~isempty(xCmp)
                    deltaMat(i,p-1) = xCmp(1) - xPre(1);
                end
            end
        end

        y = nanmean(deltaMat, 1);
        e = nanstd(deltaMat, 0, 1) ./ sqrt(sum(~isnan(deltaMat), 1));

        errorbar( ...
            1:(numel(phaseOrder)-1), y, e, ...
            '-o', ...
            'Color', modelColors(mm,:), ...
            'MarkerFaceColor', modelColors(mm,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 2.0, ...
            'CapSize', 8, ...
            'MarkerSize', 6);
    end

    yline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.1);

    xlim([0.75 numel(phaseOrder)-0.75]);
    xticks(1:(numel(phaseOrder)-1));
    xticklabels(phaseShort(2:end));

    ylabel('\Delta\eta relative to pre-ictal', 'FontWeight','bold');
    title('Phase-related \Delta\eta under control models', 'FontWeight','bold');

    legend(modelOrder, 'Location','best', 'Box','off');

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',11, ...
        'LineWidth',1.2, ...
        'Box','off', ...
        'TickDir','out', ...
        'Layer','top');

    ax = gca;
    ax.YAxis.Exponent = 0;

    exportgraphics(fig2, fullfile(outDir, 'FIG_delta_eta_control_by_model.png'), 'Resolution', 600);
    savefig(fig2, fullfile(outDir, 'FIG_delta_eta_control_by_model.fig'));
end

%% 8) 保存 MAT
save(fullfile(outDir, 'ALL_CASES_eta_density_randomization_control_results.mat'), ...
    'ALL_windowTable', ...
    'ALL_stage5Table', ...
    'CASE_density_summary', ...
    'groupDescTable', ...
    'groupPhaseStats', ...
    'groupDeltaStats', ...
    'nRand', ...
    'densityTargetMode', ...
    'densityQuantile', ...
    'doFullSpectrum');

fprintf('\n==================================================\n');
fprintf('实验 2 完成。\n');
fprintf('结果保存到:\n%s\n', outDir);
fprintf('==================================================\n');

%% =========================================================
% Local functions
%% =========================================================

function W = sanitize_W(W)
    if isempty(W)
        W = [];
        return;
    end

    W = double(W);

    if size(W,1) ~= size(W,2)
        error('W 不是方阵。');
    end

    n = size(W,1);

    W(~isfinite(W)) = 0;
    W(W < 0) = 0;

    if n > 0
        W(1:n+1:end) = 0;
    end
end

function M = compute_eta_metrics(W, doFullSpectrum)

    W = sanitize_W(W);

    M = struct();

    if isempty(W)
        M.num_hyperedges = 0;
        M.sum_W = 0;
        M.resource_R = 0;
        M.density_W = 0;
        M.binary_density = 0;
        M.lambda1 = 0;
        M.eta = 0;
        M.lambda2 = nan;
        M.lambda2_ratio = nan;
        M.N_eff_lambda = nan;
        return;
    end

    n = size(W,1);
    M.num_hyperedges = n;

    if n < 2
        M.sum_W = 0;
        M.resource_R = 0;
        M.density_W = 0;
        M.binary_density = 0;
        M.lambda1 = 0;
        M.eta = 0;
        M.lambda2 = nan;
        M.lambda2_ratio = nan;
        M.N_eff_lambda = nan;
        return;
    end

    sum_W = sum(W(:));
    M.sum_W = sum_W;
    M.resource_R = sum_W / n;
    M.density_W = sum_W / (n * (n - 1));
    M.binary_density = nnz(W > 0) / (n * (n - 1));

    if sum_W == 0
        M.lambda1 = 0;
        M.eta = 0;
        M.lambda2 = nan;
        M.lambda2_ratio = nan;
        M.N_eff_lambda = nan;
        return;
    end

    if doFullSpectrum
        eigVals = eig(W);
        eigAbs = sort(abs(eigVals), 'descend');

        lambda1 = eigAbs(1);

        if numel(eigAbs) >= 2
            lambda2 = eigAbs(2);
        else
            lambda2 = 0;
        end

        M.lambda1 = lambda1;
        M.lambda2 = lambda2;

        if lambda1 > 0
            M.eta = M.resource_R / lambda1;
            M.lambda2_ratio = lambda2 / lambda1;
        else
            M.eta = 0;
            M.lambda2_ratio = nan;
        end

        eigAbs_nonzero = eigAbs(eigAbs > 0);

        if isempty(eigAbs_nonzero)
            M.N_eff_lambda = 0;
        else
            p = eigAbs_nonzero ./ sum(eigAbs_nonzero);
            M.N_eff_lambda = 1 / sum(p.^2);
        end

    else
        lambda1 = spectral_radius_fast(W);

        M.lambda1 = lambda1;

        if lambda1 > 0
            M.eta = M.resource_R / lambda1;
        else
            M.eta = 0;
        end

        M.lambda2 = nan;
        M.lambda2_ratio = nan;
        M.N_eff_lambda = nan;
    end
end

function lambda1 = spectral_radius_fast(W)

    n = size(W,1);

    if n == 0 || sum(W(:)) == 0
        lambda1 = 0;
        return;
    end

    try
        if n >= 20
            opts.disp = 0;
            val = eigs(sparse(W), 1, 'largestabs', opts);
            lambda1 = abs(val);
        else
            lambda1 = max(abs(eig(W)));
        end
    catch
        lambda1 = max(abs(eig(W)));
    end

    if ~isfinite(lambda1)
        lambda1 = 0;
    end
end

function Wd = density_match_topK(W, targetDensity)

    W = sanitize_W(W);
    n = size(W,1);

    Wd = zeros(size(W));

    if n < 2 || targetDensity <= 0 || sum(W(:)) == 0
        return;
    end

    mask = W > 0;
    idx = find(mask);
    weights = W(idx);

    K_orig = numel(idx);
    K_target = floor(targetDensity * n * (n - 1));

    if K_target < 1 && K_orig > 0
        K_target = 1;
    end

    K_target = min(K_target, K_orig);

    if K_target <= 0
        return;
    end

    [~, order] = sort(weights, 'descend');
    keepIdx = idx(order(1:K_target));

    Wd(keepIdx) = W(keepIdx);
    Wd(1:n+1:end) = 0;
end

function Wr = weight_shuffle_keep_topology(W)

    W = sanitize_W(W);
    n = size(W,1);

    Wr = zeros(size(W));

    if n < 2 || sum(W(:)) == 0
        return;
    end

    idx = find(W > 0);
    weights = W(idx);

    if isempty(weights)
        return;
    end

    weights_shuffled = weights(randperm(numel(weights)));

    Wr(idx) = weights_shuffled;
    Wr(1:n+1:end) = 0;
end

function Wr = topology_shuffle_keep_weights(W)

    W = sanitize_W(W);
    n = size(W,1);

    Wr = zeros(size(W));

    if n < 2 || sum(W(:)) == 0
        return;
    end

    idx = find(W > 0);
    weights = W(idx);

    K = numel(weights);

    if K == 0
        return;
    end

    offdiagMask = true(n,n);
    offdiagMask(1:n+1:end) = false;
    possibleIdx = find(offdiagMask);

    if K > numel(possibleIdx)
        K = numel(possibleIdx);
        weights = weights(1:K);
    end

    newIdx = possibleIdx(randperm(numel(possibleIdx), K));
    weights_shuffled = weights(randperm(numel(weights)));

    Wr(newIdx) = weights_shuffled;
    Wr(1:n+1:end) = 0;
end

function [Mmean, Msd] = average_metric_structs(metricCells, doFullSpectrum)

    baseMetrics = { ...
        'eta', ...
        'resource_R', ...
        'lambda1', ...
        'sum_W', ...
        'density_W', ...
        'binary_density', ...
        'num_hyperedges'};

    if doFullSpectrum
        baseMetrics = [baseMetrics, {'lambda2','lambda2_ratio','N_eff_lambda'}];
    else
        baseMetrics = [baseMetrics, {'lambda2','lambda2_ratio','N_eff_lambda'}];
    end

    Mmean = struct();
    Msd = struct();

    for k = 1:numel(baseMetrics)

        metric = baseMetrics{k};
        x = nan(numel(metricCells),1);

        for i = 1:numel(metricCells)
            if isfield(metricCells{i}, metric)
                x(i) = metricCells{i}.(metric);
            end
        end

        Mmean.(metric) = mean(x, 'omitnan');
        Msd.(metric) = std(x, 'omitnan');
    end
end

function row = make_window_row(case_id, window_idx, phase5, model, targetDensity, nRand, isRandomModel, M, Msd)

    if isempty(Msd)
        Msd = struct();
        Msd.eta = nan;
        Msd.resource_R = nan;
        Msd.lambda1 = nan;
        Msd.sum_W = nan;
        Msd.density_W = nan;
        Msd.binary_density = nan;
        Msd.num_hyperedges = nan;
        Msd.lambda2 = nan;
        Msd.lambda2_ratio = nan;
        Msd.N_eff_lambda = nan;
    end

    row = table( ...
        {case_id}, ...
        window_idx, ...
        {phase5}, ...
        {model}, ...
        targetDensity, ...
        nRand, ...
        isRandomModel, ...
        M.eta, ...
        M.resource_R, ...
        M.lambda1, ...
        M.sum_W, ...
        M.density_W, ...
        M.binary_density, ...
        M.num_hyperedges, ...
        M.lambda2, ...
        M.lambda2_ratio, ...
        M.N_eff_lambda, ...
        Msd.eta, ...
        Msd.resource_R, ...
        Msd.lambda1, ...
        Msd.sum_W, ...
        Msd.density_W, ...
        Msd.binary_density, ...
        Msd.num_hyperedges, ...
        Msd.lambda2, ...
        Msd.lambda2_ratio, ...
        Msd.N_eff_lambda, ...
        'VariableNames', { ...
        'case_id', ...
        'window_idx', ...
        'phase5', ...
        'model', ...
        'target_binary_density', ...
        'nRand', ...
        'isRandomModel', ...
        'eta', ...
        'resource_R', ...
        'lambda1', ...
        'sum_W', ...
        'density_W', ...
        'binary_density', ...
        'num_hyperedges', ...
        'lambda2', ...
        'lambda2_ratio', ...
        'N_eff_lambda', ...
        'eta_rand_sd', ...
        'resource_R_rand_sd', ...
        'lambda1_rand_sd', ...
        'sum_W_rand_sd', ...
        'density_W_rand_sd', ...
        'binary_density_rand_sd', ...
        'num_hyperedges_rand_sd', ...
        'lambda2_rand_sd', ...
        'lambda2_ratio_rand_sd', ...
        'N_eff_lambda_rand_sd'});
end