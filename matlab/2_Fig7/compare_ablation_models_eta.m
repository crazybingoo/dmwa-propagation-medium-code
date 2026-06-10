%%
clc; clear; close all;

%% =========================================================
% 批量运行 24 次发作：计算 eta (内部机制消融实验 Ablation Study)
%
% 对比三种 eta 值：
%   1. eta_orig: 原版机制 (Coverage + DegreeBias)
%   2. eta_cov: 仅保留覆盖率 Coverage，去除 DegreeBias
%   3. eta_deg: 仅保留度偏置 DegreeBias，无视超边阶数差异
%
% 阶段：pre-ictal, early, mid, late, post-ictal
%
% 输出内容：
%   包含 orig, cov, deg 三种指标的 CSV 和 3 张箱线图
%
% eta 定义：
%   eta = ((1/N) * sum_{i,j} W_{ij}) / lambda_1
%   其中 N 为网络节点数（超边数），lambda_1 为 W 的主特征值（谱半径）
% =========================================================

%% 0) 输出目录
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig7';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% 1) 基本参数
fs = 1024;
winLenSec = 3;

%% 2) 24 个病例信息
caseInfo = struct([]);

caseInfo(1).case_id   = 'lhs_cut07';
caseInfo(1).file_path = 'E:\wcldematlab\keep\new_idea\lhs_cut07\lihongsen_cut07_Gamma.mat';
caseInfo(1).pre_idx   = 1:200; caseInfo(1).ictal_idx = 201:400; caseInfo(1).post_idx  = 401:598;

caseInfo(2).case_id   = 'lhs_cut06';
caseInfo(2).file_path = 'E:\wcldematlab\keep\new_idea\lhs_cut06\lihongsen_cut06_Gamma.mat';
caseInfo(2).pre_idx   = 1:122; caseInfo(2).ictal_idx = 123:233; caseInfo(2).post_idx  = 234:338;

caseInfo(3).case_id   = 'gzs_cut06';
caseInfo(3).file_path = 'E:\wcldematlab\keep\new_idea\gzs\gongzishu_cut06_Gamma.mat';
caseInfo(3).pre_idx   = 1:195; caseInfo(3).ictal_idx = 196:330; caseInfo(3).post_idx  = 331:477;

caseInfo(4).case_id   = 'gzs_cut07';
caseInfo(4).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut07\gongzishu_cut07_Gamma.mat';
caseInfo(4).pre_idx   = 1:137; caseInfo(4).ictal_idx = 138:297; caseInfo(4).post_idx  = 298:398;

caseInfo(5).case_id   = 'gzs_cut08';
caseInfo(5).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut08\gongzishu_cut08_Gamma.mat';
caseInfo(5).pre_idx   = 1:252; caseInfo(5).ictal_idx = 253:413; caseInfo(5).post_idx  = 414:498;

caseInfo(6).case_id   = 'wc_cut03';
caseInfo(6).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut03\wangchun_cut03_Gamma.mat';
caseInfo(6).pre_idx   = 1:199; caseInfo(6).ictal_idx = 200:285; caseInfo(6).post_idx  = 286:398;

caseInfo(7).case_id   = 'wc_cut06';
caseInfo(7).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut06\wangchun_cut06_Gamma.mat';
caseInfo(7).pre_idx   = 1:239; caseInfo(7).ictal_idx = 240:331; caseInfo(7).post_idx  = 332:498;

caseInfo(8).case_id   = 'wc_cut08';
caseInfo(8).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut08\wangchun_cut08_Gamma.mat';
caseInfo(8).pre_idx   = 1:236; caseInfo(8).ictal_idx = 237:313; caseInfo(8).post_idx  = 314:498;

caseInfo(9).case_id   = 'cxm_cut104';
caseInfo(9).file_path = 'E:\wcldematlab\keep\new_idea\cxm_cut104\cxm_cut104_Gamma.mat';
caseInfo(9).pre_idx   = 1:59; caseInfo(9).ictal_idx = 60:112; caseInfo(9).post_idx  = 113:385;

caseInfo(10).case_id   = 'bsp_cut84';
caseInfo(10).file_path = 'E:\wcldematlab\keep\new_idea\bsp_cut84\bsp_cut84_Gamma.mat';
caseInfo(10).pre_idx   = 1:244; caseInfo(10).ictal_idx = 245:470; caseInfo(10).post_idx  = 471:598;

caseInfo(11).case_id   = 'zzy_cut151';
caseInfo(11).file_path = 'E:\wcldematlab\keep\new_idea\zzy\zzy_cut151_Gamma.mat';
caseInfo(11).pre_idx   = 1:150; caseInfo(11).ictal_idx = 151:271; caseInfo(11).post_idx  = 272:404;

caseInfo(12).case_id   = 'zzy_cut147';
caseInfo(12).file_path = 'E:\wcldematlab\keep\new_idea\zzy_cut147\zzy_cut147_Gamma.mat';
caseInfo(12).pre_idx   = 1:180; caseInfo(12).ictal_idx = 181:261; caseInfo(12).post_idx  = 262:498;

caseInfo(13).case_id   = 'lzq_cut73';
caseInfo(13).file_path = 'E:\wcldematlab\keep\new_idea\lzq\lzq_cut73_Gamma.mat';
caseInfo(13).pre_idx   = 1:69; caseInfo(13).ictal_idx = 70:130; caseInfo(13).post_idx  = 131:198;

caseInfo(14).case_id   = 'lzq_cut77';
caseInfo(14).file_path = 'E:\wcldematlab\keep\new_idea\lzq_cut77\lzq_cut77_Gamma.mat';
caseInfo(14).pre_idx   = 1:252; caseInfo(14).ictal_idx = 253:346; caseInfo(14).post_idx  = 347:448;

caseInfo(15).case_id   = 'll_process';
caseInfo(15).file_path = 'E:\wcldematlab\keep\new_idea\ll\lilei_process_Gamma.mat';
caseInfo(15).pre_idx   = 1:154; caseInfo(15).ictal_idx = 155:230; caseInfo(15).post_idx  = 231:324;

caseInfo(16).case_id   = 'gzw_cut29';
caseInfo(16).file_path = 'E:\wcldematlab\keep\new_idea\gzw\gzw_cut29_Gamma.mat';
caseInfo(16).pre_idx   = 1:253; caseInfo(16).ictal_idx = 254:331; caseInfo(16).post_idx  = 332:489;

caseInfo(17).case_id   = 'hds_process';
caseInfo(17).file_path = 'E:\wcldematlab\keep\new_idea\hds\haidongsheng_process_Gamma.mat';
caseInfo(17).pre_idx   = 1:135; caseInfo(17).ictal_idx = 136:270; caseInfo(17).post_idx  = 271:309;

caseInfo(18).case_id   = 'ssh_cut109';
caseInfo(18).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut109\ssh_cut109_Gamma.mat';
caseInfo(18).pre_idx   = 1:117; caseInfo(18).ictal_idx = 118:300; caseInfo(18).post_idx  = 301:598;

caseInfo(19).case_id   = 'ssh_cut110';
caseInfo(19).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut110\ssh_cut110_Gamma.mat';
caseInfo(19).pre_idx   = 1:255; caseInfo(19).ictal_idx = 256:482; caseInfo(19).post_idx  = 483:598;

caseInfo(20).case_id   = 'ssh_cut111';
caseInfo(20).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut111\ssh_cut111_Gamma.mat';
caseInfo(20).pre_idx   = 1:311; caseInfo(20).ictal_idx = 312:540; caseInfo(20).post_idx  = 541:598;

caseInfo(21).case_id   = 'fys_cut93';
caseInfo(21).file_path = 'E:\wcldematlab\keep\new_idea\fys_cut93\fys_cut93_Gamma.mat';
caseInfo(21).pre_idx   = 1:94; caseInfo(21).ictal_idx = 95:170; caseInfo(21).post_idx  = 171:318;

caseInfo(22).case_id   = 'fys_cut99';
caseInfo(22).file_path = 'E:\wcldematlab\keep\new_idea\fys_cut99\fys_cut99_Gamma.mat';
caseInfo(22).pre_idx   = 1:95; caseInfo(22).ictal_idx = 96:195; caseInfo(22).post_idx  = 196:350;

caseInfo(23).case_id   = 'gs_cut58';
caseInfo(23).file_path = 'E:\wcldematlab\keep\new_idea\gs_cut58\gushuai_cut58_Gamma.mat';
caseInfo(23).pre_idx   = 1:162; caseInfo(23).ictal_idx = 163:228; caseInfo(23).post_idx  = 229:310;

caseInfo(24).case_id   = 'gwh_cut104';
caseInfo(24).file_path = 'E:\wcldematlab\keep\new_idea\gwh\guan_cut104_Gamma.mat';
caseInfo(24).pre_idx   = 1:29; caseInfo(24).ictal_idx = 30:95; caseInfo(24).post_idx  = 96:148;

nCase = numel(caseInfo);

%% 3) 汇总容器
allWindowTables = cell(nCase,1);
allStage5Tables = cell(nCase,1);

%% 主循环
for c = 1:nCase
    case_id   = caseInfo(c).case_id;
    FileName  = caseInfo(c).file_path;
    pre_idx   = caseInfo(c).pre_idx(:)';
    ictal_idx = caseInfo(c).ictal_idx(:)';
    post_idx  = caseInfo(c).post_idx(:)';

    fprintf('\n==================================================\n');
    fprintf('(%d/%d) 开始处理: %s\n', c, nCase, case_id);
    fprintf('==================================================\n');

    if ~exist(FileName, 'file')
        warning('文件不存在，跳过: %s', FileName);
        continue;
    end

    S0 = load(FileName);
    if ~isfield(S0, 'X1')
        warning('%s 中没有变量 X1，跳过。', case_id);
        continue;
    end
    X1 = S0.X1;

    totalIterations = floor(size(X1, 2) / fs) - winLenSec;
    if totalIterations < 0
        continue;
    end
    T = totalIterations + 1;
    time_axis = (1:T)';

    pre_idx   = unique(pre_idx(pre_idx >= 1 & pre_idx <= T), 'stable');
    ictal_idx = unique(ictal_idx(ictal_idx >= 1 & ictal_idx <= T), 'stable');
    post_idx  = unique(post_idx(post_idx >= 1 & post_idx <= T), 'stable');

    ictal_idx = setdiff(ictal_idx, pre_idx, 'stable');
    post_idx  = setdiff(post_idx, union(pre_idx, ictal_idx), 'stable');

    if isempty(ictal_idx)
        continue;
    end

    nIctal = numel(ictal_idx);
    splitPts = round(linspace(0, nIctal, 4));
    early_idx = ictal_idx(1:splitPts(2));
    mid_idx   = ictal_idx(splitPts(2)+1:splitPts(3));
    late_idx  = ictal_idx(splitPts(3)+1:end);

    phase5 = repmat({''}, T, 1);
    phase5(pre_idx)   = {'pre-ictal'};
    phase5(early_idx) = {'early'};
    phase5(mid_idx)   = {'mid'};
    phase5(late_idx)  = {'late'};
    phase5(post_idx)  = {'post-ictal'};
    phase5(cellfun(@isempty, phase5)) = {'unused'};

    all_hyperEdges_all = cell(T, 1);
    for time = 0:totalIterations
        start_point = fs * time + 1;
        end_point   = fs * (time + winLenSec);
        datanew = X1(:, start_point:end_point);
        all_hyperEdges_all{time + 1} = gain_hyperEdges_23(datanew);
    end

    eta_orig_all = nan(T,1);
    eta_cov_all  = nan(T,1);
    eta_deg_all  = nan(T,1);
    num_qi = nan(T,1);

    for t = 1:T
        HE = all_hyperEdges_all{t};
        num_HE = length(HE);
        num_qi(t) = num_HE;

        if num_HE < 2
            eta_orig_all(t) = 0;
            eta_cov_all(t)  = 0;
            eta_deg_all(t)  = 0;
            continue;
        end

        sizes = cellfun(@length, HE);
        sizes = sizes(:);

        O = zeros(num_HE, num_HE);
        for i = 1:num_HE
            for j = i+1:num_HE
                overlap_len = length(intersect(HE{i}, HE{j}));
                if overlap_len > 0
                    O(i, j) = overlap_len;
                    O(j, i) = overlap_len;
                end
            end
        end

        D = sum(O > 0, 2);
        Sizes_i = repmat(sizes, 1, num_HE);
        Sizes_j = repmat(sizes', num_HE, 1);
        D_i     = repmat(D, 1, num_HE);
        D_j     = repmat(D', num_HE, 1);

        Coverage = O ./ Sizes_i;
        DegreeSum  = D_i + D_j;
        DegreeBias = D_j ./ DegreeSum;
        DegreeBias(DegreeSum == 0) = 0.5;

        W_orig = Coverage;
        EqualSizeMask = (Sizes_i == Sizes_j) & (O > 0);
        W_orig(EqualSizeMask) = W_orig(EqualSizeMask) .* DegreeBias(EqualSizeMask);
        W_orig(1:num_HE+1:end) = 0;

        W_cov = Coverage;
        W_cov(1:num_HE+1:end) = 0;

        W_deg = O .* DegreeBias;
        W_deg(1:num_HE+1:end) = 0;

        eta_orig_all(t) = calculate_eta_from_W(W_orig);
        eta_cov_all(t)  = calculate_eta_from_W(W_cov);
        eta_deg_all(t)  = calculate_eta_from_W(W_deg);
    end

    windowTable = table( ...
        repmat({case_id}, T, 1), time_axis, phase5, ...
        eta_orig_all, eta_cov_all, eta_deg_all, num_qi, ...
        'VariableNames', {'case_id', 'window_idx', 'phase5', 'eta_orig', 'eta_cov', 'eta_deg', 'num_hyperedges'});
    allWindowTables{c} = windowTable;

    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    nPhase = numel(phaseOrder);

    stage5Table = table(repmat({case_id}, nPhase, 1), phaseOrder', ...
        'VariableNames', {'case_id', 'phase5'});

    for p = 1:nPhase
        idx = strcmp(windowTable.phase5, phaseOrder{p});
        stage5Table.n_windows(p) = sum(idx);
        stage5Table.eta_orig_mean(p) = nanmean(windowTable.eta_orig(idx));
        stage5Table.eta_cov_mean(p)  = nanmean(windowTable.eta_cov(idx));
        stage5Table.eta_deg_mean(p)  = nanmean(windowTable.eta_deg(idx));
    end
    allStage5Tables{c} = stage5Table;

    writetable(windowTable, fullfile(outDir, [case_id '_window_eta_ablation.csv']));
    writetable(stage5Table, fullfile(outDir, [case_id '_stage5_eta_ablation.csv']));

    save(fullfile(outDir, [case_id '_eta_ablation_results.mat']), ...
        'case_id', 'FileName', 'fs', 'winLenSec', ...
        'time_axis', 'phase5', 'num_qi', ...
        'eta_orig_all', 'eta_cov_all', 'eta_deg_all', ...
        'windowTable', 'stage5Table');

    fprintf('%s 处理完毕。\n', case_id);
end

validWindow = ~cellfun(@isempty, allWindowTables);
validStage5 = ~cellfun(@isempty, allStage5Tables);

ALL_windowTable = vertcat(allWindowTables{validWindow});
ALL_stage5Table = vertcat(allStage5Tables{validStage5});

writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_eta_ablation.csv'));
writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_eta_ablation.csv'));

if ~isempty(ALL_stage5Table)
    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    metricNames = {'eta_orig_mean', 'eta_cov_mean', 'eta_deg_mean'};
    metricTitles = {'Original (Coverage + DegreeBias)', 'Coverage Only', 'DegreeBias Only'};
    fileSuffix = {'orig', 'cov', 'deg'};

    for m = 1:3
        plot_ablation_boxplot(ALL_stage5Table, phaseOrder, metricNames{m}, metricTitles{m}, ...
            fullfile(outDir, ['GROUP_stage5_eta_', fileSuffix{m}, '_boxplot.png']));
    end
end

fprintf('\n==================================================\n');
fprintf('eta 消融实验处理完成。结果已保存到:\n%s\n', outDir);
fprintf('==================================================\n');

function eta_val = calculate_eta_from_W(W)
    if isempty(W)
        eta_val = 0;
        return;
    end
    W = double(W);
    n = size(W, 1);
    if size(W, 2) ~= n
        error('Input W must be a square matrix.');
    end
    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    W(1:n+1:end) = 0;
    if n < 2
        eta_val = 0;
        return;
    end
    sum_W = sum(W(:));
    if sum_W <= 0
        eta_val = 0;
        return;
    end
    lambda_1 = max(abs(eig(W)));
    if ~isfinite(lambda_1) || lambda_1 <= eps
        eta_val = 0;
    else
        eta_val = (sum_W / n) / lambda_1;
    end
end

function plot_ablation_boxplot(T, phaseOrder, metricName, metricTitle, outFile)
    phaseColors = [76,120,168;
                   89,161,79;
                   242,142,43;
                   225,87,89;
                   128,115,172] / 255;

    fig = figure('Color','w', 'Position',[100 100 520 420]);
    hold on;

    Y = cell(numel(phaseOrder), 1);
    maxLen = 0;
    for p = 1:numel(phaseOrder)
        idx = strcmp(T.phase5, phaseOrder{p});
        y = T.(metricName)(idx);
        y = y(~isnan(y));
        Y{p} = y;
        maxLen = max(maxLen, numel(y));
    end

    Mplot = nan(maxLen, numel(phaseOrder));
    for p = 1:numel(phaseOrder)
        if ~isempty(Y{p})
            Mplot(1:numel(Y{p}), p) = Y{p};
        end
    end

    boxplot(Mplot, ...
        'Colors', [0.15 0.15 0.15], ...
        'Symbol', '', ...
        'Widths', 0.55, ...
        'MedianStyle', 'line');

    hBox = findobj(gca, 'Tag', 'Box');
    set(hBox, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);

    hMedian = findobj(gca, 'Tag', 'Median');
    set(hMedian, 'Color', [0 0 0], 'LineWidth', 1.4);

    hWhisker = findobj(gca, 'Tag', 'Whisker');
    set(hWhisker, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

    hAdj = findobj(gca, 'Tag', 'Adjacent Value');
    set(hAdj, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

    rng(2025);
    for p = 1:numel(phaseOrder)
        y = Y{p};
        if isempty(y)
            continue;
        end
        jitter = (rand(numel(y),1) - 0.5) * 0.18;
        scatter( ...
            p + jitter, y, 26, ...
            'MarkerFaceColor', phaseColors(p,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.45, ...
            'MarkerFaceAlpha', 0.9, ...
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
    xticklabels({'pre-ictal','early','mid','late','post-ictal'});
    xtickangle(25);
    ylabel('\eta', 'FontWeight', 'bold');
    title(metricTitle, 'FontWeight', 'normal');

    Xall = [];
    for p = 1:numel(phaseOrder)
        Xall = [Xall; Y{p}(:)]; %#ok<AGROW>
    end
    Xall = Xall(~isnan(Xall));

    if ~isempty(Xall)
        y_min = min(Xall);
        y_max = max(Xall);
        if y_min == y_max
            pad = max(abs(y_min) * 0.05, 0.01);
        else
            pad = 0.08 * (y_max - y_min);
        end
        ylim([y_min - pad, y_max + pad]);
        ax = gca;
        ax.YAxis.Exponent = 0;
    end

    exportgraphics(fig, outFile, 'Resolution', 600);
    [outDir0, outName0] = fileparts(outFile);
    savefig(fig, fullfile(outDir0, [outName0 '.fig']));
end
