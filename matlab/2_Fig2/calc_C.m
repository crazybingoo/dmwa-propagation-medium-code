%%
clc; clear; close all;

%% =========================================================
% 批量运行 24 次发作：计算高阶网络循环度 C
%
% C 定义：
%   C = sum(|Im(lambda_k)|) / sum(|lambda_k|)
%   其中 lambda_k 是矩阵 W 的特征值
%
% 阶段：
%   pre-ictal, early, mid, late, post-ictal
%   其中 early / mid / late 为 ictal_idx 三等分
%
% 输出内容：
%   1) 每个病人的窗级结果 CSV
%   2) 每个病人的 5 阶段结果 CSV
%   3) 所有病人的窗级总表 / 阶段级总表 CSV
%   4) 病人级组统计 CSV（5阶段）
%   5) C 的 Nature 风格箱线图（PNG + FIG）
%
% 依赖：
%   gain_hyperEdges_23.m 需在 MATLAB 路径中
% =========================================================

%% 0) 输出目录
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig2';
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
caseInfo(1).pre_idx   = 1:200;
caseInfo(1).ictal_idx = 201:400;
caseInfo(1).post_idx  = 401:598;

caseInfo(2).case_id   = 'lhs_cut06';
caseInfo(2).file_path = 'E:\wcldematlab\keep\new_idea\lhs_cut06\lihongsen_cut06_Gamma.mat';
caseInfo(2).pre_idx   = 1:122;
caseInfo(2).ictal_idx = 123:233;
caseInfo(2).post_idx  = 234:338;

caseInfo(3).case_id   = 'gzs_cut06';
caseInfo(3).file_path = 'E:\wcldematlab\keep\new_idea\gzs\gongzishu_cut06_Gamma.mat';
caseInfo(3).pre_idx   = 1:195;
caseInfo(3).ictal_idx = 196:330;
caseInfo(3).post_idx  = 331:477;

caseInfo(4).case_id   = 'gzs_cut07';
caseInfo(4).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut07\gongzishu_cut07_Gamma.mat';
caseInfo(4).pre_idx   = 1:137;
caseInfo(4).ictal_idx = 138:297;
caseInfo(4).post_idx  = 298:398;

caseInfo(5).case_id   = 'gzs_cut08';
caseInfo(5).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut08\gongzishu_cut08_Gamma.mat';
caseInfo(5).pre_idx   = 1:252;
caseInfo(5).ictal_idx = 253:413;
caseInfo(5).post_idx  = 414:498;

caseInfo(6).case_id   = 'wc_cut03';
caseInfo(6).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut03\wangchun_cut03_Gamma.mat';
caseInfo(6).pre_idx   = 1:199;
caseInfo(6).ictal_idx = 200:285;
caseInfo(6).post_idx  = 286:398;

caseInfo(7).case_id   = 'wc_cut06';
caseInfo(7).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut06\wangchun_cut06_Gamma.mat';
caseInfo(7).pre_idx   = 1:239;
caseInfo(7).ictal_idx = 240:331;
caseInfo(7).post_idx  = 332:498;

caseInfo(8).case_id   = 'wc_cut08';
caseInfo(8).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut08\wangchun_cut08_Gamma.mat';
caseInfo(8).pre_idx   = 1:236;
caseInfo(8).ictal_idx = 237:313;
caseInfo(8).post_idx  = 314:498;

caseInfo(9).case_id   = 'cxm_cut104';
caseInfo(9).file_path = 'E:\wcldematlab\keep\new_idea\cxm_cut104\cxm_cut104_Gamma.mat';
caseInfo(9).pre_idx   = 1:59;
caseInfo(9).ictal_idx = 60:112;
caseInfo(9).post_idx  = 113:385;

caseInfo(10).case_id   = 'bsp_cut84';
caseInfo(10).file_path = 'E:\wcldematlab\keep\new_idea\bsp_cut84\bsp_cut84_Gamma.mat';
caseInfo(10).pre_idx   = 1:244;
caseInfo(10).ictal_idx = 245:470;
caseInfo(10).post_idx  = 471:598;

caseInfo(11).case_id   = 'zzy_cut151';
caseInfo(11).file_path = 'E:\wcldematlab\keep\new_idea\zzy\zzy_cut151_Gamma.mat';
caseInfo(11).pre_idx   = 1:150;
caseInfo(11).ictal_idx = 151:271;
caseInfo(11).post_idx  = 272:404;

caseInfo(12).case_id   = 'zzy_cut147';
caseInfo(12).file_path = 'E:\wcldematlab\keep\new_idea\zzy_cut147\zzy_cut147_Gamma.mat';
caseInfo(12).pre_idx   = 1:180;
caseInfo(12).ictal_idx = 181:261;
caseInfo(12).post_idx  = 262:498;

caseInfo(13).case_id   = 'lzq_cut73';
caseInfo(13).file_path = 'E:\wcldematlab\keep\new_idea\lzq\lzq_cut73_Gamma.mat';
caseInfo(13).pre_idx   = 1:69;
caseInfo(13).ictal_idx = 70:130;
caseInfo(13).post_idx  = 131:198;

caseInfo(14).case_id   = 'lzq_cut77';
caseInfo(14).file_path = 'E:\wcldematlab\keep\new_idea\lzq_cut77\lzq_cut77_Gamma.mat';
caseInfo(14).pre_idx   = 1:252;
caseInfo(14).ictal_idx = 253:346;
caseInfo(14).post_idx  = 347:448;

caseInfo(15).case_id   = 'll_process';
caseInfo(15).file_path = 'E:\wcldematlab\keep\new_idea\ll\lilei_process_Gamma.mat';
caseInfo(15).pre_idx   = 1:154;
caseInfo(15).ictal_idx = 155:230;
caseInfo(15).post_idx  = 231:324;

caseInfo(16).case_id   = 'gzw_cut29';
caseInfo(16).file_path = 'E:\wcldematlab\keep\new_idea\gzw\gzw_cut29_Gamma.mat';
caseInfo(16).pre_idx   = 1:253;
caseInfo(16).ictal_idx = 254:331;
caseInfo(16).post_idx  = 332:489;

caseInfo(17).case_id   = 'hds_process';
caseInfo(17).file_path = 'E:\wcldematlab\keep\new_idea\hds\haidongsheng_process_Gamma.mat';
caseInfo(17).pre_idx   = 1:135;
caseInfo(17).ictal_idx = 136:270;
caseInfo(17).post_idx  = 271:309;

caseInfo(18).case_id   = 'ssh_cut109';
caseInfo(18).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut109\ssh_cut109_Gamma.mat';
caseInfo(18).pre_idx   = 1:117;
caseInfo(18).ictal_idx = 118:300;
caseInfo(18).post_idx  = 301:598;

caseInfo(19).case_id   = 'ssh_cut110';
caseInfo(19).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut110\ssh_cut110_Gamma.mat';
caseInfo(19).pre_idx   = 1:255;
caseInfo(19).ictal_idx = 256:482;
caseInfo(19).post_idx  = 483:598;

caseInfo(20).case_id   = 'ssh_cut111';
caseInfo(20).file_path = 'E:\wcldematlab\keep\new_idea\ssh_cut111\ssh_cut111_Gamma.mat';
caseInfo(20).pre_idx   = 1:311;
caseInfo(20).ictal_idx = 312:540;
caseInfo(20).post_idx  = 541:598;

caseInfo(21).case_id   = 'fys_cut93';
caseInfo(21).file_path = 'E:\wcldematlab\keep\new_idea\fys_cut93\fys_cut93_Gamma.mat';
caseInfo(21).pre_idx   = 1:94;
caseInfo(21).ictal_idx = 95:170;
caseInfo(21).post_idx  = 171:318;

caseInfo(22).case_id   = 'fys_cut99';
caseInfo(22).file_path = 'E:\wcldematlab\keep\new_idea\fys_cut99\fys_cut99_Gamma.mat';
caseInfo(22).pre_idx   = 1:95;
caseInfo(22).ictal_idx = 96:195;
caseInfo(22).post_idx  = 196:350;

caseInfo(23).case_id   = 'gs_cut58';
caseInfo(23).file_path = 'E:\wcldematlab\keep\new_idea\gs_cut58\gushuai_cut58_Gamma.mat';
caseInfo(23).pre_idx   = 1:162;
caseInfo(23).ictal_idx = 163:228;
caseInfo(23).post_idx  = 229:310;

caseInfo(24).case_id   = 'gwh_cut104';
caseInfo(24).file_path = 'E:\wcldematlab\keep\new_idea\gwh\guan_cut104_Gamma.mat';
caseInfo(24).pre_idx   = 1:29;
caseInfo(24).ictal_idx = 30:95;
caseInfo(24).post_idx  = 96:148;

nCase = numel(caseInfo);

%% 3) 汇总容器
allWindowTables = cell(nCase,1);
allStage5Tables = cell(nCase,1);

%% =========================================================
% 主循环
%% =========================================================
for c = 1:nCase

    case_id   = caseInfo(c).case_id;
    FileName  = caseInfo(c).file_path;
    pre_idx   = caseInfo(c).pre_idx(:)';
    ictal_idx = caseInfo(c).ictal_idx(:)';
    post_idx  = caseInfo(c).post_idx(:)';

    fprintf('\n==================================================\n');
    fprintf('(%d/%d) 开始处理: %s\n', c, nCase, case_id);
    fprintf('文件: %s\n', FileName);
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

    %% 3.1 总时间窗数
    totalIterations = floor(size(X1, 2) / fs) - winLenSec;
    if totalIterations < 0
        warning('%s 数据长度不足，无法形成 3 秒窗。', case_id);
        continue;
    end
    T = totalIterations + 1;
    time_axis = (1:T)';

    %% 3.2 修正阶段索引（去越界、去重、去重叠）
    pre_idx   = pre_idx(pre_idx >= 1 & pre_idx <= T);
    ictal_idx = ictal_idx(ictal_idx >= 1 & ictal_idx <= T);
    post_idx  = post_idx(post_idx >= 1 & post_idx <= T);

    pre_idx   = unique(pre_idx, 'stable');
    ictal_idx = unique(ictal_idx, 'stable');
    post_idx  = unique(post_idx, 'stable');

    overlap_pre_ictal  = intersect(pre_idx, ictal_idx);
    overlap_pre_post   = intersect(pre_idx, post_idx);
    overlap_ictal_post = intersect(ictal_idx, post_idx);

    if ~isempty(overlap_pre_ictal) || ~isempty(overlap_pre_post) || ~isempty(overlap_ictal_post)
        warning('%s 的阶段索引存在重叠。代码将按 pre > ictal > post 的优先级自动去重。', case_id);
    end

    ictal_idx = setdiff(ictal_idx, pre_idx, 'stable');
    post_idx  = setdiff(post_idx, union(pre_idx, ictal_idx), 'stable');

    if isempty(ictal_idx)
        warning('%s 的 ictal_idx 为空，跳过。', case_id);
        continue;
    end

    %% 3.3 ictal 三等分 -> early / mid / late
    nIctal = numel(ictal_idx);
    splitPts = round(linspace(0, nIctal, 4));
    i1 = splitPts(2);
    i2 = splitPts(3);

    early_idx = ictal_idx(1:i1);
    mid_idx   = ictal_idx(i1+1:i2);
    late_idx  = ictal_idx(i2+1:end);

    %% 3.4 提取超边
    all_hyperEdges_all = cell(T, 1);

    tic_extract = tic;
    for time = 0:totalIterations
        start_point = fs * time + 1;
        end_point   = fs * (time + winLenSec);
        datanew = X1(:, start_point:end_point);
        all_hyperEdges_all{time + 1} = gain_hyperEdges_23(datanew);
    end
    fprintf('%s: 超边提取完成，用时 %.3f 秒。\n', case_id, toc(tic_extract));

    %% 3.5 计算 W
    W_all = cell(T, 1);

    tic_w = tic;
    for time = 0:totalIterations
        HE = all_hyperEdges_all{time + 1};
        num_HE = length(HE);

        if num_HE < 2
            W_all{time + 1} = zeros(num_HE);
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

        W = Coverage;
        EqualSizeMask = (Sizes_i == Sizes_j) & (O > 0);
        W(EqualSizeMask) = W(EqualSizeMask) .* DegreeBias(EqualSizeMask);

        W(1:num_HE+1:end) = 0;
        W_all{time + 1} = W;
    end
    fprintf('%s: W 计算完成，用时 %.3f 秒。\n', case_id, toc(tic_w));

    %% 3.6 计算每窗 C
    C_all = nan(T,1);
    num_qi  = nan(T,1);   

    tic_C = tic;
    for t = 1:T
        W = W_all{t};

        if isempty(W)
            continue;
        end

        W = double(W);
        n = size(W, 1);

        if size(W,2) ~= n
            error('病例 %s，第 %d 个时间窗的 W 不是方阵。', case_id, t);
        end

        W(~isfinite(W)) = 0;
        W(W < 0) = 0;
        W(1:n+1:end) = 0; % 对角线清零

        num_qi(t) = n;

        if n >= 2
            % 计算特征值
            E = eig(W);
            % C = sum(|Im(lambda_k)|) / sum(|lambda_k|)
            sum_abs_imag_E = sum(abs(imag(E)));
            sum_abs_E = sum(abs(E));
            
            if sum_abs_E == 0
                C_all(t) = 0;
            else
                C_all(t) = sum_abs_imag_E / sum_abs_E;
%                 C_all(t) = sum_abs_imag_E;
            end
        else
            C_all(t) = 0;
        end
    end
    fprintf('%s: C 计算完成，用时 %.3f 秒。\n', case_id, toc(tic_C));

    %% 3.7 每窗打 5 阶段标签
    phase5 = repmat({''}, T, 1);

    phase5(pre_idx)   = {'pre-ictal'};
    phase5(early_idx) = {'early'};
    phase5(mid_idx)   = {'mid'};
    phase5(late_idx)  = {'late'};
    phase5(post_idx)  = {'post-ictal'};

    emptyMask = cellfun(@isempty, phase5);
    phase5(emptyMask) = {'unused'};

    %% 3.8 窗级结果表
    windowTable = table( ...
        repmat({case_id}, T, 1), ...
        repmat({FileName}, T, 1), ...
        time_axis, ...
        phase5, ...
        C_all, ...
        num_qi, ...
        'VariableNames', { ...
        'case_id', 'file_name', 'window_idx', 'phase5', ...
        'C', 'num_hyperedges'});

    allWindowTables{c} = windowTable;
    writetable(windowTable, fullfile(outDir, [case_id '_window_level_C.csv']));

    %% 3.9 5阶段结果表（每个病人）
    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    nPhase = numel(phaseOrder);

    phase_col = cell(nPhase,1);
    n_window_col = nan(nPhase,1);
    C_mean = nan(nPhase,1);
    C_std  = nan(nPhase,1);
    C_median = nan(nPhase,1);

    for p = 1:nPhase
        idx = strcmp(windowTable.phase5, phaseOrder{p});

        phase_col{p} = phaseOrder{p};
        n_window_col(p) = sum(idx);

        x = windowTable.C(idx);
        x = x(~isnan(x));

        if ~isempty(x)
            C_mean(p)   = mean(x);
            C_std(p)    = std(x);
            C_median(p) = median(x);
        end
    end

    stage5Table = table( ...
        repmat({case_id}, nPhase, 1), ...
        phase_col, ...
        n_window_col, ...
        C_mean, ...
        C_std, ...
        C_median, ...
        'VariableNames', { ...
        'case_id', 'phase5', 'n_windows', ...
        'C_mean', 'C_std', 'C_median'});

    allStage5Tables{c} = stage5Table;
    writetable(stage5Table, fullfile(outDir, [case_id '_stage5_C.csv']));

    %% 3.10 保存每个病人的 MAT
    save(fullfile(outDir, [case_id '_C_results.mat']), ...
        'case_id', 'FileName', 'fs', 'winLenSec', ...
        'pre_idx', 'ictal_idx', 'post_idx', ...
        'early_idx', 'mid_idx', 'late_idx', ...
        'all_hyperEdges_all', 'W_all', ...
        'C_all', 'windowTable', 'stage5Table');

    fprintf('%s: 结果保存完成。\n', case_id);
end

%% 4) 合并所有病例
validWindow = ~cellfun(@isempty, allWindowTables);
validStage5 = ~cellfun(@isempty, allStage5Tables);

if any(validWindow)
    ALL_windowTable = vertcat(allWindowTables{validWindow});
    writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_level_C.csv'));
else
    ALL_windowTable = table();
end

if any(validStage5)
    ALL_stage5Table = vertcat(allStage5Tables{validStage5});
    writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_C.csv'));
else
    ALL_stage5Table = table();
end

%% 5) 病人级组统计（5阶段）
groupDescTable = table();
groupPairwiseTable = table();
groupOmnibusTable = table();

if ~isempty(ALL_stage5Table)

    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    caseIDs = unique(ALL_stage5Table.case_id, 'stable');

    %% 5.1 描述统计
    desc_phase  = {};
    desc_n      = [];
    desc_mean   = [];
    desc_std    = [];
    desc_sem    = [];
    desc_median = [];

    for p = 1:numel(phaseOrder)
        idx = strcmp(ALL_stage5Table.phase5, phaseOrder{p});
        x = ALL_stage5Table.C_mean(idx);
        x = x(~isnan(x));

        desc_phase{end+1,1}  = phaseOrder{p};
        desc_n(end+1,1)      = numel(x);

        if isempty(x)
            desc_mean(end+1,1)   = nan;
            desc_std(end+1,1)    = nan;
            desc_sem(end+1,1)    = nan;
            desc_median(end+1,1) = nan;
        else
            desc_mean(end+1,1)   = mean(x);
            desc_std(end+1,1)    = std(x);
            desc_sem(end+1,1)    = std(x) / sqrt(numel(x));
            desc_median(end+1,1) = median(x);
        end
    end

    groupDescTable = table( ...
        desc_phase, desc_n, desc_mean, desc_std, desc_sem, desc_median, ...
        'VariableNames', {'phase5','n_case','mean','std','sem','median'});

    writetable(groupDescTable, fullfile(outDir, 'GROUP_stage5_C_descriptive_stats.csv'));

    %% 5.2 Friedman + 两两 signrank
    M = nan(numel(caseIDs), numel(phaseOrder));
    for i = 1:numel(caseIDs)
        rows_i = ALL_stage5Table(strcmp(ALL_stage5Table.case_id, caseIDs{i}), :);
        for p = 1:numel(phaseOrder)
            tmp = rows_i.C_mean(strcmp(rows_i.phase5, phaseOrder{p}));
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

    groupOmnibusTable = table( ...
        size(Mvalid,1), p_friedman, ...
        'VariableNames', {'n_case_complete','p_friedman'});

    writetable(groupOmnibusTable, fullfile(outDir, 'GROUP_stage5_C_friedman_stats.csv'));

    pair_comp   = {};
    pair_n      = [];
    pair_p      = [];
    pair_mean_diff = [];
    pair_median_diff = [];

    for p1 = 1:numel(phaseOrder)-1
        for p2 = p1+1:numel(phaseOrder)
            x1 = M(:, p1);
            x2 = M(:, p2);
            valid = ~isnan(x1) & ~isnan(x2);

            pair_comp{end+1,1} = [phaseOrder{p1} '_vs_' phaseOrder{p2}];
            pair_n(end+1,1) = sum(valid);

            if sum(valid) >= 2 && exist('signrank', 'file') == 2
                pair_p(end+1,1) = signrank(x1(valid), x2(valid));
                d = x2(valid) - x1(valid);
                pair_mean_diff(end+1,1) = mean(d);
                pair_median_diff(end+1,1) = median(d);
            else
                pair_p(end+1,1) = nan;
                d = x2(valid) - x1(valid);
                if isempty(d)
                    pair_mean_diff(end+1,1) = nan;
                    pair_median_diff(end+1,1) = nan;
                else
                    pair_mean_diff(end+1,1) = mean(d);
                    pair_median_diff(end+1,1) = median(d);
                end
            end
        end
    end

    groupPairwiseTable = table( ...
        pair_comp, pair_n, pair_p, pair_mean_diff, pair_median_diff, ...
        'VariableNames', {'comparison','n_pair','p_signrank','mean_diff','median_diff'});

    writetable(groupPairwiseTable, fullfile(outDir, 'GROUP_stage5_C_pairwise_signrank.csv'));

    %% 5.3 Nature 风格箱线图
    phaseColors = [
        76,120,168;
        89,161,79;
        242,142,43;
        225,87,89;
        128,115,172
        ] / 255;

    fig = figure('Color','w', 'Position',[100 100 520 420]);
    hold on;

    Y = cell(numel(phaseOrder), 1);
    maxLen = 0;
    for p = 1:numel(phaseOrder)
        idx = strcmp(ALL_stage5Table.phase5, phaseOrder{p});
        y = ALL_stage5Table.C_mean(idx);
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
    xticklabels(phaseOrder);
    xtickangle(25);
    ylabel('C', 'FontWeight', 'bold'); % 更新为 C

    Xall = [];
    for p = 1:numel(phaseOrder)
        Xall = [Xall; Y{p}(:)];
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

    exportgraphics(fig, fullfile(outDir, 'GROUP_stage5_C_boxplot.png'), 'Resolution', 600);
    savefig(fig, fullfile(outDir, 'GROUP_stage5_C_boxplot.fig'));
end

%% 6) 总 MAT 保存
save(fullfile(outDir, 'ALL_CASES_C_group_results.mat'), ...
    'caseInfo', ...
    'ALL_windowTable', ...
    'ALL_stage5Table', ...
    'groupDescTable', ...
    'groupOmnibusTable', ...
    'groupPairwiseTable');

fprintf('\n==================================================\n');
fprintf('全部处理完成。\n');
fprintf('结果已保存到:\n%s\n', outDir);
fprintf('==================================================\n');


%% 个例展示，lhs07
T = numel(C_all);
time_axis = 0:(T-1);

fig = figure('Color', 'w', 'Position', [100, 100, 820, 300]);
ax = axes(fig);
hold(ax, 'on');

plot(time_axis, C_all, 'k-', 'LineWidth', 1.6);

xline(201, '--', 'LineWidth', 1.1, 'Color', [0.45 0.45 0.45]);
xline(400, '--', 'LineWidth', 1.1, 'Color', [0.45 0.45 0.45]);

xlabel('Time windows', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('C', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold'); % 更新为 C

set(ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.2, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [0.015 0.015], ...
    'Layer', 'top');

xlim([time_axis(1), time_axis(end)]);

valid_C = C_all(~isnan(C_all) & isfinite(C_all));
if ~isempty(valid_C)
    ymin = min(valid_C);
    ymax = max(valid_C);
    if ymin == ymax
        pad = max(abs(ymin)*0.05, 0.01);
    else
        pad = 0.08 * (ymax - ymin);
    end
    ylim([ymin - pad, ymax + pad]);
end

ax.YAxis.Exponent = 0;

% 顶部阶段标注
yl = ylim;
y_text = yl(2) - 0.04 * (yl(2) - yl(1));
text(100, y_text, 'pre-ictal',  'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text(300, y_text, 'ictal',      'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text((400 + T - 1)/2, y_text, 'post-ictal', 'HorizontalAlignment', 'center', ...
    'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(outDir, 'lhs_cut07_C_all_timecourse_nature_labeled.png'), 'Resolution', 600);
savefig(fig, fullfile(outDir, 'lhs_cut07_C_all_timecourse_nature_labeled.fig'));