%%
clc; clear; close all;

%% =========================================================
% 实验 1：eta 分解实验
%
% 目标：
%   在原 calc_eta.m 的基础上，将 eta 分解为：
%
%       eta = R / lambda1
%
%   其中：
%       R       = sum(W(:)) / N
%       lambda1 = W 的谱半径
%
% 同时计算：
%       lambda2
%       spectral_gap = lambda1 - lambda2
%       lambda2_ratio = lambda2 / lambda1
%       spectral_PR = (sum(lambda_k)^2) / sum(lambda_k^2)
%       density_W = sum(W(:)) / (N*(N-1))
%
% 用于回答：
%   eta 的变化源于什么？
%   是传播资源 R 变化，还是主导谱模态 lambda1 变化？
%
% 依赖：
%   gain_hyperEdges_23.m 需在 MATLAB 路径中
% =========================================================

%% 0) 输出目录
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig3';
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

    %% 3.2 修正阶段索引
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

    %% 3.6 计算 eta 及分解指标
    eta_all           = nan(T,1);
    resource_R_all    = nan(T,1);
    sum_W_all         = nan(T,1);
    lambda1_all       = nan(T,1);
    lambda2_all       = nan(T,1);
    spectral_gap_all  = nan(T,1);
    lambda2_ratio_all = nan(T,1);
    spectral_PR_all   = nan(T,1);
    density_W_all     = nan(T,1);
    num_hyperedges_all = nan(T,1);

    tic_eta = tic;

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
        W(1:n+1:end) = 0;

        num_hyperedges_all(t) = n;

        if n < 2
            eta_all(t) = 0;
            resource_R_all(t) = 0;
            sum_W_all(t) = 0;
            lambda1_all(t) = 0;
            lambda2_all(t) = 0;
            spectral_gap_all(t) = 0;
            lambda2_ratio_all(t) = 0;
            spectral_PR_all(t) = 0;
            density_W_all(t) = 0;
            continue;
        end

        sum_W = sum(W(:));
        sum_W_all(t) = sum_W;

        % R = sum(W) / N
        resource_R = sum_W / n;
        resource_R_all(t) = resource_R;

        % 广义加权密度
        density_W_all(t) = sum_W / (n * (n - 1));

        if sum_W == 0
            eta_all(t) = 0;
            lambda1_all(t) = 0;
            lambda2_all(t) = 0;
            spectral_gap_all(t) = 0;
            lambda2_ratio_all(t) = 0;
            spectral_PR_all(t) = 0;
            continue;
        end

        % W 非对称，eig 可能为复数；使用 abs(eig(W)) 计算谱半径
        eigVals = eig(W);
        eigAbs = sort(abs(eigVals), 'descend');

        lambda1 = eigAbs(1);

        if numel(eigAbs) >= 2
            lambda2 = eigAbs(2);
        else
            lambda2 = 0;
        end

        lambda1_all(t) = lambda1;
        lambda2_all(t) = lambda2;

        if lambda1 > 0
            eta_all(t) = resource_R / lambda1;
            spectral_gap_all(t) = lambda1 - lambda2;
            lambda2_ratio_all(t) = lambda2 / lambda1;
        else
            eta_all(t) = 0;
            spectral_gap_all(t) = 0;
            lambda2_ratio_all(t) = 0;
        end

        % 谱参与比：越大说明越多谱模态参与
        eigAbs_nonzero = eigAbs(eigAbs > 0);
        if isempty(eigAbs_nonzero)
            spectral_PR_all(t) = 0;
        else
            spectral_PR_all(t) = (sum(eigAbs_nonzero)^2) / sum(eigAbs_nonzero.^2);
        end
    end

    fprintf('%s: eta 分解指标计算完成，用时 %.3f 秒。\n', case_id, toc(tic_eta));

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
        eta_all, ...
        resource_R_all, ...
        lambda1_all, ...
        lambda2_all, ...
        spectral_gap_all, ...
        lambda2_ratio_all, ...
        spectral_PR_all, ...
        sum_W_all, ...
        density_W_all, ...
        num_hyperedges_all, ...
        'VariableNames', { ...
        'case_id', ...
        'file_name', ...
        'window_idx', ...
        'phase5', ...
        'eta', ...
        'resource_R', ...
        'lambda1', ...
        'lambda2', ...
        'spectral_gap', ...
        'lambda2_ratio', ...
        'spectral_PR', ...
        'sum_W', ...
        'density_W', ...
        'num_hyperedges'});

    allWindowTables{c} = windowTable;
    writetable(windowTable, fullfile(outDir, [case_id '_window_level_eta_decomposition.csv']));

    %% 3.9 5 阶段结果表
    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    nPhase = numel(phaseOrder);

    metricList = { ...
        'eta', ...
        'resource_R', ...
        'lambda1', ...
        'lambda2', ...
        'spectral_gap', ...
        'lambda2_ratio', ...
        'spectral_PR', ...
        'sum_W', ...
        'density_W', ...
        'num_hyperedges'};

    phase_col = cell(nPhase,1);
    n_window_col = nan(nPhase,1);

    stageStats = struct();
    for m = 1:numel(metricList)
        metric = metricList{m};
        stageStats.([metric '_mean']) = nan(nPhase,1);
        stageStats.([metric '_std']) = nan(nPhase,1);
        stageStats.([metric '_median']) = nan(nPhase,1);
    end

    for p = 1:nPhase
        idx = strcmp(windowTable.phase5, phaseOrder{p});

        phase_col{p} = phaseOrder{p};
        n_window_col(p) = sum(idx);

        for m = 1:numel(metricList)
            metric = metricList{m};
            x = windowTable.(metric)(idx);
            x = x(~isnan(x) & isfinite(x));

            if ~isempty(x)
                stageStats.([metric '_mean'])(p)   = mean(x);
                stageStats.([metric '_std'])(p)    = std(x);
                stageStats.([metric '_median'])(p) = median(x);
            end
        end
    end

    stage5Table = table( ...
        repmat({case_id}, nPhase, 1), ...
        phase_col, ...
        n_window_col, ...
        'VariableNames', { ...
        'case_id', ...
        'phase5', ...
        'n_windows'});

    for m = 1:numel(metricList)
        metric = metricList{m};
        stage5Table.([metric '_mean'])   = stageStats.([metric '_mean']);
        stage5Table.([metric '_std'])    = stageStats.([metric '_std']);
        stage5Table.([metric '_median']) = stageStats.([metric '_median']);
    end

    allStage5Tables{c} = stage5Table;
    writetable(stage5Table, fullfile(outDir, [case_id '_stage5_eta_decomposition.csv']));

    %% 3.10 保存每个病人的 MAT
    save(fullfile(outDir, [case_id '_eta_decomposition_results.mat']), ...
        'case_id', ...
        'FileName', ...
        'fs', ...
        'winLenSec', ...
        'pre_idx', ...
        'ictal_idx', ...
        'post_idx', ...
        'early_idx', ...
        'mid_idx', ...
        'late_idx', ...
        'all_hyperEdges_all', ...
        'W_all', ...
        'eta_all', ...
        'resource_R_all', ...
        'lambda1_all', ...
        'lambda2_all', ...
        'spectral_gap_all', ...
        'lambda2_ratio_all', ...
        'spectral_PR_all', ...
        'sum_W_all', ...
        'density_W_all', ...
        'num_hyperedges_all', ...
        'windowTable', ...
        'stage5Table');

    fprintf('%s: 结果保存完成。\n', case_id);
end

%% 4) 合并所有病例
validWindow = ~cellfun(@isempty, allWindowTables);
validStage5 = ~cellfun(@isempty, allStage5Tables);

if any(validWindow)
    ALL_windowTable = vertcat(allWindowTables{validWindow});
    writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_level_eta_decomposition.csv'));
else
    ALL_windowTable = table();
end

if any(validStage5)
    ALL_stage5Table = vertcat(allStage5Tables{validStage5});
    writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_eta_decomposition.csv'));
else
    ALL_stage5Table = table();
end

%% 5) 组水平统计
groupDescTable = table();
groupOmnibusTable = table();
groupPairwiseTable = table();
ALL_decompStats = table();

if ~isempty(ALL_stage5Table)

    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    caseIDs = unique(ALL_stage5Table.case_id, 'stable');

    %% 5.1 描述统计
    decompMetricList = { ...
        'eta_mean', ...
        'resource_R_mean', ...
        'lambda1_mean', ...
        'lambda2_mean', ...
        'spectral_gap_mean', ...
        'lambda2_ratio_mean', ...
        'spectral_PR_mean', ...
        'sum_W_mean', ...
        'density_W_mean', ...
        'num_hyperedges_mean'};

    descRows = table();

    for mm = 1:numel(decompMetricList)

        metricName = decompMetricList{mm};

        if ~ismember(metricName, ALL_stage5Table.Properties.VariableNames)
            continue;
        end

        for p = 1:numel(phaseOrder)
            idx = strcmp(ALL_stage5Table.phase5, phaseOrder{p});
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

            tmpRow = table( ...
                {metricName}, ...
                {phaseOrder{p}}, ...
                n_case, ...
                mean_x, ...
                std_x, ...
                sem_x, ...
                median_x, ...
                'VariableNames', { ...
                'metric', ...
                'phase5', ...
                'n_case', ...
                'mean', ...
                'std', ...
                'sem', ...
                'median'});

            descRows = [descRows; tmpRow];
        end
    end

    groupDescTable = descRows;
    writetable(groupDescTable, fullfile(outDir, 'GROUP_eta_decomposition_descriptive_stats.csv'));

    %% 5.2 Friedman + pre 对其他阶段 signrank
    for mm = 1:numel(decompMetricList)

        metricName = decompMetricList{mm};

        if ~ismember(metricName, ALL_stage5Table.Properties.VariableNames)
            warning('指标 %s 不存在，跳过。', metricName);
            continue;
        end

        Mmetric = nan(numel(caseIDs), numel(phaseOrder));

        for i = 1:numel(caseIDs)
            rows_i = ALL_stage5Table(strcmp(ALL_stage5Table.case_id, caseIDs{i}), :);

            for p = 1:numel(phaseOrder)
                tmp = rows_i.(metricName)(strcmp(rows_i.phase5, phaseOrder{p}));
                if ~isempty(tmp)
                    Mmetric(i,p) = tmp(1);
                end
            end
        end

        validRows = all(~isnan(Mmetric), 2);
        Mvalid = Mmetric(validRows, :);

        if size(Mvalid,1) >= 2 && exist('friedman', 'file') == 2
            p_friedman = friedman(Mvalid, 1, 'off');
        else
            p_friedman = nan;
        end

        tmpOmni = table( ...
            {metricName}, ...
            size(Mvalid,1), ...
            p_friedman, ...
            'VariableNames', { ...
            'metric', ...
            'n_case_complete', ...
            'p_friedman_all_phase'});

        groupOmnibusTable = [groupOmnibusTable; tmpOmni];

        preCol = 1;

        for p = 2:numel(phaseOrder)

            x_pre = Mmetric(:, preCol);
            x_cmp = Mmetric(:, p);

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

            tmpPair = table( ...
                {metricName}, ...
                {['pre-ictal_vs_' phaseOrder{p}]}, ...
                sum(valid), ...
                p_signrank, ...
                mean_diff, ...
                median_diff, ...
                'VariableNames', { ...
                'metric', ...
                'comparison', ...
                'n_pair', ...
                'p_signrank', ...
                'mean_diff', ...
                'median_diff'});

            groupPairwiseTable = [groupPairwiseTable; tmpPair];

            tmpAll = table( ...
                {metricName}, ...
                {['pre-ictal_vs_' phaseOrder{p}]}, ...
                sum(valid), ...
                p_friedman, ...
                p_signrank, ...
                mean_diff, ...
                median_diff, ...
                'VariableNames', { ...
                'metric', ...
                'comparison', ...
                'n_pair', ...
                'p_friedman_all_phase', ...
                'p_signrank', ...
                'mean_diff', ...
                'median_diff'});

            ALL_decompStats = [ALL_decompStats; tmpAll];
        end
    end

    writetable(groupOmnibusTable, fullfile(outDir, 'GROUP_eta_decomposition_friedman_stats.csv'));
    writetable(groupPairwiseTable, fullfile(outDir, 'GROUP_eta_decomposition_pairwise_signrank.csv'));
    writetable(ALL_decompStats, fullfile(outDir, 'GROUP_eta_decomposition_all_stats.csv'));

    %% 5.3 画图：eta / R / lambda1 / spectral_gap / lambda2_ratio / spectral_PR
    plotMetricList = { ...
        'eta_mean', ...
        'resource_R_mean', ...
        'lambda1_mean', ...
        'spectral_gap_mean', ...
        'lambda2_ratio_mean', ...
        'spectral_PR_mean'};

    plotYLabels = { ...
        '\eta', ...
        'R = \Sigma W / N', ...
        '\lambda_1', ...
        '\lambda_1 - \lambda_2', ...
        '\lambda_2 / \lambda_1', ...
        'Spectral PR'};

    phaseColors = [
        76,120,168;
        89,161,79;
        242,142,43;
        225,87,89;
        128,115,172
        ] / 255;

    for mm = 1:numel(plotMetricList)

        metricName = plotMetricList{mm};

        if ~ismember(metricName, ALL_stage5Table.Properties.VariableNames)
            warning('绘图指标 %s 不存在，跳过。', metricName);
            continue;
        end

        fig = figure('Color','w', 'Position',[100 100 520 420]);
        hold on;

        Y = cell(numel(phaseOrder), 1);
        maxLen = 0;

        for p = 1:numel(phaseOrder)
            idx = strcmp(ALL_stage5Table.phase5, phaseOrder{p});
            y = ALL_stage5Table.(metricName)(idx);
            y = y(~isnan(y) & isfinite(y));

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

        ylabel(plotYLabels{mm}, 'FontWeight', 'bold');

        Xall = [];
        for p = 1:numel(phaseOrder)
            Xall = [Xall; Y{p}(:)];
        end
        Xall = Xall(~isnan(Xall) & isfinite(Xall));

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

        outName = ['GROUP_stage5_' metricName '_boxplot'];
        exportgraphics(fig, fullfile(outDir, [outName '.png']), 'Resolution', 600);
        savefig(fig, fullfile(outDir, [outName '.fig']));
    end

    %% 5.4 一张 3 联图：eta / R / lambda1
    tripleMetrics = {'eta_mean', 'resource_R_mean', 'lambda1_mean'};
    tripleLabels = {'\eta', 'R = \Sigma W / N', '\lambda_1'};

    fig = figure('Color','w', 'Position',[100 100 1200 360]);

    for mm = 1:3
        metricName = tripleMetrics{mm};

        subplot(1,3,mm);
        hold on;

        Y = cell(numel(phaseOrder), 1);
        maxLen = 0;

        for p = 1:numel(phaseOrder)
            idx = strcmp(ALL_stage5Table.phase5, phaseOrder{p});
            y = ALL_stage5Table.(metricName)(idx);
            y = y(~isnan(y) & isfinite(y));

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

        for p = 1:numel(phaseOrder)
            y = Y{p};
            if isempty(y)
                continue;
            end

            jitter = (rand(numel(y),1) - 0.5) * 0.18;
            scatter( ...
                p + jitter, y, 22, ...
                'MarkerFaceColor', phaseColors(p,:), ...
                'MarkerEdgeColor', 'w', ...
                'LineWidth', 0.45, ...
                'MarkerFaceAlpha', 0.9, ...
                'MarkerEdgeAlpha', 0.95);
        end

        set(gca, ...
            'FontName', 'Arial', ...
            'FontSize', 10, ...
            'LineWidth', 1.1, ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'Layer', 'top');

        xticks(1:numel(phaseOrder));
        xticklabels(phaseOrder);
        xtickangle(25);
        ylabel(tripleLabels{mm}, 'FontWeight', 'bold');

        title(strrep(metricName, '_mean', ''), ...
            'FontName', 'Arial', ...
            'FontSize', 11, ...
            'FontWeight', 'bold');
    end

    exportgraphics(fig, fullfile(outDir, 'GROUP_eta_R_lambda1_triple_boxplot.png'), 'Resolution', 600);
    savefig(fig, fullfile(outDir, 'GROUP_eta_R_lambda1_triple_boxplot.fig'));
end

%% 6) 保存总 MAT
save(fullfile(outDir, 'ALL_CASES_eta_decomposition_group_results.mat'), ...
    'caseInfo', ...
    'ALL_windowTable', ...
    'ALL_stage5Table', ...
    'groupDescTable', ...
    'groupOmnibusTable', ...
    'groupPairwiseTable', ...
    'ALL_decompStats');

fprintf('\n==================================================\n');
fprintf('全部处理完成。\n');
fprintf('eta 分解实验结果已保存到:\n%s\n', outDir);
fprintf('==================================================\n');

