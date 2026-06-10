%%
clc; clear; close all;

%% =========================================================
% 批量运行 24 次发作：对比实验 (Comparison Study)
% 比较三种网络架构下的 eta 指标：
%   1. DMW-HLG (多阶混合：2节点 + 3节点)
%   2. Low-Order Baseline (低阶基线：仅保留 2节点)
%   3. High-Order (纯高阶对照：仅保留 3节点)
%
% 依赖：
%   gain_hyperEdges_23.m 需在 MATLAB 路径中
%
% eta 定义：
%   eta = ((1/N) * sum_{i,j} W_{ij}) / lambda_1
%   其中 N 为超边数，lambda_1 为 W 的主特征值（谱半径）
% =========================================================

%% 0) 输出目录
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig8';
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
caseInfo(1).pre_idx   = 1:200;  caseInfo(1).ictal_idx = 201:400; caseInfo(1).post_idx  = 401:598;

caseInfo(2).case_id   = 'lhs_cut06';
caseInfo(2).file_path = 'E:\wcldematlab\keep\new_idea\lhs_cut06\lihongsen_cut06_Gamma.mat';
caseInfo(2).pre_idx   = 1:122;  caseInfo(2).ictal_idx = 123:233; caseInfo(2).post_idx  = 234:338;

caseInfo(3).case_id   = 'gzs_cut06';
caseInfo(3).file_path = 'E:\wcldematlab\keep\new_idea\gzs\gongzishu_cut06_Gamma.mat';
caseInfo(3).pre_idx   = 1:195;  caseInfo(3).ictal_idx = 196:330; caseInfo(3).post_idx  = 331:477;

caseInfo(4).case_id   = 'gzs_cut07';
caseInfo(4).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut07\gongzishu_cut07_Gamma.mat';
caseInfo(4).pre_idx   = 1:137;  caseInfo(4).ictal_idx = 138:297; caseInfo(4).post_idx  = 298:398;

caseInfo(5).case_id   = 'gzs_cut08';
caseInfo(5).file_path = 'E:\wcldematlab\keep\new_idea\gzs_cut08\gongzishu_cut08_Gamma.mat';
caseInfo(5).pre_idx   = 1:252;  caseInfo(5).ictal_idx = 253:413; caseInfo(5).post_idx  = 414:498;

caseInfo(6).case_id   = 'wc_cut03';
caseInfo(6).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut03\wangchun_cut03_Gamma.mat';
caseInfo(6).pre_idx   = 1:199;  caseInfo(6).ictal_idx = 200:285; caseInfo(6).post_idx  = 286:398;

caseInfo(7).case_id   = 'wc_cut06';
caseInfo(7).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut06\wangchun_cut06_Gamma.mat';
caseInfo(7).pre_idx   = 1:239;  caseInfo(7).ictal_idx = 240:331; caseInfo(7).post_idx  = 332:498;

caseInfo(8).case_id   = 'wc_cut08';
caseInfo(8).file_path = 'E:\wcldematlab\keep\new_idea\wc_cut08\wangchun_cut08_Gamma.mat';
caseInfo(8).pre_idx   = 1:236;  caseInfo(8).ictal_idx = 237:313; caseInfo(8).post_idx  = 314:498;

caseInfo(9).case_id   = 'cxm_cut104';
caseInfo(9).file_path = 'E:\wcldematlab\keep\new_idea\cxm_cut104\cxm_cut104_Gamma.mat';
caseInfo(9).pre_idx   = 1:59;   caseInfo(9).ictal_idx = 60:112;  caseInfo(9).post_idx  = 113:385;

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
caseInfo(13).pre_idx   = 1:69;  caseInfo(13).ictal_idx = 70:130;  caseInfo(13).post_idx  = 131:198;

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
caseInfo(21).pre_idx   = 1:94;  caseInfo(21).ictal_idx = 95:170;  caseInfo(21).post_idx  = 171:318;

caseInfo(22).case_id   = 'fys_cut99';
caseInfo(22).file_path = 'E:\wcldematlab\keep\new_idea\fys_cut99\fys_cut99_Gamma.mat';
caseInfo(22).pre_idx   = 1:95;  caseInfo(22).ictal_idx = 96:195;  caseInfo(22).post_idx  = 196:350;

caseInfo(23).case_id   = 'gs_cut58';
caseInfo(23).file_path = 'E:\wcldematlab\keep\new_idea\gs_cut58\gushuai_cut58_Gamma.mat';
caseInfo(23).pre_idx   = 1:162;caseInfo(23).ictal_idx = 163:228; caseInfo(23).post_idx  = 229:310;

caseInfo(24).case_id   = 'gwh_cut104';
caseInfo(24).file_path = 'E:\wcldematlab\keep\new_idea\gwh\guan_cut104_Gamma.mat';
caseInfo(24).pre_idx   = 1:29;  caseInfo(24).ictal_idx = 30:95;   caseInfo(24).post_idx  = 96:148;

nCase = numel(caseInfo);

%% 3) 汇总容器
allWindowTables = cell(nCase, 1);
allStage5Tables = cell(nCase, 1);

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
    fprintf('(%d/%d) 开始处理对比实验: %s\n', c, nCase, case_id);
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
        continue;
    end
    T = totalIterations + 1;
    time_axis = (1:T)';

    %% 3.2 修正阶段索引
    pre_idx   = intersect(unique(pre_idx, 'stable'), 1:T);
    ictal_idx = setdiff(intersect(unique(ictal_idx, 'stable'), 1:T), pre_idx, 'stable');
    post_idx  = setdiff(intersect(unique(post_idx, 'stable'), 1:T), union(pre_idx, ictal_idx), 'stable');

    if isempty(ictal_idx)
        continue;
    end

    %% 3.3 ictal 三等分
    nIctal = numel(ictal_idx);
    splitPts = round(linspace(0, nIctal, 4));
    early_idx = ictal_idx(1:splitPts(2));
    mid_idx   = ictal_idx(splitPts(2)+1:splitPts(3));
    late_idx  = ictal_idx(splitPts(3)+1:end);

    %% 3.4 提取超边
    all_hyperEdges_all = cell(T, 1);
    tic_extract = tic;
    for time = 0:totalIterations
        datanew = X1(:, (fs * time + 1) : (fs * (time + winLenSec)));
        all_hyperEdges_all{time + 1} = gain_hyperEdges_23(datanew);
    end
    fprintf('%s: 超边提取完成，用时 %.3f 秒。\n', case_id, toc(tic_extract));

    %% 3.5 分别计算 3 种模式的 eta 指标
    eta_DMW   = nan(T,1); num_DMW   = nan(T,1);
    eta_Low2  = nan(T,1); num_Low2  = nan(T,1);
    eta_High3 = nan(T,1); num_High3 = nan(T,1);

    tic_s = tic;
    for t = 1:T
        HE = all_hyperEdges_all{t};
        if isempty(HE)
            continue;
        end
        sizes = cellfun(@length, HE);

        % 模式 1：DMW-HLG（多阶混合）
        [eta_DMW(t), num_DMW(t)] = calc_eta_from_HE(HE);

        % 模式 2：Low-Order（仅2节点）
        [eta_Low2(t), num_Low2(t)] = calc_eta_from_HE(HE(sizes == 2));

        % 模式 3：High-Order（仅3节点）
        [eta_High3(t), num_High3(t)] = calc_eta_from_HE(HE(sizes == 3));
    end
    fprintf('%s: 对比 eta 计算完成，用时 %.3f 秒。\n', case_id, toc(tic_s));

    %% 3.6 打标签
    phase5 = repmat({'unused'}, T, 1);
    phase5(pre_idx)   = {'pre-ictal'};
    phase5(early_idx) = {'early'};
    phase5(mid_idx)   = {'mid'};
    phase5(late_idx)  = {'late'};
    phase5(post_idx)  = {'post-ictal'};

    %% 3.7 构建窗级总表
    T_len = length(time_axis);

    tab_DMW = table(repmat({case_id}, T_len, 1), time_axis, phase5, ...
        repmat({'DMW-HLG'}, T_len, 1), eta_DMW, num_DMW, ...
        'VariableNames', {'case_id', 'window_idx', 'phase5', 'Model', 'eta', 'num_hyperedges'});

    tab_Low2 = table(repmat({case_id}, T_len, 1), time_axis, phase5, ...
        repmat({'Low-Order(2-nodes)'}, T_len, 1), eta_Low2, num_Low2, ...
        'VariableNames', {'case_id', 'window_idx', 'phase5', 'Model', 'eta', 'num_hyperedges'});

    tab_High3 = table(repmat({case_id}, T_len, 1), time_axis, phase5, ...
        repmat({'High-Order(3-nodes)'}, T_len, 1), eta_High3, num_High3, ...
        'VariableNames', {'case_id', 'window_idx', 'phase5', 'Model', 'eta', 'num_hyperedges'});

    windowTable = vertcat(tab_DMW, tab_Low2, tab_High3);
    allWindowTables{c} = windowTable;
    writetable(windowTable, fullfile(outDir, [case_id '_window_level_eta_Comparison.csv']));

    %% 3.8 构建阶段统计表
    phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
    models = {'DMW-HLG', 'Low-Order(2-nodes)', 'High-Order(3-nodes)'};

    stage_rows = [];
    for m = 1:length(models)
        for p = 1:length(phaseOrder)
            idx = strcmp(windowTable.phase5, phaseOrder{p}) & strcmp(windowTable.Model, models{m});
            x = windowTable.eta(idx);
            x = x(~isnan(x));

            if ~isempty(x)
                stage_rows = [stage_rows; {case_id, models{m}, phaseOrder{p}, sum(idx), mean(x), std(x), median(x)}]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(stage_rows)
        stage5Table = cell2table(stage_rows, 'VariableNames', ...
            {'case_id', 'Model', 'phase5', 'n_windows', 'eta_mean', 'eta_std', 'eta_median'});
        allStage5Tables{c} = stage5Table;
        writetable(stage5Table, fullfile(outDir, [case_id '_stage5_eta_Comparison.csv']));
    end
end

%% 4) 合并所有病例
validWindow = ~cellfun(@isempty, allWindowTables);
validStage5 = ~cellfun(@isempty, allStage5Tables);

ALL_windowTable = vertcat(allWindowTables{validWindow});
ALL_stage5Table = vertcat(allStage5Tables{validStage5});
writetable(ALL_windowTable, fullfile(outDir, 'ALL_CASES_window_level_eta_Comparison.csv'));
writetable(ALL_stage5Table, fullfile(outDir, 'ALL_CASES_stage5_eta_Comparison.csv'));

fprintf('\n对比实验批量计算完毕！查看文件夹：%s\n', outDir);

%% =========================================================
% 额外输出：方差统计表（供绘图脚本直接读取）
%% =========================================================
phaseOrder = {'pre-ictal', 'early', 'mid', 'late', 'post-ictal'};
modelOrder = {'Low-Order(2-nodes)', 'High-Order(3-nodes)', 'DMW-HLG'};

statsRows = [];
for m = 1:numel(modelOrder)
    for p = 1:numel(phaseOrder)
        idx = strcmp(ALL_stage5Table.Model, modelOrder{m}) & strcmp(ALL_stage5Table.phase5, phaseOrder{p});
        x = ALL_stage5Table.eta_mean(idx);
        x = x(~isnan(x));
        if ~isempty(x)
            statsRows = [statsRows; {modelOrder{m}, phaseOrder{p}, var(x, 'omitnan')}]; %#ok<AGROW>
        end
    end
end

if ~isempty(statsRows)
    StatsVarTable = cell2table(statsRows, 'VariableNames', {'Model', 'Phase', 'Variance'});
    writetable(StatsVarTable, fullfile(outDir, 'Stats_Variance_Comparison_eta.csv'));
end

%% =========================================================
% 局部函数：从超边集合计算 eta
%% =========================================================
function [eta_val, num_nodes] = calc_eta_from_HE(HE)
    num_HE = length(HE);
    num_nodes = num_HE;

    if num_HE < 2
        eta_val = 0;
        return;
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
    D_i = repmat(D, 1, num_HE);
    D_j = repmat(D', num_HE, 1);

    Coverage = O ./ Sizes_i;
    DegreeSum  = D_i + D_j;
    DegreeBias = D_j ./ DegreeSum;
    DegreeBias(DegreeSum == 0) = 0.5;

    W = Coverage;
    EqualSizeMask = (Sizes_i == Sizes_j) & (O > 0);
    W(EqualSizeMask) = W(EqualSizeMask) .* DegreeBias(EqualSizeMask);
    W(1:num_HE+1:end) = 0;

    W = double(W);
    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    W(1:num_HE+1:end) = 0;

    sum_W = sum(W(:));
    if sum_W <= 0
        eta_val = 0;
        return;
    end

    lambda_1 = max(abs(eig(W)));
    if ~isfinite(lambda_1) || lambda_1 <= eps
        eta_val = 0;
    else
        eta_val = (sum_W / num_HE) / lambda_1;
    end
end
