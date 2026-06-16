%% calc_eta_exp4_region_hyperedge_burden_size_normalized.m
% =========================================================
% 实验 4：区域与超边贡献分析（扣除区域大小差异）
%
% 核心改动：
%   1) 保留原始 raw 区域贡献；
%   2) 新增按区域总通道数归一化的 size-adjusted 指标（per-node）；
%   3) 新增相对区域基线占比的 enrichment 指标：share / region_size_fraction；
%   4) 组图 A 改为显示 size-adjusted 的区域 eta burden（per-node）。
%
% 问题：eta 的变化由谁承担？
%   1) SOZ / PZ / NIZ 区域是否共同承担 eta 变化？
%   2) 宏观状态组合，如 SOZ_PZ、SOZ_NIZ、PZ_NIZ、SOZ_PZ_NIZ 是否承担更大贡献？
%   3) 少数 top contribution hyperedges 是否承担关键的 eta 重分配作用？
%   4) source-like / sink-like / balanced 超边在不同阶段如何改变？
%
% 输入：calc_eta.m 输出的每个病例 MAT：
%       *_eta_results.mat
% 需要变量：
%       W_all, all_hyperEdges_all, eta_all, case_id,
%       pre_idx, early_idx, mid_idx, late_idx, post_idx
%
% 区域编码：
%       1 = SOZ
%       2 = PZ
%       3 = NIZ
%
% 主要输出：
%   1) CASE_region_map_SOZ_PZ_NIZ_embedded.csv
%   2) ALL_CASES_exp4_hyperedge_contribution.csv
%   3) ALL_CASES_exp4_window_region_contribution_size_adjusted.csv
%   4) ALL_CASES_exp4_stage_region_contribution_size_adjusted.csv
%   5) ALL_CASES_exp4_window_macro_role_contribution.csv
%   6) ALL_CASES_exp4_stage_macro_role_contribution.csv
%   7) TOP_exp4_window_hyperedges.csv
%   8) TOP_exp4_stage_hyperedges.csv
%   9) GROUP_exp4_region_phase_stats_size_adjusted.csv
%  10) GROUP_exp4_macro_role_phase_stats.csv
%  11) FIG_exp4_eta_burden_summary_size_adjusted.png / .pdf / .fig
%
% 说明：
%   本脚本不重新计算 W，也不重新计算 eta 主流程；它直接读取 calc_eta.m
%   已保存的 W_all 与 all_hyperEdges_all，对 eta 的区域/超边承担机制做解释性分解。
%
%   区域层面同时输出两类结果：
%   A. raw（原始总承担量）
%   B. size-adjusted（扣除区域大小差异）
%      - per-node：除以该病例该区域总通道数
%      - enrichment：份额 / 区域通道占比
% =========================================================

clc; clear; close all;

%% 1) 路径设置
inDir  = 'example_project\2_Fig1';
outDir = 'example_project\2_Fig5_size_adjusted';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

matFiles = dir(fullfile(inDir, '*_eta_results.mat'));
if isempty(matFiles)
    error('在 inDir 中没有找到 *_eta_results.mat。请先运行 calc_eta.m。');
end

%% 2) 参数
phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
phaseShort = {'Pre','Early','Mid','Late','Post'};
regionNames = {'SOZ','PZ','NIZ'};
macroOrder = {'SOZ_only','PZ_only','NIZ_only', ...
              'SOZ_PZ','SOZ_NIZ','PZ_NIZ','SOZ_PZ_NIZ'};
roleOrder = {'source-like','balanced','sink-like'};

% source/sink 判定阈值：
% flow_bias = (out_strength - in_strength) / (out_strength + in_strength)
sourceThr = 0.20;
sinkThr   = -0.20;

% 每个时间窗保留 top-K 超边；每个病例-阶段再汇总 top-K 超边
topKWindow = 10;
topKStage  = 10;

% 是否做 leave-one-out eta 删除检验。
% true 更直观，但会慢很多；建议初次 false，最终补充分析可 true。
doLeaveOneOut = false;

% 小数保护
EPS0 = 1e-12;

%% 3) 内嵌 SOZ/PZ/NIZ 区域映射
% 每个向量的长度必须等于对应病例的原始通道数。
% 1 = SOZ, 2 = PZ, 3 = NIZ
regionMap = build_embedded_region_map();

% 导出一份区域映射表，方便核对
regionMapTable = export_region_map_table(regionMap);
writetable(regionMapTable, fullfile(outDir, 'CASE_region_map_SOZ_PZ_NIZ_embedded.csv'));

%% 4) 汇总容器
allHyperedgeRows = {};
allWindowRegionRows = {};
allWindowMacroRoleRows = {};
allTopWindowRows = {};
allTopStageRows = {};

fprintf('\n==================================================\n');
fprintf('实验 4：区域与超边贡献分析（size-adjusted）开始\n');
fprintf('输入目录: %s\n', inDir);
fprintf('输出目录: %s\n', outDir);
fprintf('病例数: %d\n', numel(matFiles));
fprintf('==================================================\n');

%% =========================================================
% 主循环：逐病例读取 calc_eta 输出结果
%% =========================================================
for f = 1:numel(matFiles)

    matPath = fullfile(matFiles(f).folder, matFiles(f).name);
    S = load(matPath);

    if isfield(S, 'case_id')
        case_id = char(S.case_id);
    else
        [~, case_id] = fileparts(matFiles(f).name);
        case_id = strrep(case_id, '_eta_results', '');
    end

    fprintf('\n--------------------------------------------------\n');
    fprintf('(%d/%d) 处理病例: %s\n', f, numel(matFiles), case_id);
    fprintf('文件: %s\n', matFiles(f).name);

    requiredVars = {'W_all','all_hyperEdges_all','eta_all', ...
                    'pre_idx','early_idx','mid_idx','late_idx','post_idx'};
    miss = requiredVars(~isfield_multi(S, requiredVars));
    if ~isempty(miss)
        warning('%s 缺少变量: %s，跳过。', case_id, strjoin(miss, ', '));
        continue;
    end

    if ~isKey(regionMap, case_id)
        warning('%s 没有区域映射，跳过。请在 build_embedded_region_map() 中补充。', case_id);
        continue;
    end

    node_regions = regionMap(case_id);
    node_regions = node_regions(:);

    regionNodeCount = [sum(node_regions == 1); sum(node_regions == 2); sum(node_regions == 3)];
    totalNodeCount = numel(node_regions);
    regionNodeFraction = regionNodeCount ./ max(totalNodeCount, 1);

    W_all = S.W_all;
    HE_all = S.all_hyperEdges_all;
    eta_all = S.eta_all(:);
    T = numel(W_all);

    if numel(HE_all) ~= T
        warning('%s: all_hyperEdges_all 与 W_all 长度不一致，按较短长度处理。', case_id);
        T = min(numel(HE_all), T);
    end

    %% 4.1 阶段标签
    phase5 = repmat({'unused'}, T, 1);

    pre_idx   = S.pre_idx(:)';
    early_idx = S.early_idx(:)';
    mid_idx   = S.mid_idx(:)';
    late_idx  = S.late_idx(:)';
    post_idx  = S.post_idx(:)';

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

    %% 4.2 逐时间窗计算超边贡献
    caseHyperRows = {};
    caseWindowRegionRows = {};
    caseWindowMacroRoleRows = {};
    caseTopWindowRows = {};

    tic_case = tic;

    for t = 1:T

        W = sanitize_W(W_all{t});
        HE = HE_all{t};

        if isempty(W) || isempty(HE)
            continue;
        end

        nHE = size(W, 1);
        if numel(HE) ~= nHE
            warning('%s window %d: 超边数和 W 尺寸不一致，按较短长度截断。', case_id, t);
            nUse = min(numel(HE), nHE);
            W = W(1:nUse, 1:nUse);
            HE = HE(1:nUse);
            nHE = nUse;
        end

        if nHE == 0
            continue;
        end

        eta0 = compute_eta_from_W(W);
        if t <= numel(eta_all) && isfinite(eta_all(t))
            eta_saved = eta_all(t);
        else
            eta_saved = eta0;
        end

        M = compute_hyperedge_metrics(W, eta0, doLeaveOneOut, EPS0);

        % 每个超边的区域组合、宏观状态、source/sink 角色
        HEinfo = parse_hyperedges_region(HE, node_regions);
        roleLabel = classify_role(M.flow_bias, sourceThr, sinkThr);

        % 4.2.1 超边级明细表
        for i = 1:nHE
            row = {case_id, t, phase5{t}, i, HEinfo(i).hyperedge_key, ...
                   HEinfo(i).n_nodes, HEinfo(i).n_SOZ, HEinfo(i).n_PZ, HEinfo(i).n_NIZ, ...
                   HEinfo(i).prop_SOZ, HEinfo(i).prop_PZ, HEinfo(i).prop_NIZ, ...
                   HEinfo(i).macro_state, roleLabel{i}, ...
                   M.out_strength(i), M.in_strength(i), M.total_strength(i), M.flow_bias(i), ...
                   M.resource_share(i), M.lambda_share(i), ...
                   M.eta_resource_contrib(i), M.eta_lambda_contrib(i), ...
                   M.eta_proxy_delta(i), abs(M.eta_proxy_delta(i)), ...
                   M.eta_leave_one_out_delta(i), eta0, eta_saved};
            caseHyperRows(end+1, :) = row; %#ok<SAGROW>
        end

        % 4.2.2 窗级：SOZ/PZ/NIZ 区域承担（raw + size-adjusted）
        regionProps = [[HEinfo.prop_SOZ]', [HEinfo.prop_PZ]', [HEinfo.prop_NIZ]'];
        for r = 1:3
            wr = regionProps(:, r);
            row = make_window_region_row(case_id, t, phase5{t}, regionNames{r}, ...
                eta0, eta_saved, regionNodeCount(r), regionNodeFraction(r), wr, M, EPS0);
            caseWindowRegionRows(end+1, :) = row; %#ok<SAGROW>
        end

        % 4.2.3 窗级：宏观状态组合 × source/sink/balanced
        for m = 1:numel(macroOrder)
            for rr = 1:numel(roleOrder)
                idx = strcmp({HEinfo.macro_state}', macroOrder{m}) & strcmp(roleLabel(:), roleOrder{rr});
                row = make_window_macro_role_row(case_id, t, phase5{t}, macroOrder{m}, roleOrder{rr}, eta0, eta_saved, idx, M, HEinfo, EPS0);
                caseWindowMacroRoleRows(end+1, :) = row; %#ok<SAGROW>
            end
        end

        % 4.2.4 每窗 top contribution hyperedges
        score = abs(M.eta_proxy_delta);
        score(~isfinite(score)) = -inf;
        [~, ord] = sort(score, 'descend');
        ord = ord(1:min(topKWindow, numel(ord)));
        for kk = 1:numel(ord)
            i = ord(kk);
            row = {case_id, t, phase5{t}, kk, i, HEinfo(i).hyperedge_key, ...
                   HEinfo(i).macro_state, roleLabel{i}, HEinfo(i).n_nodes, ...
                   HEinfo(i).prop_SOZ, HEinfo(i).prop_PZ, HEinfo(i).prop_NIZ, ...
                   M.eta_proxy_delta(i), abs(M.eta_proxy_delta(i)), ...
                   M.eta_resource_contrib(i), M.eta_lambda_contrib(i), ...
                   M.flow_bias(i), M.out_strength(i), M.in_strength(i), ...
                   M.eta_leave_one_out_delta(i), eta0};
            caseTopWindowRows(end+1, :) = row; %#ok<SAGROW>
        end
    end

    fprintf('%s: 超边机制分析完成，用时 %.2f 秒。\n', case_id, toc(tic_case));

    %% 4.3 病例内 top stage hyperedges：按 case × phase × hyperedge_key 汇总
    if ~isempty(caseHyperRows)
        Hcase = cell2table(caseHyperRows, 'VariableNames', hyperedge_varnames());
        TopStageCase = summarize_top_stage_hyperedges(Hcase, phaseOrder, topKStage);
        if ~isempty(TopStageCase)
            allTopStageRows = [allTopStageRows; table2cell(TopStageCase)]; %#ok<AGROW>
        end
    end

    allHyperedgeRows = [allHyperedgeRows; caseHyperRows]; %#ok<AGROW>
    allWindowRegionRows = [allWindowRegionRows; caseWindowRegionRows]; %#ok<AGROW>
    allWindowMacroRoleRows = [allWindowMacroRoleRows; caseWindowMacroRoleRows]; %#ok<AGROW>
    allTopWindowRows = [allTopWindowRows; caseTopWindowRows]; %#ok<AGROW>
end

%% 5) 明细表输出
if isempty(allHyperedgeRows)
    error('没有得到任何超边级结果。请检查输入 MAT 和区域映射。');
end

ALL_hyperedgeTable = cell2table(allHyperedgeRows, 'VariableNames', hyperedge_varnames());
ALL_windowRegionTable = cell2table(allWindowRegionRows, 'VariableNames', window_region_varnames());
ALL_windowMacroRoleTable = cell2table(allWindowMacroRoleRows, 'VariableNames', window_macro_role_varnames());
TOP_windowTable = cell2table(allTopWindowRows, 'VariableNames', top_window_varnames());

writetable(ALL_hyperedgeTable, fullfile(outDir, 'ALL_CASES_exp4_hyperedge_contribution.csv'));
writetable(ALL_windowRegionTable, fullfile(outDir, 'ALL_CASES_exp4_window_region_contribution_size_adjusted.csv'));
writetable(ALL_windowMacroRoleTable, fullfile(outDir, 'ALL_CASES_exp4_window_macro_role_contribution.csv'));
writetable(TOP_windowTable, fullfile(outDir, 'TOP_exp4_window_hyperedges.csv'));

if ~isempty(allTopStageRows)
    TOP_stageTable = cell2table(allTopStageRows, 'VariableNames', top_stage_varnames());
else
    TOP_stageTable = table();
end
writetable(TOP_stageTable, fullfile(outDir, 'TOP_exp4_stage_hyperedges.csv'));

%% 6) 阶段级汇总：case × phase × region / macro-role
ALL_stageRegionTable = summarize_stage_region(ALL_windowRegionTable, phaseOrder, regionNames);
ALL_stageMacroRoleTable = summarize_stage_macro_role(ALL_windowMacroRoleTable, phaseOrder, macroOrder, roleOrder);

writetable(ALL_stageRegionTable, fullfile(outDir, 'ALL_CASES_exp4_stage_region_contribution_size_adjusted.csv'));
writetable(ALL_stageMacroRoleTable, fullfile(outDir, 'ALL_CASES_exp4_stage_macro_role_contribution.csv'));

%% 7) 组水平统计：跨病例均值、SEM、Pre 差值、signrank
GROUP_regionStats = summarize_group_region(ALL_stageRegionTable, phaseOrder, regionNames);
GROUP_macroRoleStats = summarize_group_macro_role(ALL_stageMacroRoleTable, phaseOrder, macroOrder, roleOrder);

writetable(GROUP_regionStats, fullfile(outDir, 'GROUP_exp4_region_phase_stats_size_adjusted.csv'));
writetable(GROUP_macroRoleStats, fullfile(outDir, 'GROUP_exp4_macro_role_phase_stats.csv'));

%% 8) 主要总结图：一张图，三个机制面板
plot_exp4_summary_figure_size_adjusted(ALL_stageRegionTable, ALL_stageMacroRoleTable, ...
    phaseOrder, phaseShort, regionNames, macroOrder, roleOrder, outDir);

%% 9) 保存 MAT
save(fullfile(outDir, 'ALL_CASES_exp4_region_hyperedge_contribution_results_size_adjusted.mat'), ...
    'ALL_hyperedgeTable', ...
    'ALL_windowRegionTable', ...
    'ALL_stageRegionTable', ...
    'ALL_windowMacroRoleTable', ...
    'ALL_stageMacroRoleTable', ...
    'TOP_windowTable', ...
    'TOP_stageTable', ...
    'GROUP_regionStats', ...
    'GROUP_macroRoleStats', ...
    'regionMapTable', ...
    'sourceThr', 'sinkThr', 'topKWindow', 'topKStage', 'doLeaveOneOut');

fprintf('\n==================================================\n');
fprintf('实验 4（size-adjusted）完成。\n');
fprintf('结果保存到:\n%s\n', outDir);
fprintf('==================================================\n');

%% =========================================================
% Local functions
%% =========================================================

function regionMap = build_embedded_region_map()
    % 1 = SOZ, 2 = PZ, 3 = NIZ
    regionMap = containers.Map();

    lhs = [3;3;2;2;1;2;2;2;1;3;3;3;2;2;2;2;2;3];
    gzs = [2;3;3;3;3;3;2;3;3;3;3;2;1;2;2;3];
    wc  = [1;2;1;3;2;2;2;2;3;3;2;3];
    cxm = [1;3;1;3;2;3;3;3;3;3;2;2;2;3;3;3;2;3];
    bsp = [1;3;3;3;2;2;2;3;3;3;2;3];
    zzy = [3;2;2;1;2;3;3;3;3;1;3;1;3;2];
    lzq = [1;3;1;3;2;2;2;1;2;2;2;1];
    ll  = [3;1;2;1;2;2;1;2;1;2;1;2;3;1];
    gzw = [2;1;3;2;1;2;1;2;2;1];
    hds = [2;3;3;3;3;3;2;3;3;3;3;2;1;2;3;3];
    ssh = [1;1;2;2;2;2;3;2;2;2;3;3;2;3];
    fys = [2;2;2;2;2;3;3;3;2;2;1;2;2;3;1;3];
    gs  = [2;3;2;3;1;3;2;3;2;3;3;3;3;3];
    gwh = [1;3;1;3;2;3;2;3;3;3;2;3;3;3;3;3];

    % 同一病人多次发作使用同一套区域映射
    regionMap('seizure_01') = lhs;
    regionMap('seizure_02') = lhs;

    regionMap('seizure_03') = gzs;
    regionMap('seizure_04') = gzs;
    regionMap('seizure_05') = gzs;

    regionMap('seizure_06') = wc;
    regionMap('seizure_07') = wc;
    regionMap('seizure_08') = wc;

    regionMap('seizure_09') = cxm;
    regionMap('seizure_10')  = bsp;

    regionMap('seizure_11') = zzy;
    regionMap('seizure_12') = zzy;

    regionMap('seizure_13') = lzq;
    regionMap('seizure_14') = lzq;

    regionMap('seizure_15')  = ll;
    regionMap('seizure_16')   = gzw;
    regionMap('seizure_17') = hds;

    regionMap('seizure_18') = ssh;
    regionMap('seizure_19') = ssh;
    regionMap('seizure_20') = ssh;

    regionMap('seizure_21') = fys;
    regionMap('seizure_22') = fys;

    regionMap('seizure_23')   = gs;
    regionMap('seizure_24') = gwh;
end

function tf = isfield_multi(S, names)
    tf = false(size(names));
    for i = 1:numel(names)
        tf(i) = isfield(S, names{i});
    end
end

function T = export_region_map_table(regionMap)
    keys0 = keys(regionMap);
    case_id = {};
    channel_idx = [];
    region_code = [];
    region_label = {};

    for k = 1:numel(keys0)
        cid = keys0{k};
        r = regionMap(cid);
        for i = 1:numel(r)
            case_id{end+1,1} = cid; %#ok<AGROW>
            channel_idx(end+1,1) = i; %#ok<AGROW>
            region_code(end+1,1) = r(i); %#ok<AGROW>
            region_label{end+1,1} = code_to_region(r(i)); %#ok<AGROW>
        end
    end

    T = table(case_id, channel_idx, region_code, region_label, ...
        'VariableNames', {'case_id','channel_idx','region_code','region_label'});
end

function label = code_to_region(code)
    switch code
        case 1
            label = 'SOZ';
        case 2
            label = 'PZ';
        case 3
            label = 'NIZ';
        otherwise
            label = 'Unknown';
    end
end

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

function eta = compute_eta_from_W(W)
    W = sanitize_W(W);
    if isempty(W)
        eta = 0;
        return;
    end
    n = size(W,1);
    if n < 2
        eta = 0;
        return;
    end
    sumW = sum(W(:));
    if sumW <= 0
        eta = 0;
        return;
    end
    lambda1 = spectral_radius_fast(W);
    if lambda1 <= 0 || ~isfinite(lambda1)
        eta = 0;
    else
        eta = (sumW / n) / lambda1;
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
            val = eigs(sparse(W), 1, 'lm', opts);
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

function M = compute_hyperedge_metrics(W, eta0, doLeaveOneOut, EPS0)
    W = sanitize_W(W);
    n = size(W,1);

    M.out_strength = sum(W, 2);
    M.in_strength  = sum(W, 1)';
    M.total_strength = M.out_strength + M.in_strength;
    M.flow_bias = (M.out_strength - M.in_strength) ./ (M.total_strength + EPS0);

    sumW = sum(W(:));
    if sumW > 0
        M.resource_share = M.out_strength ./ sumW;
    else
        M.resource_share = zeros(n,1);
    end

    M.lambda_share = compute_lambda_share(W, EPS0);

    M.eta_resource_contrib = eta0 .* M.resource_share;
    M.eta_lambda_contrib   = eta0 .* M.lambda_share;

    % 一阶近似：eta 的相对承担 = 资源份额 - 谱模态份额
    % 正值：该超边更偏向增加 eta；负值：更偏向增强主导谱模态、压低 eta。
    M.eta_proxy_delta = eta0 .* (M.resource_share - M.lambda_share);

    M.eta_leave_one_out_delta = nan(n,1);
    if doLeaveOneOut && n >= 3
        for i = 1:n
            W2 = W;
            W2(i,:) = [];
            W2(:,i) = [];
            eta2 = compute_eta_from_W(W2);
            M.eta_leave_one_out_delta(i) = eta0 - eta2;
        end
    end
end

function lambdaShare = compute_lambda_share(W, EPS0)
    n = size(W,1);
    lambdaShare = zeros(n,1);
    if n == 0 || sum(W(:)) <= 0
        return;
    end

    try
        if n >= 20
            opts.disp = 0;
            [V,~] = eigs(sparse(W), 1, 'lm', opts);
            [U,~] = eigs(sparse(W'), 1, 'lm', opts);
            v = abs(real(V(:,1)));
            u = abs(real(U(:,1)));
        else
            [V,D] = eig(W);
            eigVals = diag(D);
            [~, idx] = max(abs(eigVals));
            v = abs(real(V(:,idx)));

            [U,D2] = eig(W');
            eigVals2 = diag(D2);
            [~, idx2] = max(abs(eigVals2));
            u = abs(real(U(:,idx2)));
        end
    catch
        % 若特征向量不稳定，则退化为按出强度分解
        rowSum = sum(W,2);
        if sum(rowSum) > 0
            lambdaShare = rowSum ./ sum(rowSum);
        end
        return;
    end

    if sum(v) <= EPS0 || sum(u) <= EPS0
        rowSum = sum(W,2);
        if sum(rowSum) > 0
            lambdaShare = rowSum ./ sum(rowSum);
        end
        return;
    end

    den = u' * v;
    if abs(den) <= EPS0 || ~isfinite(den)
        rowSum = sum(W,2);
        if sum(rowSum) > 0
            lambdaShare = rowSum ./ sum(rowSum);
        end
        return;
    end

    % lambda sensitivity: d lambda / d W_ij = u_i v_j / (u'v)
    sens = (u * v') ./ den;
    contribMat = W .* sens;
    rowContrib = sum(contribMat, 2);
    rowContrib(rowContrib < 0 | ~isfinite(rowContrib)) = 0;

    if sum(rowContrib) > EPS0
        lambdaShare = rowContrib ./ sum(rowContrib);
    else
        rowSum = sum(W,2);
        if sum(rowSum) > 0
            lambdaShare = rowSum ./ sum(rowSum);
        end
    end
end

function HEinfo = parse_hyperedges_region(HE, node_regions)
    nHE = numel(HE);
    HEinfo = repmat(struct( ...
        'hyperedge_key', '', ...
        'n_nodes', 0, ...
        'n_SOZ', 0, ...
        'n_PZ', 0, ...
        'n_NIZ', 0, ...
        'prop_SOZ', 0, ...
        'prop_PZ', 0, ...
        'prop_NIZ', 0, ...
        'macro_state', 'Unknown'), nHE, 1);

    for i = 1:nHE
        nodes = HE{i};
        nodes = nodes(:)';
        nodes = unique(nodes(isfinite(nodes) & nodes >= 1), 'stable');
        nodes = nodes(nodes <= numel(node_regions));

        HEinfo(i).n_nodes = numel(nodes);
        HEinfo(i).hyperedge_key = nodes_to_key(nodes);

        if isempty(nodes)
            HEinfo(i).macro_state = 'Unknown';
            continue;
        end

        regs = node_regions(nodes);
        nSOZ = sum(regs == 1);
        nPZ  = sum(regs == 2);
        nNIZ = sum(regs == 3);
        nn = numel(regs);

        HEinfo(i).n_SOZ = nSOZ;
        HEinfo(i).n_PZ  = nPZ;
        HEinfo(i).n_NIZ = nNIZ;
        HEinfo(i).prop_SOZ = nSOZ / nn;
        HEinfo(i).prop_PZ  = nPZ  / nn;
        HEinfo(i).prop_NIZ = nNIZ / nn;
        HEinfo(i).macro_state = macro_state_from_counts(nSOZ, nPZ, nNIZ);
    end
end

function key = nodes_to_key(nodes)
    if isempty(nodes)
        key = 'empty';
        return;
    end
    nodes = sort(nodes(:)');
    parts = cell(1,numel(nodes));
    for i = 1:numel(nodes)
        parts{i} = sprintf('%d', nodes(i));
    end
    key = strjoin(parts, '-');
end

function state = macro_state_from_counts(nSOZ, nPZ, nNIZ)
    hasSOZ = nSOZ > 0;
    hasPZ  = nPZ  > 0;
    hasNIZ = nNIZ > 0;

    if hasSOZ && ~hasPZ && ~hasNIZ
        state = 'SOZ_only';
    elseif ~hasSOZ && hasPZ && ~hasNIZ
        state = 'PZ_only';
    elseif ~hasSOZ && ~hasPZ && hasNIZ
        state = 'NIZ_only';
    elseif hasSOZ && hasPZ && ~hasNIZ
        state = 'SOZ_PZ';
    elseif hasSOZ && ~hasPZ && hasNIZ
        state = 'SOZ_NIZ';
    elseif ~hasSOZ && hasPZ && hasNIZ
        state = 'PZ_NIZ';
    elseif hasSOZ && hasPZ && hasNIZ
        state = 'SOZ_PZ_NIZ';
    else
        state = 'Unknown';
    end
end

function roles = classify_role(flowBias, sourceThr, sinkThr)
    n = numel(flowBias);
    roles = cell(n,1);
    for i = 1:n
        if flowBias(i) > sourceThr
            roles{i} = 'source-like';
        elseif flowBias(i) < sinkThr
            roles{i} = 'sink-like';
        else
            roles{i} = 'balanced';
        end
    end
end

function row = make_window_region_row(case_id, window_idx, phase5, region, eta0, eta_saved, region_node_count, region_node_fraction, weights, M, EPS0)
    weights = weights(:);
    nEff = sum(weights);
    denomN = max(region_node_count, 1);

    % -------- raw 指标 --------
    eta_resource_contrib = sum(weights .* M.eta_resource_contrib, 'omitnan');
    eta_lambda_contrib   = sum(weights .* M.eta_lambda_contrib, 'omitnan');
    eta_proxy_delta      = sum(weights .* M.eta_proxy_delta, 'omitnan');
    eta_abs_burden       = sum(weights .* abs(M.eta_proxy_delta), 'omitnan');

    directional_load = sum(weights .* M.total_strength, 'omitnan');
    total_directional = sum(M.total_strength, 'omitnan');
    directional_load_share = directional_load / (total_directional + EPS0);

    sourceLoadVec = max(M.out_strength - M.in_strength, 0);
    sinkLoadVec   = max(M.in_strength - M.out_strength, 0);

    source_load = sum(weights .* sourceLoadVec, 'omitnan');
    sink_load   = sum(weights .* sinkLoadVec, 'omitnan');
    source_load_share = source_load / (sum(sourceLoadVec) + EPS0);
    sink_load_share   = sink_load   / (sum(sinkLoadVec)   + EPS0);

    mean_flow_bias = weighted_nanmean(M.flow_bias, weights);

    % -------- size-adjusted: per-node --------
    n_hyperedges_effective_per_node = nEff / denomN;
    eta_resource_contrib_per_node = eta_resource_contrib / denomN;
    eta_lambda_contrib_per_node   = eta_lambda_contrib / denomN;
    eta_proxy_delta_per_node      = eta_proxy_delta / denomN;
    eta_abs_burden_per_node       = eta_abs_burden / denomN;
    directional_load_per_node     = directional_load / denomN;
    source_load_per_node          = source_load / denomN;
    sink_load_per_node            = sink_load / denomN;

    % -------- size-adjusted: enrichment vs 区域大小占比 --------
    total_abs_burden = sum(abs(M.eta_proxy_delta), 'omitnan');
    eta_abs_burden_share = eta_abs_burden / (total_abs_burden + EPS0);
    eta_abs_burden_enrichment = eta_abs_burden_share / (region_node_fraction + EPS0);
    directional_load_enrichment = directional_load_share / (region_node_fraction + EPS0);
    source_load_enrichment = source_load_share / (region_node_fraction + EPS0);
    sink_load_enrichment   = sink_load_share   / (region_node_fraction + EPS0);

    row = {case_id, window_idx, phase5, region, eta0, eta_saved, ...
           region_node_count, region_node_fraction, ...
           nEff, n_hyperedges_effective_per_node, ...
           eta_resource_contrib, eta_lambda_contrib, eta_proxy_delta, eta_abs_burden, ...
           eta_resource_contrib_per_node, eta_lambda_contrib_per_node, ...
           eta_proxy_delta_per_node, eta_abs_burden_per_node, eta_abs_burden_share, eta_abs_burden_enrichment, ...
           directional_load, directional_load_share, directional_load_per_node, directional_load_enrichment, ...
           source_load, source_load_share, source_load_per_node, source_load_enrichment, ...
           sink_load, sink_load_share, sink_load_per_node, sink_load_enrichment, ...
           mean_flow_bias};
end

function row = make_window_macro_role_row(case_id, window_idx, phase5, macro_state, role, eta0, eta_saved, idx, M, HEinfo, EPS0)
    idx = idx(:);
    n_hyperedges = sum(idx);
    n_nodes_mean = nan;
    if n_hyperedges > 0
        n_nodes_mean = mean([HEinfo(idx).n_nodes], 'omitnan');
    end

    eta_resource_contrib = sum(M.eta_resource_contrib(idx), 'omitnan');
    eta_lambda_contrib   = sum(M.eta_lambda_contrib(idx), 'omitnan');
    eta_proxy_delta      = sum(M.eta_proxy_delta(idx), 'omitnan');
    eta_abs_burden       = sum(abs(M.eta_proxy_delta(idx)), 'omitnan');

    directional_load = sum(M.total_strength(idx), 'omitnan');
    directional_load_share = directional_load / (sum(M.total_strength, 'omitnan') + EPS0);
    mean_flow_bias = mean(M.flow_bias(idx), 'omitnan');

    row = {case_id, window_idx, phase5, macro_state, role, eta0, eta_saved, ...
           n_hyperedges, n_nodes_mean, eta_resource_contrib, eta_lambda_contrib, ...
           eta_proxy_delta, eta_abs_burden, directional_load, directional_load_share, mean_flow_bias};
end

function y = weighted_nanmean(x, w)
    x = x(:);
    w = w(:);
    valid = isfinite(x) & isfinite(w) & w > 0;
    if ~any(valid)
        y = nan;
    else
        y = sum(x(valid) .* w(valid)) / sum(w(valid));
    end
end

function names = hyperedge_varnames()
    names = {'case_id','window_idx','phase5','hyperedge_idx','hyperedge_key', ...
        'n_nodes','n_SOZ','n_PZ','n_NIZ','prop_SOZ','prop_PZ','prop_NIZ', ...
        'macro_state','role','out_strength','in_strength','total_strength','flow_bias', ...
        'resource_share','lambda_share','eta_resource_contrib','eta_lambda_contrib', ...
        'eta_proxy_delta','eta_abs_proxy_delta','eta_leave_one_out_delta', ...
        'eta_recomputed','eta_saved'};
end

function names = window_region_varnames()
    names = {'case_id','window_idx','phase5','region','eta_recomputed','eta_saved', ...
        'region_node_count','region_node_fraction', ...
        'n_hyperedges_effective','n_hyperedges_effective_per_node', ...
        'eta_resource_contrib','eta_lambda_contrib','eta_proxy_delta','eta_abs_burden', ...
        'eta_resource_contrib_per_node','eta_lambda_contrib_per_node', ...
        'eta_proxy_delta_per_node','eta_abs_burden_per_node','eta_abs_burden_share','eta_abs_burden_enrichment', ...
        'directional_load','directional_load_share','directional_load_per_node','directional_load_enrichment', ...
        'source_load','source_load_share','source_load_per_node','source_load_enrichment', ...
        'sink_load','sink_load_share','sink_load_per_node','sink_load_enrichment', ...
        'mean_flow_bias'};
end

function names = window_macro_role_varnames()
    names = {'case_id','window_idx','phase5','macro_state','role','eta_recomputed','eta_saved', ...
        'n_hyperedges','n_nodes_mean','eta_resource_contrib','eta_lambda_contrib', ...
        'eta_proxy_delta','eta_abs_burden','directional_load','directional_load_share','mean_flow_bias'};
end

function names = top_window_varnames()
    names = {'case_id','window_idx','phase5','rank_in_window','hyperedge_idx','hyperedge_key', ...
        'macro_state','role','n_nodes','prop_SOZ','prop_PZ','prop_NIZ', ...
        'eta_proxy_delta','eta_abs_proxy_delta','eta_resource_contrib','eta_lambda_contrib', ...
        'flow_bias','out_strength','in_strength','eta_leave_one_out_delta','eta_recomputed'};
end

function names = top_stage_varnames()
    names = {'case_id','phase5','rank_in_stage','hyperedge_key','macro_state_mode','role_mode', ...
        'n_windows_observed','mean_eta_proxy_delta','median_eta_proxy_delta', ...
        'mean_abs_eta_proxy_delta','mean_eta_resource_contrib','mean_eta_lambda_contrib', ...
        'mean_flow_bias','mean_prop_SOZ','mean_prop_PZ','mean_prop_NIZ'};
end

function Tstage = summarize_top_stage_hyperedges(Hcase, phaseOrder, topKStage)
    rows = {};
    case_id = Hcase.case_id{1};

    for p = 1:numel(phaseOrder)
        ph = phaseOrder{p};
        Hp = Hcase(strcmp(Hcase.phase5, ph), :);
        if isempty(Hp)
            continue;
        end

        keys0 = unique(Hp.hyperedge_key, 'stable');
        tmpRows = {};
        for k = 1:numel(keys0)
            idx = strcmp(Hp.hyperedge_key, keys0{k});
            score = mean(Hp.eta_abs_proxy_delta(idx), 'omitnan');
            if ~isfinite(score)
                score = -inf;
            end
            tmpRows(end+1,:) = {keys0{k}, score}; %#ok<AGROW>
        end

        if isempty(tmpRows)
            continue;
        end

        scoreVec = cell2mat(tmpRows(:,2));
        [~, ord] = sort(scoreVec, 'descend');
        ord = ord(1:min(topKStage, numel(ord)));

        for rr = 1:numel(ord)
            key = tmpRows{ord(rr),1};
            idx = strcmp(Hp.hyperedge_key, key);
            row = {case_id, ph, rr, key, ...
                mode_cell(Hp.macro_state(idx)), mode_cell(Hp.role(idx)), ...
                sum(idx), ...
                mean(Hp.eta_proxy_delta(idx), 'omitnan'), ...
                median(Hp.eta_proxy_delta(idx), 'omitnan'), ...
                mean(Hp.eta_abs_proxy_delta(idx), 'omitnan'), ...
                mean(Hp.eta_resource_contrib(idx), 'omitnan'), ...
                mean(Hp.eta_lambda_contrib(idx), 'omitnan'), ...
                mean(Hp.flow_bias(idx), 'omitnan'), ...
                mean(Hp.prop_SOZ(idx), 'omitnan'), ...
                mean(Hp.prop_PZ(idx), 'omitnan'), ...
                mean(Hp.prop_NIZ(idx), 'omitnan')};
            rows(end+1,:) = row; %#ok<AGROW>
        end
    end

    if isempty(rows)
        Tstage = table();
    else
        Tstage = cell2table(rows, 'VariableNames', top_stage_varnames());
    end
end

function val = mode_cell(c)
    if isempty(c)
        val = '';
        return;
    end
    u = unique(c, 'stable');
    counts = zeros(numel(u),1);
    for i = 1:numel(u)
        counts(i) = sum(strcmp(c, u{i}));
    end
    [~, idx] = max(counts);
    val = u{idx};
end

function Tout = summarize_stage_region(Twin, phaseOrder, regionNames)
    caseIDs = unique(Twin.case_id, 'stable');
    metrics = {'eta_recomputed','eta_saved','n_hyperedges_effective','n_hyperedges_effective_per_node', ...
        'eta_resource_contrib','eta_lambda_contrib','eta_proxy_delta','eta_abs_burden', ...
        'eta_resource_contrib_per_node','eta_lambda_contrib_per_node', ...
        'eta_proxy_delta_per_node','eta_abs_burden_per_node','eta_abs_burden_share','eta_abs_burden_enrichment', ...
        'directional_load','directional_load_share','directional_load_per_node','directional_load_enrichment', ...
        'source_load','source_load_share','source_load_per_node','source_load_enrichment', ...
        'sink_load','sink_load_share','sink_load_per_node','sink_load_enrichment', ...
        'mean_flow_bias'};

    rows = {};
    for i = 1:numel(caseIDs)
        for r = 1:numel(regionNames)
            for p = 1:numel(phaseOrder)
                idx = strcmp(Twin.case_id, caseIDs{i}) & strcmp(Twin.region, regionNames{r}) & strcmp(Twin.phase5, phaseOrder{p});
                xCount = Twin.region_node_count(idx);
                xFrac  = Twin.region_node_fraction(idx);
                if any(isfinite(xCount))
                    region_node_count = xCount(find(isfinite(xCount), 1, 'first'));
                else
                    region_node_count = nan;
                end
                if any(isfinite(xFrac))
                    region_node_fraction = xFrac(find(isfinite(xFrac), 1, 'first'));
                else
                    region_node_fraction = nan;
                end

                row = {caseIDs{i}, regionNames{r}, phaseOrder{p}, sum(idx), region_node_count, region_node_fraction};
                for m = 1:numel(metrics)
                    x = Twin.(metrics{m})(idx);
                    row{end+1} = mean(x, 'omitnan'); %#ok<AGROW>
                    row{end+1} = std(x, 'omitnan'); %#ok<AGROW>
                    row{end+1} = median(x, 'omitnan'); %#ok<AGROW>
                end
                rows(end+1,:) = row; %#ok<AGROW>
            end
        end
    end

    names = {'case_id','region','phase5','n_windows','region_node_count','region_node_fraction'};
    for m = 1:numel(metrics)
        names{end+1} = [metrics{m} '_mean']; %#ok<AGROW>
        names{end+1} = [metrics{m} '_std']; %#ok<AGROW>
        names{end+1} = [metrics{m} '_median']; %#ok<AGROW>
    end
    Tout = cell2table(rows, 'VariableNames', names);
end

function Tout = summarize_stage_macro_role(Twin, phaseOrder, macroOrder, roleOrder)
    caseIDs = unique(Twin.case_id, 'stable');
    metrics = {'eta_recomputed','eta_saved','n_hyperedges','n_nodes_mean', ...
        'eta_resource_contrib','eta_lambda_contrib','eta_proxy_delta','eta_abs_burden', ...
        'directional_load','directional_load_share','mean_flow_bias'};

    rows = {};
    for i = 1:numel(caseIDs)
        for m0 = 1:numel(macroOrder)
            for r0 = 1:numel(roleOrder)
                for p = 1:numel(phaseOrder)
                    idx = strcmp(Twin.case_id, caseIDs{i}) & strcmp(Twin.macro_state, macroOrder{m0}) & ...
                          strcmp(Twin.role, roleOrder{r0}) & strcmp(Twin.phase5, phaseOrder{p});
                    row = {caseIDs{i}, macroOrder{m0}, roleOrder{r0}, phaseOrder{p}, sum(idx)};
                    for mm = 1:numel(metrics)
                        x = Twin.(metrics{mm})(idx);
                        row{end+1} = mean(x, 'omitnan'); %#ok<AGROW>
                        row{end+1} = std(x, 'omitnan'); %#ok<AGROW>
                        row{end+1} = median(x, 'omitnan'); %#ok<AGROW>
                    end
                    rows(end+1,:) = row; %#ok<AGROW>
                end
            end
        end
    end

    names = {'case_id','macro_state','role','phase5','n_windows'};
    for mm = 1:numel(metrics)
        names{end+1} = [metrics{mm} '_mean']; %#ok<AGROW>
        names{end+1} = [metrics{mm} '_std']; %#ok<AGROW>
        names{end+1} = [metrics{mm} '_median']; %#ok<AGROW>
    end
    Tout = cell2table(rows, 'VariableNames', names);
end

function G = summarize_group_region(Tstage, phaseOrder, regionNames)
    metricNames = {'eta_resource_contrib_mean','eta_lambda_contrib_mean','eta_proxy_delta_mean', ...
        'eta_abs_burden_mean', ...
        'eta_resource_contrib_per_node_mean','eta_lambda_contrib_per_node_mean', ...
        'eta_proxy_delta_per_node_mean','eta_abs_burden_per_node_mean', ...
        'eta_abs_burden_share_mean','eta_abs_burden_enrichment_mean', ...
        'directional_load_share_mean','directional_load_per_node_mean','directional_load_enrichment_mean', ...
        'source_load_share_mean','source_load_per_node_mean','source_load_enrichment_mean', ...
        'sink_load_share_mean','sink_load_per_node_mean','sink_load_enrichment_mean', ...
        'mean_flow_bias_mean'};
    rows = {};

    for r = 1:numel(regionNames)
        for m = 1:numel(metricNames)
            metric = metricNames{m};
            preIdx = strcmp(Tstage.region, regionNames{r}) & strcmp(Tstage.phase5, 'pre-ictal');
            preCase = Tstage.case_id(preIdx);
            preVal = Tstage.(metric)(preIdx);

            for p = 1:numel(phaseOrder)
                idx = strcmp(Tstage.region, regionNames{r}) & strcmp(Tstage.phase5, phaseOrder{p});
                x = Tstage.(metric)(idx);
                x = x(isfinite(x));
                n = numel(x);
                mean_x = mean(x, 'omitnan');
                std_x = std(x, 'omitnan');
                sem_x = std_x / sqrt(max(n,1));
                med_x = median(x, 'omitnan');

                % paired delta vs pre
                delta = [];
                p_signrank = nan;
                if ~strcmp(phaseOrder{p}, 'pre-ictal')
                    currRows = Tstage(idx, :);
                    for ii = 1:numel(preCase)
                        j = find(strcmp(currRows.case_id, preCase{ii}), 1, 'first');
                        if ~isempty(j) && isfinite(preVal(ii)) && isfinite(currRows.(metric)(j))
                            delta(end+1,1) = currRows.(metric)(j) - preVal(ii); %#ok<AGROW>
                        end
                    end
                    if numel(delta) >= 2 && exist('signrank', 'file') == 2
                        p_signrank = signrank(delta);
                    end
                end

                rows(end+1,:) = {regionNames{r}, metric, phaseOrder{p}, n, ...
                    mean_x, std_x, sem_x, med_x, mean(delta, 'omitnan'), median(delta, 'omitnan'), p_signrank}; %#ok<AGROW>
            end
        end
    end

    G = cell2table(rows, 'VariableNames', {'region','metric','phase5','n_case', ...
        'mean','std','sem','median','mean_delta_vs_pre','median_delta_vs_pre','p_signrank_delta_vs_pre'});
end

function G = summarize_group_macro_role(Tstage, phaseOrder, macroOrder, roleOrder)
    metricNames = {'eta_resource_contrib_mean','eta_lambda_contrib_mean','eta_proxy_delta_mean', ...
        'eta_abs_burden_mean','directional_load_share_mean','mean_flow_bias_mean'};
    rows = {};

    for m0 = 1:numel(macroOrder)
        for r0 = 1:numel(roleOrder)
            for mm = 1:numel(metricNames)
                metric = metricNames{mm};
                preIdx = strcmp(Tstage.macro_state, macroOrder{m0}) & strcmp(Tstage.role, roleOrder{r0}) & strcmp(Tstage.phase5, 'pre-ictal');
                preCase = Tstage.case_id(preIdx);
                preVal = Tstage.(metric)(preIdx);

                for p = 1:numel(phaseOrder)
                    idx = strcmp(Tstage.macro_state, macroOrder{m0}) & strcmp(Tstage.role, roleOrder{r0}) & strcmp(Tstage.phase5, phaseOrder{p});
                    x = Tstage.(metric)(idx);
                    x = x(isfinite(x));
                    n = numel(x);
                    mean_x = mean(x, 'omitnan');
                    std_x = std(x, 'omitnan');
                    sem_x = std_x / sqrt(max(n,1));
                    med_x = median(x, 'omitnan');

                    delta = [];
                    p_signrank = nan;
                    if ~strcmp(phaseOrder{p}, 'pre-ictal')
                        currRows = Tstage(idx, :);
                        for ii = 1:numel(preCase)
                            j = find(strcmp(currRows.case_id, preCase{ii}), 1, 'first');
                            if ~isempty(j) && isfinite(preVal(ii)) && isfinite(currRows.(metric)(j))
                                delta(end+1,1) = currRows.(metric)(j) - preVal(ii); %#ok<AGROW>
                            end
                        end
                        if numel(delta) >= 2 && exist('signrank', 'file') == 2
                            p_signrank = signrank(delta);
                        end
                    end

                    rows(end+1,:) = {macroOrder{m0}, roleOrder{r0}, metric, phaseOrder{p}, n, ...
                        mean_x, std_x, sem_x, med_x, mean(delta, 'omitnan'), median(delta, 'omitnan'), p_signrank}; %#ok<AGROW>
                end
            end
        end
    end

    G = cell2table(rows, 'VariableNames', {'macro_state','role','metric','phase5','n_case', ...
        'mean','std','sem','median','mean_delta_vs_pre','median_delta_vs_pre','p_signrank_delta_vs_pre'});
end
