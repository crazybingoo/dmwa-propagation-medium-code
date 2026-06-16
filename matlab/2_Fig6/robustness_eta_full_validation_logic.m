%% 全量 24 次发作鲁棒性分析：eta 的区分度与显著性验证
% 逻辑对齐 robustness_full_validation_2.m：
% 1) 病例 × 窗长 预计算目标窗口的 PLV
% 2) 阈值网格遍历时直接把 PLV 转成超边，再计算 W 与 eta
% 3) 只比较 pre vs ictal
% 4) 单病例实时保存中间结果
%
% eta 定义：
%   eta = ((1/N) * sum_{i,j} W_{ij}) / lambda_1
% 其中 N 为超边数，lambda_1 为 W 的主特征值（谱半径）

clc; clear; close all;

%% 1. 参数设置与路径
outDir = 'example_project\2_Fig6\3\Robustness_Results_eta';
if ~exist(outDir, 'dir'), mkdir(outDir); end

winLenRange = [1, 2, 3, 4, 5];
retainedFractionRange = 0.45 : 0.05 : 0.90;
nWin   = numel(winLenRange);
nThres = numel(retainedFractionRange);

P_val_matrix = ones(nWin, nThres);     % 配对 signrank 的 p 值
Diff_matrix  = zeros(nWin, nThres);    % Mean eta(Ictal) - Mean eta(Pre-ictal)

%% 2. 病例信息定义（24 个病例）
caseInfo = struct([]);

caseInfo(1).case_id   = 'seizure_01';
caseInfo(1).file_path = fullfile('data', 'seizure_01_Gamma.mat');
caseInfo(1).pre_idx   = 1:200;
caseInfo(1).ictal_idx = 201:400;
caseInfo(1).post_idx  = 401:598;

caseInfo(2).case_id   = 'seizure_02';
caseInfo(2).file_path = fullfile('data', 'seizure_02_Gamma.mat');
caseInfo(2).pre_idx   = 1:122;
caseInfo(2).ictal_idx = 123:233;
caseInfo(2).post_idx  = 234:338;

caseInfo(3).case_id   = 'seizure_03';
caseInfo(3).file_path = fullfile('data', 'seizure_03_Gamma.mat');
caseInfo(3).pre_idx   = 1:195;
caseInfo(3).ictal_idx = 196:330;
caseInfo(3).post_idx  = 331:477;

caseInfo(4).case_id   = 'seizure_04';
caseInfo(4).file_path = fullfile('data', 'seizure_04_Gamma.mat');
caseInfo(4).pre_idx   = 1:137;
caseInfo(4).ictal_idx = 138:297;
caseInfo(4).post_idx  = 298:398;

caseInfo(5).case_id   = 'seizure_05';
caseInfo(5).file_path = fullfile('data', 'seizure_05_Gamma.mat');
caseInfo(5).pre_idx   = 1:252;
caseInfo(5).ictal_idx = 253:413;
caseInfo(5).post_idx  = 414:498;

caseInfo(6).case_id   = 'seizure_06';
caseInfo(6).file_path = fullfile('data', 'seizure_06_Gamma.mat');
caseInfo(6).pre_idx   = 1:199;
caseInfo(6).ictal_idx = 200:285;
caseInfo(6).post_idx  = 286:398;

caseInfo(7).case_id   = 'seizure_07';
caseInfo(7).file_path = fullfile('data', 'seizure_07_Gamma.mat');
caseInfo(7).pre_idx   = 1:239;
caseInfo(7).ictal_idx = 240:331;
caseInfo(7).post_idx  = 332:498;

caseInfo(8).case_id   = 'seizure_08';
caseInfo(8).file_path = fullfile('data', 'seizure_08_Gamma.mat');
caseInfo(8).pre_idx   = 1:236;
caseInfo(8).ictal_idx = 237:313;
caseInfo(8).post_idx  = 314:498;

caseInfo(9).case_id   = 'seizure_09';
caseInfo(9).file_path = fullfile('data', 'seizure_09_Gamma.mat');
caseInfo(9).pre_idx   = 1:59;
caseInfo(9).ictal_idx = 60:112;
caseInfo(9).post_idx  = 113:385;

caseInfo(10).case_id   = 'seizure_10';
caseInfo(10).file_path = fullfile('data', 'seizure_10_Gamma.mat');
caseInfo(10).pre_idx   = 1:244;
caseInfo(10).ictal_idx = 245:470;
caseInfo(10).post_idx  = 471:598;

caseInfo(11).case_id   = 'seizure_11';
caseInfo(11).file_path = fullfile('data', 'seizure_11_Gamma.mat');
caseInfo(11).pre_idx   = 1:150;
caseInfo(11).ictal_idx = 151:271;
caseInfo(11).post_idx  = 272:404;

caseInfo(12).case_id   = 'seizure_12';
caseInfo(12).file_path = fullfile('data', 'seizure_12_Gamma.mat');
caseInfo(12).pre_idx   = 1:180;
caseInfo(12).ictal_idx = 181:261;
caseInfo(12).post_idx  = 262:498;

caseInfo(13).case_id   = 'seizure_13';
caseInfo(13).file_path = fullfile('data', 'seizure_13_Gamma.mat');
caseInfo(13).pre_idx   = 1:69;
caseInfo(13).ictal_idx = 70:130;
caseInfo(13).post_idx  = 131:198;

caseInfo(14).case_id   = 'seizure_14';
caseInfo(14).file_path = fullfile('data', 'seizure_14_Gamma.mat');
caseInfo(14).pre_idx   = 1:252;
caseInfo(14).ictal_idx = 253:346;
caseInfo(14).post_idx  = 347:448;

caseInfo(15).case_id   = 'seizure_15';
caseInfo(15).file_path = fullfile('data', 'seizure_15_Gamma.mat');
caseInfo(15).pre_idx   = 1:154;
caseInfo(15).ictal_idx = 155:230;
caseInfo(15).post_idx  = 231:324;

caseInfo(16).case_id   = 'seizure_16';
caseInfo(16).file_path = fullfile('data', 'seizure_16_Gamma.mat');
caseInfo(16).pre_idx   = 1:253;
caseInfo(16).ictal_idx = 254:331;
caseInfo(16).post_idx  = 332:489;

caseInfo(17).case_id   = 'seizure_17';
caseInfo(17).file_path = fullfile('data', 'seizure_17_Gamma.mat');
caseInfo(17).pre_idx   = 1:135;
caseInfo(17).ictal_idx = 136:270;
caseInfo(17).post_idx  = 271:309;

caseInfo(18).case_id   = 'seizure_18';
caseInfo(18).file_path = fullfile('data', 'seizure_18_Gamma.mat');
caseInfo(18).pre_idx   = 1:117;
caseInfo(18).ictal_idx = 118:300;
caseInfo(18).post_idx  = 301:598;

caseInfo(19).case_id   = 'seizure_19';
caseInfo(19).file_path = fullfile('data', 'seizure_19_Gamma.mat');
caseInfo(19).pre_idx   = 1:255;
caseInfo(19).ictal_idx = 256:482;
caseInfo(19).post_idx  = 483:598;

caseInfo(20).case_id   = 'seizure_20';
caseInfo(20).file_path = fullfile('data', 'seizure_20_Gamma.mat');
caseInfo(20).pre_idx   = 1:311;
caseInfo(20).ictal_idx = 312:540;
caseInfo(20).post_idx  = 541:598;

caseInfo(21).case_id   = 'seizure_21';
caseInfo(21).file_path = fullfile('data', 'seizure_21_Gamma.mat');
caseInfo(21).pre_idx   = 1:94;
caseInfo(21).ictal_idx = 95:170;
caseInfo(21).post_idx  = 171:318;

caseInfo(22).case_id   = 'seizure_22';
caseInfo(22).file_path = fullfile('data', 'seizure_22_Gamma.mat');
caseInfo(22).pre_idx   = 1:95;
caseInfo(22).ictal_idx = 96:195;
caseInfo(22).post_idx  = 196:350;

caseInfo(23).case_id   = 'seizure_23';
caseInfo(23).file_path = fullfile('data', 'seizure_23_Gamma.mat');
caseInfo(23).pre_idx   = 1:162;
caseInfo(23).ictal_idx = 163:228;
caseInfo(23).post_idx  = 229:310;

caseInfo(24).case_id   = 'seizure_24';
caseInfo(24).file_path = fullfile('data', 'seizure_24_Gamma.mat');
caseInfo(24).pre_idx   = 1:29;
caseInfo(24).ictal_idx = 30:95;
caseInfo(24).post_idx  = 96:148;

nCase = numel(caseInfo);
fs = 1024;

%% 3. 用于存储所有病例、所有参数组合的结果
all_eta_pre_storage   = nan(nCase, nWin, nThres);
all_eta_ictal_storage = nan(nCase, nWin, nThres);

fprintf('开始全量 %d 次发作的 eta 鲁棒性分析...\n', nCase);
totalTic = tic;

%% 4. 主循环
for c = 1:nCase
    case_id = caseInfo(c).case_id;
    fprintf('\n--------------------------------------------------\n');
    fprintf('(%d/%d) 开始处理病例: %s\n', c, nCase, case_id);
    caseTic = tic;

    if ~exist(caseInfo(c).file_path, 'file')
        warning('文件不存在: %s', caseInfo(c).file_path);
        continue;
    end

    %% 4.1 加载数据
    stepTic = tic;
    S0 = load(caseInfo(c).file_path);
    if ~isfield(S0, 'X1')
        warning('%s 中没有变量 X1，跳过。', case_id);
        continue;
    end
    X1 = S0.X1;
    fprintf('  [1/4] 加载数据完成，耗时 %.2f 秒\n', toc(stepTic));

    %% 4.2 窗长循环
    for i = 1:nWin
        currWin = winLenRange(i);
        stepTic = tic;

        T_max = floor(size(X1, 2) / fs) - currWin + 1;
        if T_max < 1
            warning('%s 在窗长 %d s 下数据不足，跳过。', case_id, currWin);
            continue;
        end

        valid_pre   = sanitize_idx(caseInfo(c).pre_idx, T_max);
        valid_ictal = sanitize_idx(caseInfo(c).ictal_idx, T_max);

        % 优先级：pre > ictal
        valid_ictal = setdiff(valid_ictal, valid_pre, 'stable');

        target_indices = [valid_pre(:)', valid_ictal(:)'];
        if isempty(target_indices)
            warning('%s 在窗长 %d s 下 pre 与 ictal 都为空，跳过。', case_id, currWin);
            continue;
        end

        %% 4.2.1 预计算目标窗口的 PLV
        plv_cache = cell(numel(target_indices), 1);
        for k = 1:numel(target_indices)
            t = target_indices(k);
            dat = X1(:, fs*(t-1)+1 : fs*(t-1+currWin));
            plv_cache{k} = fast_plv_matrix(dat);
        end
        fprintf('  [2/4] 窗长 %ds: PLV 预计算完成（%d 个窗口），耗时 %.2f 秒\n', ...
            currWin, numel(target_indices), toc(stepTic));

        %% 4.2.2 阈值网格遍历：PLV -> hyperedges -> W -> eta
        stepTic = tic;
        nPre = numel(valid_pre);

        for j = 1:nThres
            currRetainedFraction = retainedFractionRange(j);
            eta_vals = nan(numel(target_indices), 1);

            for k = 1:numel(target_indices)
                HE = plv_to_hyperedges_retained_fraction(plv_cache{k}, currRetainedFraction);
                W  = calculate_W_internal(HE);
                eta_vals(k) = calculate_eta_from_W_internal(W);
            end

            if nPre > 0
                all_eta_pre_storage(c, i, j) = mean(eta_vals(1:nPre), 'omitnan');
            else
                all_eta_pre_storage(c, i, j) = NaN;
            end

            if nPre < numel(target_indices)
                all_eta_ictal_storage(c, i, j) = mean(eta_vals(nPre+1:end), 'omitnan');
            else
                all_eta_ictal_storage(c, i, j) = NaN;
            end
        end
        fprintf('  [3/4] 窗长 %ds: 阈值网格计算完成，耗时 %.2f 秒\n', currWin, toc(stepTic));
    end

    %% 4.3 单病例实时保存
    stepTic = tic;
    save_name = fullfile(outDir, [case_id '_eta_robustness_temp.mat']);
    case_pre_eta   = squeeze(all_eta_pre_storage(c, :, :));
    case_ictal_eta = squeeze(all_eta_ictal_storage(c, :, :));

    save(save_name, 'case_id', 'case_pre_eta', 'case_ictal_eta', ...
        'winLenRange', 'retainedFractionRange');

    fprintf('  [4/4] 病例结果实时保存完成: %s，耗时 %.2f 秒\n', ...
        [case_id '_eta_robustness_temp.mat'], toc(stepTic));
    fprintf('病例 %s 总处理用时: %.2f 秒\n', case_id, toc(caseTic));
end

%% 5. 汇总统计与检验
fprintf('\n==================================================\n');
fprintf('所有病例 eta 计算完成，正在进行汇总统计检验...\n');

for i = 1:nWin
    for j = 1:nThres
        pre_group   = all_eta_pre_storage(:, i, j);
        ictal_group = all_eta_ictal_storage(:, i, j);

        v_idx = ~isnan(pre_group) & ~isnan(ictal_group);

        if sum(v_idx) >= 5
            P_val_matrix(i,j) = signrank(pre_group(v_idx), ictal_group(v_idx));
            Diff_matrix(i,j)  = mean(ictal_group(v_idx)) - mean(pre_group(v_idx));
        else
            P_val_matrix(i,j) = 1;
            Diff_matrix(i,j)  = 0;
        end
    end
end
%% 6. 绘图与最终保存
xLabels = string(retainedFractionRange);
yLabels = string(winLenRange) + " s";

% 图 1：显著性热图
fig1 = figure('Color','w','Position',[100 100 650 500]);
h1 = heatmap(xLabels, yLabels, -log10(P_val_matrix));
h1.Title = 'Discriminative Power of \eta: -log10(p-value)';
h1.XLabel = 'Retained PLV edge fraction';
h1.YLabel = 'Window Length';
h1.Colormap = parula;

exportgraphics(fig1, fullfile(outDir, 'Robustness_eta_P_Value_Heatmap.png'), 'Resolution', 600);
savefig(fig1, fullfile(outDir, 'Robustness_eta_P_Value_Heatmap.fig'));

% 图 2：效应量热图
fig2 = figure('Color','w','Position',[150 150 650 500]);
h2 = heatmap(xLabels, yLabels, Diff_matrix);
h2.Title = 'Effect Size of \eta: Mean(Ictal) - Mean(Pre-ictal)';
h2.XLabel = 'Retained PLV edge fraction';
h2.YLabel = 'Window Length';
h2.Colormap = parula;

exportgraphics(fig2, fullfile(outDir, 'Robustness_eta_Diff_Heatmap.png'), 'Resolution', 600);
savefig(fig2, fullfile(outDir, 'Robustness_eta_Diff_Heatmap.fig'));

save(fullfile(outDir, 'Robustness_eta_Data_Workspace.mat'), ...
    'P_val_matrix', 'Diff_matrix', ...
    'all_eta_pre_storage', 'all_eta_ictal_storage', ...
    'winLenRange', 'retainedFractionRange', 'caseInfo');

fprintf('全量 eta 分析总耗时: %.2f 分钟\n', toc(totalTic)/60);
fprintf('结果已保存到:\n%s\n', outDir);

%% ===================== 内部函数 =====================

function idx = sanitize_idx(idx, nMax)
    if isempty(idx)
        idx = [];
        return;
    end
    idx = idx(:).';
    idx = idx(isfinite(idx));
    idx = round(idx);
    idx = idx(idx >= 1 & idx <= nMax);
    idx = unique(idx, 'stable');
end

function plvMat = fast_plv_matrix(dat)
    [nCh, ~] = size(dat);
    plvMat = zeros(nCh, nCh);
    phases = zeros(size(dat));
    for i = 1:nCh
        phases(i,:) = angle(hilbert(dat(i,:)));
    end
    for i = 1:nCh
        for j = i:nCh
            val = abs(mean(exp(1i * (phases(i,:) - phases(j,:)))));
            plvMat(i,j) = val;
            plvMat(j,i) = val;
        end
    end
end

function HE = plv_to_hyperedges_retained_fraction(plvMat, retainedFraction)
    dc = size(plvMat, 1);

    plv_core = plvMat - eye(dc);
    vals = plv_core(triu(true(dc), 1));
    vals = vals(isfinite(vals) & vals > 0);

    if isempty(vals)
        HE = {};
        return;
    end

    vals = sort(vals, 'descend');
    idx = round(numel(vals) * retainedFraction);
    idx = max(1, min(numel(vals), idx));
    threshold = vals(idx);

    conn = double(plvMat >= threshold);
    conn(1:dc+1:end) = 0;

    [I, J] = find(triu(conn, 1));

    HE2 = arrayfun(@(ii) [I(ii), J(ii)], 1:numel(I), 'UniformOutput', false).';

    d = [];
    for m = 1:numel(I)
        u = I(m);
        v = J(m);
        common = intersect(find(conn(u, :)), find(conn(v, :)));
        common = common(common > v);
        if ~isempty(common)
            for k = 1:numel(common)
                d = [d; u, v, common(k)]; %#ok<AGROW>
            end
        end
    end

    if isempty(d)
        HE3 = {};
    else
        HE3 = mat2cell(d, ones(size(d,1),1), 3);
    end

    HE = [HE2; HE3];
    HE = unique_hyperedges_internal(HE);
end

function HE_out = unique_hyperedges_internal(HE_in)
    if isempty(HE_in)
        HE_out = {};
        return;
    end

    tmp = cell(size(HE_in));
    for ii = 1:numel(HE_in)
        tmp{ii} = sort(unique(HE_in{ii}));
    end

    keys = cellfun(@(x) sprintf('%d_', x), tmp, 'UniformOutput', false);
    [~, ia] = unique(keys, 'stable');
    HE_out = tmp(ia);
end

function W = calculate_W_internal(HE)
    num_HE = length(HE);
    if num_HE < 2
        W = zeros(num_HE);
        return;
    end

    sizes = cellfun(@length, HE);
    sizes = sizes(:);

    O = zeros(num_HE);
    for i = 1:num_HE
        for j = i+1:num_HE
            ov = length(intersect(HE{i}, HE{j}));
            if ov > 0
                O(i,j) = ov;
                O(j,i) = ov;
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
end

function eta = calculate_eta_from_W_internal(W)
    if isempty(W)
        eta = 0;
        return;
    end

    W = double(W);
    n = size(W, 1);

    if size(W,2) ~= n
        error('W must be square.');
    end

    if n < 2
        eta = 0;
        return;
    end

    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    W(1:n+1:end) = 0;

    sum_W = sum(W(:));
    if sum_W <= 0
        eta = 0;
        return;
    end

    lambda_1 = max(abs(eig(W)));
    if ~isfinite(lambda_1) || lambda_1 <= eps
        eta = 0;
    else
        eta = (sum_W / n) / lambda_1;
    end
end
