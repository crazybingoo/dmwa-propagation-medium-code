%% eta 机制仿真：二维相图/切片图
% 思路：
% c = 核心集中度（active support 收缩到局部核心的程度）
% h = 结构均匀化程度（active support 内节点参与概率从 hub-dominant 向 uniform 过渡）
%
% 生成逻辑：
% 1) 先确定 active support 的大小：c 越大，support 越小，重叠越强
% 2) 在 active support 内，根据 h 混合“偏斜权重”和“均匀权重”
% 3) 生成 2-节点与 3-节点超边集合 HE
% 4) 按原方法构造 O, Coverage, DegreeBias 与 W
% 5) 计算
%       eta = ((1/N) * sum(W(:))) / lambda_1
%       avgW = (1/N) * sum(W(:))
%       lambda1 = spectral radius of W
%
% 说明：
% - 脚本会优先读取缓存 eta_phase_mechanism_data.mat
% - 若缓存不存在，则自动进行蒙特卡洛仿真并保存

clc; clear; close all;

%% 0) 路径与开关
outDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig9_eta_mechanism';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cacheFile = fullfile(outDir, 'eta_phase_mechanism_data.mat');
forceRecompute = false;

%% 1) 读取或计算
R = load_or_run_simulation(cacheFile, forceRecompute);

c_range = R.c_range;
h_range = R.h_range;

%% 2) 图 2：Sigma W / N 相图
fig = figure('Color', 'w', 'Position', [120, 120, 700, 560]);
imagesc(c_range, h_range, R.avgW_mean);
axis xy;
colormap(parula);
cb = colorbar;
cb.Label.String = '\Sigma W / N';

xlabel('Core concentration c', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Homogenization h', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
title('\Sigma W / N phase diagram', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.1, ...
    'Box', 'off', 'TickDir', 'out');

exportgraphics(fig, fullfile(outDir, 'Fig2_avgW_phase_diagram.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'Fig2_avgW_phase_diagram.pdf'), 'ContentType', 'vector');
savefig(fig, fullfile(outDir, 'Fig2_avgW_phase_diagram.fig'));


%% ========================= 局部函数 =========================
function R = load_or_run_simulation(cacheFile, forceRecompute)
    if exist(cacheFile, 'file') && ~forceRecompute
        S = load(cacheFile);
        R = S.R;
        fprintf('Loaded cache: %s\n', cacheFile);
        return;
    end

    fprintf('Cache not found. Running Monte Carlo simulation...\n');

    params.N_nodes        = 100;
    params.num_E2         = 120;    % 二元超边数
    params.num_E3         = 80;     % 三元超边数
    params.core_min_frac  = 0.25;   % c=1 时 active support 占全网比例

%     params.c_range        = 0:0.05:1;
%     params.h_range        = 0:0.05:1;
%     params.num_trials     = 40;     % 可增大以换取更平滑图像

    params.c_range    = linspace(0, 1, 41);
    params.h_range    = linspace(0, 1, 41);
    params.num_trials = 30;

    params.skewAlpha      = 4.0;    % h=0 时 active support 内权重偏斜强度

    c_range = params.c_range;
    h_range = params.h_range;

    nC = numel(c_range);
    nH = numel(h_range);

    eta_mean      = nan(nH, nC);
    eta_std       = nan(nH, nC);
    avgW_mean     = nan(nH, nC);
    avgW_std      = nan(nH, nC);
    lambda1_mean  = nan(nH, nC);
    lambda1_std   = nan(nH, nC);
    nHE_mean      = nan(nH, nC);
    degreecv_mean = nan(nH, nC);

    tic_all = tic;
    for ih = 1:nH
        for ic = 1:nC
            c = c_range(ic);
            h = h_range(ih);

            eta_vals     = nan(params.num_trials, 1);
            avgW_vals    = nan(params.num_trials, 1);
            lambda1_vals = nan(params.num_trials, 1);
            nHE_vals     = nan(params.num_trials, 1);
            dcv_vals     = nan(params.num_trials, 1);

            for t = 1:params.num_trials
                stat = simulate_one_realization(params, c, h);
                eta_vals(t)     = stat.eta;
                avgW_vals(t)    = stat.avgW;
                lambda1_vals(t) = stat.lambda1;
                nHE_vals(t)     = stat.nHE;
                dcv_vals(t)     = stat.degreeCV;
            end

            eta_mean(ih, ic)      = mean(eta_vals, 'omitnan');
            eta_std(ih, ic)       = std(eta_vals,  'omitnan');
            avgW_mean(ih, ic)     = mean(avgW_vals, 'omitnan');
            avgW_std(ih, ic)      = std(avgW_vals,  'omitnan');
            lambda1_mean(ih, ic)  = mean(lambda1_vals, 'omitnan');
            lambda1_std(ih, ic)   = std(lambda1_vals,  'omitnan');
            nHE_mean(ih, ic)      = mean(nHE_vals, 'omitnan');
            degreecv_mean(ih, ic) = mean(dcv_vals, 'omitnan');
        end

        fprintf('Progress %2d/%2d | elapsed %.1fs\n', ih, nH, toc(tic_all));
    end

    R = struct();
    R.params        = params;
    R.c_range       = c_range;
    R.h_range       = h_range;
    R.eta_mean      = eta_mean;
    R.eta_std       = eta_std;
    R.avgW_mean     = avgW_mean;
    R.avgW_std      = avgW_std;
    R.lambda1_mean  = lambda1_mean;
    R.lambda1_std   = lambda1_std;
    R.nHE_mean      = nHE_mean;
    R.degreecv_mean = degreecv_mean;

    save(cacheFile, 'R');
    fprintf('Saved cache: %s\n', cacheFile);
end

function stat = simulate_one_realization(params, c, h)
    N = params.N_nodes;
    num_E2 = params.num_E2;
    num_E3 = params.num_E3;

    min_support = max(4, round(params.core_min_frac * N));
    active_size = round(min_support + (1 - c) * (N - min_support));
    active_size = max(4, min(N, active_size));

    active_nodes = randperm(N, active_size);

    % active support 内的“偏斜”与“均匀”混合
    skew = exp(-params.skewAlpha * linspace(0, 1, active_size));
    skew = skew(randperm(active_size));  % 打乱，以免固定节点排序产生偏差
    skew = skew / sum(skew);

    uniform_w = ones(1, active_size) / active_size;
    w_active = (1 - h) * skew + h * uniform_w;
    w_active = w_active / sum(w_active);

    p = zeros(1, N);
    p(active_nodes) = w_active;

    edge_sizes = [2 * ones(1, num_E2), 3 * ones(1, num_E3)];
    edge_sizes = edge_sizes(randperm(numel(edge_sizes)));

    HE = cell(numel(edge_sizes), 1);
    for k = 1:numel(edge_sizes)
        HE{k} = weighted_sample_without_replacement(p, edge_sizes(k));
    end

    HE = unique_hyperedges_local(HE);
    W = build_W_from_hyperedges_local(HE);

    stat.nHE = size(W, 1);

    if isempty(W) || size(W,1) < 2
        stat.avgW = 0;
        stat.lambda1 = 0;
        stat.eta = 0;
        stat.degreeCV = 0;
        return;
    end

    W = double(W);
    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    n = size(W,1);
    W(1:n+1:end) = 0;

    deg = sum(W > 0, 2);
    if mean(deg) <= eps
        degreeCV = 0;
    else
        degreeCV = std(deg) / mean(deg);
    end

    avgW = sum(W(:)) / n;
    lambda1 = max(abs(eig(W)));

    if ~isfinite(lambda1) || lambda1 <= eps
        eta = 0;
    else
        eta = avgW / lambda1;
    end

    stat.avgW = avgW;
    stat.lambda1 = lambda1;
    stat.eta = eta;
    stat.degreeCV = degreeCV;
end

function nodes = weighted_sample_without_replacement(p, k)
    p = p(:)';
    nodes = zeros(1, k);
    avail = 1:numel(p);
    weights = p;

    for ii = 1:k
        if sum(weights) <= 0
            error('No positive probability mass left for sampling.');
        end
        weights = weights / sum(weights);
        u = rand();
        cs = cumsum(weights);
        idx = find(u <= cs, 1, 'first');
        nodes(ii) = avail(idx);

        avail(idx) = [];
        weights(idx) = [];
    end

    nodes = sort(nodes);
end

function HE_out = unique_hyperedges_local(HE_in)
    if isempty(HE_in)
        HE_out = {};
        return;
    end

    tmp = cell(size(HE_in));
    for i = 1:numel(HE_in)
        tmp{i} = sort(unique(HE_in{i}));
    end

    keys = cellfun(@(x) sprintf('%d_', x), tmp, 'UniformOutput', false);
    [~, ia] = unique(keys, 'stable');
    HE_out = tmp(ia);
end

function W = build_W_from_hyperedges_local(HE)
    num_HE = numel(HE);

    if num_HE < 2
        W = zeros(num_HE);
        return;
    end

    sizes = cellfun(@numel, HE);
    sizes = sizes(:);

    O = zeros(num_HE, num_HE);
    for i = 1:num_HE
        for j = (i + 1):num_HE
            overlap_len = numel(intersect(HE{i}, HE{j}));
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
end
