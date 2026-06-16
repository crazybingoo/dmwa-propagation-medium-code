function all_hyperEdges = gain_hyperEdges_23_robust(datanew, retainedFraction)
% Extract hyperedges from PLV. Manuscript main scripts select and pass one
% seizure-level elbow threshold. If a numeric value is supplied here, it is
% used as the retained PLV edge fraction required by the Fig. 6 robustness grid.

% -------------------- Step 1: 计算 PLV --------------------
nCh = size(datanew, 1);

sig = datanew;
analytic_sig = hilbert(sig.').';           % 每行一个通道
phase_sig = angle(analytic_sig);

PLV = zeros(nCh, nCh);
for i = 1:nCh
    for j = i+1:nCh
        dphi = phase_sig(i,:) - phase_sig(j,:);
        plv_ij = abs(mean(exp(1i * dphi)));
        PLV(i,j) = plv_ij;
        PLV(j,i) = plv_ij;
    end
end

% -------------------- Step 2: threshold --------------------
PLV_U = triu(PLV, 1);
plv_vals = PLV_U(PLV_U > 0);

if isempty(plv_vals)
    all_hyperEdges = {};
    return;
end

if nargin < 2 || isempty(retainedFraction)
    thr = select_plv_elbow_threshold_local(plv_vals);
else
    plv_vals = sort(plv_vals, 'descend');
    n = numel(plv_vals);
    idx = round(n * retainedFraction);
    idx = max(1, min(n, idx));
    thr = plv_vals(idx);
end

A = double(PLV >= thr);
A(1:nCh+1:end) = 0;

% -------------------- Step 3: 从二值图提“超边” --------------------
% 这里提供一个通用占位版本：
% 用极大团近似高同步子集，作为超边候选
% 如果你原 gain_hyperEdges_23.m 有自己的超边提取逻辑，
% 强烈建议直接用你原来的逻辑替换本部分。

cliques = maximalCliques_bruteforce(A);

% 只保留大小 >= 3 的 clique 作为高阶超边
all_hyperEdges = {};
for k = 1:numel(cliques)
    he = cliques{k};
    if numel(he) >= 3
        all_hyperEdges{end+1,1} = he(:).'; %#ok<AGROW>
    end
end

% 如果没有 >=3 的超边，就退化为边对
if isempty(all_hyperEdges)
    [r,c] = find(triu(A,1) > 0);
    for k = 1:numel(r)
        all_hyperEdges{end+1,1} = [r(k), c(k)]; %#ok<AGROW>
    end
end

% 去重
all_hyperEdges = unique_hyperedges(all_hyperEdges);

end


function threshold = select_plv_elbow_threshold_local(plvVals)
plvVals = plvVals(isfinite(plvVals) & plvVals > 0);
if isempty(plvVals)
    threshold = inf;
    return;
end

candidates = 0:0.0025:1;
density = arrayfun(@(x) mean(plvVals >= x), candidates);

if numel(candidates) < 3 || (max(candidates) - min(candidates)) == 0 || ...
        (max(density) - min(density)) == 0
    threshold = median(plvVals);
    return;
end

xn = (candidates - min(candidates)) ./ (max(candidates) - min(candidates));
yn = (density - min(density)) ./ (max(density) - min(density));
p1 = [xn(1), yn(1)];
p2 = [xn(end), yn(end)];
lineVec = p2 - p1;

dist = abs(lineVec(1) .* (p1(2) - yn) - (p1(1) - xn) .* lineVec(2)) ./ norm(lineVec);
[~, idx] = max(dist);
threshold = candidates(idx);
end
