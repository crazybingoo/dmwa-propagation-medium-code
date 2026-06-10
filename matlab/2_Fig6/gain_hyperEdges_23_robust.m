function all_hyperEdges = gain_hyperEdges_23_robust(datanew, plvQuantile)

if nargin < 2 || isempty(plvQuantile)
    plvQuantile = 0.55;
end

% =========================================================
% 说明：
% 1) 这里给的是一个可运行的“占位框架”
% 2) 你最稳的做法是：把你原 gain_hyperEdges_23.m 的主体完整粘到这里
% 3) 然后仅把“固定 0.55 阈值”那一行替换成 plvQuantile 版本
% =========================================================

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

% -------------------- Step 2: 分位数阈值 --------------------
PLV_U = triu(PLV, 1);
plv_vals = PLV_U(PLV_U > 0);

if isempty(plv_vals)
    all_hyperEdges = {};
    return;
end

plv_vals = sort(plv_vals, 'ascend');
n = numel(plv_vals);
idx = round(n * plvQuantile);
idx = max(1, min(n, idx));
thr = plv_vals(idx);

A = double(PLV >= thr);
A(1:nCh+1:end) = 0;

% -------------------- Step 3: 从二值图提“超边” --------------------
% 这里提供一个通用占位版本：
% 用极大团近似高同步子集，作为超边候选
% 如果你原 gain_hyperEdges_23.m 有自己的超边提取逻辑，
% 强烈建议直接用你原来的逻辑替换本部分。

G = graph(A);
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