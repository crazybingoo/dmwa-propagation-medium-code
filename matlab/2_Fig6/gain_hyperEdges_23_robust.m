function hyperEdges = gain_hyperEdges_23_robust(datanew, retainedFraction)
% Extract 2-node and 3-node hyperedges for the Fig. 6 robustness grid.
%
% Main manuscript scripts pass a seizure-level elbow PLV threshold. Fig. 6
% instead stress-tests retained PLV edge fractions from 0.45 to 0.90 while
% keeping the same 2-node edge plus closed 3-node triangle hyperedge rule.

plvMatrix = compute_plv_matrix_local(datanew);

if nargin < 2 || isempty(retainedFraction)
    plvThreshold = select_plv_elbow_threshold_local(plvMatrix);
else
    plvThreshold = threshold_from_retained_fraction(plvMatrix, retainedFraction);
end

connection_matrix = double(plvMatrix >= plvThreshold);
nCh = size(connection_matrix, 1);
connection_matrix(1:nCh+1:end) = 0;

triangles = [];
for i = 1:nCh-2
    for j = i+1:nCh-1
        if ~connection_matrix(i, j)
            continue;
        end
        for k = j+1:nCh
            if connection_matrix(i, k) && connection_matrix(j, k)
                triangles = [triangles; i, j, k]; %#ok<AGROW>
            end
        end
    end
end

[row, col] = find(triu(connection_matrix, 1));
if isempty(row)
    hyperEdges_dim2 = cell(0, 1);
else
    hyperEdges_dim2 = arrayfun(@(ii) [row(ii), col(ii)], 1:numel(row), ...
        'UniformOutput', false).';
end

if isempty(triangles)
    hyperEdges_dim3 = cell(0, 1);
else
    hyperEdges_dim3 = mat2cell(triangles, ones(size(triangles, 1), 1), 3);
end

hyperEdges = [hyperEdges_dim2; hyperEdges_dim3(:)];
hyperEdges = unique_hyperedges(hyperEdges);
end


function plvMatrix = compute_plv_matrix_local(datanew)
nCh = size(datanew, 1);
phaseSig = angle(hilbert(datanew.')).';
plvMatrix = eye(nCh);

for i = 1:nCh-1
    for j = i+1:nCh
        dphi = phaseSig(i, :) - phaseSig(j, :);
        val = abs(mean(exp(1i * dphi)));
        plvMatrix(i, j) = val;
        plvMatrix(j, i) = val;
    end
end
end


function threshold = threshold_from_retained_fraction(plvMatrix, retainedFraction)
vals = plvMatrix(triu(true(size(plvMatrix)), 1));
vals = vals(isfinite(vals) & vals > 0);

if isempty(vals)
    threshold = inf;
    return;
end

retainedFraction = max(0, min(1, retainedFraction));
vals = sort(vals, 'descend');
idx = round(numel(vals) * retainedFraction);
idx = max(1, min(numel(vals), idx));
threshold = vals(idx);
end


function threshold = select_plv_elbow_threshold_local(plvMatrix)
vals = plvMatrix(triu(true(size(plvMatrix)), 1));
vals = vals(isfinite(vals) & vals > 0);

if isempty(vals)
    threshold = inf;
    return;
end

candidates = 0:0.0025:1;
density = arrayfun(@(x) mean(vals >= x), candidates);
threshold = elbow_from_curve_local(candidates, density);
end


function threshold = elbow_from_curve_local(x, y)
if numel(x) < 3 || (max(x) - min(x)) == 0 || (max(y) - min(y)) == 0
    threshold = median(x);
    return;
end

xn = (x - min(x)) ./ (max(x) - min(x));
yn = (y - min(y)) ./ (max(y) - min(y));
p1 = [xn(1), yn(1)];
p2 = [xn(end), yn(end)];
lineVec = p2 - p1;

dist = abs(lineVec(1) .* (p1(2) - yn) - (p1(1) - xn) .* lineVec(2)) ./ norm(lineVec);
[~, idx] = max(dist);
threshold = x(idx);
end
