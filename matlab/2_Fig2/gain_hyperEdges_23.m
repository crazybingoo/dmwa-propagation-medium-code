function hyperEdges = gain_hyperEdges_23(datanew, plvThreshold)
% Extract 2-node and 3-node hyperedges from a PLV-thresholded scaffold.
% The manuscript analysis passes a seizure-level PLV threshold selected from
% the elbow of the density-threshold curve.

plvMatrix = compute_plv_matrix_local(datanew);

if nargin < 2 || isempty(plvThreshold)
    plvThreshold = select_plv_elbow_threshold_local(plvMatrix);
end

connection_matrix = double(plvMatrix >= plvThreshold);
nCh = size(connection_matrix, 1);
connection_matrix(1:nCh+1:end) = 0;

% 3-node hyperedges are closed triangles in the binary scaffold.
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
hyperEdges = unique_hyperedges_local(hyperEdges);
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


function hyperEdges = unique_hyperedges_local(hyperEdges)
if isempty(hyperEdges)
    return;
end

keys = cellfun(@(x) sprintf('%d_', sort(x(:).')), hyperEdges, 'UniformOutput', false);
[~, ia] = unique(keys, 'stable');
hyperEdges = hyperEdges(ia);
end
