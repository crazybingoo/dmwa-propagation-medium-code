function [threshold, curve] = select_seizure_plv_elbow_threshold(X1, fs, winLenSec, windowOffsets)
% Select one PLV threshold per seizure from the density-threshold elbow.
% The selected threshold is applied to every window from the same seizure.

if nargin < 4 || isempty(windowOffsets)
    nSamples = size(X1, 2);
    totalIterations = floor(nSamples / fs) - winLenSec;
    windowOffsets = 0:totalIterations;
end

allVals = [];
winSamples = round(winLenSec * fs);

for ii = 1:numel(windowOffsets)
    startPoint = fs * windowOffsets(ii) + 1;
    endPoint = startPoint + winSamples - 1;
    if startPoint < 1 || endPoint > size(X1, 2)
        continue;
    end

    plvMatrix = compute_plv_matrix_for_threshold(X1(:, startPoint:endPoint));
    vals = plvMatrix(triu(true(size(plvMatrix)), 1));
    vals = vals(isfinite(vals) & vals > 0);
    allVals = [allVals; vals(:)]; %#ok<AGROW>
end

if isempty(allVals)
    threshold = inf;
    curve = table([], [], 'VariableNames', {'plv_threshold', 'mean_density'});
    return;
end

candidates = (0:0.0025:1).';
density = nan(numel(candidates), 1);

for i = 1:numel(candidates)
    density(i) = mean(allVals >= candidates(i));
end

threshold = elbow_from_density_curve(candidates, density);
curve = table(candidates, density, 'VariableNames', {'plv_threshold', 'mean_density'});
end


function plvMatrix = compute_plv_matrix_for_threshold(dat)
nCh = size(dat, 1);
phaseSig = angle(hilbert(dat.')).';
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


function threshold = elbow_from_density_curve(x, y)
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

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
