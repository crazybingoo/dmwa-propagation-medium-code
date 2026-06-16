function result = analyze_one_seizure_eta(seizure, fs, windowSec, stepSec, retainedFraction)

dataFile = seizure.file;
if ~isfile(dataFile)
    error('Data file not found: %s', dataFile);
end

S = load(dataFile);
if ~isfield(S, 'X1')
    error('File %s does not contain variable X1.', dataFile);
end
X1 = S.X1;

winSamples  = round(windowSec * fs);
stepSamples = round(stepSec * fs);

if winSamples <= 0 || stepSamples <= 0
    error('windowSec and stepSec must be positive.');
end

nSamples = size(X1, 2);
if nSamples < winSamples
    error('Data is shorter than the target window length.');
end

startIdx = 1 : stepSamples : (nSamples - winSamples + 1);
nWin = numel(startIdx);

eta_all   = nan(nWin, 1);
centerSec = nan(nWin, 1);
numHE_all = nan(nWin, 1);

sameAsReference = ...
    isfield(seizure, 'refWindowSec') && isfield(seizure, 'refStepSec') && ...
    abs(windowSec - seizure.refWindowSec) < 1e-12 && ...
    abs(stepSec   - seizure.refStepSec)   < 1e-12;

if sameAsReference
    stageCode = build_stage_code_from_reference_indices(seizure, nWin);
else
    stageCode = build_stage_code_by_reference_overlap(seizure, fs, startIdx, winSamples, nSamples);
end

for i = 1:nWin
    s1 = startIdx(i);
    s2 = s1 + winSamples - 1;

    datanew = X1(:, s1:s2);

    HE = gain_hyperEdges_23_robust(datanew, retainedFraction);
    numHE_all(i) = numel(HE);

    W = build_W_from_hyperedges(HE);
    eta_all(i) = compute_eta_from_W(W);

    centerSec(i) = ((s1 - 1) + (s2 - 1)) / 2 / fs;
end

stageMean.pre   = mean(eta_all(stageCode == 1), 'omitnan');
stageMean.early = mean(eta_all(stageCode == 2), 'omitnan');
stageMean.mid   = mean(eta_all(stageCode == 3), 'omitnan');
stageMean.late  = mean(eta_all(stageCode == 4), 'omitnan');
stageMean.post  = mean(eta_all(stageCode == 5), 'omitnan');
stageMean.ictal = mean(eta_all(stageCode >= 2 & stageCode <= 4), 'omitnan');

result.id                  = seizure.id;
result.file                = seizure.file;
result.windowSec           = windowSec;
result.stepSec             = stepSec;
result.retainedFraction    = retainedFraction;
result.eta_all             = eta_all;
result.centerSec           = centerSec;
result.stageCode           = stageCode;
result.stageMean           = stageMean;
result.numHE_all           = numHE_all;
result.sameAsReference     = sameAsReference;

if isfield(seizure, 'ictalStartCenterSec')
    result.ictalStartCenterSec = seizure.ictalStartCenterSec;
else
    result.ictalStartCenterSec = NaN;
end

if isfield(seizure, 'ictalEndCenterSec')
    result.ictalEndCenterSec = seizure.ictalEndCenterSec;
else
    result.ictalEndCenterSec = NaN;
end

end


function stageCode = build_stage_code_from_reference_indices(seizure, nWin)

stageCode = nan(nWin, 1);

[pre_idx, early_idx, mid_idx, late_idx, post_idx] = ...
    get_clean_stage_indices(seizure, nWin);

stageCode(pre_idx)   = 1;
stageCode(early_idx) = 2;
stageCode(mid_idx)   = 3;
stageCode(late_idx)  = 4;
stageCode(post_idx)  = 5;

end


function stageCode = build_stage_code_by_reference_overlap(seizure, fs, startIdx, winSamples, nSamples)

nWin = numel(startIdx);
stageCode = nan(nWin, 1);

if ~isfield(seizure, 'refWindowSec') || ~isfield(seizure, 'refStepSec')
    error('seizure.refWindowSec / seizure.refStepSec is required.');
end

refWinSamples  = round(seizure.refWindowSec * fs);
refStepSamples = round(seizure.refStepSec   * fs);

if refWinSamples <= 0 || refStepSamples <= 0
    error('Invalid reference window/step setting.');
end

nRefWin = floor((nSamples - refWinSamples) / refStepSamples) + 1;
if nRefWin < 1
    error('Reference grid is invalid for current data length.');
end

[pre_idx, early_idx, mid_idx, late_idx, post_idx] = ...
    get_clean_stage_indices(seizure, nRefWin);

stageIntervals = nan(5, 2);
stageIntervals(1, :) = idx_to_time_interval(pre_idx,   refStepSamples, refWinSamples, nSamples);
stageIntervals(2, :) = idx_to_time_interval(early_idx, refStepSamples, refWinSamples, nSamples);
stageIntervals(3, :) = idx_to_time_interval(mid_idx,   refStepSamples, refWinSamples, nSamples);
stageIntervals(4, :) = idx_to_time_interval(late_idx,  refStepSamples, refWinSamples, nSamples);
stageIntervals(5, :) = idx_to_time_interval(post_idx,  refStepSamples, refWinSamples, nSamples);

for i = 1:nWin
    s1 = startIdx(i);
    s2 = s1 + winSamples - 1;

    overlaps = zeros(1, 5);
    for k = 1:5
        overlaps(k) = interval_overlap_len([s1, s2], stageIntervals(k, :));
    end

    maxOv = max(overlaps);

    if maxOv <= 0
        centerSec = ((s1 - 1) + (s2 - 1)) / 2 / fs;
        stageCode(i) = fallback_stage_code_by_center(centerSec, seizure);
        continue;
    end

    tied = find(overlaps == maxOv);

    if numel(tied) == 1
        stageCode(i) = tied;
    else
        centerSec = ((s1 - 1) + (s2 - 1)) / 2 / fs;
        code0 = fallback_stage_code_by_center(centerSec, seizure);

        if any(tied == code0)
            stageCode(i) = code0;
        else
            stageCode(i) = tied(1);
        end
    end
end

end


function [pre_idx, early_idx, mid_idx, late_idx, post_idx] = get_clean_stage_indices(seizure, nMax)

pre_idx   = sanitize_idx(seizure.pre_idx,   nMax);
ictal_idx = sanitize_idx(seizure.ictal_idx, nMax);
post_idx  = sanitize_idx(seizure.post_idx,  nMax);

ictal_idx = setdiff(ictal_idx, pre_idx, 'stable');
post_idx  = setdiff(post_idx, [pre_idx, ictal_idx], 'stable');

[early_idx, mid_idx, late_idx] = split_ictal_into_three(ictal_idx);

end


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


function [early_idx, mid_idx, late_idx] = split_ictal_into_three(ictal_idx)

n = numel(ictal_idx);

if n == 0
    early_idx = [];
    mid_idx   = [];
    late_idx  = [];
    return;
end

cut1 = floor(n / 3);
cut2 = floor(2 * n / 3);

early_idx = ictal_idx(1:cut1);
mid_idx   = ictal_idx(cut1 + 1 : cut2);
late_idx  = ictal_idx(cut2 + 1 : end);

if isempty(early_idx) && isempty(mid_idx) && ~isempty(late_idx) && n == 1
    mid_idx  = late_idx;
    late_idx = [];
elseif isempty(early_idx) && isempty(mid_idx) && isempty(late_idx) && n > 0
    late_idx = ictal_idx;
end

end


function interval = idx_to_time_interval(idxVec, stepSamples, winSamples, nSamples)

if isempty(idxVec)
    interval = [NaN, NaN];
    return;
end

s1 = (idxVec(1)   - 1) * stepSamples + 1;
s2 = (idxVec(end) - 1) * stepSamples + winSamples;

s1 = max(1, s1);
s2 = min(nSamples, s2);

if s2 < s1
    interval = [NaN, NaN];
else
    interval = [s1, s2];
end

end


function ov = interval_overlap_len(a, b)

if any(~isfinite(b))
    ov = 0;
    return;
end

left  = max(a(1), b(1));
right = min(a(2), b(2));

ov = max(0, right - left + 1);

end


function code = fallback_stage_code_by_center(centerSec, seizure)

if isfield(seizure, 'ictalStartCenterSec') && isfield(seizure, 'ictalEndCenterSec')
    ictalStartCenterSec = seizure.ictalStartCenterSec;
    ictalEndCenterSec   = seizure.ictalEndCenterSec;
else
    code = 3;
    return;
end

if centerSec < ictalStartCenterSec
    code = 1;
    return;
end

if centerSec > ictalEndCenterSec
    code = 5;
    return;
end

dur = ictalEndCenterSec - ictalStartCenterSec;
if dur <= 0
    code = 3;
    return;
end

frac = (centerSec - ictalStartCenterSec) / dur;

if frac < 1/3
    code = 2;
elseif frac < 2/3
    code = 3;
else
    code = 4;
end

end
