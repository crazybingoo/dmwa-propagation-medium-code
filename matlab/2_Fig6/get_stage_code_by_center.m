function code = get_stage_code_by_center(centerSec, ictalStartCenterSec, ictalEndCenterSec)

% 1 = pre
% 2 = early
% 3 = mid
% 4 = late
% 5 = post

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