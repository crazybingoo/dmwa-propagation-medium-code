function summary = summarize_robustness_results_eta(seizureResults, windowSec, plvQuantile)

n = numel(seizureResults);

preVals   = nan(n, 1);
earlyVals = nan(n, 1);
midVals   = nan(n, 1);
lateVals  = nan(n, 1);
postVals  = nan(n, 1);
ictalVals = nan(n, 1);
trendKeep = nan(n, 1);

for i = 1:n
    sm = seizureResults(i).stageMean;
    preVals(i)   = sm.pre;
    earlyVals(i) = sm.early;
    midVals(i)   = sm.mid;
    lateVals(i)  = sm.late;
    postVals(i)  = sm.post;
    ictalVals(i) = sm.ictal;

    if ~isnan(sm.pre) && ~isnan(sm.early) && ~isnan(sm.mid) && ~isnan(sm.late)
        trendKeep(i) = double(sm.pre > sm.early && sm.pre > sm.mid && sm.pre > sm.late);
    end
end

validPair = ~isnan(preVals) & ~isnan(ictalVals);
nValidPairs = sum(validPair);

if nValidPairs >= 3
    try
        p_pre_vs_ictal = signrank(preVals(validPair), ictalVals(validPair));
    catch
        [~, p_pre_vs_ictal] = ttest(preVals(validPair), ictalVals(validPair));
    end
    dz_pre_vs_ictal = paired_effect_size_d(preVals(validPair), ictalVals(validPair));
else
    p_pre_vs_ictal = NaN;
    dz_pre_vs_ictal = NaN;
end

summary.windowSec       = windowSec;
summary.plvQuantile     = plvQuantile;
summary.nSeizures       = n;

summary.preMean         = mean(preVals, 'omitnan');
summary.earlyMean       = mean(earlyVals, 'omitnan');
summary.midMean         = mean(midVals, 'omitnan');
summary.lateMean        = mean(lateVals, 'omitnan');
summary.postMean        = mean(postVals, 'omitnan');
summary.ictalMean       = mean(ictalVals, 'omitnan');

summary.trendKeepRate   = mean(trendKeep, 'omitnan');
summary.p_pre_vs_ictal  = p_pre_vs_ictal;
summary.dz_pre_vs_ictal = dz_pre_vs_ictal;
summary.nValidPairs     = nValidPairs;

end
