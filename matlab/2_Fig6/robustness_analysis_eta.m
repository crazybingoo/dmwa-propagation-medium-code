clc; clear; close all;

%% User settings
rootDir = 'example_project\2_Fig6';
outDir  = fullfile(rootDir, 'robustness_out_eta');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fs = 1024;
stepSec = 1;

windowSecList = 1:5;
retainedFractionList = 0.45:0.05:0.90;

defaultWindowSec = 3;
defaultRetainedFraction = 0.55;

%% Seizure information
seizures = get_seizure_info(rootDir);

if isempty(seizures)
    error('No seizure information was returned by get_seizure_info.m.');
end

%% Main robustness loop
allSummary = struct([]);
cnt = 0;

for iw = 1:numel(windowSecList)
    for iq = 1:numel(retainedFractionList)
        windowSec = windowSecList(iw);
        retainedFraction = retainedFractionList(iq);

        fprintf('\n==================================================\n');
        fprintf('Running eta: window = %.1f s, retained PLV edge fraction = %.2f\n', windowSec, retainedFraction);
        fprintf('==================================================\n');

        fprintf('  [%02d/%02d] %s\n', 1, numel(seizures), seizures(1).id);
        firstResult = analyze_one_seizure_eta(seizures(1), fs, windowSec, stepSec, retainedFraction);

        seizureResults = repmat(firstResult, numel(seizures), 1);
        seizureResults(1) = firstResult;

        for k = 2:numel(seizures)
            fprintf('  [%02d/%02d] %s\n', k, numel(seizures), seizures(k).id);
            seizureResults(k) = analyze_one_seizure_eta( ...
                seizures(k), fs, windowSec, stepSec, retainedFraction);
        end

        summary = summarize_robustness_results_eta(seizureResults, windowSec, retainedFraction);

        cnt = cnt + 1;
        allSummary(cnt).windowSec        = windowSec;
        allSummary(cnt).retainedFraction = retainedFraction;
        allSummary(cnt).nSeizures        = summary.nSeizures;
        allSummary(cnt).preMean          = summary.preMean;
        allSummary(cnt).earlyMean        = summary.earlyMean;
        allSummary(cnt).midMean          = summary.midMean;
        allSummary(cnt).lateMean         = summary.lateMean;
        allSummary(cnt).postMean         = summary.postMean;
        allSummary(cnt).ictalMean        = summary.ictalMean;
        allSummary(cnt).trendKeepRate    = summary.trendKeepRate;
        allSummary(cnt).p_pre_vs_ictal   = summary.p_pre_vs_ictal;
        allSummary(cnt).dz_pre_vs_ictal  = summary.dz_pre_vs_ictal;
        allSummary(cnt).nValidPairs      = summary.nValidPairs;

        save(fullfile(outDir, sprintf('detail_eta_win%.1f_retained%.2f.mat', windowSec, retainedFraction)), ...
            'seizureResults', 'summary');
    end
end

%% Save summary
SummaryTable = struct2table(allSummary);
writetable(SummaryTable, fullfile(outDir, 'robustness_summary_eta.csv'));

save(fullfile(outDir, 'robustness_summary_eta.mat'), ...
    'allSummary', 'SummaryTable', 'seizures', ...
    'windowSecList', 'retainedFractionList', ...
    'defaultWindowSec', 'defaultRetainedFraction');

%% Plot
plot_robustness_heatmap_eta(allSummary, windowSecList, retainedFractionList, ...
    defaultWindowSec, defaultRetainedFraction, outDir);

disp('All done.');
