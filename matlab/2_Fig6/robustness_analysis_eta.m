clc; clear; close all;

%% ===================== 用户设置 =====================
rootDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig6';   % 改成你的数据根目录
outDir  = fullfile(rootDir, 'robustness_out_eta');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fs = 1024;
stepSec = 1;

windowSecList   = [2 3 4 5];
plvQuantileList = [0.45 0.50 0.55 0.60 0.65];

defaultWindowSec   = 3;
defaultPlvQuantile = 0.55;

%% ===================== 24次发作信息 =====================
seizures = get_seizure_info(rootDir);

if isempty(seizures)
    error('get_seizure_info.m 里没有填入发作信息。');
end

%% ===================== 主循环 =====================
allSummary = struct([]);
cnt = 0;

for iw = 1:numel(windowSecList)
    for iq = 1:numel(plvQuantileList)

        windowSec = windowSecList(iw);
        plvQ = plvQuantileList(iq);

        fprintf('\n==================================================\n');
        fprintf('Running eta: window = %.1f s, PLV quantile = %.2f\n', windowSec, plvQ);
        fprintf('==================================================\n');

        fprintf('  [%02d/%02d] %s\n', 1, numel(seizures), seizures(1).id);
        firstResult = analyze_one_seizure_eta(seizures(1), fs, windowSec, stepSec, plvQ);

        seizureResults = repmat(firstResult, numel(seizures), 1);
        seizureResults(1) = firstResult;

        for k = 2:numel(seizures)
            fprintf('  [%02d/%02d] %s\n', k, numel(seizures), seizures(k).id);
            seizureResults(k) = analyze_one_seizure_eta( ...
                seizures(k), fs, windowSec, stepSec, plvQ);
        end

        summary = summarize_robustness_results_eta(seizureResults, windowSec, plvQ);

        cnt = cnt + 1;
        allSummary(cnt).windowSec         = windowSec;
        allSummary(cnt).plvQuantile       = plvQ;
        allSummary(cnt).nSeizures         = summary.nSeizures;
        allSummary(cnt).preMean           = summary.preMean;
        allSummary(cnt).earlyMean         = summary.earlyMean;
        allSummary(cnt).midMean           = summary.midMean;
        allSummary(cnt).lateMean          = summary.lateMean;
        allSummary(cnt).postMean          = summary.postMean;
        allSummary(cnt).ictalMean         = summary.ictalMean;
        allSummary(cnt).trendKeepRate     = summary.trendKeepRate;
        allSummary(cnt).p_pre_vs_ictal    = summary.p_pre_vs_ictal;
        allSummary(cnt).dz_pre_vs_ictal   = summary.dz_pre_vs_ictal;
        allSummary(cnt).nValidPairs       = summary.nValidPairs;

        save(fullfile(outDir, sprintf('detail_eta_win%.1f_q%.2f.mat', windowSec, plvQ)), ...
            'seizureResults', 'summary');
    end
end

%% ===================== 保存汇总 =====================
SummaryTable = struct2table(allSummary);
writetable(SummaryTable, fullfile(outDir, 'robustness_summary_eta.csv'));

save(fullfile(outDir, 'robustness_summary_eta.mat'), ...
    'allSummary', 'SummaryTable', 'seizures', ...
    'windowSecList', 'plvQuantileList', ...
    'defaultWindowSec', 'defaultPlvQuantile');

%% ===================== 画图 =====================
plot_robustness_heatmap_eta(allSummary, windowSecList, plvQuantileList, ...
    defaultWindowSec, defaultPlvQuantile, outDir);

disp('All done.');
