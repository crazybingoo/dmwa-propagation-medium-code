%% Hindmarsh-Rose signal simulation for the DMW-HLG eta pipeline
% This script tests whether epileptiform neural dynamics can induce eta
% changes when the same SEEG -> DMW-HLG construction used in paper 2 is held
% fixed. The HR model only generates simulated X1 signals. Hyperedge
% construction and eta calculation reuse the native paper-2
% functions:
%   gain_hyperEdges_23.m
%   build_W_from_hyperedges.m
%   compute_eta_from_W.m

clc; clear; close all;

%% Paths
rootDir = 'example_data\8 - n_u_v';
outDir  = fullfile(rootDir, '2_Fig9_HR_eta_simulation');
if ~exist(outDir, 'dir'), mkdir(outDir); end

addpath(fullfile(rootDir, '2'));
addpath(fullfile(rootDir, '2_Fig6'));
addpath(outDir);

requiredFunctions = {'gain_hyperEdges_23', 'select_seizure_plv_elbow_threshold', ...
    'build_W_from_hyperedges', 'compute_eta_from_W'};
for i = 1:numel(requiredFunctions)
    if exist(requiredFunctions{i}, 'file') ~= 2
        error('Required paper-2 function not found on path: %s', requiredFunctions{i});
    end
end

%% Simulation parameters
params = struct();
params.seed              = 260519;
params.nNodes            = 18;
params.coreNodes         = 1:4;
params.relayNodes        = 5:9;
params.peripheryNodes    = 10:18;
params.betaLevels        = [0 0.25 0.50 0.75 1.00];
params.numTrials         = 6;

% HR time is used as an analysis grid for simulated SEEG-like traces.
params.dt                = 0.02;
params.analysisFs        = round(1 / params.dt);
params.durationSec       = 14;
params.burnSec           = 4;
params.winLenSec         = 3;
params.stepSec           = 1;
params.winSamples        = round(params.winLenSec * params.analysisFs);
params.stepSamples       = round(params.stepSec * params.analysisFs);

% Hindmarsh-Rose parameters.
params.a                 = 1;
params.b                 = 3;
params.c                 = 1;
params.d                 = 5;
params.r                 = 0.006;
params.s                 = 4;
params.xR                = -1.6;
params.Ibase             = 3.05;
params.noiseAmp          = 0.025;
params.commonFreq        = 0.18;
params.couplingBase      = 0.010;
params.couplingSlope     = 1.200;
params.commonBase        = 0.015;
params.commonSlope       = 0.050;
params.heterogeneityDrop = 0.900;

fprintf('HR-DMW-HLG eta simulation\n');
fprintf('Output directory: %s\n', outDir);
fprintf('Native DMW-HLG functions:\n');
for i = 1:numel(requiredFunctions)
    fprintf('  %s -> %s\n', requiredFunctions{i}, which(requiredFunctions{i}));
end

rng(params.seed);
nBeta = numel(params.betaLevels);
nRowsEstimate = params.numTrials * nBeta * 64;
rows = cell(nRowsEstimate, 1);
rowCount = 0;
trialSummary = cell(params.numTrials * nBeta, 1);
summaryCount = 0;
representative = struct();

totalTimer = tic;
for trial = 1:params.numTrials
    trialSpec = make_trial_spec(params, trial);

    for ib = 1:nBeta
        beta = params.betaLevels(ib);
        sim = simulate_hr_network(params, trialSpec, beta);
        X1 = sim.X1;

        if trial == 1 && (abs(beta - params.betaLevels(1)) < eps || abs(beta - params.betaLevels(end)) < eps)
            fieldName = sprintf('beta_%03d', round(beta * 100));
            representative.(fieldName) = sim;
        end

        nSamples = size(X1, 2);
        startIdx = 1:params.stepSamples:(nSamples - params.winSamples + 1);
        nWin = numel(startIdx);
        plvThreshold = select_seizure_plv_elbow_threshold( ...
            X1, params.analysisFs, params.winLenSec, 0:(nWin - 1));

        etaVals      = nan(nWin, 1);
        avgWVals     = nan(nWin, 1);
        lambdaVals   = nan(nWin, 1);
        pairDensVals = nan(nWin, 1);
        heVals       = nan(nWin, 1);
        he2Vals      = nan(nWin, 1);
        he3Vals      = nan(nWin, 1);
        meanSizeVals = nan(nWin, 1);
        meanPlvVals  = nan(nWin, 1);
        corePartVals = nan(nWin, 1);
        coreFracVals = nan(nWin, 1);

        for iw = 1:nWin
            idx1 = startIdx(iw);
            idx2 = idx1 + params.winSamples - 1;
            datanew = X1(:, idx1:idx2);

            % Native paper-2 signal-to-hyperedge step.
            HE = gain_hyperEdges_23(datanew, plvThreshold);

            % Native paper-2 hyperedge-to-DMWA and eta steps.
            W = build_W_from_hyperedges(HE);
            etaNative = compute_eta_from_W(W);
            [etaCheck, avgW, lambda1] = eta_parts_from_W(W);
            if abs(etaNative - etaCheck) > 1e-10
                warning('eta mismatch in trial %d beta %.2f window %d: native %.12g, check %.12g', ...
                    trial, beta, iw, etaNative, etaCheck);
            end

            diag = hyperedge_diagnostics(HE, params.coreNodes, params.nNodes);
            plvMat = compute_plv_matrix_local(datanew);
            meanPlv = mean(plvMat(triu(true(params.nNodes), 1)), 'omitnan');

            etaVals(iw)      = etaNative;
            avgWVals(iw)     = avgW;
            lambdaVals(iw)   = lambda1;
            pairDensVals(iw) = diag.pairDensity;
            heVals(iw)       = diag.numHE;
            he2Vals(iw)      = diag.numHE2;
            he3Vals(iw)      = diag.numHE3;
            meanSizeVals(iw) = diag.meanSize;
            meanPlvVals(iw)  = meanPlv;
            corePartVals(iw) = diag.coreParticipation;
            coreFracVals(iw) = diag.coreNodeFraction;

            rowCount = rowCount + 1;
            if rowCount > numel(rows)
                rows = [rows; cell(numel(rows), 1)]; %#ok<AGROW>
            end
            rows{rowCount} = {trial, beta, iw, ...
                ((idx1 - 1) + (idx2 - 1)) / 2 / params.analysisFs, ...
                etaNative, avgW, lambda1, diag.pairDensity, meanPlv, ...
                diag.numHE, diag.numHE2, diag.numHE3, diag.meanSize, ...
                diag.coreParticipation, diag.coreNodeFraction, plvThreshold};
        end

        summaryCount = summaryCount + 1;
        trialSummary{summaryCount} = {trial, beta, ...
            mean(etaVals, 'omitnan'), mean(avgWVals, 'omitnan'), mean(lambdaVals, 'omitnan'), ...
            mean(pairDensVals, 'omitnan'), mean(meanPlvVals, 'omitnan'), ...
            mean(heVals, 'omitnan'), mean(he2Vals, 'omitnan'), mean(he3Vals, 'omitnan'), ...
            mean(meanSizeVals, 'omitnan'), mean(corePartVals, 'omitnan'), ...
            mean(coreFracVals, 'omitnan'), plvThreshold, nWin};

        fprintf('trial %02d/%02d | beta %.2f | windows %d | eta %.4f | R %.4f | lambda1 %.4f | pair density %.3f | nHE %.1f\n', ...
            trial, params.numTrials, beta, nWin, ...
            mean(etaVals, 'omitnan'), mean(avgWVals, 'omitnan'), ...
            mean(lambdaVals, 'omitnan'), mean(pairDensVals, 'omitnan'), ...
            mean(heVals, 'omitnan'));
    end
end

rows = rows(1:rowCount);
trialSummary = trialSummary(1:summaryCount);

WindowTable = cell2table(vertcat(rows{:}), 'VariableNames', { ...
    'trial', 'beta', 'window_idx', 'center_time', ...
    'eta', 'avgW', 'lambda1', 'pair_density', 'mean_plv', ...
    'num_hyperedges', 'num_HE2', 'num_HE3', 'mean_hyperedge_size', ...
    'core_participation', 'core_node_fraction', 'plv_elbow_threshold'});

TrialTable = cell2table(vertcat(trialSummary{:}), 'VariableNames', { ...
    'trial', 'beta', 'eta_mean', 'avgW_mean', 'lambda1_mean', ...
    'pair_density_mean', 'mean_plv_mean', 'num_hyperedges_mean', ...
    'num_HE2_mean', 'num_HE3_mean', 'mean_hyperedge_size_mean', ...
    'core_participation_mean', 'core_node_fraction_mean', ...
    'plv_elbow_threshold', 'n_windows'});

BetaSummary = summarize_by_beta(TrialTable);
StatsTable = trend_stats(TrialTable);

writetable(WindowTable, fullfile(outDir, 'HR_eta_window_metrics.csv'));
writetable(TrialTable,  fullfile(outDir, 'HR_eta_trial_metrics.csv'));
writetable(BetaSummary, fullfile(outDir, 'HR_eta_beta_summary.csv'));
writetable(StatsTable,  fullfile(outDir, 'HR_eta_trend_stats.csv'));

save(fullfile(outDir, 'HR_eta_simulation_workspace.mat'), ...
    'params', 'WindowTable', 'TrialTable', 'BetaSummary', 'StatsTable', 'representative');

plot_hr_eta_summary(params, TrialTable, BetaSummary, representative, outDir);

fprintf('Finished in %.2f s.\n', toc(totalTimer));
fprintf('Saved outputs to: %s\n', outDir);

%% ========================= local functions =========================

function spec = make_trial_spec(params, trial)
    rng(params.seed + trial * 101);
    spec.A = build_fixed_coupling(params);
    spec.x0 = -1.4 + 0.55 * randn(params.nNodes, 1);
    spec.y0 = -7.0 + 0.55 * randn(params.nNodes, 1);
    spec.z0 =  3.0 + 0.30 * randn(params.nNodes, 1);
    spec.Ioffset = 0.70 * randn(params.nNodes, 1);

    nStepsTotal = round(params.durationSec / params.dt);
    spec.noise = randn(params.nNodes, nStepsTotal);
    spec.commonPhase = 2 * pi * rand();
end

function A = build_fixed_coupling(params)
    N = params.nNodes;
    A = zeros(N);
    core = params.coreNodes;
    relay = params.relayNodes;
    peri = params.peripheryNodes;

    A(core, core) = 1.0;
    A(relay, relay) = 0.45;
    A(peri, peri) = 0.16;
    A(core, relay) = 0.75;
    A(relay, core) = 0.55;
    A(relay, peri) = 0.30;
    A(peri, relay) = 0.20;
    A(core, peri) = 0.12;
    A(peri, core) = 0.08;
    A(1:N+1:end) = 0;

    jitter = 0.75 + 0.5 * rand(N);
    A = A .* jitter;
    A = max(A, A');
    A(1:N+1:end) = 0;

    rowSum = sum(A, 2);
    rowSum(rowSum == 0) = 1;
    A = A ./ rowSum;
end

function sim = simulate_hr_network(params, spec, beta)
    dt = params.dt;
    nStepsTotal = round(params.durationSec / dt);
    burnSteps = round(params.burnSec / dt);
    N = params.nNodes;

    x = spec.x0;
    y = spec.y0;
    z = spec.z0;
    X = nan(N, nStepsTotal - burnSteps);

    coreWeight = zeros(N, 1);
    coreWeight(params.coreNodes) = 1.00;
    coreWeight(params.relayNodes) = 0.55;
    coreWeight(params.peripheryNodes) = 0.20;

    excitability = zeros(N, 1);
    excitability(params.coreNodes) = 0.46;
    excitability(params.relayNodes) = 0.25;
    excitability(params.peripheryNodes) = 0.08;

    g = params.couplingBase + params.couplingSlope * beta;
    commonAmp = params.commonBase + params.commonSlope * beta;
    rowSum = sum(spec.A, 2);

    outIdx = 0;
    for step = 1:nStepsTotal
        t = (step - 1) * dt;
        commonDrive = commonAmp * sin(2 * pi * params.commonFreq * t + spec.commonPhase);
        heterogeneity = (1 - params.heterogeneityDrop * beta) * spec.Ioffset;
        inputCurrent = params.Ibase + heterogeneity + beta * excitability + commonDrive * coreWeight;
        inputCurrent = inputCurrent + params.noiseAmp * spec.noise(:, step);

        coupling = g * (spec.A * x - rowSum .* x);

        dx = y - params.a * x.^3 + params.b * x.^2 - z + inputCurrent + coupling;
        dy = params.c - params.d * x.^2 - y;
        dz = params.r * (params.s * (x - params.xR) - z);

        x = x + dt * dx;
        y = y + dt * dy;
        z = z + dt * dz;

        if step > burnSteps
            outIdx = outIdx + 1;
            X(:, outIdx) = x;
        end
    end

    X = X(:, 1:outIdx);
    X = zscore_channels(X);

    sim.X1 = X;
    sim.beta = beta;
    sim.time = (0:(size(X, 2)-1)) / params.analysisFs;
end

function Xz = zscore_channels(X)
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd <= eps | ~isfinite(sd)) = 1;
    Xz = (X - mu) ./ sd;
    Xz(~isfinite(Xz)) = 0;
end

function [eta, avgW, lambda1] = eta_parts_from_W(W)
    if isempty(W)
        eta = 0; avgW = 0; lambda1 = 0; return;
    end
    W = double(W);
    n = size(W, 1);
    if n < 2 || size(W, 2) ~= n
        eta = 0; avgW = 0; lambda1 = 0; return;
    end
    W(~isfinite(W)) = 0;
    W(W < 0) = 0;
    W(1:n+1:end) = 0;
    avgW = sum(W(:)) / n;
    if avgW <= 0
        eta = 0; lambda1 = 0; return;
    end
    lambda1 = max(abs(eig(W)));
    if ~isfinite(lambda1) || lambda1 <= eps
        eta = 0;
    else
        eta = avgW / lambda1;
    end
end

function diag = hyperedge_diagnostics(HE, coreNodes, nNodes)
    if isempty(HE)
        diag.numHE = 0;
        diag.numHE2 = 0;
        diag.numHE3 = 0;
        diag.meanSize = 0;
        diag.pairDensity = 0;
        diag.coreParticipation = 0;
        diag.coreNodeFraction = 0;
        return;
    end

    sizes = cellfun(@numel, HE);
    isCore = cellfun(@(he) any(ismember(he, coreNodes)), HE);
    coreFrac = cellfun(@(he) sum(ismember(he, coreNodes)) / numel(he), HE);

    diag.numHE = numel(HE);
    diag.numHE2 = sum(sizes == 2);
    diag.numHE3 = sum(sizes == 3);
    diag.meanSize = mean(sizes, 'omitnan');
    diag.pairDensity = diag.numHE2 / nchoosek(nNodes, 2);
    diag.coreParticipation = mean(isCore, 'omitnan');
    diag.coreNodeFraction = mean(coreFrac, 'omitnan');
end

function plvMat = compute_plv_matrix_local(dat)
    [nCh, ~] = size(dat);
    phases = angle(hilbert(dat.')).';
    plvMat = eye(nCh);
    for i = 1:nCh
        for j = i+1:nCh
            val = abs(mean(exp(1i * (phases(i,:) - phases(j,:)))));
            plvMat(i,j) = val;
            plvMat(j,i) = val;
        end
    end
end

function BetaSummary = summarize_by_beta(TrialTable)
    betas = unique(TrialTable.beta);
    rows = cell(numel(betas), 1);
    metricNames = {'eta_mean', 'avgW_mean', 'lambda1_mean', 'pair_density_mean', ...
        'mean_plv_mean', 'num_hyperedges_mean', 'num_HE2_mean', 'num_HE3_mean', ...
        'mean_hyperedge_size_mean', 'core_participation_mean', 'core_node_fraction_mean'};

    for ib = 1:numel(betas)
        beta = betas(ib);
        idx = TrialTable.beta == beta;
        row = cell(1, 2 + numel(metricNames) * 2);
        row{1} = beta;
        row{2} = sum(idx);
        col = 3;
        for im = 1:numel(metricNames)
            x = TrialTable.(metricNames{im})(idx);
            row{col} = mean(x, 'omitnan');
            row{col+1} = std(x, 'omitnan') / sqrt(sum(~isnan(x)));
            col = col + 2;
        end
        rows{ib} = row;
    end

    names = {'beta', 'n_trials'};
    for im = 1:numel(metricNames)
        names{end+1} = [metricNames{im}, '_avg']; %#ok<AGROW>
        names{end+1} = [metricNames{im}, '_sem']; %#ok<AGROW>
    end
    BetaSummary = cell2table(vertcat(rows{:}), 'VariableNames', names);
end

function StatsTable = trend_stats(TrialTable)
    metricNames = {'eta_mean', 'avgW_mean', 'lambda1_mean', 'pair_density_mean', ...
        'mean_plv_mean', 'num_hyperedges_mean', 'num_HE3_mean', ...
        'core_participation_mean', 'core_node_fraction_mean'};
    rows = cell(numel(metricNames), 1);
    beta = TrialTable.beta;
    for im = 1:numel(metricNames)
        x = TrialTable.(metricNames{im});
        valid = isfinite(beta) & isfinite(x);
        if sum(valid) >= 3 && numel(unique(x(valid))) > 1
            [rho, pval] = corr(beta(valid), x(valid), 'Type', 'Spearman');
        else
            rho = NaN;
            pval = NaN;
        end
        rows{im} = {metricNames{im}, rho, pval, sum(valid)};
    end
    StatsTable = cell2table(vertcat(rows{:}), ...
        'VariableNames', {'metric', 'spearman_rho_vs_beta', 'p_value', 'n'});
end

function plot_hr_eta_summary(~, TrialTable, BetaSummary, representative, outDir)
    betas = BetaSummary.beta;
    fig = figure('Color', 'w', 'Position', [80 80 1180 820]);
    tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(1);
    hold(ax1, 'on');
    keys = fieldnames(representative);
    traceColors = [0.25 0.35 0.85; 0.85 0.25 0.20];
    labels = cell(numel(keys), 1);
    for k = 1:numel(keys)
        sim = representative.(keys{k});
        t = sim.time;
        idx = t <= min(8, max(t));
        offset = (0:7)' * 5.0;
        Xplot = sim.X1(1:8, idx) + offset;
        for ch = 1:8
            plot(ax1, t(idx), Xplot(ch,:), 'Color', traceColors(k,:), 'LineWidth', 0.75);
        end
        labels{k} = sprintf('\\beta = %.2f', sim.beta);
    end
    xlabel(ax1, 'Simulated time');
    ylabel(ax1, 'HR membrane potential, offset');
    title(ax1, 'Representative simulated neural signals');
    set(ax1, 'Box', 'off', 'TickDir', 'out');
    text(ax1, 0.02, 0.95, strjoin(labels, ' / '), 'Units', 'normalized', 'FontWeight', 'bold');

    ax2 = nexttile(2);
    yyaxis(ax2, 'left');
    errorbar(ax2, betas, BetaSummary.eta_mean_avg, BetaSummary.eta_mean_sem, '-o', ...
        'LineWidth', 1.8, 'MarkerSize', 5);
    ylabel(ax2, '\eta');
    yyaxis(ax2, 'right');
    errorbar(ax2, betas, BetaSummary.pair_density_mean_avg, BetaSummary.pair_density_mean_sem, '-s', ...
        'LineWidth', 1.4, 'MarkerSize', 5);
    ylabel(ax2, 'Pair-edge density');
    xlabel(ax2, 'Epileptiform drive \beta');
    title(ax2, '\eta versus density control');
    set(ax2, 'Box', 'off', 'TickDir', 'out');

    ax3 = nexttile(3);
    hold(ax3, 'on');
    plot_metric_with_sem(ax3, betas, BetaSummary.avgW_mean_avg, BetaSummary.avgW_mean_sem, [0.20 0.55 0.65], 'R = \SigmaW/N');
    plot_metric_with_sem(ax3, betas, BetaSummary.lambda1_mean_avg, BetaSummary.lambda1_mean_sem, [0.80 0.35 0.25], '\lambda_1');
    plot_metric_with_sem(ax3, betas, BetaSummary.eta_mean_avg, BetaSummary.eta_mean_sem, [0.20 0.25 0.30], '\eta');
    xlabel(ax3, 'Epileptiform drive \beta');
    ylabel(ax3, 'Normalized to \beta = 0');
    title(ax3, 'Decomposition of \eta = R/\lambda_1');
    legend(ax3, 'Location', 'best', 'Box', 'off');
    set(ax3, 'Box', 'off', 'TickDir', 'out');

    ax4 = nexttile(4);
    scatter(ax4, TrialTable.pair_density_mean, TrialTable.eta_mean, 48, TrialTable.beta, 'filled', ...
        'MarkerEdgeColor', [0.2 0.2 0.2]);
    cb = colorbar(ax4);
    cb.Label.String = '\beta';
    xlabel(ax4, 'Pair-edge density');
    ylabel(ax4, '\eta');
    title(ax4, 'Trial-level density check');
    set(ax4, 'Box', 'off', 'TickDir', 'out');

    exportgraphics(fig, fullfile(outDir, 'HR_eta_signal_simulation_summary.png'), 'Resolution', 600);
    exportgraphics(fig, fullfile(outDir, 'HR_eta_signal_simulation_summary.pdf'), 'ContentType', 'vector');
    savefig(fig, fullfile(outDir, 'HR_eta_signal_simulation_summary.fig'));
end

function plot_metric_with_sem(ax, x, y, sem, colorVal, labelVal)
    y0 = y(1);
    if ~isfinite(y0) || abs(y0) <= eps
        yNorm = y;
        semNorm = sem;
    else
        yNorm = y / y0;
        semNorm = sem / abs(y0);
    end
    errorbar(ax, x, yNorm, semNorm, '-o', 'Color', colorVal, ...
        'MarkerFaceColor', colorVal, 'MarkerEdgeColor', 'w', ...
        'LineWidth', 1.7, 'MarkerSize', 5, 'DisplayName', labelVal);
end
