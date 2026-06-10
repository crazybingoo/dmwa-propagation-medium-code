%% Probe HR parameters before running the full DMW-HLG eta pipeline.
% This file only checks whether simulated HR signals show increasing
% phase-locking with beta. It does not replace the paper-2 DMW-HLG pipeline.

clc; clear; close all;
rng(260519);

base = struct();
base.nNodes = 18;
base.coreNodes = 1:4;
base.relayNodes = 5:9;
base.peripheryNodes = 10:18;
base.betaLevels = [0 0.25 0.50 0.75 1.00];
base.numTrials = 4;
base.analysisFs = 256;
base.durationSec = 12;
base.burnSec = 4;
base.dt = 0.02;
base.a = 1; base.b = 3; base.c = 1; base.d = 5;
base.r = 0.006; base.s = 4; base.xR = -1.6;
base.Ibase = 3.05;
base.commonFreq = 0.18;

cfg = [
    0.000 0.500 0.030 0.080 0.55
    0.000 0.900 0.030 0.080 0.55
    0.010 1.200 0.025 0.050 0.70
    0.020 1.800 0.020 0.030 0.85
    ];
% cols: coupling_base, coupling_slope, noise_amp, common_slope, heterogeneity_scale

for ic = 1:size(cfg, 1)
    params = base;
    params.couplingBase = cfg(ic, 1);
    params.couplingSlope = cfg(ic, 2);
    params.noiseAmp = cfg(ic, 3);
    params.commonSlope = cfg(ic, 4);
    params.heterogeneityScale = cfg(ic, 5);

    fprintf('\nConfig %d | g %.3f+%.3f beta | noise %.3f | common %.3f beta | hetero %.2f\n', ...
        ic, params.couplingBase, params.couplingSlope, params.noiseAmp, ...
        params.commonSlope, params.heterogeneityScale);
    for ib = 1:numel(params.betaLevels)
        beta = params.betaLevels(ib);
        vals = nan(params.numTrials, 4);
        for trial = 1:params.numTrials
            spec = make_trial_spec_probe(params, trial);
            X = simulate_hr_probe(params, spec, beta);
            PLV = compute_plv_probe(X);
            upper = PLV(triu(true(params.nNodes), 1));
            vals(trial, 1) = mean(upper, 'omitnan');
            vals(trial, 2) = median(upper, 'omitnan');
            vals(trial, 3) = mean(std(diff(X, 1, 2), 0, 2), 'omitnan');
            vals(trial, 4) = mean(sum(X > 1.5, 1) >= 4);
        end
        fprintf('  beta %.2f | meanPLV %.3f +/- %.3f | medPLV %.3f | diffstd %.4f | syncBurst %.3f\n', ...
            beta, mean(vals(:,1)), std(vals(:,1)), mean(vals(:,2)), ...
            mean(vals(:,3)), mean(vals(:,4)));
    end
end

function spec = make_trial_spec_probe(params, trial)
    rng(1000 + trial * 31);
    spec.A = build_coupling_probe(params);
    spec.x0 = -1.5 + 0.60 * randn(params.nNodes, 1);
    spec.y0 = -7.0 + 0.60 * randn(params.nNodes, 1);
    spec.z0 =  3.0 + 0.35 * randn(params.nNodes, 1);
    spec.Ioffset = params.heterogeneityScale * randn(params.nNodes, 1);
    totalSteps = round((params.durationSec + params.burnSec) * params.analysisFs);
    spec.noise = randn(params.nNodes, totalSteps);
    spec.commonPhase = 2 * pi * rand();
end

function A = build_coupling_probe(params)
    N = params.nNodes;
    A = zeros(N);
    c = params.coreNodes; r = params.relayNodes; p = params.peripheryNodes;
    A(c,c)=1; A(r,r)=0.35; A(p,p)=0.10; A(c,r)=0.85; A(r,c)=0.65;
    A(r,p)=0.35; A(p,r)=0.22; A(c,p)=0.12; A(p,c)=0.08;
    A(1:N+1:end)=0;
    A = max(A, A');
    A = A .* (0.65 + 0.70 * rand(N));
    A = max(A, A');
    A(1:N+1:end)=0;
    rowSum = sum(A, 2); rowSum(rowSum == 0) = 1;
    A = A ./ rowSum;
end

function X = simulate_hr_probe(params, spec, beta)
    totalSteps = round((params.durationSec + params.burnSec) * params.analysisFs);
    burnSteps = round(params.burnSec * params.analysisFs);
    outSteps = totalSteps - burnSteps;
    x = spec.x0; y = spec.y0; z = spec.z0;
    X = nan(params.nNodes, outSteps);
    rowSum = sum(spec.A, 2);

    coreWeight = 0.30 * ones(params.nNodes, 1);
    coreWeight(params.coreNodes) = 1.00;
    coreWeight(params.relayNodes) = 0.70;
    coreWeight(params.peripheryNodes) = 0.42;

    excit = zeros(params.nNodes, 1);
    excit(params.coreNodes) = 0.38;
    excit(params.relayNodes) = 0.20;
    excit(params.peripheryNodes) = 0.06;

    g = params.couplingBase + params.couplingSlope * beta;
    commonAmp = 0.015 + params.commonSlope * beta;
    out = 0;
    for step = 1:totalSteps
        t = (step - 1) * params.dt;
        commonDrive = commonAmp * sin(2 * pi * params.commonFreq * t + spec.commonPhase);
        hetero = (1 - 0.90 * beta) * spec.Ioffset;
        inputCurrent = params.Ibase + hetero + beta * excit + commonDrive * coreWeight;
        inputCurrent = inputCurrent + params.noiseAmp * spec.noise(:, step);
        coupling = g * (spec.A * x - rowSum .* x);

        dx = y - params.a*x.^3 + params.b*x.^2 - z + inputCurrent + coupling;
        dy = params.c - params.d*x.^2 - y;
        dz = params.r * (params.s*(x - params.xR) - z);

        x = x + params.dt * dx;
        y = y + params.dt * dy;
        z = z + params.dt * dz;
        if step > burnSteps
            out = out + 1;
            X(:, out) = x;
        end
    end
    X = zscore_probe(X);
end

function Xz = zscore_probe(X)
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd <= eps | ~isfinite(sd)) = 1;
    Xz = (X - mu) ./ sd;
    Xz(~isfinite(Xz)) = 0;
end

function PLV = compute_plv_probe(X)
    ph = angle(hilbert(X.')).';
    N = size(X, 1);
    PLV = eye(N);
    for i = 1:N
        for j = i+1:N
            v = abs(mean(exp(1i * (ph(i,:) - ph(j,:)))));
            PLV(i,j) = v;
            PLV(j,i) = v;
        end
    end
end
