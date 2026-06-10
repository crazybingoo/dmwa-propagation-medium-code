clc; clear; close all;

workDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig7';
csvFile = fullfile(workDir, 'ALL_CASES_stage5_eta_ablation.csv');

if ~exist(csvFile, 'file')
    error('找不到文件: %s', csvFile);
end

T = readtable(csvFile);

phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};
metricNames  = {'eta_orig_mean', 'eta_cov_mean', 'eta_deg_mean'};

phaseColors = [76,120,168;
               89,161,79;
               242,142,43;
               225,87,89;
               128,115,172] / 255;

cases = unique(T.case_id, 'stable');
nCase = numel(cases);

blockStarts = [1 6.5 12];
blockCenters = blockStarts + 2;

allPositions = [];
groupCounter = 0;
for b = 1:3
    for p = 1:5
        groupCounter = groupCounter + 1;
        allPositions(groupCounter) = blockStarts(b) + (p-1); %#ok<SAGROW>
    end
end

dataVec = [];
groupVec = [];
globalMin = inf;
globalMax = -inf;

for b = 1:3
    for p = 1:5
        y = T.(metricNames{b})(strcmp(T.phase5, phaseOrder{p}));
        y = y(~isnan(y));

        if ~isempty(y)
            globalMin = min(globalMin, min(y));
            globalMax = max(globalMax, max(y));
        end

        groupID = (b-1)*5 + p;
        dataVec = [dataVec; y(:)]; %#ok<AGROW>
        groupVec = [groupVec; repmat(groupID, numel(y), 1)]; %#ok<AGROW>
    end
end

if ~isfinite(globalMin) || ~isfinite(globalMax)
    error('没有可用数据。');
end

Xmethod = cell(1,3);
for m = 1:3
    X = nan(nCase, 5);
    for i = 1:nCase
        idxCase = strcmp(T.case_id, cases{i});
        for p = 1:5
            idx = idxCase & strcmp(T.phase5, phaseOrder{p});
            vals = T.(metricNames{m})(idx);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                X(i,p) = vals(1);
            end
        end
    end
    X = X(all(~isnan(X),2), :);
    Xmethod{m} = X;
end

fig = figure('Color','w','Position',[40 80 1000 600]);
ax = axes('Position',[0.05 0.14 0.93 0.76]);
hold(ax, 'on');

yr = globalMax - globalMin;
if yr <= 0
    yr = 0.05;
end

ylim(ax, [globalMin - 0.08*yr, globalMax + 0.66*yr]);
xlim(ax, [0 15]);
yl = ylim(ax);

boxplot(ax, dataVec, groupVec, ...
    'Positions', allPositions, ...
    'Symbol', '', ...
    'Widths', 0.55, ...
    'Colors', [0.18 0.18 0.18], ...
    'MedianStyle', 'line');

set(findobj(ax, 'Tag', 'Box'), ...
    'Color', [0.20 0.20 0.20], 'LineWidth', 1.1);

set(findobj(ax, 'Tag', 'Median'), ...
    'Color', [0 0 0], 'LineWidth', 1.3);

set(findobj(ax, 'Tag', 'Whisker'), ...
    'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);

set(findobj(ax, 'Tag', 'Adjacent Value'), ...
    'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);

rng(2025);
for b = 1:3
    for p = 1:5
        x0 = blockStarts(b) + (p-1);
        y = T.(metricNames{b})(strcmp(T.phase5, phaseOrder{p}));
        y = y(~isnan(y));
        if isempty(y)
            continue;
        end

        jitter = (rand(numel(y),1) - 0.5) * 0.32;

        scatter(ax, x0 + jitter, y, 22, ...
            'MarkerFaceColor', phaseColors(p,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.40, ...
            'MarkerFaceAlpha', 0.92, ...
            'MarkerEdgeAlpha', 0.95);
    end
end

set(ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', [0.012 0.012], ...
    'Layer', 'top');

xticks(ax, blockCenters);
xticklabels(ax, {'Original', 'Coverage only', 'DegreeBias only'});
ylabel(ax, '\eta', 'FontWeight', 'bold');

hPre = scatter(ax, nan, nan, 38, 'MarkerFaceColor', phaseColors(1,:), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
hEarly = scatter(ax, nan, nan, 38, 'MarkerFaceColor', phaseColors(2,:), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
hMid = scatter(ax, nan, nan, 38, 'MarkerFaceColor', phaseColors(3,:), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
hLate = scatter(ax, nan, nan, 38, 'MarkerFaceColor', phaseColors(4,:), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
hPost = scatter(ax, nan, nan, 38, 'MarkerFaceColor', phaseColors(5,:), 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);

legend(ax, [hPre hEarly hMid hLate hPost], ...
    {'Pre', 'Early', 'Mid', 'Late', 'Post'}, ...
    'Location', 'northoutside', ...
    'Orientation', 'horizontal', ...
    'Box', 'on', ...
    'FontSize', 10);

plot(ax, [0.7 5.3],  [yl(1)+0.002*yr yl(1)+0.002*yr], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
plot(ax, [7.7 12.3], [yl(1)+0.002*yr yl(1)+0.002*yr], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
plot(ax, [14.7 19.3],[yl(1)+0.002*yr yl(1)+0.002*yr], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

pairIdx_adj = [1 2; 2 3; 3 4; 4 5];
pairIdx_key = [1 3; 3 5];

baseY = globalMax + 0.05*yr;
stepY = 0.075*yr;

for b = 1:3
    X = Xmethod{b};
    xpos = blockStarts(b) + (0:4);

    for k = 1:size(pairIdx_adj,1)
        i1 = pairIdx_adj(k,1);
        i2 = pairIdx_adj(k,2);
        [pval, dval] = paired_stats_pd(X(:,i1), X(:,i2));
        label = sprintf('%s\n%s', format_p_value(pval), format_d_value(dval));
        yLine = baseY + (k-1)*stepY;
        add_sig_line_with_label(ax, xpos(i1), xpos(i2), yLine, label, 0.9, 8.0);
    end

    for k = 1:size(pairIdx_key,1)
        i1 = pairIdx_key(k,1);
        i2 = pairIdx_key(k,2);
        [pval, dval] = paired_stats_pd(X(:,i1), X(:,i2));
        label = sprintf('%s\n%s', format_p_value(pval), format_d_value(dval));
        yLine = baseY + (size(pairIdx_adj,1)+k-1)*stepY;
        add_sig_line_with_label(ax, xpos(i1), xpos(i2), yLine, label, 1.1, 8.2);
    end
end

pngFile = fullfile(workDir, 'GROUP_stage5_eta_ablation_oneaxis_3blocks.png');
pdfFile = fullfile(workDir, 'GROUP_stage5_eta_ablation_oneaxis_3blocks.pdf');

exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

fprintf('已保存:\n%s\n%s\n', pngFile, pdfFile);

function [pval, dval] = paired_stats_pd(x1, x2)
    x1 = x1(:);
    x2 = x2(:);
    mask = ~(isnan(x1) | isnan(x2));
    x1 = x1(mask);
    x2 = x2(mask);

    if numel(x1) < 2
        pval = NaN;
        dval = NaN;
        return;
    end

    pval = signrank(x1, x2);
    d = x2 - x1;
    sd_d = std(d, 'omitnan');

    if sd_d == 0 || isnan(sd_d)
        dval = 0;
    else
        dval = mean(d, 'omitnan') / sd_d;
    end
end

function add_sig_line_with_label(ax, x1, x2, y, label, lw, fs)
    yr = diff(ylim(ax));
    h = 0.016 * yr;

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], ...
        'Color', [0.18 0.18 0.18], ...
        'LineWidth', lw);

    text(ax, mean([x1 x2]), y + 1.10*h, label, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', ...
        'FontSize', fs);
end

function txt = format_p_value(p)
    if isnan(p)
        txt = 'p = NaN';
    elseif p < 0.001
        txt = 'p < 0.001';
    else
        txt = sprintf('p = %.3f', p);
    end
end

function txt = format_d_value(d)
    if isnan(d)
        txt = 'd = NaN';
    else
        txt = sprintf('d = %.2f', d);
    end
end
