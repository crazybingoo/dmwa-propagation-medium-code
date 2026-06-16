function plot_robustness_heatmap_eta(allSummary, windowSecList, retainedFractionList, ...
    defaultWindowSec, defaultRetainedFraction, outDir)

nW = numel(windowSecList);
nQ = numel(retainedFractionList);

DZ    = nan(nW, nQ);
PVAL  = nan(nW, nQ);
TREND = nan(nW, nQ);

PRE   = nan(nW, nQ);
EARLY = nan(nW, nQ);
MID   = nan(nW, nQ);
LATE  = nan(nW, nQ);
POST  = nan(nW, nQ);

for i = 1:numel(allSummary)
    w = allSummary(i).windowSec;
    q = allSummary(i).retainedFraction;

    iw = find(abs(windowSecList - w) < 1e-12, 1);
    iq = find(abs(retainedFractionList - q) < 1e-12, 1);

    DZ(iw, iq)    = allSummary(i).dz_pre_vs_ictal;
    PVAL(iw, iq)  = allSummary(i).p_pre_vs_ictal;
    TREND(iw, iq) = allSummary(i).trendKeepRate;

    PRE(iw, iq)   = allSummary(i).preMean;
    EARLY(iw, iq) = allSummary(i).earlyMean;
    MID(iw, iq)   = allSummary(i).midMean;
    LATE(iw, iq)  = allSummary(i).lateMean;
    POST(iw, iq)  = allSummary(i).postMean;
end

iw0 = find(abs(windowSecList - defaultWindowSec) < 1e-12, 1);
iq0 = find(abs(retainedFractionList - defaultRetainedFraction) < 1e-12, 1);

baseDZ = DZ(iw0, iq0);

if isempty(baseDZ) || isnan(baseDZ) || abs(baseDZ) < eps
    RET = DZ;
else
    RET = DZ / baseDZ;
end

repTrend = [PRE(iw0,iq0), EARLY(iw0,iq0), MID(iw0,iq0), LATE(iw0,iq0), POST(iw0,iq0)];

fontName   = 'Arial';
fontSizeA  = 9;
fontSizeB  = 8;
markerSize = 5;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18.4 8.2]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t, 1);

imagesc(RET, 'Parent', ax1);
hold(ax1, 'on');

cmap = parula(256);
colormap(ax1, cmap);
caxis(ax1, [0.4 1.1]);

set(ax1, 'XTick', 1:nQ, 'XTickLabel', compose('%.2f', retainedFractionList), ...
    'YTick', 1:nW, 'YTickLabel', compose('%.0f', windowSecList), ...
    'FontName', fontName, 'FontSize', fontSizeA, 'LineWidth', 0.8, ...
    'TickDir', 'out', 'Box', 'off');

xlabel(ax1, 'Retained PLV edge fraction', 'FontName', fontName, 'FontSize', fontSizeA);
ylabel(ax1, 'Window length (s)', 'FontName', fontName, 'FontSize', fontSizeA);
title(ax1, 'Robustness of \eta across parameter settings', ...
    'FontName', fontName, 'FontSize', fontSizeA + 1, 'FontWeight', 'normal');

xlim(ax1, [0.5, nQ + 0.5]);
ylim(ax1, [0.5, nW + 0.5]);

for x = 0.5 : 1 : nQ + 0.5
    plot(ax1, [x x], [0.5 nW+0.5], '-', 'Color', [1 1 1]*0.92, 'LineWidth', 0.6);
end
for y = 0.5 : 1 : nW + 0.5
    plot(ax1, [0.5 nQ+0.5], [y y], '-', 'Color', [1 1 1]*0.92, 'LineWidth', 0.6);
end

for i = 1:nW
    for j = 1:nQ
        if ~isnan(DZ(i,j))

            v = RET(i,j);
            if isnan(v)
                txtColor = [1 1 1];
            elseif v > 0.72
                txtColor = [1 1 1];
            else
                txtColor = [0.08 0.08 0.08];
            end

            txt1 = sprintf('dz = %.2f', DZ(i,j));
            txt2 = sprintf('%d%%', round(TREND(i,j)*100));

            if isnan(PVAL(i,j))
                txt3 = 'p = NA';
            elseif PVAL(i,j) < 0.001
                txt3 = 'p < 0.001';
            else
                txt3 = sprintf('p = %.3f', PVAL(i,j));
            end

            text(j, i-0.20, txt1, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontName', fontName, 'FontSize', 6.8, ...
                'FontWeight', 'bold', ...
                'Color', txtColor, 'Parent', ax1);

            text(j, i+0.02, txt2, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontName', fontName, 'FontSize', 6.8, ...
                'Color', txtColor, 'Parent', ax1);

            text(j, i+0.24, txt3, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontName', fontName, 'FontSize', 6.5, ...
                'Color', txtColor, 'Parent', ax1);
        end
    end
end

rectangle(ax1, 'Position', [iq0-0.5, iw0-0.5, 1, 1], ...
    'EdgeColor', [0 0 0], 'LineWidth', 1.8, 'LineStyle', '-');

cb = colorbar(ax1);
cb.Label.String = 'Effect retention (relative to default)';
cb.Label.FontName = fontName;
cb.Label.FontSize = fontSizeA;
cb.FontName = fontName;
cb.FontSize = fontSizeB;
cb.Box = 'off';

ax2 = nexttile(t, 2);
hold(ax2, 'on');

x = 1:5;
plot(ax2, x, repTrend, '-o', ...
    'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.5, ...
    'MarkerSize', markerSize, ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0.10 0.10 0.10]);

patch(ax2, [1.5 4.5 4.5 1.5], ...
    [min(repTrend)-0.02 min(repTrend)-0.02 max(repTrend)+0.02 max(repTrend)+0.02], ...
    [0.93 0.93 0.93], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

plot(ax2, x, repTrend, '-o', ...
    'Color', [0.10 0.10 0.10], ...
    'LineWidth', 1.5, ...
    'MarkerSize', markerSize, ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0.10 0.10 0.10]);

set(ax2, 'XLim', [0.7 5.3], ...
    'XTick', 1:5, ...
    'XTickLabel', {'pre','early','mid','late','post'}, ...
    'FontName', fontName, 'FontSize', fontSizeA, ...
    'LineWidth', 0.8, 'TickDir', 'out', 'Box', 'off');

ylabel(ax2, 'Mean \eta', 'FontName', fontName, 'FontSize', fontSizeA);
title(ax2, sprintf('Representative stage profile (%.0f s, retained = %.2f)', ...
    defaultWindowSec, defaultRetainedFraction), ...
    'FontName', fontName, 'FontSize', fontSizeA + 1, 'FontWeight', 'normal');

ymin = min(repTrend) - 0.03;
ymax = max(repTrend) + 0.03;
ylim(ax2, [ymin ymax]);

text(ax2, 3.0, ymax - 0.01, 'ictal', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontName', fontName, 'FontSize', fontSizeB, ...
    'Color', [0.25 0.25 0.25]);

if isnan(PVAL(iw0, iq0))
    pTxt0 = 'p = NA';
elseif PVAL(iw0, iq0) < 0.001
    pTxt0 = 'p < 0.001';
else
    pTxt0 = sprintf('p = %.3f', PVAL(iw0, iq0));
end

txtInfo = sprintf('default: %.0f s, retained = %.2f\ndz = %.2f\ntrend = %d%%\n%s', ...
    defaultWindowSec, defaultRetainedFraction, DZ(iw0, iq0), round(TREND(iw0, iq0)*100), pTxt0);

text(ax2, 0.98, 0.05, txtInfo, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'FontName', fontName, 'FontSize', 7.2, ...
    'BackgroundColor', 'w', ...
    'Margin', 4);

set(fig, 'Renderer', 'painters');

exportgraphics(fig, fullfile(outDir, 'Figure_Robustness_eta_NatureStyle.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'Figure_Robustness_eta_NatureStyle.pdf'), 'ContentType', 'vector');

disp('Figure saved:');
disp(fullfile(outDir, 'Figure_Robustness_eta_NatureStyle.png'));
disp(fullfile(outDir, 'Figure_Robustness_eta_NatureStyle.pdf'));

end
