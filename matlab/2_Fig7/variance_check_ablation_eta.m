clc; clear; close all;

workDir = 'example_project\2_Fig7';
csvFile = fullfile(workDir, 'ALL_CASES_stage5_eta_ablation.csv');

T = readtable(csvFile);

phaseOrder = {'pre-ictal','early','mid','late','post-ictal'};

VarTable = table('Size',[5 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Phase','Var_Orig','Var_Cov','Var_Deg'});

for p = 1:5
    idx = strcmp(T.phase5, phaseOrder{p});
    x1 = T.eta_orig_mean(idx);
    x2 = T.eta_cov_mean(idx);
    x3 = T.eta_deg_mean(idx);

    VarTable.Phase(p)    = string(phaseOrder{p});
    VarTable.Var_Orig(p) = var(x1, 'omitnan');
    VarTable.Var_Cov(p)  = var(x2, 'omitnan');
    VarTable.Var_Deg(p)  = var(x3, 'omitnan');
end

allVar_orig = var(T.eta_orig_mean, 'omitnan');
allVar_cov  = var(T.eta_cov_mean,  'omitnan');
allVar_deg  = var(T.eta_deg_mean,  'omitnan');

fig = figure('Color','w','Position',[100 100 860 420]);

c3 = [0, 160, 135]/255;
c2 = [60,  84, 136]/255;
c1 = [230,  75,  53]/255;
colorSet = [c1; c2; c3];

ax1 = axes('Position',[0.10 0.18 0.68 0.66]);
Y1 = [VarTable.Var_Orig, VarTable.Var_Cov, VarTable.Var_Deg];
b1 = bar(ax1, Y1, 'grouped', 'BarWidth', 0.78);
hold(ax1, 'on');

for i = 1:3
    b1(i).FaceColor = colorSet(i,:);
    b1(i).EdgeColor = [0.15 0.15 0.15];
    b1(i).LineWidth = 0.8;
end

set(ax1, ...
    'FontName','Arial', ...
    'FontSize',10.5, ...
    'LineWidth',1.0, ...
    'Box','off', ...
    'TickDir','out', ...
    'TickLength',[0.015 0.015], ...
    'Layer','top');

xticks(ax1, 1:5);
xticklabels(ax1, {'Pre','Early','Mid','Late','Post'});
ylabel(ax1, 'Variance', 'FontWeight','bold');
grid(ax1, 'off');
ax1.XRuler.Axle.LineStyle = 'solid';
ax1.YRuler.Axle.LineStyle = 'solid';

legend(ax1, {'Original','Coverage only','DegreeBias only'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontSize',10);

ax2 = axes('Position',[0.57 0.50 0.18 0.22]);
Y2 = [allVar_orig, allVar_cov, allVar_deg];
x2 = 1:3;
hold(ax2, 'on');

plot(ax2, x2, Y2, '-o', ...
    'Color', [0.20 0.20 0.20], ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5.5, ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0.15 0.15 0.15]);

for i = 1:3
    scatter(ax2, x2(i), Y2(i), 34, ...
        'MarkerFaceColor', colorSet(i,:), ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.8);
end

set(ax2, ...
    'FontName','Arial', ...
    'FontSize',8.5, ...
    'LineWidth',0.9, ...
    'Box','off', ...
    'TickDir','out', ...
    'TickLength',[0.020 0.020], ...
    'Layer','top');

xticks(ax2, 1:3);
xticklabels(ax2, {'Ori','Cov','Deg'});
title(ax2, 'Overall', 'FontSize',9.5, 'FontWeight','normal');
ylabel(ax2, 'Var', 'FontSize',8.5);
xlim(ax2, [0.7 3.3]);
grid(ax2, 'off');

pngFile = fullfile(workDir, 'variance_check_eta_ablation.png');
pdfFile = fullfile(workDir, 'variance_check_eta_ablation.pdf');
exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

writetable(VarTable, fullfile(workDir, 'variance_table_eta_ablation.csv'));
