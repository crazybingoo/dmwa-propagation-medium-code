clc; clear; close all;

%% 1. 数据准备 (保持不变)
resultDir = 'E:\wcldematlab\keep\new_idea\8 - n_u_v\2_Fig1'; 
inFile    = fullfile(resultDir, 'ALL_CASES_stage5_eta.csv');

T = readtable(inFile);
if iscell(T.case_id), T.case_id = string(T.case_id); end
if iscell(T.phase5),  T.phase5  = string(T.phase5);  end

phaseOrder  = ["pre-ictal","early","mid","late","post-ictal"];
phaseLabels = {'pre-ictal','early','mid','late','post-ictal'};
caseIDs = unique(T.case_id, 'stable');

M = nan(numel(caseIDs), numel(phaseOrder));
for i = 1:numel(caseIDs)
    rows_i = T(T.case_id == caseIDs(i), :);
    for p = 1:numel(phaseOrder)
        idx = (rows_i.phase5 == phaseOrder(p));
        tmp = rows_i.eta_mean(idx);
        if ~isempty(tmp), M(i,p) = tmp(1); end
    end
end

%% 2. 统计计算 (计算每一组相对于 pre-ictal 的差异)
pvals = nan(numel(phaseOrder), 1);
dvals = nan(numel(phaseOrder), 1);
for p = 2:numel(phaseOrder)
    x1 = M(:, 1); x2 = M(:, p);
    valid = ~isnan(x1) & ~isnan(x2);
    if sum(valid) >= 2
        pvals(p) = signrank(x1(valid), x2(valid));
        diffV = x2(valid) - x1(valid);
        dvals(p) = mean(diffV) / std(diffV);
    end
end

%% 3. 绘图参数 (优化垂直偏移以防重叠)
phaseColors = [76,120,168; 89,161,79; 242,142,43; 225,87,89; 128,115,172] / 255;
cloudH  = 0.45;   % 密度图高度系数
boxH    = 0.08;   % 箱子宽度（变窄一点更精致）
cloudGap = 0.05;  % 云与基准线的间隙
rainGap  = 0.25;  % 雨与基准线的间隙

fig = figure('Color','w', 'Position', [100, 100, 800, 750]);
hold on;

% 垂直基准线
refVal = median(M(:,1), 'omitnan');
xline(refVal, '--', 'Color', [0.2 0.4 0.7], 'LineWidth', 1.5, 'Alpha', 0.8);

for p = 1:numel(phaseOrder)
    data = M(:,p); data = data(~isnan(data));
    if isempty(data), continue; end
    color = phaseColors(p,:);
    
    % --- (A) 绘制 "云"：向上移动 cloudGap ---
    [f, xi] = ksdensity(data);
    f = f / max(f) * cloudH;
    % p - cloudGap 是起点，再向上减去 f
    patch(xi, (p - cloudGap) - f, color, 'FaceAlpha', 0.6, 'EdgeColor', color, 'LineWidth', 1);
    
    % --- (B) 绘制标准箱线图：正好压在基准线 p 上 ---
    q1 = prctile(data, 25); q3 = prctile(data, 75); med = median(data);
    whisL = min(data(data >= q1 - 1.5*(q3-q1)));
    whisR = max(data(data <= q3 + 1.5*(q3-q1)));
    
    % 须线
    line([whisL whisR], [p p], 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    % 箱体 (实体填充)
    rectangle('Position', [q1, p-boxH/2, q3-q1, boxH], ...
        'FaceColor', color*0.9, 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1);
    % 中位数竖线 (白色更醒目)
    line([med med], [p-boxH/2 p+boxH/2], 'Color', 'w', 'LineWidth', 2);
    
    % --- (C) 绘制 "雨"：向下移动 rainGap ---
    jitter = (rand(numel(data),1) - 0.5) * 0.15;
    scatter(data, p + rainGap + jitter, 30, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4);
    
    % --- (D) 显著性星号：放在云的顶部边缘 ---
    if p > 1 && ~isnan(pvals(p))
        stars = get_stars(pvals(p));
        % 放在该组分布最大值附近，垂直位置在云的上方
        text(max(xi), p - cloudGap - 0.1, stars, 'FontSize', 14, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end
%% 4. 界面润色 (关键格式调整)
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.5, ...
    'Box', 'on', 'TickDir', 'out', ...
    'YTick', 1:numel(phaseOrder), 'YTickLabel', phaseLabels, ...
    'YDir', 'reverse', ...
    'YAxisLocation', 'right', 'FontWeight', 'bold'); % 重点：坐标标签移到右侧

xlabel('\eta (mean value)', 'FontWeight', 'bold', 'FontSize', 14);

% 坐标范围微调
Xall = M(:); Xall = Xall(~isnan(Xall));
xlim([min(Xall)-0.02, 1.03]); 
ylim([0.5, numel(phaseOrder) + 0.6]);

function stars = get_stars(p)
    if p < 0.001, stars = '***';
    elseif p < 0.01, stars = '**';
    elseif p < 0.05, stars = '*';
    else, stars = 'ns';
    end
end