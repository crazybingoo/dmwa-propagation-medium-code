%% 原生逻辑仿真 2.0：基于过度同步(Hypersynchrony)机制的 DMW-HLG 验证（eta版本）
% 目标：考察真实发作期(Ictal)对应的“高阶网络折射率” eta 随同步强度 beta 的变化
clc; clear; close all;

% 开始全局计时
tic_total = tic;
fprintf('==================================================\n');
fprintf('开始基于【过度同步(Hypersynchrony)】机制的大规模仿真（eta版本）\n');
fprintf('==================================================\n\n');

%% 0. 设定输出路径并创建文件夹
outDir = 'example_project\2_Fig9';
if ~exist(outDir, 'dir')
    mkdir(outDir);
    fprintf('已创建输出目录: %s\n\n', outDir);
else
    fprintf('输出目录已存在: %s\n\n', outDir);
end

%% 1. 基础参数
N_nodes = 100;           % 底层节点数量
num_E1 = 150;            % 1阶超边数量 (Size=2)
num_E2 = 100;            % 2阶超边数量 (Size=3)

% beta (sync_level) 模拟癫痫发作的“过度同步化”强度
% beta = 0: 正常态(Pre-ictal)，全脑稀疏连接，异质性高，有方向偏置
% beta = 1: 发作态(Ictal)，网络高度致密同质化，度数趋同
sync_level_range = 0:0.05:1;
num_lambdas = length(sync_level_range);
num_trials = 100;        % 每个参数下的蒙特卡洛次数

% 预分配存储空间
eta_all  = zeros(num_trials, num_lambdas);
eta_mean = zeros(num_lambdas, 1);
eta_std  = zeros(num_lambdas, 1);

% 定义“同步核心池”：发作时，超边将极度密集地在这个局部池中爆发
sync_core_size = 25; % 假设全脑25%的节点陷入核心共振
sync_nodes = 1:sync_core_size;

%% 2. 开始主循环
for idx = 1:num_lambdas
    beta = sync_level_range(idx);

    % 初始化当前 beta 下各步骤的耗时统计
    time_StepA = 0;
    time_StepB = 0;
    time_StepC = 0;

    for trial = 1:num_trials

        %% [步骤A: 根据同步强度生成超边拓扑]
        t_A = tic;
        HE = cell(num_E1 + num_E2, 1);

        for k = 1:(num_E1 + num_E2)
            if k <= num_E1
                edge_size = 2;
            else
                edge_size = 3;
            end

            % 核心发病机制：过度致密同步
            if rand() < beta
                % 同步态：超边被强制完全在紧凑的"同步池"内部生成
                tmp_nodes = randsample(sync_nodes, edge_size);
            else
                % 正常态：超边在全脑范围内稀疏随机生成
                tmp_nodes = randsample(N_nodes, edge_size);
            end

            % 强制转化为行向量，避免维度拼接报错
            HE{k} = tmp_nodes(:)';
        end
        time_StepA = time_StepA + toc(t_A);

        %% [步骤B: 重叠度 O、Coverage 与 Degree Bias 计算]
        t_B = tic;
        num_HE = length(HE);
        sizes = cellfun(@length, HE);
        sizes = sizes(:);

        % 计算重叠矩阵 O
        O = zeros(num_HE, num_HE);
        for i = 1:num_HE
            for j = i+1:num_HE
                overlap_len = length(intersect(HE{i}, HE{j}));
                if overlap_len > 0
                    O(i, j) = overlap_len;
                    O(j, i) = overlap_len;
                end
            end
        end

        % 1. 计算 Coverage
        Sizes_i = repmat(sizes, 1, num_HE);
        Sizes_j = repmat(sizes', num_HE, 1);
        Coverage = O ./ Sizes_i;

        % 2. 计算 Degree Bias
        D = sum(O > 0, 2);
        D_i = repmat(D, 1, num_HE);
        D_j = repmat(D', num_HE, 1);

        DegreeSum = D_i + D_j;
        DegreeBias = D_j ./ DegreeSum;
        DegreeBias(DegreeSum == 0) = 0.5;

        % 3. 生成 W 矩阵
        W = Coverage;
        EqualSizeMask = (Sizes_i == Sizes_j) & (O > 0);

        % 同阶应用 DegreeBias 产生不对称
        W(EqualSizeMask) = W(EqualSizeMask) .* DegreeBias(EqualSizeMask);

        % 清理自身环
        W(1:num_HE+1:end) = 0;
        time_StepB = time_StepB + toc(t_B);

        %% [步骤C: 计算 eta 值]
        t_C = tic;
        W_double = double(W);
        W_double(~isfinite(W_double)) = 0;
        W_double(W_double < 0) = 0;
        W_double(1:num_HE+1:end) = 0;

        if num_HE < 2
            eta_all(trial, idx) = 0;
        else
            sum_W = sum(W_double(:));

            if sum_W <= 0
                eta_all(trial, idx) = 0;
            else
                lambda_1 = max(abs(eig(W_double)));

                if ~isfinite(lambda_1) || lambda_1 <= eps
                    eta_all(trial, idx) = 0;
                else
                    eta_all(trial, idx) = (sum_W / num_HE) / lambda_1;
                end
            end
        end
        time_StepC = time_StepC + toc(t_C);

    end % 结束当前 beta 下的所有 trials

    % 记录统计量
    eta_mean(idx) = mean(eta_all(:, idx));
    eta_std(idx)  = std(eta_all(:, idx));

    % 打印该组 (beta) 的各项用时
    total_idx_time = time_StepA + time_StepB + time_StepC;
    fprintf('[进度 %2d/%2d] Beta = %.2f | 拓扑生成: %.3fs | W矩阵计算: %.3fs | eta值计算: %.3fs | 该组总计: %.3fs\n', ...
        idx, num_lambdas, beta, time_StepA, time_StepB, time_StepC, total_idx_time);
end

%% 3. 可视化绘图 (Nature Style)
fig = figure('Color', 'w', 'Position', [200, 200, 600, 480]);
hold on;

% 绘制误差带 (半透明)
patch([sync_level_range, fliplr(sync_level_range)], ...
      [eta_mean' + eta_std', fliplr(eta_mean' - eta_std')], ...
      [0.2 0.6 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.25);

% 绘制均值主线
plot(sync_level_range, eta_mean, '-o', 'LineWidth', 2.5, ...
    'Color', [0.15 0.45 0.75], 'MarkerFaceColor', 'w', 'MarkerSize', 7);

% 图形修饰
xlabel('\beta', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('\eta', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');

set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2, ...
         'Box', 'off', 'TickDir', 'out');

% 动态调整Y轴，保证视觉美观
y_min = min(eta_mean - eta_std);
y_max = max(eta_mean + eta_std);
y_pad = (y_max - y_min) * 0.1;
if y_pad == 0
    y_pad = max(abs(y_max), 0.01) * 0.1;
end
ylim([y_min - y_pad, y_max + y_pad]);
grid on;

% 添加辅助解释文本
text(0.05, y_min, 'Pre-ictal (Sparse & Directional)', 'Color', [0.4 0.4 0.4], 'FontSize', 10, 'FontWeight', 'bold');
text(0.70, y_max - y_pad*0.5, 'Ictal (Dense & Symmetric)', 'Color', [0.8 0.2 0.2], 'FontSize', 10, 'FontWeight', 'bold');

%% 4. 自动保存所有结果到指定路径
fprintf('\n正在保存结果至: %s ...\n', outDir);

% 保存图片
exportgraphics(fig, fullfile(outDir, 'Validation_Hypersynchrony_eta_vs_Beta.png'), 'Resolution', 600);
savefig(fig, fullfile(outDir, 'Validation_Hypersynchrony_eta_vs_Beta.fig'));

% 保存核心数据矩阵
save(fullfile(outDir, 'Validation_Hypersynchrony_eta_Results.mat'), ...
    'sync_level_range', 'eta_all', 'eta_mean', 'eta_std', ...
    'N_nodes', 'num_E1', 'num_E2', 'num_trials', 'sync_core_size');

%% 5. 结束全局计时
total_time = toc(tic_total);
fprintf('\n==================================================\n');
fprintf('全部仿真及文件保存已完成！\n');
fprintf('仿真总耗时: %.2f 秒 (约 %.2f 分钟)\n', total_time, total_time/60);
fprintf('==================================================\n');
