function plot_exp4_summary_figure_size_adjusted(Tregion, TmacroRole, phaseOrder, phaseShort, regionNames, macroOrder, roleOrder, outDir)
    % 一张图：A 区域承担（size-adjusted）；B 宏观状态组合承担；C source/sink/balanced 承担
    fig = figure('Color','w', 'Position',[100 100 1280 380]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    regionColors = [213 94 0; 0 114 178; 0 158 115] / 255;
    roleColors = [204 121 167; 80 80 80; 86 180 233] / 255;
    macroColors = [0.10 0.10 0.10; 0.35 0.35 0.35; 0.60 0.60 0.60; ...
                   0.80 0.35 0.20; 0.20 0.45 0.75; 0.35 0.60 0.30; 0.50 0.30 0.70];

    %% A: SOZ/PZ/NIZ eta_abs_burden_per_node（扣除区域大小差异）
    nexttile; hold on;
    for r = 1:numel(regionNames)
        [y,e] = group_mean_sem(Tregion, 'region', regionNames{r}, 'eta_abs_burden_per_node_mean', phaseOrder);
        errorbar(1:numel(phaseOrder), y, e, '-o', ...
            'Color', regionColors(r,:), 'MarkerFaceColor', regionColors(r,:), ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.8, 'CapSize', 7, 'MarkerSize', 5.5);
    end
    format_phase_axis(phaseShort);
    ylabel('eta burden per node', 'FontWeight','bold');
    title('Regional eta burden (size-adjusted)', 'FontWeight','bold');
    legend(regionNames, 'Box','off', 'Location','best');

    %% B: 宏观状态组合 eta_abs_burden，合并 role 后统计
    nexttile; hold on;
    for m = 1:numel(macroOrder)
        [y,e] = group_mean_sem_macro_all_roles(TmacroRole, macroOrder{m}, 'eta_abs_burden_mean', phaseOrder, roleOrder);
        errorbar(1:numel(phaseOrder), y, e, '-o', ...
            'Color', macroColors(m,:), 'MarkerFaceColor', macroColors(m,:), ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.5, 'CapSize', 6, 'MarkerSize', 5.0);
    end
    format_phase_axis(phaseShort);
    ylabel('eta burden, abs proxy', 'FontWeight','bold');
    title('Macro-state combination', 'FontWeight','bold');
    legend(macroOrder, 'Box','off', 'Location','eastoutside', 'FontSize',8);

    %% C: source-like / balanced / sink-like eta_abs_burden，合并 macro 后统计
    nexttile; hold on;
    for rr = 1:numel(roleOrder)
        [y,e] = group_mean_sem_role_all_macro(TmacroRole, roleOrder{rr}, 'eta_abs_burden_mean', phaseOrder, macroOrder);
        errorbar(1:numel(phaseOrder), y, e, '-o', ...
            'Color', roleColors(rr,:), 'MarkerFaceColor', roleColors(rr,:), ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.8, 'CapSize', 7, 'MarkerSize', 5.5);
    end
    format_phase_axis(phaseShort);
    ylabel('eta burden, abs proxy', 'FontWeight','bold');
    title('Source/sink role', 'FontWeight','bold');
    legend(roleOrder, 'Box','off', 'Location','best');

%     title(tl, 'Experiment 4: regional and hyperedge-level contributors of eta (size-adjusted)', ...
%         'FontName','Arial', 'FontSize',13, 'FontWeight','bold');

    exportgraphics(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted.png'), 'Resolution', 600);
    exportgraphics(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted.pdf'), 'ContentType','vector');
    savefig(fig, fullfile(outDir, 'FIG_exp4_eta_burden_summary_size_adjusted.fig'));

end

function [y,e] = group_mean_sem(T, groupVar, groupValue, metricName, phaseOrder)
    y = nan(1,numel(phaseOrder));
    e = nan(1,numel(phaseOrder));
    for p = 1:numel(phaseOrder)
        idx = strcmp(T.(groupVar), groupValue) & strcmp(T.phase5, phaseOrder{p});
        x = T.(metricName)(idx);
        x = x(isfinite(x));
        y(p) = mean(x, 'omitnan');
        e(p) = std(x, 'omitnan') / sqrt(max(numel(x),1));
    end
end

function [y,e] = group_mean_sem_macro_all_roles(T, macroState, metricName, phaseOrder, roleOrder)
    caseIDs = unique(T.case_id, 'stable');
    Y = nan(numel(caseIDs), numel(phaseOrder));
    for i = 1:numel(caseIDs)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for r = 1:numel(roleOrder)
                idx = strcmp(T.case_id, caseIDs{i}) & strcmp(T.macro_state, macroState) & ...
                      strcmp(T.role, roleOrder{r}) & strcmp(T.phase5, phaseOrder{p});
                x = T.(metricName)(idx);
                if ~isempty(x) && isfinite(x(1))
                    val = val + x(1);
                    hasAny = true;
                end
            end
            if hasAny
                Y(i,p) = val;
            end
        end
    end
    y = mean(Y, 1, 'omitnan');
    e = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
end

function [y,e] = group_mean_sem_role_all_macro(T, roleName, metricName, phaseOrder, macroOrder)
    caseIDs = unique(T.case_id, 'stable');
    Y = nan(numel(caseIDs), numel(phaseOrder));
    for i = 1:numel(caseIDs)
        for p = 1:numel(phaseOrder)
            val = 0;
            hasAny = false;
            for m = 1:numel(macroOrder)
                idx = strcmp(T.case_id, caseIDs{i}) & strcmp(T.role, roleName) & ...
                      strcmp(T.macro_state, macroOrder{m}) & strcmp(T.phase5, phaseOrder{p});
                x = T.(metricName)(idx);
                if ~isempty(x) && isfinite(x(1))
                    val = val + x(1);
                    hasAny = true;
                end
            end
            if hasAny
                Y(i,p) = val;
            end
        end
    end
    y = mean(Y, 1, 'omitnan');
    e = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(Y), 1));
end

function format_phase_axis(phaseShort)
    xlim([0.75 numel(phaseShort)+0.25]);
    xticks(1:numel(phaseShort));
    xticklabels(phaseShort);
    xtickangle(25);
    set(gca, 'FontName','Arial', 'FontSize',10.5, 'LineWidth',1.15, ...
        'Box','off', 'TickDir','out', 'Layer','top');
    ax = gca;
    ax.YAxis.Exponent = 0;
end

