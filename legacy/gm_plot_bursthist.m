function n = gm_plot_bursthist(ax, R, xmax, showY)
% GM_PLOT_BURSTHIST  Burst size distribution of the granules (red bars) and of their
% burst-free copies (black outline), as in MS Figs. 6-9. R from GelMetricsReport.m.
% Writes the median burst sizes and the turnover beyond noise in the corner.
% Returns the number of granule bursts.
if nargin < 3, xmax = 60; end
if nargin < 4, showY = true; end
r = 100 * vertcat(R.items.r0); c = 100 * vertcat(R.items.c0);
edges = 0:2.5:xmax;
histogram(ax, min(r, xmax - 1e-6), edges, 'Normalization', 'pdf', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
hold(ax, 'on');
histogram(ax, min(c, xmax - 1e-6), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
xlim(ax, [0 xmax]); xlabel(ax, 'Burst size (% of granule)');
if showY, ylabel(ax, 'Frequency'); end
text(ax, 0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
    R.sizeGranules, R.sizeCopies, R.turnover, R.turnover_lo, R.turnover_hi), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 6.5);
gm_style(ax);
n = numel(r);
end
