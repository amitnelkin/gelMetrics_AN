function gm_plot_bars(ax, Rs, field, ylab, colors, x, tickLabels)
% GM_PLOT_BARS  Bars with 95% confidence intervals for one measure of several
% conditions, as in MS Figs. 6(d-f) and 10(b-d).
%   Rs     cell array of GelMetricsReport results
%   field  'meanBurst', 'turnover', 'net', 'sizeBeyond', 'medianBurst', ...
%   colors n x 3 matrix (one row per bar) or one row for all
%   x, tickLabels  (optional) bar positions and labels (default 1..n, R.name)
n = numel(Rs);
if nargin < 6 || isempty(x), x = 1:n; end
if nargin < 7 || isempty(tickLabels), tickLabels = cellfun(@(R) R.name, Rs, 'UniformOutput', false); end
if size(colors, 1) == 1, colors = repmat(colors, n, 1); end
y = cellfun(@(R) R.(field), Rs); lo = cellfun(@(R) R.([field '_lo']), Rs); hi = cellfun(@(R) R.([field '_hi']), Rs);
hold(ax, 'on');
for k = 1:n
    bar(ax, x(k), y(k), 0.6, 'FaceColor', colors(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
end
errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
yline(ax, 0, 'k', 'LineWidth', 0.5);
set(ax, 'XTick', x, 'XTickLabel', tickLabels, 'TickLabelInterpreter', 'tex');
xlim(ax, [min(x) - 0.6, max(x) + 0.6]); ylabel(ax, ylab);
gm_style(ax);
end
