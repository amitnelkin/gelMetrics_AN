% Local copy (lab PC): the two helpers added to the figure composer (compose_ms_figs.m, not in the repo)
% for Figs. 6(d-f), 9(d, e), 10(b-e) with one dot per granule and no error bars (RESULTS section 81).
% dotPanel: grey bar at the condition value + one dot per granule (in vitro: coloured by movie, side by side).
% drawBrackets: the p brackets of compareConditions.m, now clearing the dots, also on log axes.

function drawBrackets(ax, file, measure, br)
% the p brackets of compareConditions.m drawn again from its own results (file: the saved
% results), stacked so that no two brackets or p values overlap: solid for a two-sided
% p < 0.05, dashed and marked one-sided when only the one-sided p < 0.05 (as the script decides).
% br (optional): the brackets to draw instead, rows [a b p oneSided]. The brackets clear the
% error bars, the bars and the granule dots; on a log axis the steps are equal in log space.
if nargin < 4
    V = load(file, 'results', 'conditions', 'pairs');
    br = zeros(0, 4);
    for q = 1:size(V.pairs, 1)
        r = V.results{q}; i = strcmp(r.measure, measure);
        p2 = r.('2-side-p')(i); p1 = r.('1-side-p')(i);
        a = find(strcmp(V.conditions, V.pairs{q, 1})); b = find(strcmp(V.conditions, V.pairs{q, 2}));
        if p2 < 0.05, br(end + 1, :) = [min(a, b), max(a, b), p2, 0]; %#ok<AGROW>
        elseif p1 < 0.05, br(end + 1, :) = [min(a, b), max(a, b), p1, 1]; %#ok<AGROW>
        end
    end
end
br(:, 5) = 0;                                  % a, b, p, one-sided?, level
delete(findobj(ax, 'Type', 'line', '-function', @(h) numel(h.XData) == 4));
t = findobj(ax, 'Type', 'text');
for k = 1:numel(t)
    if ~isempty(regexp(char(t(k).String), '^p [=<]', 'once')), delete(t(k)); end
end
isLog = strcmp(ax.YScale, 'log');
ex = []; etop = []; elo = [];
e = findobj(ax, 'Type', 'errorbar');
for k = 1:numel(e)
    ex = [ex, e(k).XData]; etop = [etop, e(k).YData + e(k).YPositiveDelta]; elo = [elo, e(k).YData - e(k).YNegativeDelta]; %#ok<AGROW>
end
d = [findobj(ax, 'Type', 'scatter'); findobj(ax, 'Type', 'bar')];
for k = 1:numel(d)
    ex = [ex, d(k).XData(:)']; etop = [etop, d(k).YData(:)']; elo = [elo, d(k).YData(:)']; %#ok<AGROW>
end
if isLog
    ok = etop > 0; ex = ex(ok); etop = log10(etop(ok)); elo = log10(elo(ok));
    step = 0.07 * (max(etop) - min(elo));
else
    step = 0.07 * (max([etop 0]) - min([elo 0]));
end
[~, o] = sort(br(:, 2) - br(:, 1)); br = br(o, :);
for k = 1:size(br, 1)
    a = br(k, 1); b = br(k, 2);
    base = max(etop(ex >= a - 0.5 & ex <= b + 0.5));
    for j = 1:k - 1                               % above every bracket already drawn over these bars
        if br(j, 1) <= b && br(j, 2) >= a, base = max(base, br(j, 5) + 1.6 * step); end
    end
    br(k, 5) = base + step;
    lvl = br(k, 5);
    style = '-'; if br(k, 4), style = '--'; end
    yy = [lvl - step / 3, lvl, lvl, lvl - step / 3]; ty = lvl;
    if isLog, yy = 10 .^ yy; ty = 10 ^ ty; end
    plot(ax, [a a b b], yy, style, 'Color', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    if br(k, 3) < 0.001, s = 'p < 0.001'; else, s = sprintf('p = %.3f', br(k, 3)); end
    if br(k, 4), s = [s ' (one-sided)']; end %#ok<AGROW>
    text(ax, (a + b) / 2, ty, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', 'FontSize', 6.5, 'Interpreter', 'none');
end
top = max([etop, br(:, 5)' + 1.6 * step]);
if isLog
    ylim(ax, 10 .^ [min(elo) - 0.5 * step, top + 0.3 * step]);
else
    yl = ylim(ax);
    lo = min([elo 0]) - 0.5 * step;
    if isempty(findobj(ax, 'Type', 'scatter')), lo = min(yl(1), lo); end   % with dots: just below the lowest dot
    ylim(ax, [lo, top + 0.3 * step]);
end
end

function dotPanel(ax, bars, vals, cols, groups)
% the bars at the condition values (grey, as before) and one dot per granule on top; no error bars.
% vals{i}: the granule values of condition i (NaN: no dot); cols{i}: one colour per dot (rows), or
% cols = {one colour} for all; groups{i} (optional): the movie of every dot, drawn side by side
delete(findobj(ax, 'Type', 'errorbar')); delete(findobj(ax, 'Type', 'bar')); delete(findobj(ax, 'Type', 'scatter'));
hold(ax, 'on');
n = numel(bars); w = 0.6;
base = 0;
if strcmp(ax.YScale, 'log')
    v = vertcat(vals{:}); base = 10 ^ floor(log10(min(v(v > 0))));
end
bar(ax, 1:n, bars, w, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'BaseValue', base, 'ShowBaseLine', 'off');
rng(1);                                        % the same jitter every time
for i = 1:n
    v = vals{i}(:); ok = ~isnan(v);
    if strcmp(ax.YScale, 'log'), ok = ok & v > 0; end
    if isempty(groups)
        x = i + (rand(size(v)) - 0.5) * 0.8 * w;
    else
        g = groups{i}(:); u = unique(g(ok)); slot = 0.9 * w / numel(u); x = zeros(size(v));
        for j = 1:numel(u)
            mm = g == u(j);
            x(mm) = i - 0.45 * w + (j - 0.5) * slot + (rand(nnz(mm), 1) - 0.5) * 0.6 * slot;
        end
    end
    c = cols{min(i, numel(cols))};
    if size(c, 1) > 1, c = c(ok, :); end
    scatter(ax, x(ok), v(ok), 5, c, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
end
xlim(ax, [0.4, n + 0.6]); ax.XTick = 1:n;
uistack(findobj(ax, 'Type', 'constantline'), 'top');
end
