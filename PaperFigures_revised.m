%% PaperFigures_revised.m - Figures 5-9 of the revised manuscript
% Run one section at a time (Ctrl+Enter) or the whole script. Each section is
% self-contained: it starts from the repository data, runs the same functions a
% user runs on their own data (GelMetricsReport.m, GelMetricsCompare.m) and draws
% the figure in the manuscript style. Figures are saved in figures/ (PNG, 300 dpi, and PDF).
%
% The noise copies use the seeds of the paper (one copy per granule), so the
% numbers are those of the manuscript. Confidence intervals and permutation p values
% are recomputed and can differ in the last digit. Set nCopies = 20 for values that
% do not depend on the random draw of the copies (RESULTS.md section 45).
%
% Run time: Figure 5 about 10 min (simulations); Figures 6-9 about 1-2 min each.

%% Settings
clearvars; close all; clc;
repo = fileparts(mfilename('fullpath'));
if isempty(repo) || ~exist(fullfile(repo, 'analysis'), 'dir'), repo = pwd; end   % run from the repository folder
addpath(repo, fullfile(repo, 'analysis'), fullfile(repo, 'legacy'));
outDir = fullfile(repo, 'figures'); if ~exist(outDir, 'dir'), mkdir(outDir); end
makePlots = true;          % false: compute and print the numbers only
nCopies = 1;               % 1 = as in the paper; 20 = stable noise reference
RED = [1 0.4 0.4]; BLUE = [0 0.447 0.741]; ORANGE = [0.835 0.369 0];
CPP = [0 0.447 0.698]; CFUS = [0 0.620 0.451]; CQB = [0.835 0.369 0];
report = @(name) GelMetricsReport(gm_paper_dataset(name), struct('nCopies', nCopies, ...
    'seed', gm_paper_dataset(name).paperSeed, 'verbose', false));

%% Figure 5 - estimating the signal parameters on simulated data
% (a-c) burst amplitudes of simulated bursty granules (lambda 1, 2, 4; k0 = 60; noise 40),
% (d) lambda estimated by the weighted score (WS) and by maximum likelihood vs the truth,
% (e) mean burst size lambda*k0 (shape fit) vs the truth,
% (f) share of true bursts found and share of false detections vs the stringency parameter,
%     at the original simulation noise and at the PP7 noise level.
tVec = (0:360) * 10/60; lams = [0.5 1 2 3 4 5];
% (a-d): MS section III simulations, 1000 granules, SP 15, five repeats per lambda
amps = cell(6, 5); lamWS = zeros(6, 5); lamNLL = zeros(6, 5);
for li = 1:6
    for rep = 1:5
        rng(900 + (li - 1) * 5 + rep);
        [S, B] = generateSim(1000, 361, 3, 0.0026, 40, lams(li), 60);
        ws = warning('off', 'all'); [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, 15, 0, 0); warning(ws);
        amps{li, rep} = abs(seg(seg(:, 2) < 0, 5));
        [lamWS(li, rep), ~, lamNLL(li, rep)] = gm_eq19_fit(amps{li, rep});
    end
    fprintf('Fig 5(d) lambda %g: WS %s | likelihood %s\n', lams(li), mat2str(round(lamWS(li, :), 2)), mat2str(round(lamNLL(li, :), 2)));
end
% (e): 200 granules, SP 12, shape fit on all bursts, 95% CI from resampling granules
shp = zeros(6, 3);
for li = 1:6
    rng(700 + (li - 1) * 3 + 2);
    [S, B] = generateSim(200, 361, 3, 0.0026, 40, lams(li), 60);
    ws = warning('off', 'all'); [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, 12, 0, 0); warning(ws);
    neg = seg(seg(:, 2) < 0, :);
    rng(27); [l, k, lo, hi] = gm_shape_fit(abs(neg(:, 5)), neg(:, 1), 200);
    shp(li, :) = [l * k, lo(3), hi(3)];
    fprintf('Fig 5(e) true lambda*k0 %3d: estimated %.0f [%.0f, %.0f]\n', 60 * lams(li), shp(li, :));
end
% (f): bursty granules (lambda 2) and burst-free twins, 500 each, at two noise levels
SPs = [10 12 14 16 18 20 22 25]; sets = {'original noise', 40, 60; 'PP7 noise', 1100, 300};
recall = zeros(2, numel(SPs)); falseShare = recall;
for j = 1:2
    nz = sets{j, 2}; k0 = sets{j, 3}; seed = 3400 + j;
    rng(seed); [S, B] = generateSim(500, 361, 3, 0.0026, nz, 2, k0);
    rng(seed); S0 = generateSim(500, 361, 3, 0.0026, 0, 2, k0);
    rng(seed + 50); [F, BF] = generateSim(500, 361, 1, 0.0026, nz, 2, k0);
    u = S0 ./ exp(-0.0026 * linspace(0, 60, 361)); du = diff(u, 1, 2);
    [ti, tj] = find(du < -1e-9); tj = tj + 1;                     % granule and frame of every true loss
    for q = 1:numel(SPs)
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, SPs(q), 0, 0);
        [~, ~, ~, ~, ~, ~, segF] = analysis(num2cell(F, 2), num2cell(BF, 2), tVec, 13, SPs(q), 0, 0);
        warning(ws);
        neg = seg(seg(:, 2) < 0, :); nF = sum(segF(:, 2) < 0);
        found = false(numel(ti), 1);
        for b = 1:size(neg, 1)
            found(ti == neg(b, 1) & tj >= neg(b, 3) - 6 & tj <= neg(b, 4) + 6) = true;
        end
        recall(j, q) = 100 * mean(found); falseShare(j, q) = 100 * min(1, nF / max(size(neg, 1), 1));
    end
    fprintf('Fig 5(f) %s: bursts found %s %%, false share %s %%\n', sets{j, 1}, mat2str(round(recall(j, :))), mat2str(round(falseShare(j, :))));
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 3.9], 'Color', 'w');
    t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:3
        ax = nexttile(t); show = [1 2 4]; li = find(lams == show(j));
        histogram(ax, amps{li, 1}, 0:10:1000, 'Normalization', 'pdf', 'FaceColor', RED, 'EdgeColor', 'none', 'FaceAlpha', 1);
        xlim(ax, [0 1000]); ylim(ax, [0 0.012]); xlabel(ax, 'Amplitude (a.u.)');
        if j == 1, ylabel(ax, 'Frequency'); end
        text(ax, 0.95, 0.85, sprintf('\\lambda_{sim}=%d', lams(li)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 9);
        gm_style(ax, char('a' + j - 1));
    end
    ax = nexttile(t); hold(ax, 'on');
    m = median(lamWS, 2); errorbar(ax, lams, m, m - min(lamWS, [], 2), max(lamWS, [], 2) - m, 'o', 'Color', BLUE, 'MarkerSize', 4, 'CapSize', 3);
    m = median(lamNLL, 2); errorbar(ax, lams, m, m - min(lamNLL, [], 2), max(lamNLL, [], 2) - m, 's', 'Color', ORANGE, 'MarkerSize', 4, 'CapSize', 3);
    plot(ax, [0 5.5], [0 5.5], 'r--', 'LineWidth', 1);
    xlim(ax, [0 5.5]); ylim(ax, [0 5.5]); xlabel(ax, '\lambda_{sim}'); ylabel(ax, '\lambda_{est}');
    legend(ax, {'WS score', 'likelihood', 'truth'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 6.5); gm_style(ax, 'd');
    ax = nexttile(t); hold(ax, 'on');
    errorbar(ax, 60 * lams, shp(:, 1), shp(:, 1) - shp(:, 2), shp(:, 3) - shp(:, 1), 'o', 'Color', BLUE, 'MarkerSize', 4, 'CapSize', 3);
    plot(ax, [0 330], [0 330], 'r--', 'LineWidth', 1); xlim(ax, [0 330]); ylim(ax, [0 330]);
    xlabel(ax, 'true \lambda\cdotk_0 (a.u.)'); ylabel(ax, 'estimated \lambda\cdotk_0 (a.u.)'); gm_style(ax, 'e');
    ax = nexttile(t); hold(ax, 'on');
    h1 = plot(ax, SPs, recall(1, :), '-o', 'Color', BLUE, 'MarkerSize', 4); plot(ax, SPs, falseShare(1, :), ':', 'Color', BLUE, 'LineWidth', 1);
    h2 = plot(ax, SPs, recall(2, :), '-s', 'Color', ORANGE, 'MarkerSize', 4); plot(ax, SPs, falseShare(2, :), ':', 'Color', ORANGE, 'LineWidth', 1);
    h3 = plot(ax, NaN, NaN, 'k-o', 'MarkerSize', 3); h4 = plot(ax, NaN, NaN, 'k:', 'LineWidth', 1);
    legend(ax, [h1 h2 h3 h4], {'original noise', 'PP7 noise', 'bursts found', 'false share'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 5.5);
    xlabel(ax, 'Stringency parameter j'); ylabel(ax, '%'); ylim(ax, [0 105]); gm_style(ax, 'f');
    saveFig(f, outDir, 'Fig5');
end

%% Figure 6 - in vitro granules
V = {'vitro_1_1', 'vitro_10_1', 'vitro_100_1'}; ratio = {'1:1', '10:1', '100:1'};
RV = cellfun(report, V, 'UniformOutput', false);
for k = 1:3
    R = RV{k};
    fprintf('Fig 6 %s: median %.1f vs %.1f %%, mean burst %.0f [%.0f, %.0f], turnover %.1f [%.1f, %.1f], net %.1f [%.1f, %.1f]\n', ratio{k}, ...
        R.sizeGranules, R.sizeCopies, R.meanBurst, R.meanBurst_lo, R.meanBurst_hi, R.turnover, R.turnover_lo, R.turnover_hi, R.net, R.net_lo, R.net_hi);
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 3.7], 'Color', 'w');
    t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:3
        ax = nexttile(t); n = gm_plot_bursthist(ax, RV{k}, 60, k == 1);
        title(ax, sprintf('protein:RNA %s  (n = %d)', ratio{k}, n), 'FontWeight', 'normal', 'FontSize', 7.5); gm_style(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6); end
    end
    yl = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
    fl = {'meanBurst', 'turnover', 'net'};
    for k = 1:3
        ax = nexttile(t); gm_plot_bars(ax, RV, fl{k}, yl{k}, RED, 1:3, ratio); xlabel(ax, 'protein:RNA'); gm_style(ax, char('d' + k - 1));
    end
    saveFig(f, outDir, 'Fig6');
end

%% Figure 7 - PP7 granules in vivo
R7 = {report('4PP7'), report('30PP7')}; lab7 = {'PP7 + QCP-5x/PCP-4x', 'PP7 + PCP-30x'};
for k = 1:2, fprintf('Fig 7 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', lab7{k}, R7{k}.sizeGranules, R7{k}.sizeCopies, R7{k}.turnover, R7{k}.turnover_lo, R7{k}.turnover_hi); end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        ax = nexttile(t); n = gm_plot_bursthist(ax, R7{k}, 60, k == 1);
        title(ax, sprintf('%s  (n = %d)', lab7{k}, n), 'FontWeight', 'normal', 'FontSize', 7.5); gm_style(ax, char('a' + k - 1));
    end
    legend(nexttile(t, 1), {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5);
    saveFig(f, outDir, 'Fig7');
end

%% Figure 8 - Qbeta granules in vivo
R8 = {report('5Qb'), report('10Qb')}; lab8 = {'Q\beta + QCP-5x/PCP-4x', 'Q\beta + QCP-10x'};
for k = 1:2, fprintf('Fig 8 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', lab8{k}, R8{k}.sizeGranules, R8{k}.sizeCopies, R8{k}.turnover, R8{k}.turnover_lo, R8{k}.turnover_hi); end
T = GelMetricsCompare(R8{1}, R7{1}); T = GelMetricsCompare(R8{2}, R7{2});   % Qbeta vs PP7 (text of Sec. V.B)
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        ax = nexttile(t); n = gm_plot_bursthist(ax, R8{k}, 60, k == 1);
        title(ax, sprintf('%s  (n = %d)', lab8{k}, n), 'FontWeight', 'normal', 'FontSize', 7.5); gm_style(ax, char('a' + k - 1));
    end
    legend(nexttile(t, 1), {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5);
    saveFig(f, outDir, 'Fig8');
end

%% Figure 9 - PP7-FUS granules in vivo
R9 = {report('4PP7-FUS'), report('30PP7-FUS')}; lab9 = {'PP7-FUS + QCP-5x/PCP-4x', 'PP7-FUS + PCP-30x'};
P = {report('4PP7'), report('30PP7')};
T4 = GelMetricsCompare(R9{1}, P{1}); T30 = GelMetricsCompare(R9{2}, P{2});
pv = @(T, m, c) T.(c)(strcmp(T.measure, m));
% (d) PP7-FUS minus PP7: total gap in median burst size, the part reproduced by the copies, the part beyond noise
gap = zeros(2, 3); gapLo = gap; gapHi = gap;
for i = 1:2
    A = R9{i}; Bp = P{i}; fld = {'sizeGranules', 'sizeCopies', 'sizeBeyond'};
    for k = 1:3
        d = A.boot(:, strcmp(A.bootKeys, fld{k})) - Bp.boot(:, strcmp(Bp.bootKeys, fld{k}));
        gap(i, k) = A.(fld{k}) - Bp.(fld{k}); q = prctile(d, [2.5 97.5]); gapLo(i, k) = q(1); gapHi(i, k) = q(2);
    end
    fprintf('Fig 9(d) %s: total %.1f [%.1f, %.1f], noise %.1f [%.1f, %.1f], beyond %.1f [%.1f, %.1f]\n', lab9{i}, [gap(i, :); gapLo(i, :); gapHi(i, :)]);
end
fprintf('Fig 9(c) brightness p = %.4f (4x), %.2f (30x); 9(d) beyond noise 4x: p = %.3f one-sided, %.3f two-sided; 30x: p = %.2f\n', ...
    pv(T4, 'brightnessPooled', 'p'), pv(T30, 'brightnessPooled', 'p'), pv(T4, 'sizeBeyond', 'p_AgtB'), pv(T4, 'sizeBeyond', 'p'), pv(T30, 'sizeBeyond', 'p'));
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 4.4], 'Color', 'w');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        ax = nexttile(t); n = gm_plot_bursthist(ax, R9{k}, 60, k == 1);
        title(ax, sprintf('%s  (n = %d)', lab9{k}, n), 'FontWeight', 'normal', 'FontSize', 7.5); gm_style(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    % (c) granule brightness (median over granules, 95% CI)
    ax = nexttile(t); Rs = {P{1}, R9{1}, P{2}, R9{2}}; x = [1 2 3.6 4.6];
    y = cellfun(@(R) R.brightnessPooled, Rs) / 1000; lo = cellfun(@(R) R.brightnessPooled_lo, Rs) / 1000; hi = cellfun(@(R) R.brightnessPooled_hi, Rs) / 1000;
    hold(ax, 'on'); cl = [CPP; CFUS; CPP; CFUS];
    for k = 1:4, bar(ax, x(k), y(k), 0.7, 'FaceColor', cl(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85); end
    errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    set(ax, 'XTick', x, 'XTickLabel', {'PP7', 'PP7-FUS', 'PP7', 'PP7-FUS'}); ylabel(ax, 'Granule brightness (a.u./1000)');
    yt = max(hi) * 1.08; ylim(ax, [0 max(hi) * 1.2]);
    text(ax, 1.5, yt, sprintf('p = %.3f', pv(T4, 'brightnessPooled', 'p')), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 7);
    text(ax, 4.1, yt, sprintf('p = %.2f', pv(T30, 'brightnessPooled', 'p')), 'HorizontalAlignment', 'center', 'FontSize', 7);
    xlabel(ax, 'QCP-5x/PCP-4x                    PCP-30x'); gm_style(ax, 'c');
    % (d) burst size gap
    ax = nexttile(t); hold(ax, 'on'); cl = [0.55 0.55 0.55; 0.85 0.85 0.85; CFUS]; hb = gobjects(1, 3);
    for i = 1:2
        for k = 1:3
            xx = (i - 1) * 1.1 + (k - 2) * 0.28;
            hb(k) = bar(ax, xx, gap(i, k), 0.26, 'FaceColor', cl(k, :), 'EdgeColor', 'k', 'LineWidth', 0.4);
            errorbar(ax, xx, gap(i, k), gap(i, k) - gapLo(i, k), gapHi(i, k) - gap(i, k), 'k', 'LineWidth', 0.8, 'CapSize', 3);
        end
    end
    text(ax, 0, 13.5, sprintf('p = %.3f (one-sided)\np = %.3f (two-sided)', pv(T4, 'sizeBeyond', 'p_AgtB'), pv(T4, 'sizeBeyond', 'p')), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 6.5);
    text(ax, 1.1, 13.5, sprintf('p = %.2f', pv(T30, 'sizeBeyond', 'p')), 'HorizontalAlignment', 'center', 'FontSize', 6.5);
    yline(ax, 0, 'k', 'LineWidth', 0.5); ylim(ax, [-4 17]);
    set(ax, 'XTick', [0 1.1], 'XTickLabel', {'QCP-5x/PCP-4x', 'PCP-30x'});
    ylabel(ax, sprintf('PP7-FUS - PP7, median burst\nsize (%% points)'));
    legend(ax, hb, {'total', 'explained by noise', 'beyond noise'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 6);
    gm_style(ax, 'd');
    saveFig(f, outDir, 'Fig9');
end


function saveFig(f, outDir, name)
% saves a figure as PNG (300 dpi) and PDF
exportgraphics(f, fullfile(outDir, [name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile(outDir, [name '.pdf']), 'ContentType', 'vector');
end
