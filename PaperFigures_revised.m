%% PaperFigures_revised.m - Figures 5-9 of the revised manuscript
% Run one section at a time (Ctrl+Enter) or the whole script, from the repository
% folder. Each section starts from the repository data, runs the same code a user
% runs on their own data (burstMeasures.m, compareConditions.m) and draws the figure
% in the manuscript style. Figures are saved in figures/ (PNG, 300 dpi, and PDF).
%
% The noise copies use the seeds of the paper (one copy per granule), so the numbers
% are those of the manuscript. Set nCopies = 20 for values that do not depend on the
% random draw of the copies (RESULTS.md section 45).
% In vitro uses the spot sum, as in vivo (the user's decision, RESULTS.md section 53).
%
% Run time: Figure 5 about 8 min (simulations); Figures 6-9 about 1 min each.

%% Settings
clearvars; close all; clc;
addpath('legacy');                % the original lambda/k0 code (only for reference)
if ~exist('figures', 'dir'), mkdir('figures'); end
makePlots = true;                 % false: compute and print the numbers only
nCopies = 1;                      % 1 = as in the paper; 20 = stable noise reference
% the random seeds of the copies used for the paper
seeds = containers.Map({'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'}, ...
                       {18, 19, 20, 21, 22, 23, 18, 19, 20});
RED = [1 0.4 0.4]; BLUE = [0 0.447 0.741]; ORANGE = [0.835 0.369 0];
CPP = [0 0.447 0.698]; CFUS = [0 0.620 0.451];

%% Figure 5 - estimating the signal parameters on simulated data
% (a-c) burst amplitudes of simulated bursty granules (lambda 1, 2, 4; k0 = 60; noise 40)
% (d)   lambda estimated by the weighted score (WS) and by maximum likelihood, vs the truth
% (e)   mean burst size lambda*k0 (shape fit), vs the truth
% (f)   share of the true bursts found, and share of false detections, vs the stringency
%       parameter, at the original simulation noise and at the PP7 noise level
tVec = (0:360) * 10/60;
lams = [0.5 1 2 3 4 5];

% (a-d): simulations as in MS section III: 1000 granules, SP 15, five repeats per lambda
amps = cell(6, 5); lamWS = zeros(6, 5); lamNLL = zeros(6, 5);
for li = 1:6
    for rep = 1:5
        rng(900 + (li - 1) * 5 + rep);
        [S, B] = generateSim(1000, 361, 3, 0.0026, 40, lams(li), 60);
        oldWarnings = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, 15, 0, 0);
        warning(oldWarnings);
        amps{li, rep} = abs(seg(seg(:, 2) < 0, 5));                % negative bursts
        [lamWS(li, rep), lamNLL(li, rep)] = eq19Fit(amps{li, rep});
    end
    fprintf('Fig 5(d) lambda %g: WS %s | likelihood %s\n', lams(li), mat2str(round(lamWS(li, :), 2)), mat2str(round(lamNLL(li, :), 2)));
end

% (e): 200 granules, SP 12, shape fit on all bursts, 95% CI from resampling granules
lamk0 = zeros(6, 3);                                             % estimate, CI low, CI high
for li = 1:6
    rng(700 + (li - 1) * 3 + 2);
    [S, B] = generateSim(200, 361, 3, 0.0026, 40, lams(li), 60);
    oldWarnings = warning('off', 'all');
    [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, 12, 0, 0);
    warning(oldWarnings);
    neg = seg(seg(:, 2) < 0, :);
    a = abs(neg(:, 5)); granule = neg(:, 1);
    [lam, k0] = shapeFit(a);
    rng(27); boot = zeros(200, 1);
    withBursts = unique(granule); nG = numel(withBursts);       % granules that have bursts
    for b = 1:200                                                % resample those granules
        chosen = withBursts(randi(nG, nG, 1)); aBoot = [];
        for g = 1:nG, aBoot = [aBoot; a(granule == chosen(g))]; end %#ok<AGROW>
        [lb, kb] = shapeFit(aBoot); boot(b) = lb * kb;
    end
    lamk0(li, :) = [lam * k0, prctile(boot, 2.5), prctile(boot, 97.5)];
    fprintf('Fig 5(e) true lambda*k0 %3d: estimated %.0f [%.0f, %.0f]\n', 60 * lams(li), lamk0(li, :));
end

% (f): bursty granules (lambda 2) and burst-free twins, 500 each, at two noise levels
SPs = [10 12 14 16 18 20 22 25];
noiseLevels = [40 1100]; unitSizes = [60 300];                  % original noise; PP7-like noise
recall = zeros(2, numel(SPs)); falseShare = zeros(2, numel(SPs));
for j = 1:2
    seed = 3400 + j;
    rng(seed); [S, B] = generateSim(500, 361, 3, 0.0026, noiseLevels(j), 2, unitSizes(j));
    rng(seed); S0 = generateSim(500, 361, 3, 0.0026, 0, 2, unitSizes(j));          % the same granules without noise
    rng(seed + 50); [F, BF] = generateSim(500, 361, 1, 0.0026, noiseLevels(j), 2, unitSizes(j));   % burst-free twins
    % every true loss of material: granule and frame
    u = S0 ./ exp(-0.0026 * linspace(0, 60, 361)); du = diff(u, 1, 2);
    [lossGranule, lossFrame] = find(du < -1e-9); lossFrame = lossFrame + 1;
    for q = 1:numel(SPs)
        oldWarnings = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, 13, SPs(q), 0, 0);
        [~, ~, ~, ~, ~, ~, segF] = analysis(num2cell(F, 2), num2cell(BF, 2), tVec, 13, SPs(q), 0, 0);
        warning(oldWarnings);
        neg = seg(seg(:, 2) < 0, :); nFalse = sum(segF(:, 2) < 0);
        % a true loss is found if it lies within a detected negative burst (+- 6 frames)
        found = false(numel(lossGranule), 1);
        for b = 1:size(neg, 1)
            found(lossGranule == neg(b, 1) & lossFrame >= neg(b, 3) - 6 & lossFrame <= neg(b, 4) + 6) = true;
        end
        recall(j, q) = 100 * mean(found);
        falseShare(j, q) = 100 * min(1, nFalse / max(size(neg, 1), 1));
    end
    fprintf('Fig 5(f) noise %d: bursts found %s %%, false share %s %%\n', noiseLevels(j), mat2str(round(recall(j, :))), mat2str(round(falseShare(j, :))));
end

if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 3.9], 'Color', 'w');
    t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    show = [1 2 4];
    for j = 1:3
        ax = nexttile(t); li = find(lams == show(j));
        histogram(ax, amps{li, 1}, 0:10:1000, 'Normalization', 'pdf', 'FaceColor', RED, 'EdgeColor', 'none', 'FaceAlpha', 1);
        xlim(ax, [0 1000]); ylim(ax, [0 0.012]); set(ax, 'XTick', 0:250:1000, 'YTick', 0:0.005:0.01); ytickformat(ax, '%.3f');
        xlabel(ax, 'Amplitude (a.u.)');
        if j == 1, ylabel(ax, 'Frequency'); end
        text(ax, 0.95, 0.85, sprintf('\\lambda_{sim}=%d', show(j)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 9);
        styleAxes(ax, char('a' + j - 1));
    end
    ax = nexttile(t); hold(ax, 'on');
    m = median(lamWS, 2); errorbar(ax, lams, m, m - min(lamWS, [], 2), max(lamWS, [], 2) - m, 'o', 'Color', BLUE, 'MarkerSize', 4, 'CapSize', 3);
    m = median(lamNLL, 2); errorbar(ax, lams, m, m - min(lamNLL, [], 2), max(lamNLL, [], 2) - m, 's', 'Color', ORANGE, 'MarkerSize', 4, 'CapSize', 3);
    plot(ax, [0 5.5], [0 5.5], 'r--', 'LineWidth', 1);
    xlim(ax, [0 5.5]); ylim(ax, [0 5.5]); xlabel(ax, '\lambda_{sim}'); ylabel(ax, '\lambda_{est}');
    legend(ax, {'WS score', 'likelihood', 'truth'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 6.5);
    styleAxes(ax, 'd');
    ax = nexttile(t); hold(ax, 'on');
    errorbar(ax, 60 * lams, lamk0(:, 1), lamk0(:, 1) - lamk0(:, 2), lamk0(:, 3) - lamk0(:, 1), 'o', 'Color', BLUE, 'MarkerSize', 4, 'CapSize', 3);
    plot(ax, [0 330], [0 330], 'r--', 'LineWidth', 1); xlim(ax, [0 330]); ylim(ax, [0 330]);
    xlabel(ax, 'true \lambda\cdotk_0 (a.u.)'); ylabel(ax, 'estimated \lambda\cdotk_0 (a.u.)');
    styleAxes(ax, 'e');
    ax = nexttile(t); hold(ax, 'on');
    h1 = plot(ax, SPs, recall(1, :), '-o', 'Color', BLUE, 'MarkerSize', 4); plot(ax, SPs, falseShare(1, :), ':', 'Color', BLUE, 'LineWidth', 1);
    h2 = plot(ax, SPs, recall(2, :), '-s', 'Color', ORANGE, 'MarkerSize', 4); plot(ax, SPs, falseShare(2, :), ':', 'Color', ORANGE, 'LineWidth', 1);
    h3 = plot(ax, NaN, NaN, 'k-o', 'MarkerSize', 3); h4 = plot(ax, NaN, NaN, 'k:', 'LineWidth', 1);
    legend(ax, [h1 h2 h3 h4], {'original noise', 'PP7 noise', 'bursts found', 'false share'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 5.5);
    xlabel(ax, 'Stringency parameter j'); ylabel(ax, '%'); ylim(ax, [0 105]);
    styleAxes(ax, 'f');
    saveFig(f, 'Fig5');
end

%% Figure 6 - in vitro granules
names = {'vitro_1_1', 'vitro_10_1', 'vitro_100_1'}; ratio = {'1:1', '10:1', '100:1'};
M6 = cell(1, 3);
for k = 1:3
    d = loadPaperDataset(names{k});
    M6{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 6 %s: median %.1f vs %.1f %%, mean burst %.0f [%.0f, %.0f], turnover %.1f [%.1f, %.1f], net %.1f [%.1f, %.1f]\n', ratio{k}, ...
        M6{k}.sizeGranules, M6{k}.sizeCopies, M6{k}.meanBurst, M6{k}.meanBurst_lo, M6{k}.meanBurst_hi, ...
        M6{k}.turnover, M6{k}.turnover_lo, M6{k}.turnover_hi, M6{k}.net, M6{k}.net_lo, M6{k}.net_hi);
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 3.7], 'Color', 'w');
    t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:3
        ax = nexttile(t); plotBurstHist(ax, M6{k});
        title(ax, sprintf('protein:RNA %s  (n = %d)', ratio{k}, numel(vertcat(M6{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6); end
    end
    measures = {'meanBurst', 'turnover', 'net'};
    labels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
    for j = 1:3
        ax = nexttile(t); hold(ax, 'on');
        y = zeros(1, 3); lo = y; hi = y;
        for k = 1:3
            y(k) = M6{k}.(measures{j}); lo(k) = M6{k}.([measures{j} '_lo']); hi(k) = M6{k}.([measures{j} '_hi']);
            bar(ax, k, y(k), 0.6, 'FaceColor', RED, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
        end
        errorbar(ax, 1:3, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
        yline(ax, 0, 'k', 'LineWidth', 0.5);
        set(ax, 'XTick', 1:3, 'XTickLabel', ratio); xlim(ax, [0.4 3.6]);
        xlabel(ax, 'protein:RNA'); ylabel(ax, labels{j});
        styleAxes(ax, char('d' + j - 1));
    end
    saveFig(f, 'Fig6');
end

%% Figure 7 - PP7 granules in vivo
names = {'4PP7', '30PP7'}; labels = {'PP7 + QCP-5x/PCP-4x', 'PP7 + PCP-30x'};
M7 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M7{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 7 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', labels{k}, ...
        M7{k}.sizeGranules, M7{k}.sizeCopies, M7{k}.turnover, M7{k}.turnover_lo, M7{k}.turnover_hi);
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M7{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M7{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    saveFig(f, 'Fig7');
end

%% Figure 8 - Qbeta granules in vivo (uses M7 from Figure 7 for the tests)
names = {'5Qb', '10Qb'}; labels = {['Q' char(946) ' + QCP-5x/PCP-4x'], ['Q' char(946) ' + QCP-10x']};
M8 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M8{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 8 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', labels{k}, ...
        M8{k}.sizeGranules, M8{k}.sizeCopies, M8{k}.turnover, M8{k}.turnover_lo, M8{k}.turnover_hi);
end
compareConditions(M8{1}, M7{1}, 'Qbeta 5x', 'PP7 4x');          % Qbeta vs PP7 (text of Sec. V.B)
compareConditions(M8{2}, M7{2}, 'Qbeta 10x', 'PP7 30x');
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M8{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M8{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    saveFig(f, 'Fig8');
end

%% Figure 9 - PP7-FUS granules in vivo (uses M7 from Figure 7 as the PP7 reference)
names = {'4PP7-FUS', '30PP7-FUS'}; labels = {'PP7-FUS + QCP-5x/PCP-4x', 'PP7-FUS + PCP-30x'};
M9 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M9{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
end
T4 = compareConditions(M9{1}, M7{1}, 'PP7-FUS 4x', 'PP7 4x');
T30 = compareConditions(M9{2}, M7{2}, 'PP7-FUS 30x', 'PP7 30x');
pBright4 = T4.p(strcmp(T4.measure, 'brightnessPooled')); pBright30 = T30.p(strcmp(T30.measure, 'brightnessPooled'));
pBeyond4 = T4.p(strcmp(T4.measure, 'sizeBeyond')); pBeyond4one = T4.p_AgtB(strcmp(T4.measure, 'sizeBeyond'));
pBeyond30 = T30.p(strcmp(T30.measure, 'sizeBeyond'));
% (d) PP7-FUS minus PP7: total gap in median burst size, the part reproduced by the copies,
% and the part beyond noise; 95% CI from the two conditions' bootstrap samples
parts = {'sizeGranules', 'sizeCopies', 'sizeBeyond'};
gap = zeros(2, 3); gapLo = gap; gapHi = gap;
for i = 1:2
    for k = 1:3
        gap(i, k) = M9{i}.(parts{k}) - M7{i}.(parts{k});
        bootDiff = M9{i}.boot.(parts{k}) - M7{i}.boot.(parts{k});
        gapLo(i, k) = prctile(bootDiff, 2.5); gapHi(i, k) = prctile(bootDiff, 97.5);
    end
    fprintf('Fig 9(d) %s: total %.1f [%.1f, %.1f], noise %.1f [%.1f, %.1f], beyond %.1f [%.1f, %.1f]\n', labels{i}, [gap(i, :); gapLo(i, :); gapHi(i, :)]);
end
fprintf('Fig 9(c) brightness p = %.4f (4x), %.2f (30x); 9(d) beyond noise 4x: p = %.3f one-sided, %.3f two-sided; 30x: p = %.2f\n', ...
    pBright4, pBright30, pBeyond4one, pBeyond4, pBeyond30);
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 4.6], 'Color', 'w');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M9{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M9{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    % (c) granule brightness (median over granules, 95% CI)
    ax = nexttile(t); hold(ax, 'on');
    Ms = {M7{1}, M9{1}, M7{2}, M9{2}}; x = [1 2 3.6 4.6]; colors = [CPP; CFUS; CPP; CFUS];
    y = zeros(1, 4); lo = y; hi = y;
    for k = 1:4
        y(k) = Ms{k}.brightnessPooled / 1000; lo(k) = Ms{k}.brightnessPooled_lo / 1000; hi(k) = Ms{k}.brightnessPooled_hi / 1000;
        bar(ax, x(k), y(k), 0.7, 'FaceColor', colors(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    set(ax, 'XTick', x, 'XTickLabel', {'PP7', 'PP7-FUS', 'PP7', 'PP7-FUS'}); ylabel(ax, 'Granule brightness (a.u./1000)');
    yt = max(hi) * 1.06; ylim(ax, [0 max(hi) * 1.3]);         % group names and p values above the bars
    text(ax, 1.5, yt, {'QCP-5x/PCP-4x', sprintf('p = %.3f', pBright4)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 7);
    text(ax, 4.1, yt, {'PCP-30x', sprintf('p = %.2f', pBright30)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    styleAxes(ax, 'c');
    % (d) burst size gap
    ax = nexttile(t); hold(ax, 'on');
    colors = [0.55 0.55 0.55; 0.85 0.85 0.85; CFUS]; hb = gobjects(1, 3);
    for i = 1:2
        for k = 1:3
            xx = (i - 1) * 1.1 + (k - 2) * 0.28;
            hb(k) = bar(ax, xx, gap(i, k), 0.26, 'FaceColor', colors(k, :), 'EdgeColor', 'k', 'LineWidth', 0.4);
            errorbar(ax, xx, gap(i, k), gap(i, k) - gapLo(i, k), gapHi(i, k) - gap(i, k), 'k', 'LineWidth', 0.8, 'CapSize', 3);
        end
    end
    text(ax, 0, 11.8, sprintf('p = %.3f (one-sided)\np = %.3f (two-sided)', pBeyond4one, pBeyond4), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 6.5);
    text(ax, 1.1, 11.8, sprintf('p = %.2f', pBeyond30), 'HorizontalAlignment', 'center', 'FontSize', 6.5);
    yline(ax, 0, 'k', 'LineWidth', 0.5); ylim(ax, [-4 17]);
    set(ax, 'XTick', [0 1.1], 'XTickLabel', {'QCP-5x/PCP-4x', 'PCP-30x'}); xlim(ax, [-0.55 1.65]);
    ylabel(ax, sprintf('PP7-FUS - PP7, median burst\nsize (%% points)'));
    lg = legend(ax, hb, {'total', 'explained by noise', 'beyond noise'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 6);
    lg.ItemTokenSize = [10 6];
    styleAxes(ax, 'd');
    saveFig(f, 'Fig9');
end


%% ------------------------------------------------------------------ functions
function d = loadPaperDataset(name)
% One dataset of the paper: spot mean, background, spot sum and the movie of every trace.
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/';
    file = strrep(name, 'vitro', 'same');                  % e.g. vitro_1_1 -> same_1_1
    info = readtable([folder 'extraction_info_same.csv'], 'Delimiter', ',', 'TextType', 'string');
    info = info(info.dataset == file, :);
    d.spotMean = loadTraces([folder file '_data.mat']);
    d.background = loadTraces([folder file '_background.mat']);
    d.spotSum = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    info = [readtable([folder 'extraction_info.csv'], 'TextType', 'string'); ...
            readtable([folder 'extraction_info_5Qb_10Qb.csv'], 'TextType', 'string')];
    info = info(info.dataset == name, :);
    d.spotMean = loadTraces([folder name '_data.mat']);
    d.background = loadTraces([folder name '_background.mat']);
    d.spotSum = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
d.movie = repelem((1:height(info))', info.nKept);         % the traces are stored movie by movie
end


function traces = loadTraces(file)
% The extraction stores a cell array of traces whose first cell is empty.
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end


function plotBurstHist(ax, M)
% Burst sizes (% of the granule) of the granules (red bars) and of their copies (black outline).
fracReal = 100 * vertcat(M.fracReal{:}); fracCopy = 100 * vertcat(M.fracCopy{:});
edges = 0:2.5:60;
histogram(ax, min(fracReal, 59.99), edges, 'Normalization', 'pdf', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
hold(ax, 'on');
histogram(ax, min(fracCopy, 59.99), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
xlim(ax, [0 60]); xlabel(ax, 'Burst size (% of granule)');
text(ax, 0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
    M.sizeGranules, M.sizeCopies, M.turnover, M.turnover_lo, M.turnover_hi), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 6.5);
end


function styleAxes(ax, panel)
% Manuscript style: boxed axes, ticks inward, Arial 8, bold panel label "(a)".
set(ax, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.6, 'Layer', 'top', ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'TitleFontWeight', 'normal');
ax.Title.FontSize = 7.5;
if isempty(ax.Title.String), title(ax, ' '); end            % leaves room for the panel label
text(ax, -0.2, 1.04, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
    'VerticalAlignment', 'bottom', 'FontName', 'Arial');
end


function saveFig(f, name)
% Saves a figure in figures/ as PNG (300 dpi) and PDF.
exportgraphics(f, fullfile('figures', [name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile('figures', [name '.pdf']), 'ContentType', 'vector');
end


function [lamWS, lamNLL] = eq19Fit(a)
% Fits the manuscript's Eq. (1.9) to burst amplitudes a (MS Fig. 5(d)):
%   p(I) = (1/k0) lambda^(I/k0) exp(-lambda) / Gamma(I/k0 + 1)
% on a grid of lambda (0.05-20, plus 1..5) and k0 (2-600), with two scores:
%   WS  = the manuscript's weighted rank score: 0.5 rank(MSE) + 0.25 rank(QQ) + 0.25 rank(KS)
%   NLL = maximum likelihood
% Returns the lambda of the best WS and of the best likelihood.
lamGrid = unique([exp(linspace(log(0.05), log(20), 80)), 1:5]);
k0Grid = exp(linspace(log(2), log(600), 140));
a = sort(a(:))'; n = numel(a);
edges = linspace(a(1), a(end), 11); counts = histcounts(a, edges);         % 10-bin histogram, as the MS
centers = (edges(1:end-1) + edges(2:end)) / 2;
NLL = zeros(numel(lamGrid), numel(k0Grid)); KS = NLL; QQ = NLL; MSE = NLL;
for i = 1:numel(lamGrid)
    lam = lamGrid(i);
    % the model (in units of k0), normalised, and its cumulative distribution
    x = linspace(0, lam + 15 * sqrt(lam) + 30, 6000);
    pdf = exp(x * log(lam) - lam - gammaln(x + 1)); pdf = pdf / trapz(x, pdf);
    cdf = [0, cumsum((pdf(2:end) + pdf(1:end-1)) / 2 .* diff(x))];
    [cdfUnique, iu] = unique(cdf);
    for j = 1:numel(k0Grid)
        k0 = k0Grid(j);
        F = interp1(x, cdf, a / k0, 'linear', 1);                % model CDF at the data
        f = interp1(x, pdf, a / k0, 'linear', 0) / k0;           % model density at the data
        NLL(i, j) = -sum(log(f + 1e-300));
        KS(i, j) = max(max(abs(F - (1:n) / n)), max(abs(F - (0:n-1) / n)));
        quantiles = interp1(cdfUnique, x(iu), min(max(((1:n) - 0.5) / n, cdfUnique(1)), cdfUnique(end)));
        QQ(i, j) = sqrt(mean((a - k0 * quantiles) .^ 2));
        m = interp1(x, pdf, centers / k0, 'linear', 0) / k0;      % model at the histogram bins,
        m = m / max(max(m), 1e-300) * max(counts);               % scaled to the histogram peak
        MSE(i, j) = mean((m - counts) .^ 2);
    end
end
WS = 0.5 * rankAll(MSE) + 0.25 * rankAll(QQ) + 0.25 * rankAll(KS);
lamWS = lamGrid(bestRow(WS)); lamNLL = lamGrid(bestRow(NLL));
end


function R = rankAll(A)
% rank of every grid cell among all cells, scaled to (0, 1]
[~, order] = sort(A(:)); r = zeros(numel(A), 1); r(order) = 1:numel(A);
R = reshape(r, size(A)) / numel(A);
end


function i = bestRow(A)
% row (lambda) of the smallest value; ties go to the first lambda, then the first k0
At = A'; [~, q] = min(At(:)); [~, i] = ind2sub(size(At), q);
end


function [lam, k0] = shapeFit(a)
% Shape-only fit of lambda and k0 to burst amplitudes a (MS Fig. 5(e)), as the original
% model curve: the Poisson probabilities at whole numbers of units (zero units excluded)
% joined by straight lines, at A/k0, normalised over [0, max(A)]; maximum likelihood.
a = a(:); top = max(a) * 1.001; units = 0:399;
best = inf;
for lam0 = [0.3 1 3 10]                                         % several starting points
    for factor = [0.7 1.5]
        start = [log(lam0), log(factor * mean(a) / max(lam0, 1))];
        [p, value] = fminsearch(@(p) shapeNLL(p, a, top, units), start, ...
            optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 3000, 'MaxFunEvals', 6000, 'Display', 'off'));
        if value < best, best = value; lam = exp(p(1)); k0 = exp(p(2)); end
    end
end
end


function v = shapeNLL(p, a, top, units)
lam = exp(p(1)); k0 = exp(p(2));
if lam < 0.05 || lam > 60 || k0 < top / 300 || k0 > 3 * top, v = 1e12; return; end   % bounds
prob = exp(units * log(lam) - lam - gammaln(units + 1)); prob(1) = 0;
grid = linspace(0, top / k0, 2000);
area = trapz(grid, interp1(units, prob, grid, 'linear', 0)) * k0;
v = -sum(log(interp1(units, prob, a / k0, 'linear', 0) + 1e-300) - log(area));
end
