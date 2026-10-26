% SIMULATIONRUNWITHANALYSIS  The GelMetrics analysis on simulated granules with known
% dynamics (generateSim.m), through the same code as the experimental data
% (burstMeasures.m, as in ExperimentalAnalysis.m). Makes MS Figure 5.
%
% Simulated granules: material leaves and enters in bursts of a Poisson-distributed
% number of units (lambda units of k0 each, so the mean burst size is lambda*k0).
% Every combination of lambdaList and snrList is one simulated condition; the noise
% is set by the signal-to-noise ratio of one unit: noise = k0 / SNR.
%
% For every condition the script simulates nFits "experiments" (numTracks granules
% each, a new random draw every time) and runs on each:
%   - the manuscript's original fit of Eq. (1.9) with the weighted score (WS),
%     which gives lambda and k0;
%   - burstMeasures.m, which gives the mean burst size (lambda*k0) and, for the first
%     experiment, all the measures of MS Table I with 95% confidence intervals.
% The true values are taken from the same granules simulated without noise.
%
% Figure (MS Fig. 5):
%   (a) lambda from the WS fit vs the true lambda, one dot per experiment
%   (b) k0 vs lambda from the same fits: they move together along lambda*k0 = constant
%   (c) lambda*k0: the WS product and the mean burst size vs the true value
%   (d-f) burst sizes of the granules and of their burst-free copies, at SNR histSNR
%   (g-i) mean burst size, turnover beyond noise and net change of every condition,
%         as ExperimentalAnalysis.m (black dots: true values)
%
% Run time with the default settings: about 15 min.
close all; clearvars; clc;

% Simulation parameters
lambdaList = [1 2 4];      % mean number of units per burst
snrList = [5 0.5];         % size of one unit / frame noise
k0_sim = 60;               % intensity of one unit (a.u.)
numTracks = 200;           % granules per simulated experiment
trackLength = 361;         % time points (60 min at 10 s)
expConst = 0.0026;         % photobleaching rate (1/min), median of the measured backgrounds
nFits = 20;                % simulated experiments per condition
nCopies = 20;              % burst-free copies per granule (first experiment of each condition)
histSNR = 0.5;             % SNR shown in panels (d-f)
saveFigures = false;      % true: save the figure as figures/Fig5.png and .pdf

nLam = numel(lambdaList); nSNR = numel(snrList);
lamWS = zeros(nLam, nSNR, nFits);      % lambda from the WS fit
k0WS = zeros(nLam, nSNR, nFits);       % k0 from the WS fit
meanBurst = zeros(nLam, nSNR, nFits);  % mean burst size from burstMeasures (a.u.)
M = cell(nLam, nSNR);                  % all measures of the first experiment
truth = cell(nLam, nSNR);              % true values of the first experiment
for li = 1:nLam
    for si = 1:nSNR
        noiseStd = k0_sim / snrList(si);
        for it = 1:nFits
            seed = 1000 * li + 100 * si + it;                   % a fixed seed for every experiment
            rng(seed);
            [signals, backgrounds] = generateSim(numTracks, trackLength, 3, expConst, noiseStd, lambdaList(li), k0_sim);
            spotMean = cell(numTracks, 1);
            background = cell(numTracks, 1);
            for i = 1:numTracks
                spotMean{i} = signals(i, :)';
                background{i} = backgrounds(i, :)';
            end
            % the copies only matter for the noise measures; the fit and the mean burst
            % size use the granule bursts only, so the later experiments use 1 copy (faster)
            if it == 1, copies = nCopies; else, copies = 1; end
            % simulations have no movies (every granule is its own replicate) and no spot sum
            Mi = burstMeasures(spotMean, background, [], [], copies, seed);
            amps = vertcat(Mi.ampAU{:});                        % negative bursts, a.u.
            [lamWS(li, si, it), k0WS(li, si, it)] = eq19Fit(amps);
            meanBurst(li, si, it) = Mi.meanBurst;
            if it == 1
                M{li, si} = Mi;
                truth{li, si} = trueValues(seed, numTracks, trackLength, expConst, lambdaList(li), k0_sim);
            end
        end
        fprintf('\nlambda %g, k0 %g, SNR %g (noise %.0f a.u.)\n', lambdaList(li), k0_sim, snrList(si), noiseStd);
        fprintf('  WS fit, lambda over %d experiments:   median %.1f, range %.1f - %.1f\n', nFits, ...
            median(lamWS(li, si, :)), min(lamWS(li, si, :)), max(lamWS(li, si, :)));
        fprintf('  WS fit, k0:                          median %.0f, range %.0f - %.0f\n', ...
            median(k0WS(li, si, :)), min(k0WS(li, si, :)), max(k0WS(li, si, :)));
        fprintf('  WS fit, lambda*k0:                   median %.0f, range %.0f - %.0f  (true %g)\n', ...
            median(lamWS(li, si, :) .* k0WS(li, si, :)), min(lamWS(li, si, :) .* k0WS(li, si, :)), ...
            max(lamWS(li, si, :) .* k0WS(li, si, :)), lambdaList(li) * k0_sim);
        fprintf('  mean burst size (burstMeasures):     median %.0f, range %.0f - %.0f\n', ...
            median(meanBurst(li, si, :)), min(meanBurst(li, si, :)), max(meanBurst(li, si, :)));
        Mc = M{li, si}; T = truth{li, si};
        fprintf('  first experiment, %d copies per granule (true values in brackets):\n', nCopies);
        fprintf('    mean burst size (a.u.)               %.0f [%.0f, %.0f]   (true %.0f)\n', Mc.meanBurst, Mc.meanBurst_lo, Mc.meanBurst_hi, T.meanBurst);
        fprintf('    burst size, granules / copies (%%)    %.1f / %.1f   (true %.1f)\n', Mc.sizeGranules, Mc.sizeCopies, T.size);
        fprintf('    burst size beyond noise (%% points)   %.1f [%.1f, %.1f]\n', Mc.sizeBeyond, Mc.sizeBeyond_lo, Mc.sizeBeyond_hi);
        fprintf('    turnover beyond noise (%%/h)          %.1f [%.1f, %.1f]   (true %.1f)\n', Mc.turnover, Mc.turnover_lo, Mc.turnover_hi, T.turnover);
        fprintf('    net change beyond bleaching (%%/h)    %.1f [%.1f, %.1f]   (true %.1f)\n', Mc.net, Mc.net_lo, Mc.net_hi, T.net);
    end
end

%% Figure
BLUE = [0 0.447 0.741]; ORANGE = [0.835 0.369 0];
snrColors = [BLUE; ORANGE; 0.47 0.67 0.19; 0.49 0.18 0.56];   % one color per SNR
f = figure('Units', 'inches', 'Position', [1 1 6.5 6.2], 'Color', 'w');
t = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
lam = char(955); dot = char(183);   % the letter lambda and a middle dot, typed directly (TeX codes stretch the spaces in the labels)

% (a) lambda from the WS fit, one dot per experiment
ax = nexttile(t); hold(ax, 'on');
h = gobjects(1, nSNR);
for si = 1:nSNR
    for li = 1:nLam
        x = lambdaList(li) + 0.15 * (si - (nSNR + 1) / 2) + linspace(-0.05, 0.05, nFits)';   % small spread so the dots do not overlap
        h(si) = plot(ax, x, squeeze(lamWS(li, si, :)), 'o', 'Color', snrColors(si, :), 'MarkerSize', 2.5);
    end
end
top = max(max(lambdaList) + 1, max(lamWS(:)) + 0.5);
plot(ax, [0 top], [0 top], 'k--', 'LineWidth', 0.8);
xlim(ax, [0 top]); ylim(ax, [0 top]);
xlabel(ax, [lam ' sim']); ylabel(ax, [lam ' from the WS fit']);
legendText = cell(1, nSNR);
for si = 1:nSNR, legendText{si} = sprintf('SNR %g', snrList(si)); end
legend(ax, h, legendText, 'Location', 'north', 'Box', 'off', 'FontSize', 6);   % top middle, clear of the dots
styleAxes(ax, 'a');

% (b) k0 vs lambda from the same fits; dashed lines: lambda*k0 = the true value
ax = nexttile(t); hold(ax, 'on');
markers = 'osd^v';
for li = 1:nLam
    lamLine = linspace(0.2, top, 100);
    plot(ax, lamLine, lambdaList(li) * k0_sim ./ lamLine, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    for si = 1:nSNR
        plot(ax, squeeze(lamWS(li, si, :)), squeeze(k0WS(li, si, :)), markers(li), 'Color', snrColors(si, :), 'MarkerSize', 3);
    end
end
xlim(ax, [0 top]); ylim(ax, [0 1.1 * max(k0WS(:))]);
xlabel(ax, [lam ' from the WS fit']); ylabel(ax, 'k0 from the WS fit (a.u.)');
markerText = cell(1, nLam); hm = gobjects(1, nLam);
for li = 1:nLam
    hm(li) = plot(ax, NaN, NaN, markers(li), 'Color', 'k', 'MarkerSize', 3);
    markerText{li} = sprintf('%s sim = %g', lam, lambdaList(li));
end
legend(ax, hm, markerText, 'Location', 'northeast', 'Box', 'off', 'FontSize', 6);
styleAxes(ax, 'b');

% (c) lambda*k0: WS product and mean burst size (median and range over the experiments)
ax = nexttile(t); hold(ax, 'on');
trueProduct = lambdaList * k0_sim;
for si = 1:nSNR
    for li = 1:nLam
        p = squeeze(lamWS(li, si, :) .* k0WS(li, si, :)); mb = squeeze(meanBurst(li, si, :));
        x = trueProduct(li) + 8 * (si - (nSNR + 1) / 2);
        errorbar(ax, x - 3, median(p), median(p) - min(p), max(p) - median(p), 'o', 'Color', snrColors(si, :), 'MarkerSize', 3, 'CapSize', 2);
        errorbar(ax, x + 3, median(mb), median(mb) - min(mb), max(mb) - median(mb), 's', 'Color', snrColors(si, :), ...
            'MarkerFaceColor', snrColors(si, :), 'MarkerSize', 3, 'CapSize', 2);
    end
end
topP = 1.6 * max(trueProduct);
plot(ax, [0 topP], [0 topP], 'k--', 'LineWidth', 0.8);
xlim(ax, [0 topP]); ylim(ax, [0 topP]);
h1 = plot(ax, NaN, NaN, 'ko', 'MarkerSize', 3); h2 = plot(ax, NaN, NaN, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
legend(ax, [h1 h2], {[lam dot 'k0 from the WS fit'], 'mean burst size'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 6);
xlabel(ax, ['true ' lam dot 'k0 (a.u.)']); ylabel(ax, ['estimated ' lam dot 'k0 (a.u.)']);
styleAxes(ax, 'c');

% (d-f) burst sizes of the granules and of their copies, at SNR histSNR
si = find(snrList == histSNR, 1);
for li = 1:min(nLam, 3)
    ax = nexttile(t);
    plotBurstHist(ax, M{li, si});
    title(ax, sprintf('%s sim %g, SNR %g, n = %d', lam, lambdaList(li), histSNR, numel(vertcat(M{li, si}.fracReal{:}))));
    if li == 1
        ylabel(ax, 'Frequency');
        lg = legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6);
        lg.ItemTokenSize = [10 6];
    end
    styleAxes(ax, char('c' + li));
end
for li = nLam + 1:3, nexttile(t); axis off; end      % fewer than three lambdas: empty tiles

% (g-i) all measures of every condition, as ExperimentalAnalysis.m; black dots: true values
measures = {'meanBurst', 'turnover', 'net'};
labels = {['Mean burst size, ' lam dot 'k0 (a.u.)'], 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
trueField = {'meanBurst', 'turnover', 'net'};
for j = 1:3
    ax = nexttile(t); hold(ax, 'on');
    xTicks = zeros(1, nLam); hb = gobjects(1, nSNR);
    for li = 1:nLam
        for si = 1:nSNR
            x = (li - 1) * (nSNR + 1) + si;                      % conditions grouped by lambda
            y = M{li, si}.(measures{j}); lo = M{li, si}.([measures{j} '_lo']); hi = M{li, si}.([measures{j} '_hi']);
            hb(si) = bar(ax, x, y, 0.8, 'FaceColor', snrColors(si, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
            errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
            plot(ax, x, truth{li, si}.(trueField{j}), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
        end
        xTicks(li) = (li - 1) * (nSNR + 1) + (nSNR + 1) / 2;
    end
    yline(ax, 0, 'k', 'LineWidth', 0.5);
    tickText = cell(1, nLam);
    for li = 1:nLam, tickText{li} = sprintf('%s sim = %g', lam, lambdaList(li)); end
    set(ax, 'XTick', xTicks, 'XTickLabel', tickText, 'XTickLabelRotation', 0);
    xlim(ax, [0.2, nLam * (nSNR + 1) - 0.2]);
    ylabel(ax, labels{j});
    if j == 1
        lg = legend(ax, [hb plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3)], [legendText {'true value'}], ...
            'Location', 'northwest', 'Box', 'off', 'FontSize', 6);
        lg.ItemTokenSize = [10 6];   % short legend symbols, so the true-value dot sits next to its label
    end
    styleAxes(ax, char('f' + j));
    ax.XAxis.FontSize = 6.5;   % smaller tick labels (after the style, which sets 8), so they do not run into each other
end

if saveFigures
    if ~exist('figures', 'dir'), mkdir('figures'); end
    exportgraphics(f, fullfile('figures', 'Fig5.png'), 'Resolution', 300);
    exportgraphics(f, fullfile('figures', 'Fig5.pdf'), 'ContentType', 'vector');
end


%% ------------------------------------------------------------------ functions
function T = trueValues(seed, numTracks, trackLength, expConst, lambda, k0)
% The true values of one simulated experiment: the same granules (same seed) simulated
% without noise, before photobleaching. Burst size as burstMeasures.m measures it: the
% amplitude relative to the granule above its background, at the start of the burst.
rng(seed);
[signals, backgrounds] = generateSim(numTracks, trackLength, 3, expConst, 0, lambda, k0);
t = linspace(0, 60, trackLength);
amps = []; fractions = []; nets = zeros(numTracks, 1);
for i = 1:numTracks
    % net change as burstMeasures.m measures it: robust slope of ln(granule / background)
    % (as there: only frames where the granule is above its background; NaN if fewer than 90%)
    ratio = (signals(i, :) - backgrounds(i, :)) ./ backgrounds(i, :);
    frames = find(ratio > 0); v = log(ratio(frames));
    slopes = [];
    for j = 1:numel(frames) - 1
        slopes = [slopes, (v(j+1:end) - v(j)) ./ (frames(j+1:end) - frames(j))]; %#ok<AGROW>
    end
    if numel(frames) < 0.9 * trackLength
        nets(i) = NaN;
    else
        nets(i) = median(slopes) * 360 * 100;            % %/h
    end
    % bursts
    u = signals(i, :) ./ exp(-expConst * t);             % granule without bleaching
    b = backgrounds(i, :) ./ exp(-expConst * t);         % its background without bleaching
    for j = 2:trackLength
        step = u(j) - u(j - 1);
        if step < -1e-9                                  % material leaves the granule
            amps(end+1, 1) = -step; %#ok<AGROW>
            fractions(end+1, 1) = -step / (u(j - 1) - b(j - 1)); %#ok<AGROW>
        end
    end
end
hours = numTracks * trackLength / 360;
T.meanBurst = mean(amps);                                % a.u.
T.size = 100 * median(fractions);                        % % of the granule
T.turnover = 100 * sum(fractions) / hours;               % %/h
T.net = mean(nets, 'omitnan');                                    % %/h (about 0: bursts in and out are equally likely)
end


function [lamWS, k0WS] = eq19Fit(a)
% The manuscript's original fit (MS section III): Eq. (1.9),
%   p(I) = (1/k0) lambda^(I/k0) exp(-lambda) / Gamma(I/k0 + 1),
% on a grid of lambda (0.05-20, plus 1..5) and k0 (2-600), scored by the weighted score
%   WS = 0.5 rank(MSE) + 0.25 rank(QQ) + 0.25 rank(KS)
% (MSE on a 10-bin histogram, QQ and KS on the amplitudes). Returns the best lambda and k0.
lamGrid = unique([exp(linspace(log(0.05), log(20), 80)), 1:5]);
k0Grid = exp(linspace(log(2), log(600), 140));
a = sort(a(:))'; n = numel(a);
edges = linspace(a(1), a(end), 11); counts = histcounts(a, edges);         % 10-bin histogram, as the MS
centers = (edges(1:end-1) + edges(2:end)) / 2;
KS = zeros(numel(lamGrid), numel(k0Grid)); QQ = KS; MSE = KS;
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
        KS(i, j) = max(max(abs(F - (1:n) / n)), max(abs(F - (0:n-1) / n)));
        quantiles = interp1(cdfUnique, x(iu), min(max(((1:n) - 0.5) / n, cdfUnique(1)), cdfUnique(end)));
        QQ(i, j) = sqrt(mean((a - k0 * quantiles) .^ 2));
        m = interp1(x, pdf, centers / k0, 'linear', 0) / k0;      % model at the histogram bins,
        m = m / max(max(m), 1e-300) * max(counts);               % scaled to the histogram peak
        MSE(i, j) = mean((m - counts) .^ 2);
    end
end
WS = 0.5 * rankAll(MSE) + 0.25 * rankAll(QQ) + 0.25 * rankAll(KS);
% the best grid cell; ties go to the first lambda, then the first k0
WSt = WS'; [~, q] = min(WSt(:)); [j, i] = ind2sub(size(WSt), q);
lamWS = lamGrid(i); k0WS = k0Grid(j);
end


function R = rankAll(A)
% rank of every grid cell among all cells, scaled to (0, 1]
[~, order] = sort(A(:)); r = zeros(numel(A), 1); r(order) = 1:numel(A);
R = reshape(r, size(A)) / numel(A);
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
if isempty(ax.Title.String), title(ax, ' '); end           % leaves room for the panel label
text(ax, -0.2, 1.04, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
    'VerticalAlignment', 'bottom', 'FontName', 'Arial');
end
