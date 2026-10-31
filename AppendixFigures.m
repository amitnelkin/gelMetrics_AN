%% AppendixFigures.m - the figures of the appendices (algorithm debugging), Figs. 11-18
% Run the first section ("Settings and data") once, then any appendix section on its own
% (Ctrl+Enter). Every appendix section computes its figure from the repository data with
% the same functions as ExperimentalAnalysis.m (analysis.m, burstMeasures.m, generateSim.m,
% PoissonFit2WS.m). Two panels need the raw movies, which are not in the repository
% (Fig. 14(b) camera noise, Fig. 14(c) static synthetic granules); they are drawn from the
% committed result tables in analysis_outputs/.
% The run time of each section is given in its title.
%
%   APPENDIX C  Full analysis results         Fig. 11 (WS fit per lambda), Fig. 12 (stringency)
%   APPENDIX D  Two-channel movies            Fig. 13
%   APPENDIX E  Noise floor and its source    Fig. 14
%   APPENDIX F  Detection settings            Fig. 15
%   APPENDIX G  Estimation of lambda and k0   Fig. 16
%   APPENDIX H  Large slow changes and time   Fig. 17 (large slow changes), Fig. 18 (halves)


%% Settings and data - run this section first (a few seconds)
close all; clearvars; clc;

% Algorithm Parameters (the same as in ExperimentalAnalysis.m)
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 12;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 361);      % Time steps vector: 361 frames, 60 min at 10 s per frame
saveFigures = false;                      % true: save every figure as PNG in figures/
moviesAsReplicates = false;

vivoFolder = './in_vivo_experiments/';
vitroFolder = './in_vitro_experiments/';
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'};
labels = {'PP7 4x', 'PP7-FUS 4x', 'PP7 30x', 'PP7-FUS 30x', 'Q\beta 5x', 'Q\beta 10x', 'vitro 1:1', 'vitro 10:1', 'vitro 100:1'};
csvNames = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'};   % names in analysis_outputs/
nVivo = 6;                                % the first six datasets are in vivo
colors = [0 0.447 0.698; 0 0.620 0.451; 0 0.447 0.698; 0 0.620 0.451; 0.835 0.369 0; 0.835 0.369 0; ...
          0.4 0.4 0.4; 0.4 0.4 0.4; 0.4 0.4 0.4];    % PP7 blue, PP7-FUS green, Qbeta orange, in vitro grey

% Load every dataset; keep the traces of at least 180 frames (30 min), as the analysis does
D = struct('name', names, 'S', [], 'B', [], 'G', [], 'X', [], 'replicate', []);
for c = 1:numel(names)
    if c <= nVivo
        folder = [vivoFolder names{c} '/'];
    else
        folder = [vitroFolder names{c} '/'];
    end
    signals = struct2cell(load([folder 'data.mat']));
    backgrounds = struct2cell(load([folder 'background.mat']));
    spotSums = struct2cell(load([folder 'spotsum.mat']));
    movies = load([folder 'movies.mat']);
    signals = signals{1}(2:end);          % the first cell is the empty placeholder of the extraction
    backgrounds = backgrounds{1}(2:end);
    spotSums = spotSums{1}(2:end);
    for k = 1:numel(signals)
        if numel(signals{k}) >= 180
            D(c).S{end+1, 1} = double(signals{k}(:));                          % spot mean
            D(c).B{end+1, 1} = double(backgrounds{k}(:));                      % background
            D(c).G{end+1, 1} = double(spotSums{k}(:));                         % spot sum (granule above background)
            D(c).X{end+1, 1} = double(spotSums{k}(:)) + double(backgrounds{k}(:));   % spot sum + background, as analysis.m expects
            if moviesAsReplicates
                D(c).replicate(end+1, 1) = movies.movieId(k);
            else
                D(c).replicate(end+1, 1) = numel(D(c).S);
            end
        end
    end
    fprintf('%s: %d granules, %d replicates\n', names{c}, numel(D(c).S), numel(unique(D(c).replicate)));
end
if saveFigures && ~exist('figures', 'dir')
    mkdir('figures');
end


%% =====================================================================================
%  APPENDIX C - FULL ANALYSIS RESULTS
%  ======================================================================================

%% Appendix C, Fig. 11 - WS fit of Eq. (1.9) for every lambda, in vivo (about 2 min)
% Negative burst amplitudes (a.u.) of every in vivo dataset, the top 10% removed as in the
% original ExperimentalAnalysis.m, fitted on a grid of lambda and k0. For every lambda: the
% lowest WS over k0 (a) and the k0 where it is reached (b).
lambdaGrid = 1:5;
K0Vector = 1:500;
wsBest = zeros(nVivo, numel(lambdaGrid));
k0Best = zeros(nVivo, numel(lambdaGrid));
for c = 1:nVivo
    [~, ~, ~, ampNeg] = analysis(D(c).S, D(c).B, timeVecTotal, movingAverageWindowSpan, stringencyParameter, 0, 0);
    amps = rmoutliers(abs(ampNeg(2:end)), 'percentiles', [0 90]);
    [~, ~, WS] = PoissonFit2WS(amps, lambdaGrid, K0Vector, 10);
    for l = 1:numel(lambdaGrid)
        [wsBest(c, l), ix] = min(WS(l, :));
        k0Best(c, l) = K0Vector(ix);
    end
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
ax1 = nexttile; hold on;
ax2 = nexttile; hold on;
for c = 1:nVivo
    if any(c == [1 2 5])
        style = '-o';       % low-valency RNA: solid
    else
        style = '--o';      % high-valency RNA: dashed
    end
    plot(ax1, lambdaGrid, wsBest(c, :), style, 'Color', colors(c, :), 'MarkerSize', 4, 'DisplayName', labels{c});
    plot(ax2, lambdaGrid, k0Best(c, :), style, 'Color', colors(c, :), 'MarkerSize', 4);
end
xlabel(ax1, '\lambda_{est}'); ylabel(ax1, 'WS');
xlabel(ax2, '\lambda_{est}'); ylabel(ax2, 'k_0 (a.u.)');
legend(ax1, 'Location', 'best', 'Box', 'off', 'FontSize', 7, 'NumColumns', 2);
title(ax1, '(a)'); title(ax2, '(b)');
if saveFigures
    exportgraphics(f, 'figures/AppendixFig11.png', 'Resolution', 300);
end

%% Appendix C, Fig. 12 - PP7 vs PP7-FUS at stringency 10, 12 and 14 (about 1 min)
% Negative burst amplitudes (a.u.) of PP7 (red) and PP7-FUS (blue) granules, for the
% low-valency (top) and high-valency (bottom) RNA.
spList = [10 12 14];
rowPairs = [1 2; 3 4];          % 4PP7 vs 4PP7-FUS, 30PP7 vs 30PP7-FUS
rowNames = {'QCP-5x/PCP-4x', 'PCP-30x'};
amps = cell(2, numel(spList), 2);
for r = 1:2
    for s = 1:numel(spList)
        for k = 1:2
            c = rowPairs(r, k);
            [~, ~, ~, ampNeg] = analysis(D(c).S, D(c).B, timeVecTotal, movingAverageWindowSpan, spList(s), 0, 0);
            amps{r, s, k} = abs(ampNeg(2:end));
        end
    end
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 3.6], 'Color', 'w');
tiledlayout(2, 3, 'TileSpacing', 'compact');
barColors = {[1 0 0], [0 0 1]};
for r = 1:2
    for s = 1:numel(spList)
        nexttile; hold on;
        for k = 1:2
            c = rowPairs(r, k);
            histogram(min(amps{r, s, k}, 1199), 0:50:1200, 'Normalization', 'pdf', 'FaceColor', barColors{k}, 'FaceAlpha', 0.45, ...
                'EdgeColor', 'none', 'DisplayName', sprintf('%s (n=%d)', labels{c}, numel(amps{r, s, k})));
        end
        xlim([0 1200]);
        legend('Location', 'northeast', 'Box', 'off', 'FontSize', 6);
        if r == 1
            title(sprintf('SP = %d', spList(s)));
        else
            xlabel('Amplitude (a.u.)');
        end
        if s == 1
            ylabel({rowNames{r}, 'Frequency'});
        end
    end
end
if saveFigures
    exportgraphics(f, 'figures/AppendixFig12.png', 'Resolution', 300);
end


%% =====================================================================================
%  APPENDIX D - IMAGE ANALYSIS OF TWO-CHANNEL MOVIES
%  ======================================================================================

%% Appendix D, Fig. 13 - the zigzag of the original extraction (seconds)
% (a) A trace of the original data (a two-channel movie: fluorescence and DIC frames read
% one after the other) and the longest trace of the corrected extraction. The corrected
% traces have the camera offset (200) removed; it is added back here for the same scale.
% (b) The first 16 time points of the original trace.
old = struct2cell(load('./in_vivo_experiments-old/4PP7_data.mat'));
old = old{1};
oldTrace = double(old{50}(:));           % trace 49 (the first cell is empty), movie 007
lengths = zeros(numel(D(1).S), 1);
for k = 1:numel(D(1).S)
    lengths(k) = numel(D(1).S{k});
end
[~, longest] = max(lengths);
newTrace = D(1).S{longest} + 200;

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 2.2], 'Color', 'w');
tiledlayout(1, 4, 'TileSpacing', 'compact');
nexttile([1 3]); hold on;
plot((0:numel(oldTrace) - 1) / 6, oldTrace, 'Color', [1 0.4 0.4], 'LineWidth', 0.5, 'DisplayName', 'original extraction (two-channel movie)');
plot((0:numel(newTrace) - 1) / 6, newTrace, 'Color', [0 0.447 0.741], 'LineWidth', 0.6, 'DisplayName', 'corrected extraction (fluorescence only)');
xlim([0 60]); xlabel('Time (min)'); ylabel('Spot intensity (a.u.)');
legend('Location', 'east', 'FontSize', 7); title('(a)');
nexttile;
plot(1:16, oldTrace(1:16), 'o-', 'Color', [1 0.4 0.4], 'MarkerSize', 3);
xlabel('Time point'); ylabel('Spot intensity (a.u.)'); title('(b)');
if saveFigures
    exportgraphics(f, 'figures/AppendixFig13.png', 'Resolution', 300);
end


%% =====================================================================================
%  APPENDIX E - NOISE FLOOR AND ITS SOURCE
%  ======================================================================================

%% Appendix E, Fig. 14 - noise floor (about 5 min)
% (a) Negative bursts in the burst-free copies (5 per granule) divided by those in the
% granules, for stringency 10-20 (window 13); in vivo solid, in vitro dotted.
% (b) Measured noise / noise predicted from camera and photon statistics, and
% (c) turnover beyond noise of real and of static synthetic granules: both need the raw
% movies, so they are read from analysis_outputs/ (checks/noise_from_background.m and
% checks/static_control.m, run on the lab PC).
spList = 10:2:20;
nCopies = 5;
ratio = zeros(numel(names), numel(spList));
for c = 1:numel(names)
    rng(c);
    [C, BC] = burstFreeCopies(D(c).X, D(c).B, nCopies);
    for s = 1:numel(spList)
        [~, ~, ~, ~, ~, ~, segReal] = analysis(D(c).X, D(c).B, timeVecTotal, movingAverageWindowSpan, spList(s), 0, 0);
        [~, ~, ~, ~, ~, ~, segCopy] = analysis(C, BC, timeVecTotal, movingAverageWindowSpan, spList(s), 0, 0);
        nReal = sum(segReal(:, 2) < 0 & segReal(:, 6) > 0);          % negative bursts (as burstMeasures.m counts them)
        nCopy = sum(segCopy(:, 2) < 0 & segCopy(:, 6) > 0) / nCopies;
        ratio(c, s) = nCopy / nReal;
    end
end
noiseVivo = readtable('analysis_outputs/noise_background/summary.csv', 'TextType', 'string');
noiseVitro = readtable('analysis_outputs/noise_background_vitro/summary.csv', 'TextType', 'string');
noiseTable = [noiseVivo(:, {'dataset', 'ratioSum', 'ratioCell'}); noiseVitro(:, {'dataset', 'ratioSum', 'ratioCell'})];
staticVivo = readtable('analysis_outputs/steady/static_control_summary.csv', 'TextType', 'string');
staticVitro = readtable('analysis_outputs/steady/static_control_vitro_sum_summary.csv', 'TextType', 'string');
staticColumns = {'dataset', 'version', 'real_excess_w13', 'static_excess_w13', 'static_excess_se_w13'};
staticTable = [staticVivo(:, staticColumns); staticVitro(:, staticColumns)];
staticTable = staticTable(startsWith(staticTable.version, 'per-frame'), :);
noiseGranule = zeros(1, numel(names)); noiseCell = zeros(1, numel(names));
realExcess = zeros(1, numel(names)); staticExcess = zeros(1, numel(names)); staticSE = zeros(1, numel(names));
for c = 1:numel(names)
    row = noiseTable(noiseTable.dataset == csvNames{c}, :);
    noiseGranule(c) = str2double(strtok(row.ratioSum));                 % "4.14 [3.90-4.52]" -> 4.14
    noiseCell(c) = str2double(strtok(row.ratioCell));
    row = staticTable(staticTable.dataset == csvNames{c}, :);
    if isempty(row)                      % the static control was run for 4PP7, 4PP7-FUS and in vitro only
        realExcess(c) = NaN; staticExcess(c) = NaN; staticSE(c) = NaN;
    else
        realExcess(c) = row.real_excess_w13;
        staticExcess(c) = row.static_excess_w13;
        staticSE(c) = row.static_excess_se_w13;
    end
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 7.2 2.6], 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact');
nexttile; hold on;
for c = 1:numel(names)
    if c <= nVivo
        style = '-o';
    else
        style = ':o';
    end
    plot(spList, ratio(c, :), style, 'Color', colors(c, :), 'MarkerSize', 3, 'DisplayName', labels{c});
end
yline(1, 'k--');
xlabel('Stringency parameter j'); ylabel('Bursts in copies / bursts in granules');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 5, 'NumColumns', 2); title('(a)');
nexttile; hold on;
for c = 1:numel(names)
    bar(c - 0.18, noiseGranule(c), 0.36, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
    bar(c + 0.18, noiseCell(c), 0.36, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
end
yline(1, 'k--');
set(gca, 'XTick', 1:numel(names), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'FontSize', 7);
ylabel('Measured / camera-predicted noise');
text(0.03, 0.95, 'coloured: granule (spot sum); grey: background', 'Units', 'normalized', 'FontSize', 6);
title('(b)');
nexttile; hold on;
for c = 1:numel(names)
    bar(c - 0.18, realExcess(c), 0.36, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
    bar(c + 0.18, staticExcess(c), 0.36, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
    errorbar(c + 0.18, staticExcess(c), 1.96 * staticSE(c), 'k', 'LineStyle', 'none', 'CapSize', 2);
end
yline(0, 'k');
set(gca, 'XTick', 1:numel(names), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'FontSize', 7);
ylabel('Turnover beyond noise (%/h)');
text(0.03, 0.95, 'coloured: granules; grey: static synthetic granule', 'Units', 'normalized', 'FontSize', 6);
title('(c)');
if saveFigures
    exportgraphics(f, 'figures/AppendixFig14.png', 'Resolution', 300);
end


%% =====================================================================================
%  APPENDIX F - DETECTION SETTINGS
%  ======================================================================================

%% Appendix F, Fig. 15 - window, stringency and the detection limit (about 20 min)
% (a) Bursts in the burst-free copies (1 per granule) / bursts in the granules, median over
% the in vivo datasets, for every window and stringency.
% (b) Negative bursts in the granules (solid) and in their copies (5 per granule, dotted) at
% high stringency (window 13).
% (c) Simulations with known bursts (lambda 2): % of the true bursts found (solid) and the
% false share of the detections (bursts found in a burst-free twin / bursts found; dotted),
% at the noise of the original simulations (0.7% of the signal) and of PP7 granules (18%).
windowList = [5 7 9 13 17 21];
spGrid = [3 4 5 6 8 10 12];
gridRatio = zeros(numel(spGrid), numel(windowList), nVivo);
for c = 1:nVivo
    rng(c);
    [C, BC] = burstFreeCopies(D(c).X, D(c).B, 1);
    for w = 1:numel(windowList)
        for s = 1:numel(spGrid)
            [~, ~, ~, ~, ~, ~, segReal] = analysis(D(c).X, D(c).B, timeVecTotal, windowList(w), spGrid(s), 0, 0);
            [~, ~, ~, ~, ~, ~, segCopy] = analysis(C, BC, timeVecTotal, windowList(w), spGrid(s), 0, 0);
            gridRatio(s, w, c) = sum(segCopy(:, 2) < 0 & segCopy(:, 6) > 0) / sum(segReal(:, 2) < 0 & segReal(:, 6) > 0);
        end
    end
    fprintf('Fig. 15(a): %s done\n', names{c});
end
gridMedian = median(gridRatio, 3);

spHigh = [10 12 14 16 18 20 22 25];
nCopies = 5;
realHigh = zeros(nVivo, numel(spHigh));
copyHigh = zeros(nVivo, numel(spHigh));
for c = 1:nVivo
    rng(c);
    [C, BC] = burstFreeCopies(D(c).X, D(c).B, nCopies);
    for s = 1:numel(spHigh)
        [~, ~, ~, ~, ~, ~, segReal] = analysis(D(c).X, D(c).B, timeVecTotal, movingAverageWindowSpan, spHigh(s), 0, 0);
        [~, ~, ~, ~, ~, ~, segCopy] = analysis(C, BC, timeVecTotal, movingAverageWindowSpan, spHigh(s), 0, 0);
        realHigh(c, s) = sum(segReal(:, 2) < 0 & segReal(:, 6) > 0);
        copyHigh(c, s) = sum(segCopy(:, 2) < 0 & segCopy(:, 6) > 0) / nCopies;
    end
end

% (c) simulations: setting 1 as the original simulations, setting 2 with PP7-like noise
simNoise = [40 1100];          % frame noise (a.u.): ~0.7% and ~18% of the 2000-10000 level
simK0 = [60 300];              % one unit (a.u.); with lambda 2, bursts of ~11% of the level for setting 2
expConst = 0.0026;             % photobleaching (1/min)
nTracks = 500;
recall = zeros(2, numel(spHigh));
falseShare = zeros(2, numel(spHigh));
for j = 1:2
    rng(3400 + j);
    [signals, backgrounds] = generateSim(nTracks, 361, 3, expConst, simNoise(j), 2, simK0(j));   % with bursts
    rng(3400 + j);
    noiseFree = generateSim(nTracks, 361, 3, expConst, 0, 2, simK0(j));                      % the same granules without noise
    rng(3450 + j);
    [flat, flatBackgrounds] = generateSim(nTracks, 361, 1, expConst, simNoise(j), 2, simK0(j));  % burst-free twin, same noise
    u = noiseFree ./ exp(-expConst * timeVecTotal);                     % without photobleaching
    [trueTrace, trueFrame] = find(diff(u, 1, 2) < -1e-9);               % every true loss: its trace and frame
    trueFrame = trueFrame + 1;
    for s = 1:numel(spHigh)
        [~, ~, ~, ~, ~, ~, seg] = analysis(signals, backgrounds, timeVecTotal, movingAverageWindowSpan, spHigh(s), 0, 0);
        [~, ~, ~, ~, ~, ~, segFlat] = analysis(flat, flatBackgrounds, timeVecTotal, movingAverageWindowSpan, spHigh(s), 0, 0);
        neg = seg(seg(:, 2) < 0, :);
        found = false(numel(trueTrace), 1);
        for q = 1:size(neg, 1)       % a true loss is found if it lies within a detected burst (+- 6 frames)
            found(trueTrace == neg(q, 1) & trueFrame >= neg(q, 3) - 6 & trueFrame <= neg(q, 4) + 6) = true;
        end
        recall(j, s) = 100 * mean(found);
        falseShare(j, s) = 100 * min(1, sum(segFlat(:, 2) < 0) / max(size(neg, 1), 1));
    end
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 7.2 2.6], 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact');
nexttile;
imagesc(gridMedian, [0.6 1]);
set(gca, 'YDir', 'normal', 'XTick', 1:numel(windowList), 'XTickLabel', windowList, 'YTick', 1:numel(spGrid), 'YTickLabel', spGrid);
for s = 1:numel(spGrid)
    for w = 1:numel(windowList)
        text(w, s, sprintf('%.2f', gridMedian(s, w)), 'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 6);
    end
end
colorbar;
xlabel('Moving-average window (frames)'); ylabel('Stringency parameter j'); title('(a) copies / granules');
nexttile; hold on;
for c = 1:nVivo
    plot(spHigh, realHigh(c, :), '-o', 'Color', colors(c, :), 'MarkerSize', 3, 'DisplayName', labels{c});
    plot(spHigh, copyHigh(c, :), ':', 'Color', colors(c, :), 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log');
xlabel('Stringency parameter j'); ylabel('Negative bursts detected');
legend('Location', 'southwest', 'Box', 'off', 'FontSize', 5);
text(0.97, 0.95, {'solid: granules', 'dotted: burst-free copies'}, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 6);
title('(b)');
nexttile; hold on;
simColors = [0 0.447 0.741; 0.835 0.369 0];
simNames = {'original noise', 'PP7 noise'};
for j = 1:2
    plot(spHigh, recall(j, :), '-o', 'Color', simColors(j, :), 'MarkerSize', 3, 'DisplayName', [simNames{j} ', bursts found']);
    plot(spHigh, falseShare(j, :), ':', 'Color', simColors(j, :), 'LineWidth', 1, 'DisplayName', [simNames{j} ', false share']);
end
ylim([0 105]);
xlabel('Stringency parameter j'); ylabel('Simulated bursts (%)');
legend('Location', 'northeast', 'Box', 'off', 'FontSize', 5); title('(c)');
if saveFigures
    exportgraphics(f, 'figures/AppendixFig15.png', 'Resolution', 300);
end


%% =====================================================================================
%  APPENDIX G - ESTIMATION OF LAMBDA AND K0
%  ======================================================================================

%% Appendix G, Fig. 16 - which lambda the fit returns on simulated data (about 10 min)
% Simulated bursty granules as in MS Sec. III (k0 = 60, noise 40 a.u., i.e. SNR = 1.5;
% window 13, SP 15), 20 simulated sets of 200 granules for each true lambda. Every set is
% fitted with Eq. (1.9) on the grid lambda 1-5, k0 1-500 and scored five ways. (The MS also used
% the likelihood; with the zero-truncated model of PoissonFit2WS.m it is not defined for bursts
% smaller than k0/2, so it is left out.)
% (a) Fraction of the sets returning each lambda with the WS. (b) Mean returned lambda for
% each score; the dashed line marks perfect recovery.
lambdaSim = [0.5 1 2 3 4 5];
scoreNames = {'WS', 'MSE', 'QQ', 'KS', 'CvM'};
nSets = 20;
lambdaGrid = 1:5;
K0Vector = 1:500;
expConst = 0.0026;
chosen = zeros(numel(lambdaSim), numel(scoreNames), nSets);     % returned lambda
for l = 1:numel(lambdaSim)
    for r = 1:nSets
        rng(5000 + 100 * l + r);
        [signals, backgrounds] = generateSim(200, 361, 3, expConst, 40, lambdaSim(l), 60);
        [~, ~, ~, ampNeg] = analysis(signals, backgrounds, timeVecTotal, 13, 15, 0, 0);
        [lamWS, ~, ~, scores] = PoissonFit2WS(abs(ampNeg(2:end)), lambdaGrid, K0Vector, 10);
        chosen(l, 1, r) = lamWS;
        for k = 2:numel(scoreNames)
            [~, best] = min(scores.(scoreNames{k})(:));
            [lamId, ~] = ind2sub(size(scores.(scoreNames{k})), best);
            chosen(l, k, r) = lambdaGrid(lamId);
        end
    end
    fprintf('Fig. 16: lambda %g done\n', lambdaSim(l));
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
fraction = zeros(numel(lambdaSim), numel(lambdaGrid));
for l = 1:numel(lambdaSim)
    for g = 1:numel(lambdaGrid)
        fraction(l, g) = mean(chosen(l, 1, :) == lambdaGrid(g));
    end
end
bar(1:numel(lambdaSim), fraction, 0.7, 'stacked');
set(gca, 'XTick', 1:numel(lambdaSim), 'XTickLabel', lambdaSim);
ylim([0 1]); xlabel('\lambda_{sim}'); ylabel('Fraction of fits (WS)');
legendText = cell(1, numel(lambdaGrid));
for g = 1:numel(lambdaGrid)
    legendText{g} = sprintf('\\lambda_{est} = %d', lambdaGrid(g));
end
legend(legendText, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 6); title('(a)');
nexttile; hold on;
for k = 1:numel(scoreNames)
    plot(lambdaSim, mean(chosen(:, k, :), 3), '-o', 'MarkerSize', 3, 'DisplayName', scoreNames{k});
end
plot([0 5.5], [0 5.5], 'k--', 'HandleVisibility', 'off');
xlim([0 5.5]); ylim([0 5.5]);
xlabel('\lambda_{sim}'); ylabel('Mean \lambda_{est}');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 6, 'NumColumns', 2); title('(b)');
if saveFigures
    exportgraphics(f, 'figures/AppendixFig16.png', 'Resolution', 300);
end


%% =====================================================================================
%  APPENDIX H - LARGE SLOW CHANGES AND TIME UNIFORMITY
%  ======================================================================================

%% Appendix H, Fig. 17 - large slow changes (about 3 min)
% Rises and drops of at least 30% within 15 min (90 frames), counted without the run rule of
% the burst detector, in the granules and in their burst-free copies (5 per granule). Signal
% as the detector sees it: moving average (13 frames) of spot sum + background, divided by
% the cubic fit of the background. Excess = granules minus copies, per granule-hour, with a
% 95% CI from resampling movies, then granules. (a) Relative to each granule's own level.
% (b) The same amount of light for every granule: 30% of the median PP7 granule level.
threshold = 0.3;
maxFrames = 90;
nCopies = 5;
nBoot = 1000;
% the absolute unit: median bleaching-corrected level of the PP7 granules (4 PP7 datasets)
levels = [];
for c = 1:4
    for k = 1:numel(D(c).X)
        x = (1:numel(D(c).B{k}))';
        backgroundFit = polyval(polyfit(x, D(c).B{k}, 3), x);
        levels(end+1, 1) = median(D(c).X{k} ./ backgroundFit * backgroundFit(1)); %#ok<SAGROW>
    end
end
unitAbsolute = median(levels);

modes = {'relative', 'absolute'};
excess = zeros(2, nVivo); excessLo = zeros(2, nVivo); excessHi = zeros(2, nVivo);
perGranule = cell(2, nVivo);                % per granule: [real events, copy events, hours]
for m = 1:2
    for c = 1:nVivo
        rng(35 + c);
        [C, ~, owner] = burstFreeCopies(D(c).X, D(c).B, nCopies);
        n = numel(D(c).X);
        counts = zeros(n, 3);
        for k = 1:n
            x = (1:numel(D(c).B{k}))';
            backgroundFit = polyval(polyfit(x, D(c).B{k}, 3), x);
            if m == 1
                scale = 1; unit = [];
            else
                scale = backgroundFit(1); unit = unitAbsolute;
            end
            [up, down] = countLargeChanges(movmean(D(c).X{k}, 13) ./ backgroundFit * scale, threshold, maxFrames, unit);
            counts(k, 1) = up + down;
            copiesOfGranule = find(owner == k);
            for q = 1:numel(copiesOfGranule)
                [up, down] = countLargeChanges(movmean(C{copiesOfGranule(q)}, 13) ./ backgroundFit * scale, threshold, maxFrames, unit);
                counts(k, 2) = counts(k, 2) + (up + down) / nCopies;
            end
            counts(k, 3) = numel(D(c).X{k}) / 360;
        end
        perGranule{m, c} = counts;
        excess(m, c) = (sum(counts(:, 1)) - sum(counts(:, 2))) / sum(counts(:, 3));
        % 95% CI: resample replicates, then granules within each chosen replicate
        replicateList = unique(D(c).replicate);
        boot = zeros(nBoot, 1);
        for b = 1:nBoot
            picks = [];
            chosenReplicates = replicateList(randi(numel(replicateList), numel(replicateList), 1));
            for q = 1:numel(chosenReplicates)
                granulesOfReplicate = find(D(c).replicate == chosenReplicates(q));
                picks = [picks; granulesOfReplicate(randi(numel(granulesOfReplicate), numel(granulesOfReplicate), 1))]; %#ok<AGROW>
            end
            boot(b) = (sum(counts(picks, 1)) - sum(counts(picks, 2))) / sum(counts(picks, 3));
        end
        excessLo(m, c) = prctile(boot, 2.5);
        excessHi(m, c) = prctile(boot, 97.5);
    end
end
% PP7-FUS vs PP7 (4x): shuffle whole replicates between the two conditions (two-sided p)
pFus = zeros(1, 2);
nPerm = 2000;
for m = 1:2
    countsBoth = [perGranule{m, 2}; perGranule{m, 1}];
    [~, ~, replicateA] = unique(D(2).replicate);
    [~, ~, replicateB] = unique(D(1).replicate);
    replicateBoth = [replicateA; replicateB + max(replicateA)];
    nA = max(replicateA);
    observed = excess(m, 2) - excess(m, 1);
    rng(1);
    countExtreme = 0;
    for p = 1:nPerm
        order = randperm(max(replicateBoth));
        inA = ismember(replicateBoth, order(1:nA));
        dA = (sum(countsBoth(inA, 1)) - sum(countsBoth(inA, 2))) / sum(countsBoth(inA, 3));
        dB = (sum(countsBoth(~inA, 1)) - sum(countsBoth(~inA, 2))) / sum(countsBoth(~inA, 3));
        if abs(dA - dB) >= abs(observed) - 1e-12
            countExtreme = countExtreme + 1;
        end
    end
    pFus(m) = (countExtreme + 1) / (nPerm + 1);
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
panelTitles = {'(a) relative to the granule', '(b) same amount of light'};
for m = 1:2
    nexttile; hold on;
    for c = 1:nVivo
        bar(c, excess(m, c), 0.65, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
        errorbar(c, excess(m, c), excess(m, c) - excessLo(m, c), excessHi(m, c) - excess(m, c), 'k', 'LineStyle', 'none', 'CapSize', 2);
    end
    yline(0, 'k');
    set(gca, 'XTick', 1:nVivo, 'XTickLabel', labels(1:nVivo), 'XTickLabelRotation', 45, 'FontSize', 7);
    ylabel({'Excess changes >= 30%', 'per granule-hour'});
    text(0.03, 0.95, sprintf('PP7-FUS vs PP7 (4x): p = %.3f', pFus(m)), 'Units', 'normalized', 'FontSize', 7);
    title(panelTitles{m});
end
if saveFigures
    exportgraphics(f, 'figures/AppendixFig17.png', 'Resolution', 300);
end

%% Appendix H, Fig. 18 - first vs second 30 min (about 3 min)
% Turnover beyond noise (a) and net change beyond photobleaching (b) in the first and second
% half of the measurements, from burstMeasures.m (20 copies per granule); * p < 0.05
% (bootstrap test of the difference between the halves).
first = zeros(2, numel(names)); second = zeros(2, numel(names)); pHalves = zeros(2, numel(names));
for c = 1:numel(names)
    M = burstMeasures(D(c).S, D(c).B, D(c).G, D(c).replicate, 20, 1, movingAverageWindowSpan, stringencyParameter, numel(timeVecTotal));
    first(:, c) = [M.turnover1; M.net1];
    second(:, c) = [M.turnover2; M.net2];
    pHalves(:, c) = [M.pHalvesTurnover; M.pHalvesNet];
end

% --- figure
f = figure('Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
yLabels = {'Turnover beyond noise (%/h)', 'Net change beyond bleaching (%/h)'};
for j = 1:2
    nexttile; hold on;
    for c = 1:numel(names)
        plot([c - 0.12, c + 0.12], [first(j, c), second(j, c)], '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
        if pHalves(j, c) < 0.05
            text(c, max(first(j, c), second(j, c)) + 4, '*', 'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    end
    plot((1:numel(names)) - 0.12, first(j, :), 'o', 'Color', [0 0.447 0.741], 'MarkerSize', 4, 'DisplayName', 'first 30 min');
    plot((1:numel(names)) + 0.12, second(j, :), 's', 'Color', [0.835 0.369 0], 'MarkerSize', 4, 'DisplayName', 'second 30 min');
    yline(0, 'k', 'HandleVisibility', 'off');
    set(gca, 'XTick', 1:numel(names), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'FontSize', 7);
    ylabel(yLabels{j}); title(['(' char('a' + j - 1) ')']);
    if j == 1
        legend('Location', 'northwest', 'Box', 'off', 'FontSize', 6);
    end
end
if saveFigures
    exportgraphics(f, 'figures/AppendixFig18.png', 'Resolution', 300);
end


%% ------------------------------------------------------------------ functions
function [C, BC, owner] = burstFreeCopies(X, B, nCopies)
% Burst-free copies of every trace, made as in burstMeasures.m: the median bleaching shape of
% the condition (cubic fit of each trace, divided by its mean) scaled to the trace, plus white
% noise at the trace's own frame-to-frame level. BC: the real background of each copy;
% owner: the granule each copy belongs to. Set the random seed (rng) before calling.
n = numel(X);
L = zeros(n, 1);
for k = 1:n
    L(k) = numel(X{k});
end
minLength = min(L);
shapes = zeros(n, minLength);
sigma = zeros(n, 1);
for k = 1:n
    x = (1:L(k))';
    trend = polyval(polyfit(x, X{k}, 3), x);
    shapes(k, :) = trend(1:minLength)' / mean(trend(1:minLength));
    d = diff(X{k});
    sigma(k) = 1.4826 * median(abs(d - median(d))) / sqrt(2);
end
shape = median(shapes, 1);
C = cell(n * nCopies, 1);
BC = cell(n * nCopies, 1);
owner = zeros(n * nCopies, 1);
q = 0;
for c = 1:nCopies
    for k = 1:n
        sh = [shape, repmat(shape(end), 1, L(k) - minLength)]';
        q = q + 1;
        C{q} = mean(X{k}) / mean(sh) * sh + sigma(k) * randn(L(k), 1);
        BC{q} = B{k};
        owner(q) = k;
    end
end
end


function [nRise, nDrop] = countLargeChanges(x, threshold, maxFrames, unit)
% Rises and drops of at least `threshold` within `maxFrames` frames, counted without overlap:
% scanning forward, the first frame that reaches the threshold ends the event and the scan
% goes on from there. unit = []: the change relative to the starting value; otherwise the
% change in units of `unit` (the same amount of light for every granule).
nRise = 0;
nDrop = 0;
t = 1;
n = numel(x);
while t < n
    segment = x(t + 1:min(n, t + maxFrames));
    if isempty(unit)
        change = segment / x(t) - 1;
    else
        change = (segment - x(t)) / unit;
    end
    iRise = find(change >= threshold, 1);
    iDrop = find(change <= -threshold, 1);
    if isempty(iRise) && isempty(iDrop)
        t = t + 1;
    elseif isempty(iDrop) || (~isempty(iRise) && iRise < iDrop)
        nRise = nRise + 1;
        t = t + iRise;
    else
        nDrop = nDrop + 1;
        t = t + iDrop;
    end
end
end
