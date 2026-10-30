close all;   % Close all open figure windows
clearvars;       % Clear workspace
clc;

% Algorithm Parameters
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 12;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 361);      % Time steps vector: 361 frames, 60 min at 10 s per frame
nCopies = 20;                             % Burst-free copies per granule (the noise floor)
seed = 1;                                 % Random seed of the copies

% Load signal and background data
% Every condition is a folder with the four files written by Collect_intensity_data.m:
% data.mat, background.mat, spotsum.mat and movies.mat.
% in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
% in vitro: 'protein_RNA_1_1', 'vitro_10_1', 'vitro_100_1'
% For your own data, set dataFolder to the folder that holds one folder per condition.
% This script only reports each condition; to test conditions against each other use compareConditions.m
dataFolder = './in_vivo_experiments/';
% dataFolder = './in_vitro_experiments/';
conditions = {'4PP7', '4PP7-FUS', '30PP7'};


% Signal plotting controls
% To plot a few sample signals you can change these variables
% plotMin is the first signal plotted
% plotMax is the last signal plotted
% Example: to plot the first 5 signals change plotMax to 5
plotMin = 0;
plotMax = 0;

nCond = numel(conditions);
M = cell(1, nCond);
for c = 1:nCond
    folder = [dataFolder conditions{c} '/'];
    signals = struct2cell(load([folder 'data.mat']));
    backgrounds = struct2cell(load([folder 'background.mat']));
    spotSums = struct2cell(load([folder 'spotsum.mat']));
    movies = load([folder 'movies.mat']);

    signals = signals{1};
    backgrounds = backgrounds{1};
    spotSums = spotSums{1};

    % Remove the empty first cell left by the extraction
    if ~iscolumn(signals)
        signals = signals(2:end)';
    else
        signals = signals(2:end);
    end

    if ~iscolumn(backgrounds)
        backgrounds = backgrounds(2:end)';
    else
        backgrounds = backgrounds(2:end);
    end

    if ~iscolumn(spotSums)
        spotSums = spotSums(2:end)';
    else
        spotSums = spotSums(2:end);
    end

    % Perform analysis
    [timeBetweenPos, timeBetweenNeg, ampPos, ampNeg, ampVar, eventsPerSignal] = ...
            analysis(signals, backgrounds, timeVecTotal, ...
            movingAverageWindowSpan, stringencyParameter, plotMin, plotMax);

    % Remove the -Inf placeholder used to initialize the variable
    ampPos = ampPos(2:end);
    ampNeg = ampNeg(2:end);
    ampVar = ampVar(2:end);

    fprintf('\n%s\n', conditions{c});
    fprintf('  • Found %d positive burst events: \n', length(ampPos));
    fprintf('  • Found %d unclassified events: \n', length(ampVar));
    fprintf('  • Found %d negative burst events: \n', length(ampNeg));

    timeBetweenPos = timeBetweenPos(2:end);
    timeBetweenNeg = timeBetweenNeg(2:end);

    figure('Name', conditions{c});
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    nexttile;
    histogram(ampVar(abs(ampVar)<1500),50);
    ylabel('Counts');
    xlabel('Intesity [A.U]');
    title('Quiescent Amplitudes');

    % Plot histogram of all amplitudes
    nexttile;
    histogram(ampVar(abs(ampVar)<3000),'FaceColor','b','BinWidth',100);
    hold on;
    histogram(ampPos(ampPos<3000),'FaceColor','g','BinWidth',100);
    histogram(ampNeg(abs(ampNeg)<3000),'FaceColor','r','BinWidth',100);
    ylabel('Counts'); xlim([-3000,3000])
    title('Segment amplitudes','FontSize',12);
    xlabel('Amplitudes [A.U]','FontSize',12);

    % Plot boxplots of durations between burst events
    nexttile;
    plot_var = [timeBetweenPos;timeBetweenNeg];
    group_var = zeros(length(plot_var),1);
    group_var(1:length(timeBetweenPos))=2;
    group_var(length(timeBetweenPos)+1:end)=3;
    boxplot(plot_var,group_var);
    xticklabels({'Positive','Negative'});
    ylabel('Time [min]');
    title('Time between burst events');
    sgtitle(conditions{c}, 'Interpreter', 'none');

    % Noise-controlled measures: the bursts of the granules against the bursts of their
    % burst-free copies, with 95% confidence intervals over movies (see burstMeasures.m).
    % The mean burst size is lambda*k0; lambda and k0 cannot be estimated separately (MS Fig. 5).
    M{c} = burstMeasures(signals, backgrounds, spotSums, movies.movieId, nCopies, seed, ...
        movingAverageWindowSpan, stringencyParameter, numel(timeVecTotal));
    fprintf('  %d granules (at least 30 min), %d movies, %d negative bursts\n', M{c}.nGranules, M{c}.nMovies, M{c}.nBursts);
    fprintf('  brightness (median over movies)      %.0f\n', M{c}.brightness);
    fprintf('  frame noise (%% of granule)           %.1f\n', M{c}.frameNoise);
    fprintf('  mean burst size, lambda*k0 (a.u.)    %.0f [%.0f, %.0f]\n', M{c}.meanBurst, M{c}.meanBurst_lo, M{c}.meanBurst_hi);
    fprintf('  burst size, granules / copies (%%)    %.1f / %.1f\n', M{c}.sizeGranules, M{c}.sizeCopies);
    fprintf('  burst size beyond noise (%% points)   %.1f [%.1f, %.1f]\n', M{c}.sizeBeyond, M{c}.sizeBeyond_lo, M{c}.sizeBeyond_hi);
    fprintf('  turnover beyond noise (%%/h)          %.1f [%.1f, %.1f]\n', M{c}.turnover, M{c}.turnover_lo, M{c}.turnover_hi);
    fprintf('  net change beyond bleaching (%%/h)    %.1f [%.1f, %.1f]\n', M{c}.net, M{c}.net_lo, M{c}.net_hi);
    fprintf('  first vs second 30 min: turnover %.1f -> %.1f (p = %.3f), net %.1f -> %.1f (p = %.3f)\n', ...
        M{c}.turnover1, M{c}.turnover2, M{c}.pHalvesTurnover, M{c}.net1, M{c}.net2, M{c}.pHalvesNet);
end

% Burst size distributions: granules (red bars) and their burst-free copies (black outline)
figure('Name', 'Burst sizes', 'Units', 'inches', 'Position', [1 1 2.4 * min(nCond, 3) 2.3 * ceil(nCond / 3)], 'Color', 'w');
tiledlayout(ceil(nCond / 3), min(nCond, 3), 'TileSpacing', 'compact', 'Padding', 'compact');
edges = 0:2.5:60;
for c = 1:nCond
    nexttile; hold on;
    fracReal = 100 * vertcat(M{c}.fracReal{:});
    fracCopy = 100 * vertcat(M{c}.fracCopy{:});
    histogram(min(fracReal, 59.99), edges, 'Normalization', 'probability', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
    histogram(min(fracCopy, 59.99), edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
    xlim([0 60]);
    xlabel('Burst size (% of granule)');
    if c == 1
        ylabel('Fraction of bursts');
        legend({'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 7);
    end
    title(sprintf('%s  (n = %d)', conditions{c}, numel(fracReal)), 'Interpreter', 'none', 'FontSize', 8);
    text(0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
        M{c}.sizeGranules, M{c}.sizeCopies, M{c}.turnover, M{c}.turnover_lo, M{c}.turnover_hi), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 7);
    set(gca, 'Box', 'on', 'TickDir', 'in', 'FontSize', 8);
end

% Summary of all conditions (95% CI over movies, then granules)
figure('Name', 'Summary', 'Units', 'inches', 'Position', [1 1 7 2.6], 'Color', 'w');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
measures = {'meanBurst', 'turnover', 'net'};
labels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', 'Net change beyond bleaching (%/h)'};
colors = lines(nCond);
for j = 1:3
    nexttile; hold on;
    for c = 1:nCond
        y = M{c}.(measures{j});
        lo = M{c}.([measures{j} '_lo']);
        hi = M{c}.([measures{j} '_hi']);
        bar(c, y, 0.6, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
        errorbar(c, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    end
    yline(0, 'k', 'LineWidth', 0.5);
    set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'TickLabelInterpreter', 'none', 'Box', 'on', 'TickDir', 'in', 'FontSize', 8);
    xlim([0.4, nCond + 0.6]);
    ylabel(labels{j});
end
