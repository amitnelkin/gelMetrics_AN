close all; clearvars; clc;

% Is condition A different from condition B?
% For every pair: the measures of A and B, their difference, and the p value from shuffling
% whole replicates between the two conditions: two-sided test and one-sided test, for A > B, only when. 
% The figure shows every condition (95% CI) and a bracket with the p value over every tested pair with p < 0.05.

% Algorithm Parameters (the same as in ExperimentalAnalysis.m)
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 12;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 361);      % 361 frames, 60 min at 10 s per frame
nCopies = 20;                             % Burst-free copies per granule (the noise floor)
seed = 1;                                 % Random seed of the copies and of the shuffles
nPerm = 10000;                            % Number of shuffles
moviesAsReplicates = false;               % false: every granule is its own replicate (true: shuffles movies, not granules)

% Okabe-Ito pallette
% ORANGE = '#E69F00', GREEN = '#009E73', BLUE = '#0072B2', YELLOW = '#F0E442', RED = '#D55E00', PURPLE = '#CC79A7'
OKI = ["#E69F00","#009E73","#0072B2","#F0E442","#D55E00","#CC79A7"];

% Pairs to compare, one per row: A (first column) against B (second column).
% Every condition is a folder in dataFolder with data.mat, background.mat, spotsum.mat and movies.mat.
% in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'   (in ./in_vivo_experiments/)
% in vitro: 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'   (in ./in_vitro_experiments/)
dataFolder = './in_vivo_experiments/';
pairs = {'4PP7-FUS', '4PP7'};
conditions = unique(pairs(:), 'stable');  % the order of the bars; any order that holds every name in pairs

% sizeGranules: median burst as a fraction of granule level
% sizeCopies:   median burst as a fraction of copy level
% sizeBeyond:   sizeGranules - sizeCopies (what isn't explained by noise)
% turnover:     (sum of bursts) - (sum of copy bursts) per granule per time [total material lost/gained per hour]
% net:          granule fade rate relative to the background 
% medianBurst:  median of burst amplitudes
% meanBurst:    mean of burst amplitudes [lambda * k0]
% brightness:   median of brightness (intensity)  over replicates
measures = {'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'net', 'medianBurst', 'meanBurst', 'brightness'};

%% Measures of every condition
allM = cell(1, numel(conditions));
for c = 1:numel(conditions)
    folder = [dataFolder conditions{c} '/'];
    signals = struct2cell(load([folder 'data.mat']));       
    signals = signals{1}(2:end);       % the first cell is empty
    backgrounds = struct2cell(load([folder 'background.mat'])); 
    backgrounds = backgrounds{1}(2:end);
    spotSums = struct2cell(load([folder 'spotsum.mat']));    
    spotSums = spotSums{1}(2:end);
    movies = load([folder 'movies.mat']);
    if moviesAsReplicates
        replicateOfTrace = movies.movieId;
    else
        replicateOfTrace = [];
    end
    allM{c} = burstMeasures(signals, backgrounds, spotSums, replicateOfTrace, nCopies, seed, ...
        movingAverageWindowSpan, stringencyParameter, numel(timeVecTotal));
end

%% Compare every pair: shuffle whole replicates between A and B
results = cell(size(pairs, 1), 1);
perGranule = {'level', 'noise', 'hours', 'hours1', 'hours2', 'granuleNet', 'granuleNet1', 'granuleNet2', ...
              'ampAU', 'fracReal', 'startReal', 'fracCopy', 'startCopy'}; % field required by summarizeGranules.m
for q = 1:size(pairs, 1)
    MA = allM{strcmp(conditions, pairs{q, 1})};
    MB = allM{strcmp(conditions, pairs{q, 2})};
    strA = pairs{q, 1};
    strB = pairs{q, 2};
    G.nCopies = MA.nCopies;   % This variable is passed to summarizeGranules.m and requires nCopies (assuming MA.nCopies = MB.nCopies) 
    for j = 1:numel(perGranule) % one list of granules
        G.(perGranule{j}) = [MA.(perGranule{j}); MB.(perGranule{j})];
    end
    [~, ~, repA] = unique(MA.replicate);
    [~, ~, repB] = unique(MB.replicate);
    nA = max(repA);
    G.replicate = [repA; repB + nA]; % replicates renumbered 1..nA (A) and nA+1.. (B)
    
    VA = summarizeGranules(MA, []);
    VB = summarizeGranules(MB, []);
    
    A = zeros(numel(measures), 1); 
    B = zeros(numel(measures), 1); 
    countOneSide = zeros(numel(measures), 1);
    countTwoSide = zeros(numel(measures), 1); 
    
    for j = 1:numel(measures)
        A(j) = VA.(measures{j});
        B(j) = VB.(measures{j});
    end

    difference = A - B;
    rng(seed);
    
    for i = 1:nPerm
        order = randperm(max(G.replicate));
        V1 = summarizeGranules(G, find(ismember(G.replicate, order(1:nA))));      % random split of A condition
        V2 = summarizeGranules(G, find(ismember(G.replicate, order(nA+1:end))));  % random split of B condition
        for j = 1:numel(measures)
            d = V1.(measures{j}) - V2.(measures{j});
            countTwoSide(j) = countTwoSide(j) + (abs(d) >= abs(difference(j)) - 1e-12);  % is the random split difference larger/smaller than the original difference [two-sided test]
            countOneSide(j) = countOneSide(j) + (d >= difference(j) - 1e-12);            % is the random split difference larger than the original difference [one-sided test of A larger than B]
        end
    end
    p = (countTwoSide + 1) / (nPerm + 1);      % two-sided p_value
    p_A_B = (countOneSide + 1) / (nPerm + 1);  % one-sided p_value
    results{q} = table(measures', A, B, difference, p, p_A_B, 'VariableNames', {'measure', strA, strB, 'difference', '2-side-p', '1-side-p'});
    fprintf('\n%s (n=%d) vs %s (n=%d), %d shuffles\n', strA, nA, strB, max(repB), nPerm);
    disp(results{q});
end

%% Figure: every condition (95% CI) and the significant pairs
plotMeasures = {'meanBurst', 'sizeBeyond', 'turnover', 'net'};
plotLabels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Median burst beyond noise (a.u.)', 'Turnover beyond noise (%/h)', 'Net change beyond bleaching (%/h)'};
figure('Units', 'inches', 'Position', [1 1 7 2.8], 'Color', 'w');
tiledlayout(1, numel(plotMeasures), 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(plotMeasures)
    nexttile; hold on;
    top = 0; bottom = 0;
    for c = 1:numel(conditions)
        y = allM{c}.(plotMeasures{j}); 
        lo = allM{c}.([plotMeasures{j} '_lo']); 
        hi = allM{c}.([plotMeasures{j} '_hi']);
        bar(c, y, 0.6, 'FaceColor', OKI(c), 'EdgeColor', 'none');
        errorbar(c, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'CapSize', 3);
        top = max(top, hi); bottom = min(bottom, lo);
    end
    yline(0, 'k');
    step = 0.08 * (top - bottom);         % space between brackets
    level = top + step;
    for q = 1:size(pairs, 1)
        pValueTwo = results{q}.('2-side-p')(strcmp(results{q}.measure, plotMeasures{j}));
        pValueOne = results{q}.('1-side-p')(strcmp(results{q}.measure, plotMeasures{j}));
        a = find(strcmp(conditions, pairs{q, 1}));
        b = find(strcmp(conditions, pairs{q, 2}));
        if pValueTwo < 0.05
            plot([a a b b], [level - step / 3, level, level, level - step / 3], 'k');
            text((a + b) / 2, level, sprintf('p = %.3f', pValueTwo), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
            level = level + step;         % the next bracket goes higher
        elseif pValueOne < 0.05
            plot([a a b b], [level - step / 3, level, level, level - step / 3], '--k');
            text((a + b) / 2, level, sprintf('p = %.3f', pValueOne), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        end
    end
    ylim([min(1.1 * bottom, 0), level + step]);
    xlim([0.4, numel(conditions) + 0.6]);
    set(gca, 'XTick', 1:numel(conditions), 'XTickLabel', conditions, 'TickLabelInterpreter', 'none', 'Box', 'on', 'FontSize', 8);
    ylabel(plotLabels{j});
end

%% Figure: 
% median burst size of A and B
% the difference split into noise and beyond noise

partColors = [204/255, 121/255, 167/255; 0/255, 158/255, 115/255];               % noise purple, beyond noise green
figure('Units', 'inches', 'Position', [1 1 7 2.8], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axA = nexttile; hold on;
axB = nexttile; hold on;
tickPos = zeros(1, size(pairs, 1));
for q = 1:size(pairs, 1)
    MA = allM{strcmp(conditions, pairs{q, 1})};
    MB = allM{strcmp(conditions, pairs{q, 2})};
    x = (q - 1) * 3 + [1 2];                               % two bars per pair, a gap between pairs
    tickPos(q) = mean(x);

    % median bursts
    levels = [MA.sizeGranules, MB.sizeGranules];
    lo = [MA.sizeGranules_lo, MB.sizeGranules_lo];
    hi = [MA.sizeGranules_hi, MB.sizeGranules_hi];
        hA = bar(axA, x(1), levels(1), 0.8, 'FaceColor', [0 114/255 178/255]);      % A blue
    hB = bar(axA, x(2), levels(2), 0.8, 'FaceColor', [230/255 159/255 0]);      % B orange
    errorbar(axA, x, levels, levels - lo, hi - levels, 'k', 'LineStyle', 'none', 'CapSize', 3);

    % the difference 
    total = MA.sizeGranules - MB.sizeGranules;
    noise = MA.sizeCopies - MB.sizeCopies;
    beyond = MA.sizeBeyond - MB.sizeBeyond;                % = total - noise
    ciTotal = prctile(MA.boot.sizeGranules - MB.boot.sizeGranules, [2.5 97.5]);   % each condition resampled on its own
    ciNoise = prctile(MA.boot.sizeCopies - MB.boot.sizeCopies, [2.5 97.5]);
    ciBeyond = prctile(MA.boot.sizeBeyond - MB.boot.sizeBeyond, [2.5 97.5]);
    hTotal = bar(axB, x(1), total, 0.8, 'FaceColor', [213/255, 94/255, 0]);
    errorbar(axB, x(1), total, total - ciTotal(1), ciTotal(2) - total, 'k', 'LineStyle', 'none', 'CapSize', 3);
    hParts = bar(axB, x(2), [noise beyond], 0.8, 'stacked');
    hParts(1).FaceColor = partColors(1, :);
    hParts(2).FaceColor = partColors(2, :);
    errorbar(axB, x(2) - 0.15, noise, noise - ciNoise(1), ciNoise(2) - noise, 'k', 'LineStyle', 'none', 'CapSize', 3);                  % noise uncertainty
    errorbar(axB, x(2) + 0.15, noise + beyond, beyond - ciBeyond(1), ciBeyond(2) - beyond, 'k', 'LineStyle', 'none', 'CapSize', 3);     % beyond-noise uncertainty
    isBeyond = strcmp(results{q}.measure, 'sizeBeyond');
    text(axB, tickPos(q), max(ciTotal(2), noise + beyond) + 0.5, ...
        sprintf('beyond noise %.1f [%.1f, %.1f]\np = %.3f (one-sided)\np = %.3f (two-sided)', beyond, ciBeyond(1), ciBeyond(2), ...
        results{q}.('1-side-p')(isBeyond), results{q}.('2-side-p')(isBeyond)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
end
set(axA, 'XTick', tickPos, 'Box', 'on', 'FontSize', 8);
ylabel(axA, 'Median burst size (% of granule)');
title(axA, '(a)');
yline(axB, 0, 'k');
set(axB, 'XTick', tickPos, 'Box', 'on', 'FontSize', 8);
legend(axA, [hA hB], {strjoin(pairs(:, 1)', ', '), strjoin(pairs(:, 2)', ', ')}, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7, 'Interpreter', 'none');
ylabel(axB, {'A - B, median burst', 'size (% points)'});
legend(axB, [hTotal hParts], {'total', 'explained by noise', 'beyond noise'}, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 7);
title(axB, '(b)');
ylim(axA, [0, max(ylim(axA)) + 3]);
ylim(axB, [min(0, min(ylim(axB))), max(ylim(axB)) + 4]);