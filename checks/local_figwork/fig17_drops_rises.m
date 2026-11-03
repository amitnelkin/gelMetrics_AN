% Fig. 17 split into drops and rises (section 76). Local copy of AppendixFigures.m: the settings section
% (movies as replicates) and the Fig. 17 section counting drops and rises separately; saves the staged
% .fig files and numbers that compose_ms_figs (local) turns into Fig17_drops / Fig17_rises. Run from the repo root.
%% Settings and data - run this section first (a few seconds)
close all; clearvars; clc;

% Algorithm Parameters (the same as in ExperimentalAnalysis.m)
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 12;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 361);      % Time steps vector: 361 frames, 60 min at 10 s per frame
saveFigures = false;                      % true: save every figure as PNG in figures/
moviesAsReplicates = true;

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
directions = {'drops', 'rises'};
for d = 1:2
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
            if d == 1, counts(k, 1) = down; else, counts(k, 1) = up; end
            copiesOfGranule = find(owner == k);
            for q = 1:numel(copiesOfGranule)
                [up, down] = countLargeChanges(movmean(C{copiesOfGranule(q)}, 13) ./ backgroundFit * scale, threshold, maxFrames, unit);
                if d == 1, counts(k, 2) = counts(k, 2) + down / nCopies; else, counts(k, 2) = counts(k, 2) + up / nCopies; end
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
R17.excess{d} = excess; R17.lo{d} = excessLo; R17.hi{d} = excessHi; R17.p{d} = pFus; R17.perGranule{d} = perGranule;
f = figure(d); set(f, 'Units', 'inches', 'Position', [1 1 6.5 2.6], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
panelTitles = {[directions{d} ', relative to the granule'], [directions{d} ', absolute']};
for m = 1:2
    nexttile; hold on;
    for c = 1:nVivo
        bar(c, excess(m, c), 0.65, 'FaceColor', colors(c, :), 'EdgeColor', 'none');
        errorbar(c, excess(m, c), excess(m, c) - excessLo(m, c), excessHi(m, c) - excess(m, c), 'k', 'LineStyle', 'none', 'CapSize', 2);
    end
    yline(0, 'k');
    set(gca, 'XTick', 1:nVivo, 'XTickLabel', labels(1:nVivo), 'XTickLabelRotation', 45, 'FontSize', 7);
    ylabel({['Excess ' directions{d} ' >= 30%'], 'per granule-hour'});
    if pFus(m) < 0.001, pText = 'p < 0.001'; else, pText = sprintf('p = %.3f', pFus(m)); end
    text(0.03, 0.95, ['PP7-FUS vs PP7 (4x): ' pText], 'Units', 'normalized', 'FontSize', 7);
    title(panelTitles{m});
end
end
for d = 1:2, savefig(figure(d), fullfile('C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN\figures_revision_local\stage', sprintf('loc_af17_%d.fig', d))); end
save(fullfile('C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN\figures_revision_local\stage', 'loc_af17_vars.mat'), 'R17', 'unitAbsolute', 'names');

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
