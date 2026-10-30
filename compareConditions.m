close all;   % Close all open figure windows
clearvars;       % Clear workspace
clc;

% Is condition A different from condition B?
% For every pair: the measures of A and of B, the difference A - B, and p values from
% shuffling whole movies between the two conditions (movies are the independent replicates):
%   p        two-sided
%   p_AgtB   one-sided, for A > B (use it only when that direction was stated in advance)
% The measures of each condition alone (with confidence intervals) are in ExperimentalAnalysis.m.

% Algorithm Parameters (the same as in ExperimentalAnalysis.m)
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 12;                 % Parameter for burst events detection (the value used in the paper)
timeVecTotal = linspace(0, 60, 361);      % Time steps vector: 361 frames, 60 min at 10 s per frame
nCopies = 20;                             % Burst-free copies per granule (the noise floor)
seed = 1;                                 % Random seed of the copies and of the shuffles
nPerm = 10000;                            % Number of movie shuffles

% Conditions to compare, one pair per row: A (first column) against B (second column).
% Every condition is a folder in dataFolder with the files written by Collect_intensity_data.m.
% The datasets of the paper are in ./paper_data/:
%   in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%   in vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'  (protein:RNA)
dataFolder = './paper_data/';
pairs = {'4PP7-FUS', '4PP7';
         '30PP7-FUS', '30PP7'};

measures = {'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'net', 'medianBurst', 'meanBurst', 'brightness', 'brightnessPooled'};


% Measures of every condition that appears in a pair
conditions = unique(pairs(:), 'stable');
M = cell(1, numel(conditions));
for c = 1:numel(conditions)
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

    M{c} = burstMeasures(signals, backgrounds, spotSums, movies.movieId, nCopies, seed, ...
        movingAverageWindowSpan, stringencyParameter, numel(timeVecTotal));
    fprintf('%s: %d granules, %d movies\n', conditions{c}, M{c}.nGranules, M{c}.nMovies);
end

% Compare every pair
nM = numel(measures);
results = cell(size(pairs, 1), 1);
for q = 1:size(pairs, 1)
    MA = M{strcmp(conditions, pairs{q, 1})};
    MB = M{strcmp(conditions, pairs{q, 2})};

    % Put both conditions into one list of granules; movies renumbered 1..nA (A) and nA+1.. (B)
    [~, ~, movieA] = unique(MA.movie);
    [~, ~, movieB] = unique(MB.movie);
    nA = max(movieA);
    nTotal = nA + max(movieB);
    J = struct();
    J.movie = [movieA(:); movieB(:) + nA];
    J.nGranules = MA.nGranules + MB.nGranules;
    J.nCopies = MA.nCopies;
    perGranule = {'level', 'noise', 'hours', 'hours1', 'hours2', 'granuleNet', 'granuleNet1', 'granuleNet2', ...
                  'ampAU', 'fracReal', 'startReal', 'fracCopy', 'startCopy'};
    for j = 1:numel(perGranule)
        J.(perGranule{j}) = [MA.(perGranule{j}); MB.(perGranule{j})];
    end

    % Observed values
    VA = summarizeGranules(MA, 1:MA.nGranules);
    VB = summarizeGranules(MB, 1:MB.nGranules);
    valueA = zeros(nM, 1);
    valueB = zeros(nM, 1);
    for j = 1:nM
        valueA(j) = VA.(measures{j});
        valueB(j) = VB.(measures{j});
    end
    observed = valueA - valueB;

    % Shuffle whole movies between the two conditions
    rng(seed);
    countTwo = zeros(nM, 1);
    countOne = zeros(nM, 1);
    for p = 1:nPerm
        order = randperm(nTotal);
        picksA = find(ismember(J.movie, order(1:nA)));          % granules of the movies now called "A"
        picksB = find(ismember(J.movie, order(nA+1:end)));
        V1 = summarizeGranules(J, picksA);
        V2 = summarizeGranules(J, picksB);
        for j = 1:nM
            d = V1.(measures{j}) - V2.(measures{j});
            if abs(d) >= abs(observed(j)) - 1e-12
                countTwo(j) = countTwo(j) + 1;
            end
            if d >= observed(j) - 1e-12
                countOne(j) = countOne(j) + 1;
            end
        end
    end
    p = (countTwo + 1) / (nPerm + 1);
    p_AgtB = (countOne + 1) / (nPerm + 1);

    results{q} = table(measures', valueA, valueB, observed, p, p_AgtB, ...
        'VariableNames', {'measure', 'A', 'B', 'difference', 'p', 'p_AgtB'});
    fprintf('\n%s (A) vs %s (B): %d vs %d movies, p from %d movie shuffles\n', pairs{q, 1}, pairs{q, 2}, nA, nTotal - nA, nPerm);
    disp(results{q});
end
