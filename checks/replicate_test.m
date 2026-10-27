% REPLICATE_TEST  Movies or granules as the replicates? (RESULTS.md section 63)
% 1. Between-movie share of the variation (intraclass correlation, ICC) of the per-granule
%    values: brightness, frame noise, net change, median burst size, turnover beyond noise.
% 2. 95% CIs and comparison p values with movies as replicates (as in the paper) and with
%    every granule as its own replicate. Same data, same copies (20 per granule, seed 1);
%    only the resampling unit differs. Uses checks/burstMeasures_ref.m (copy of the committed
%    burstMeasures.m). Writes analysis_outputs/replicates/{icc,ci,p}.csv. No plotting.
cd(fileparts(fileparts(mfilename('fullpath')))); addpath('checks');
outDir = fullfile('analysis_outputs', 'replicates');
if ~exist(outDir, 'dir'), mkdir(outDir); end
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'};
ciMeasures = {'sizeBeyond', 'turnover', 'net', 'medianBurst', 'meanBurst'};
Mmov = cell(1, numel(names)); Mgra = cell(1, numel(names));
ICC = table(); CI = table();
for d = 1:numel(names)
    [S, B, G, mv] = loadData(names{d});
    Mmov{d} = burstMeasures_ref(S, B, G, mv, 20, 1);      % movies as replicates
    Mgra{d} = burstMeasures_ref(S, B, G, [], 20, 1);      % every granule its own replicate
    Mgra{d}.movieReal = Mmov{d}.movie;                    % (kept for reference)
    M = Mmov{d}; n = M.nGranules;
    % per-granule values
    medBurst = nan(n, 1); turn = zeros(n, 1);
    for k = 1:n
        if ~isempty(M.fracReal{k}), medBurst(k) = 100 * median(M.fracReal{k}); end
        turn(k) = 100 * (sum(M.fracReal{k}) - sum(M.fracCopy{k}) / M.nCopies) / M.hours(k);
    end
    values = {M.level, 100 * M.noise, M.granuleNet, medBurst, turn};
    labels = {'brightness', 'frame noise', 'net change', 'median burst size', 'turnover beyond noise'};
    for j = 1:numel(values)
        icc = iccOneWay(values{j}, M.movie);
        ICC = [ICC; table(string(names{d}), string(labels{j}), icc, M.nMovies, n, ...
            'VariableNames', {'dataset', 'measure', 'ICC', 'movies', 'granules'})]; %#ok<AGROW>
    end
    for j = 1:numel(ciMeasures)
        f = ciMeasures{j};
        CI = [CI; table(string(names{d}), string(f), M.(f), M.([f '_lo']), M.([f '_hi']), Mgra{d}.([f '_lo']), Mgra{d}.([f '_hi']), ...
            'VariableNames', {'dataset', 'measure', 'value', 'movie_lo', 'movie_hi', 'granule_lo', 'granule_hi'})]; %#ok<AGROW>
    end
    fprintf('%s done\n', names{d});
end
writetable(ICC, fullfile(outDir, 'icc.csv'));
writetable(CI, fullfile(outDir, 'ci.csv'));

% comparisons of the paper: A vs B, movie-level and granule-level shuffling
pairs = {'4PP7-FUS', '4PP7'; '30PP7-FUS', '30PP7'; '5Qb', '4PP7'; '10Qb', '30PP7'; 'vitro_10_1', 'vitro_1_1'; 'vitro_100_1', 'vitro_1_1'};
P = table();
for q = 1:size(pairs, 1)
    a = find(strcmp(names, pairs{q, 1})); b = find(strcmp(names, pairs{q, 2}));
    Tm = compareConditions(Mmov{a}, Mmov{b}, pairs{q, 1}, pairs{q, 2}, 10000, 1);
    Tg = compareConditions(Mgra{a}, Mgra{b}, pairs{q, 1}, pairs{q, 2}, 10000, 1);
    for j = 1:height(Tm)
        P = [P; table(string(pairs{q, 1}) + " vs " + string(pairs{q, 2}), string(Tm{j, 1}), Tm.p(j), Tg.p(j), Tm.p_AgtB(j), Tg.p_AgtB(j), ...
            'VariableNames', {'comparison', 'measure', 'p_movies', 'p_granules', 'p1_movies', 'p1_granules'})]; %#ok<AGROW>
    end
end
writetable(P, fullfile(outDir, 'p.csv'));


function icc = iccOneWay(y, group)
% between-group share of the variation, one-way ANOVA estimate (ICC(1)); NaN values left out
ok = ~isnan(y); y = y(ok); group = group(ok);
g = unique(group); k = numel(g); N = numel(y);
grand = mean(y); ssb = 0; ssw = 0; sizes = zeros(k, 1);
for i = 1:k
    yi = y(group == g(i)); sizes(i) = numel(yi);
    ssb = ssb + sizes(i) * (mean(yi) - grand) ^ 2;
    ssw = ssw + sum((yi - mean(yi)) .^ 2);
end
msb = ssb / (k - 1); msw = ssw / (N - k);
n0 = (N - sum(sizes .^ 2) / N) / (k - 1);               % average group size (unequal groups)
icc = (msb - msw) / (msb + (n0 - 1) * msw);
end


function [S, B, G, mv] = loadData(name)
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/'; file = strrep(name, 'vitro', 'same');
    info = readtable([folder 'extraction_info_same.csv'], 'Delimiter', ',', 'TextType', 'string'); info = info(info.dataset == file, :);
    S = loadTraces([folder file '_data.mat']); B = loadTraces([folder file '_background.mat']); G = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    info = [readtable([folder 'extraction_info.csv'], 'TextType', 'string'); readtable([folder 'extraction_info_5Qb_10Qb.csv'], 'TextType', 'string')];
    info = info(info.dataset == name, :);
    S = loadTraces([folder name '_data.mat']); B = loadTraces([folder name '_background.mat']);
    G = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
mv = repelem((1:height(info))', info.nKept);
end


function traces = loadTraces(file)
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end
