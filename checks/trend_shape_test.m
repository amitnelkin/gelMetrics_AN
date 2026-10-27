% TREND_SHAPE_TEST  Does the shared bleaching shape of the burst-free copies change if each
% trace's slow trend comes from trenddecomp (default settings) instead of a cubic fit?
% (RESULTS.md section 60.) Same steps as burstMeasures.m: traces >= 180 frames, X = spot sum
% + background, trend per trace, normalised to its mean over the common length, median over
% traces. Writes analysis_outputs/trend_shape/shapes.csv and summary.csv. No plotting.
cd(fileparts(fileparts(mfilename('fullpath'))));
outDir = fullfile('analysis_outputs', 'trend_shape');
if ~exist(outDir, 'dir'), mkdir(outDir); end
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'};
summary = table();
allShapes = table();
for d = 1:numel(names)
    [S, B, G] = loadData(names{d});
    keep = [];
    for i = 1:numel(S)
        if numel(S{i}) >= 180, keep(end+1) = i; end %#ok<AGROW>
    end
    n = numel(keep); L = zeros(n, 1); X = cell(n, 1);
    for k = 1:n
        i = keep(k);
        X{k} = double(G{i}(:)) + double(B{i}(:));
        L(k) = numel(X{k});
    end
    Lmin = min(L);
    shapesPoly = zeros(n, Lmin); shapesTD = zeros(n, Lmin); perTraceDiff = zeros(n, 1);
    for k = 1:n
        x = (1:L(k))';
        trendPoly = polyval(polyfit(x, X{k}, 3), x);
        trendTD = trenddecomp(X{k});
        shapesPoly(k, :) = trendPoly(1:Lmin)' / mean(trendPoly(1:Lmin));
        shapesTD(k, :) = trendTD(1:Lmin)' / mean(trendTD(1:Lmin));
        perTraceDiff(k) = 100 * max(abs(shapesPoly(k, :) - shapesTD(k, :)));
    end
    shapePoly = median(shapesPoly, 1); shapeTD = median(shapesTD, 1);
    diffShared = 100 * max(abs(shapePoly - shapeTD));
    fprintf('%-12s %3d traces: shared shapes differ by at most %.2f %% (per trace: median %.1f %%, 90th pct %.1f %%)\n', ...
        names{d}, n, diffShared, median(perTraceDiff), prctile(perTraceDiff, 90));
    summary = [summary; table(string(names{d}), n, diffShared, median(perTraceDiff), prctile(perTraceDiff, 90), ...
        100 * (shapePoly(end) - shapePoly(1)), 100 * (shapeTD(end) - shapeTD(1)), ...
        'VariableNames', {'dataset', 'traces', 'maxDiffShared_pct', 'medianDiffPerTrace_pct', 'p90DiffPerTrace_pct', 'dropPoly_pct', 'dropTD_pct'})]; %#ok<AGROW>
    allShapes = [allShapes; table(repmat(string(names{d}), Lmin, 1), (1:Lmin)', shapePoly', shapeTD', ...
        'VariableNames', {'dataset', 'frame', 'poly3', 'trenddecomp'})]; %#ok<AGROW>
end
writetable(summary, fullfile(outDir, 'summary.csv'));
writetable(allShapes, fullfile(outDir, 'shapes.csv'));


function [S, B, G] = loadData(name)
% spot mean, background and spot sum of one paper dataset (as ExperimentalAnalysis.m)
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/';
    file = strrep(name, 'vitro', 'same');
    S = loadTraces([folder file '_data.mat']);
    B = loadTraces([folder file '_background.mat']);
    G = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    S = loadTraces([folder name '_data.mat']);
    B = loadTraces([folder name '_background.mat']);
    G = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
end


function traces = loadTraces(file)
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end
