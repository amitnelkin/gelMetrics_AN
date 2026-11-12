% Per-granule values behind the dot previews of Figs. 9(d), 10(b-e) and 13 (section 80).
% Run from the repo root with OUTDIR set (e.g. setenv('OUTDIR', 'analysis_outputs/granule_points')).
% Writes vivo_granules.csv (brightness, mean burst, burst size beyond noise, turnover, net change per granule;
% same settings as ExperimentalAnalysis.m: 20 copies, seed 1, window 13, SP 12) and drops_granules.csv
% (Fig. 13: drops >= 30% in 15 min, granule and copy counts, hours; mode 1 relative, 2 absolute).
% Needs burstFreeCopies.m and countLargeChanges.m from AppendixFigures.m on the path (checks/granule_points_helpers/).
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

OUT = getenv('OUTDIR');
% --- per granule, in vivo: brightness, mean burst size, burst size beyond noise, turnover beyond noise, net change
fid = fopen(fullfile(OUT, 'vivo_granules.csv'), 'w');
fprintf(fid, 'condition,granule,level,bursts,meanBurst,sizeBeyond,turnover,net\n');
for c = 1:nVivo
    M = burstMeasures(D(c).S, D(c).B, D(c).G, [], 20, 1, 13, 12, 361);
    for k = 1:M.nGranules
        a = M.ampAU{k}; if isempty(a), mb = NaN; else, mb = mean(a); end
        fr = M.fracReal{k}; fc = M.fracCopy{k};
        if isempty(fr) || isempty(fc), sb = NaN; else, sb = 100 * (median(fr) - median(fc)); end
        V = summarizeGranules(M, k);
        fprintf(fid, '%s,%d,%.2f,%d,%.4f,%.4f,%.4f,%.4f\n', names{c}, k, M.level(k), numel(a), mb, sb, V.turnover, M.granuleNet(k));
    end
    fprintf('%s: %d granules; pooled turnover %.2f net %.2f\n', names{c}, M.nGranules, M.turnover, M.net);
end
fclose(fid);
% --- per granule, Fig. 13: excess drops >= 30% within 15 min per granule-hour (relative and absolute)
threshold = 0.3; maxFrames = 90; nCopies = 5;
levels = [];
for c = 1:4
    for k = 1:numel(D(c).X)
        x = (1:numel(D(c).B{k}))'; backgroundFit = polyval(polyfit(x, D(c).B{k}, 3), x);
        levels(end+1, 1) = median(D(c).X{k} ./ backgroundFit * backgroundFit(1)); %#ok<SAGROW>
    end
end
unitAbsolute = median(levels);
fid = fopen(fullfile(OUT, 'drops_granules.csv'), 'w');
fprintf(fid, 'condition,granule,mode,real,copies,hours\n');
for m = 1:2
    for c = 1:nVivo
        rng(35 + c);
        [C, ~, owner] = burstFreeCopies(D(c).X, D(c).B, nCopies);
        for k = 1:numel(D(c).X)
            x = (1:numel(D(c).B{k}))'; backgroundFit = polyval(polyfit(x, D(c).B{k}, 3), x);
            if m == 1, scale = 1; unit = []; else, scale = backgroundFit(1); unit = unitAbsolute; end
            [~, down] = countLargeChanges(movmean(D(c).X{k}, 13) ./ backgroundFit * scale, threshold, maxFrames, unit);
            cp = 0; q = find(owner == k);
            for j = 1:numel(q)
                [~, dd] = countLargeChanges(movmean(C{q(j)}, 13) ./ backgroundFit * scale, threshold, maxFrames, unit); cp = cp + dd / nCopies;
            end
            fprintf(fid, '%s,%d,%d,%d,%.4f,%.4f\n', names{c}, k, m, down, cp, numel(D(c).X{k}) / 360);
        end
    end
end
fclose(fid);
