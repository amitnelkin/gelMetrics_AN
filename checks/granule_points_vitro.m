% Fig. 6(d-f) values per granule (same analysis as ExperimentalAnalysis.m: 20 copies, seed 1, window 13, SP 12).
conditions = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'};
out = fopen(fullfile(getenv('OUTDIR'), 'fig6_granules.csv'), 'w');
fprintf(out, 'condition,granule,movie,bursts,meanBurst,turnover,net\n');
for c = 1:3
    folder = ['./in_vitro_experiments/' conditions{c} '/'];
    signals = struct2cell(load([folder 'data.mat'])); signals = signals{1};
    backgrounds = struct2cell(load([folder 'background.mat'])); backgrounds = backgrounds{1};
    spotSums = struct2cell(load([folder 'spotsum.mat'])); spotSums = spotSums{1};
    movies = load([folder 'movies.mat']);
    signals = signals(2:end); backgrounds = backgrounds(2:end); spotSums = spotSums(2:end);
    signals = signals(:); backgrounds = backgrounds(:); spotSums = spotSums(:);
    M = burstMeasures(signals, backgrounds, spotSums, movies.movieId, 20, 1, 13, 12, 361);
    for k = 1:M.nGranules
        a = M.ampAU{k};
        if isempty(a), mb = NaN; else, mb = mean(a); end
        V = summarizeGranules(M, k);
        fprintf(out, '%s,%d,%d,%d,%.4f,%.4f,%.4f\n', conditions{c}, k, M.replicate(k), numel(a), mb, V.turnover, M.granuleNet(k));
    end
end
fclose(out);
