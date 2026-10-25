function vitro_orig_verify()
% VITRO_ORIG_VERIFY  Are the titration movies the source of the stored in vitro
% data? (RESULTS.md section 40.) The ORIGINAL pipeline (f7eded0 read_tables.m +
% Collect_intensity_data.m + generate_movArray.m, unchanged; a uigetdir shim
% gives the folder, as checks for 4PP7 and Qb in sections 4 and 12) is run on the
% re-tracked movies (..\zigzag_scratch\orig_vitro\<folder>: the movie TIFF and the
% Mosaic trajectory CSV from BatchParticleTracker.ijm). Then:
%  1. traces: each stored trace (in_vitro_experiments) is searched among the
%     new traces of every folder: bit-exact (same length, max |diff| < 1e-9);
%  2. bursts: the ORIGINAL analysis.m (f7eded0; window 13, SP 14, time vector
%     linspace(0, 60, 360), as ExperimentalAnalysis.m) is run on the stored
%     traces and on their bit-exact twins in the same order, and the burst
%     amplitudes are compared; also run on all new traces of the matching folder.
% Writes analysis_outputs/vitro_orig_verify/{traces.csv, bursts.csv, per_movie.csv}.
Z = 'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\zigzag_scratch';
repo = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(repo, 'analysis_outputs', 'vitro_orig_verify'); if ~exist(out, 'dir'), mkdir(out); end
NEW = {'vitro_1000ngul_0p1mgml', 'vitro_1000ng_10ug', 'vitro_1000ngul_1mgml', 'vitro_100ngul_1mgml'};
STORED = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'};
for k = 1:numel(NEW)
    w = fullfile(Z, 'orig_vitro', NEW{k});
    if ~exist(fullfile(w, 'data.mat'), 'file')
        run_orig(Z, w);
    end
end
addpath(fullfile(Z, 'orig_analysis'));   % original analysis.m first on the path
N = struct(); Mv = struct();
for k = 1:numel(NEW)
    w = fullfile(Z, 'orig_vitro', NEW{k});
    d = load(fullfile(w, 'data.mat')); b = load(fullfile(w, 'background.mat'));
    N.(NEW{k}) = {d.intensity(2:end), b.background(2:end)};
    % movie of every kept trace (same rule as Collect_intensity_data: > = 24 points after dropping the first)
    xs = dir(fullfile(w, '*_xy.mat')); mv = [];
    for f = 1:numel(xs)
        xy = load(fullfile(w, xs(f).name));
        mv = [mv, repmat(f, 1, sum(cellfun(@numel, xy.x) - 1 >= 24))]; %#ok<AGROW>
    end
    assert(numel(mv) == numel(N.(NEW{k}){1}));
    Mv.(NEW{k}) = mv;
end
trows = cell(0, 6); brows = cell(0, 9); prow = cell(0, 5);
for s = 1:numel(STORED)
    S = struct2cell(load(fullfile(repo, 'in_vitro_experiments', [STORED{s} '_data.mat']))); S = S{1}(2:end);
    B = struct2cell(load(fullfile(repo, 'in_vitro_experiments', [STORED{s} '_background.mat']))); B = B{1}(2:end);
    for k = 1:numel(NEW)
        T = N.(NEW{k}){1}; hit = zeros(numel(S), 1);
        for i = 1:numel(S)
            for j = 1:numel(T)
                if numel(T{j}) == numel(S{i}) && max(abs(double(T{j}(:)) - double(S{i}(:)))) < 1e-9, hit(i) = j; break, end
            end
        end
        trows(end+1, :) = {STORED{s}, NEW{k}, numel(S), sum(hit > 0), numel(T), mat2str(unique(Mv.(NEW{k})(hit(hit > 0))))}; %#ok<AGROW>
        fprintf('%s vs %s: %d of %d stored traces found bit-exactly (%d new traces)\n', STORED{s}, NEW{k}, sum(hit > 0), numel(S), numel(T));
        if sum(hit > 0) == 0, continue, end
        for f = unique(Mv.(NEW{k}))
            prow(end+1, :) = {STORED{s}, NEW{k}, f, sum(Mv.(NEW{k}) == f), sum(Mv.(NEW{k})(hit(hit > 0)) == f)}; %#ok<AGROW>
        end
        % bursts: stored traces vs their twins (only stored traces that have a twin)
        ok = hit > 0; Tn = N.(NEW{k});
        [a1p, a1n] = orig_bursts(S(ok), B(ok)); [a2p, a2n] = orig_bursts(Tn{1}(hit(ok)), Tn{2}(hit(ok)));
        [a0p, a0n] = orig_bursts(S, B); [a3p, a3n] = orig_bursts(Tn{1}, Tn{2});
        same = numel(a1n) == numel(a2n) && numel(a1p) == numel(a2p) && max([abs(a1n(:) - a2n(:)); abs(a1p(:) - a2p(:)); 0]) < 1e-9;
        brows(end+1, :) = {STORED{s}, NEW{k}, sum(ok), numel(a1n), numel(a2n), same, numel(a0n), numel(a3n), ...
            sprintf('stored all: median |neg| %.1f; new all: %.1f', median(abs(a0n)), median(abs(a3n)))}; %#ok<AGROW>
        fprintf('   bursts on matched traces: stored %d neg / %d pos, new %d neg / %d pos, identical %d\n', numel(a1n), numel(a1p), numel(a2n), numel(a2p), same);
        fprintf('   all traces: stored %d neg bursts (median |amp| %.1f), new %d (%.1f)\n', numel(a0n), median(abs(a0n)), numel(a3n), median(abs(a3n)));
    end
end
writetable(cell2table(trows, 'VariableNames', {'stored','folder','storedTraces','bitExact','newTraces','moviesWithTwins'}), fullfile(out, 'traces.csv'));
writetable(cell2table(brows, 'VariableNames', {'stored','folder','matchedTraces','negBurstsStored','negBurstsNew','identical','negBurstsStoredAll','negBurstsNewAll','note'}), fullfile(out, 'bursts.csv'));
writetable(cell2table(prow, 'VariableNames', {'stored','folder','movie','newTraces','twins'}), fullfile(out, 'per_movie.csv'));
disp('VITRO_ORIG_DONE');
end

function run_orig(Z, w)
% the original scripts call clear/clearvars, so they run in a separate function workspace
addpath(fullfile(Z, 'shim_env')); addpath(fullfile(Z, 'orig_code'));
setenv('QB_DIR', w); setenv('VITRO_RET', pwd); cd(w);
read_tables; Collect_intensity_data;   % both clear the workspace
cd(getenv('VITRO_RET'));
end

function [ap, an] = orig_bursts(S, B)
ws = warning('off', 'all');
[~, ~, ap, an] = analysis(S(:), B(:), linspace(0, 60, 360), 13, 14, 0, 0);
warning(ws);
ap = ap(2:end); an = an(2:end);
end
