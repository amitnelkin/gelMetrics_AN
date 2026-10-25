function detection_grid(suFolder, outName, SPs, wins, doVitro, vitroFolder, vitroJobs, doVivo)
% DETECTION_GRID  Burst detection across stringency SP and smoothing window,
% on the real traces and on burst-free surrogates (one shared bleaching shape
% per dataset scaled to each trace + white noise at the trace's robust noise
% level, as score_export2 'bleach'), to see whether any setting separates real
% bursts from noise better (RESULTS.md section 20). Reviewer points R2 6, R2 2b,
% R2 4a. Main measure per dataset (spot sum in vivo, mean in vitro), traces >= 180
% frames, same detection code as the runners (analysis.m).
% SP in {3 4 5 6 8 10 12}, window in {5 7 9 13 17 21} frames.
% Writes analysis_outputs/detection_grid/bursts.csv, traces.csv.
% Optional (section 23): suFolder (spot-sum data, default after_R2_sums),
% outName (output folder under analysis_outputs), SPs, wins, doVitro (default true).
if nargin < 1, suFolder = 'in_vivo_reextracted/after_R2_sums'; end
if nargin < 2, outName = 'detection_grid'; end
if nargin < 3, SPs = [3 4 5 6 8 10 12]; end
if nargin < 4, wins = [5 7 9 13 17 21]; end
if nargin < 5, doVitro = true; end
% Optional (section 38): vitroFolder (default 'in_vitro_experiments'), vitroJobs = {label,
% file stem; ...} (default the three stored datasets), doVivo (default true). For a
% re-extracted folder its extraction_info*.csv gives the movie of each trace.
if nargin < 6, vitroFolder = 'in_vitro_experiments'; end
if nargin < 7, vitroJobs = {'protein_RNA_1_1','protein_RNA_1_1'; 'protein_RNA_10_1','protein_RNA_10_1'; 'protein_RNA_100_1','protein_RNA_100_1'}; end
if nargin < 8, doVivo = true; end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', outName);
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset'; SU = suFolder;
info = [readtable(fullfile(repoDir, R2, 'extraction_info.csv')); readtable(fullfile(repoDir, R2, 'extraction_info_5Qb_10Qb.csv'))];
vi = dir(fullfile(repoDir, vitroFolder, 'extraction_info*.csv'));
for q = 1:numel(vi), info = [info; readtable(fullfile(vi(q).folder, vi(q).name))]; end %#ok<AGROW>
jobs = {};
for d = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    if doVivo, jobs(end+1,:) = {d{1}, SU, [d{1} '_spotsum'], true, d{1}}; end %#ok<AGROW>
end
for q = 1:size(vitroJobs, 1)
    if doVitro, jobs(end+1,:) = {vitroJobs{q,1}, vitroFolder, vitroJobs{q,2}, false, vitroJobs{q,2}}; end %#ok<AGROW>
end
tVec = (0:360) * 10/60; bu = cell(0, 9); tr = cell(0, 4);
for j = 1:size(jobs,1)
    [S, B] = gm_loadDataset(repoDir, jobs{j,2}, jobs{j,3});
    I = info(strcmp(info.dataset, jobs{j,5}), :);
    if isempty(I), mov = zeros(numel(S), 1); else, mov = repelem((1:height(I))', I.nKept); end
    keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep); mov = mov(keep);
    if jobs{j,4}, S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false); end
    % shared-bleaching surrogate
    rng(17 + j);
    Lm = min(cellfun(@numel, S)); shp = zeros(numel(S), Lm);
    for i = 1:numel(S)
        s = S{i}(:); x = (1:numel(s))'; trd = polyval(polyfit(x, s, 3), x); shp(i,:) = trd(1:Lm)' / mean(trd(1:Lm));
    end
    shape = median(shp, 1); Ssur = S;
    for i = 1:numel(S)
        s = S{i}(:); n = numel(s); sh = [shape, repmat(shape(end), 1, max(n - Lm, 0))]'; sh = sh(1:n);
        sgm = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
        Ssur{i} = mean(s) / mean(sh) * sh + sgm * randn(n, 1);
    end
    for i = 1:numel(S), tr(end+1,:) = {jobs{j,1}, i, numel(S{i}), mov(i)}; end %#ok<AGROW>
    for w = wins
        for kind = {'real', 'surr'}
            X = S; if strcmp(kind{1}, 'surr'), X = Ssur; end
            for SP = SPs
                ws = warning('off', 'all');
                [~, ~, ~, ~, ~, ~, seg] = analysis(X, B, tVec, w, SP, 0, 0);
                warning(ws);
                seg = seg(seg(:,2) < 0, :);
                for q = 1:size(seg,1)
                    bu(end+1,:) = {jobs{j,1}, w, SP, kind{1}, seg(q,1), seg(q,5), seg(q,6), seg(q,3), seg(q,4) - seg(q,3)}; %#ok<AGROW>
                end
            end
        end
    end
    fprintf('%s done\n', jobs{j,1});
end
writetable(cell2table(bu, 'VariableNames', {'dataset','window','SP','kind','trace','amp','levelBefore','startFrame','durationFrames'}), ...
    fullfile(outDir, 'bursts.csv'));
writetable(cell2table(tr, 'VariableNames', {'dataset','trace','lengthFrames','movie'}), fullfile(outDir, 'traces.csv'));
disp('GRID_DONE');
end
