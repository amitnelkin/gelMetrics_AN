function example_traces()
% EXAMPLE_TRACES  One real trace with detected bursts and its burst-free
% surrogate, through exactly the detector of analysis.m (manuscript defaults:
% window 13, SP 12), for 4PP7-FUS (spot sum) and in vitro 10:1 (spot mean).
% Surrogate as in checks/detection_grid.m: the dataset's shared bleaching
% shape scaled to the trace + white noise at the trace's own noise level
% (seeded). Trace choice, fixed in advance: among traces >= 300 frames, the one
% whose number of detected negative bursts is closest to the 75th percentile
% (first such trace). RESULTS.md section 21.
% Writes analysis_outputs/noise_background/example_traces.csv, example_bursts.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_background');
if ~exist(outDir, 'dir'), mkdir(outDir); end
jobs = {'4PP7-FUS', 'in_vivo_reextracted/after_R2_sums', '4PP7-FUS_spotsum', true; ...
        'protein_RNA_10_1', 'in_vitro_experiments-old', 'protein_RNA_10_1', false};
w = 13; SP = 12; tVec = (0:360) * 10/60;
rows = cell(0, 6); bur = cell(0, 8);
for j = 1:size(jobs, 1)
    [S, B] = gm_loadDataset(repoDir, jobs{j,2}, jobs{j,3});
    keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep);
    if jobs{j,4}, S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false); end
    rng(101 + j);
    Lm = min(cellfun(@numel, S)); shp = zeros(numel(S), Lm);
    for i = 1:numel(S)
        s = S{i}(:); x = (1:numel(s))'; trd = polyval(polyfit(x, s, 3), x); shp(i,:) = trd(1:Lm)' / mean(trd(1:Lm));
    end
    shape = median(shp, 1);
    ws = warning('off', 'all');
    [~, ~, ~, ~, ~, ~, seg] = analysis(S, B, tVec, w, SP, 0, 0);
    warning(ws);
    nNeg = accumarray(seg(seg(:,2) < 0, 1), 1, [numel(S) 1]);
    long = find(cellfun(@numel, S) >= 300);
    target = prctile(nNeg(long), 75);
    [~, p] = min(abs(nNeg(long) - target)); i = long(p);
    s = S{i}(:); n = numel(s); sh = [shape, repmat(shape(end), 1, max(n - Lm, 0))]'; sh = sh(1:n);
    sgm = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
    sur = mean(s) / mean(sh) * sh + sgm * randn(n, 1);
    for kind = {'real', 'surrogate'}
        x = s; if strcmp(kind{1}, 'surrogate'), x = sur; end
        [proc, sg] = detect(x, B{i}(:), tVec, w, SP);
        for t = 1:n
            rows(end+1,:) = {jobs{j,1}, kind{1}, i, t, x(t), proc(t)}; %#ok<AGROW>
        end
        for q = 1:size(sg, 1)
            bur(end+1,:) = {jobs{j,1}, kind{1}, i, sg(q,2), sg(q,3), sg(q,4), sg(q,6), sg(q,7)}; %#ok<AGROW>
        end
    end
    fprintf('%s: trace %d, %d frames, noise sd %.1f\n', jobs{j,1}, i, n, sgm);
end
writetable(cell2table(rows, 'VariableNames', {'dataset','kind','trace','frame','raw','processed'}), fullfile(outDir, 'example_traces.csv'));
writetable(cell2table(bur, 'VariableNames', {'dataset','kind','trace','sign','firstFrame','lastFrame','levelStart','levelEnd'}), fullfile(outDir, 'example_bursts.csv'));
disp('EXAMPLES_DONE');
end

function [data, seg] = detect(s, b, tVec, w, SP)
% the processed signal exactly as analysis.m builds it, plus its bursts
ws = warning('off', 'all');
tv = tVec(1:numel(s));
f = fit(tv', b, 'poly3'); bf = feval(f, tv);
x = filloutliers(s', 'nearest', 'movmedian', 3 * w);
x = smoothdata(x, 'movmean', w);
data = ((x ./ bf' - 1)) * bf(1);
[~, ~, ~, ~, ~, ~, seg] = analysis({s'}, {b'}, tVec, w, SP, 0, 0);
warning(ws);
end
