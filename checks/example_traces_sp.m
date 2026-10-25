function example_traces_sp(SPs)
% EXAMPLE_TRACES_SP  The example of checks/example_traces.m (one real trace and
% its burst-free copy through the detector) repeated at higher stringency SP, for
% the same 4PP7-FUS trace and copy as before (same selection rule at SP 12, same
% seed), plus a 5Qb and a 4PP7 example chosen by the same rule (RESULTS.md section 34).
% SPs default [12 15 18 21 25]; window 13.
% Writes analysis_outputs/noise_background/example_traces_sp.csv, example_bursts_sp.csv.
if nargin < 1, SPs = [12 15 18 21 25]; end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_background');
if ~exist(outDir, 'dir'), mkdir(outDir); end
jobs = {'4PP7-FUS', 'in_vivo_reextracted/after_R2_sums', '4PP7-FUS_spotsum', true; ...
        '5Qb', 'in_vivo_reextracted/after_R2_sums', '5Qb_spotsum', true; ...
        '4PP7', 'in_vivo_reextracted/after_R2_sums', '4PP7_spotsum', true};
w = 13; SP = 12; tVec = (0:360) * 10/60;
rows = cell(0, 6); bur = cell(0, 9);
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
        for spq = SPs
            [proc, sg] = detect(x, B{i}(:), tVec, w, spq);
            if spq == SPs(1)
                for t = 1:n
                    rows(end+1,:) = {jobs{j,1}, kind{1}, i, t, x(t), proc(t)}; %#ok<AGROW>
                end
            end
            for q = 1:size(sg, 1)
                bur(end+1,:) = {jobs{j,1}, kind{1}, i, spq, sg(q,2), sg(q,3), sg(q,4), sg(q,6), sg(q,7)}; %#ok<AGROW>
            end
        end
    end
    fprintf('%s: trace %d, %d frames, noise sd %.1f\n', jobs{j,1}, i, n, sgm);
end
writetable(cell2table(rows, 'VariableNames', {'dataset','kind','trace','frame','raw','processed'}), fullfile(outDir, 'example_traces_sp.csv'));
writetable(cell2table(bur, 'VariableNames', {'dataset','kind','trace','SP','sign','firstFrame','lastFrame','levelStart','levelEnd'}), fullfile(outDir, 'example_bursts_sp.csv'));
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
