function fp1000(nCopies, SPs, datasets)
% FP1000  Calibrate burst detection so that a burst-free granule shows a false
% burst with a chance of about 1 in 1000 (per granule, whole trace, either sign),
% then count bursts in the real granules at that setting (RESULTS.md section 37).
% Burst-free copies are built exactly as in detection_grid.m (one shared bleaching
% shape per dataset scaled to each trace + white noise at the trace's robust noise
% level), but nCopies independent copies per granule so a 1/1000 rate can be
% estimated. Same detector as the runners (analysis.m, window 13).
% Writes analysis_outputs/fp1000/segs.csv (one row per burst) and traces.csv.
if nargin < 1, nCopies = 100; end
if nargin < 2, SPs = [12 16 20 24 28 32 36 40]; end
if nargin < 3, datasets = {'4PP7', '4PP7-FUS'}; end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'fp1000');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset'; SU = 'in_vivo_reextracted/after_R2_sums';
info = readtable(fullfile(repoDir, R2, 'extraction_info.csv'));
tVec = (0:360) * 10/60; win = 13;
segs = zeros(0, 10); tr = zeros(0, 5);
for j = 1:numel(datasets)
    [S, B] = gm_loadDataset(repoDir, SU, [datasets{j} '_spotsum']);
    I = info(strcmp(info.dataset, datasets{j}), :);
    mov = repelem((1:height(I))', I.nKept);
    keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep); mov = mov(keep);
    S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false);
    n = numel(S);
    % shared bleaching shape and per-trace noise level (as detection_grid.m)
    Lm = min(cellfun(@numel, S)); shp = zeros(n, Lm); sgm = zeros(n, 1);
    for i = 1:n
        s = S{i}(:); x = (1:numel(s))'; trd = polyval(polyfit(x, s, 3), x);
        shp(i,:) = trd(1:Lm)' / mean(trd(1:Lm));
        sgm(i) = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
    end
    shape = median(shp, 1);
    for i = 1:n, tr(end+1,:) = [j, i, numel(S{i}), mov(i), mean(S{i})]; end %#ok<AGROW>
    % copies: nCopies per granule, seeded
    C = cell(n * nCopies, 1); BC = C; cid = zeros(n * nCopies, 2); q = 0;
    for k = 1:nCopies
        rng(1000 * j + k);
        for i = 1:n
            s = S{i}(:); m = numel(s); sh = [shape, repmat(shape(end), 1, max(m - Lm, 0))]'; sh = sh(1:m);
            q = q + 1; C{q} = mean(s) / mean(sh) * sh + sgm(i) * randn(m, 1); BC{q} = B{i}; cid(q,:) = [k, i];
        end
    end
    for SP = SPs
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, sr] = analysis(S, B, tVec, win, SP, 0, 0);
        [~, ~, ~, ~, ~, ~, sc] = analysis(C, BC, tVec, win, SP, 0, 0);
        warning(ws);
        % [dataset, copy (0 = real), trace, SP, sign, first, last, amp, levelBefore, levelAfter]
        segs = [segs; repmat([j 0], size(sr,1), 1), sr(:,1), repmat(SP, size(sr,1), 1), sr(:,2:7); ...
                cid(sc(:,1), 1) * 0 + j, cid(sc(:,1), :), repmat(SP, size(sc,1), 1), sc(:,2:7)]; %#ok<AGROW>
        fprintf('%s SP %d: real %d, copies %d\n', datasets{j}, SP, size(sr,1), size(sc,1));
    end
end
writetable(array2table(segs, 'VariableNames', {'dataset','copy','trace','SP','sign','first','last','amp','levelBefore','levelAfter'}), fullfile(outDir, 'segs.csv'));
writetable(array2table(tr, 'VariableNames', {'dataset','trace','lengthFrames','movie','meanLevel'}), fullfile(outDir, 'traces.csv'));
writecell(datasets(:), fullfile(outDir, 'datasets.csv'));
disp('FP1000_DONE');
end
