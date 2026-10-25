% SHED_FRACTION  Each burst as a fraction of the granule intensity just before
% it ("fraction shed" for negative bursts, "fraction gained" for positive),
% instead of lambda. Unitless, so conditions can be compared without knowing
% the single-unit brightness k0. Reviewer points R2 4b/7 (interpretability),
% Qbeta vs PP7 and FUS claims.
% Measures (same burst detection as the runners, analysis.m, window 13):
%   mean     current trace (mean of the multithresh spot pixels above the
%            cell level); PP7/Qb R2 re-extraction, in vitro stored data
%   spotsum  sum of the spot pixels above the cell mean (after_R2_sums)
%   cropsum  sum of the 15x15 crop above its border median (after_R2_sums)
% For the sums the analysis gets signal = sum + background, so its
% normalisation (S/B - 1)*B(1) returns the bleach-corrected sum itself.
% Fraction = |amplitude| / fitted level at the burst start (bursts with a
% level <= 0 are dropped and counted). SP 10..15; traces >= 180 frames.
% Writes analysis_outputs/shed_fraction/fractions.csv and dropped.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'shed_fraction');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset'; SU = 'in_vivo_reextracted/after_R2_sums';
vivo = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'};
vitro = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'};
jobs = {};
for d = vivo
    jobs(end+1,:) = {d{1}, 'mean', R2, d{1}, false}; %#ok<SAGROW>
    jobs(end+1,:) = {d{1}, 'spotsum', SU, [d{1} '_spotsum'], true}; %#ok<SAGROW>
    jobs(end+1,:) = {d{1}, 'cropsum', SU, [d{1} '_cropsum'], true}; %#ok<SAGROW>
end
for d = vitro, jobs(end+1,:) = {d{1}, 'mean', 'in_vitro_experiments', d{1}, false}; end %#ok<SAGROW>
tVec = (0:360) * 10/60; rows = cell(0, 9); drop = {};
for j = 1:size(jobs,1)
    [S, B] = gm_loadDataset(repoDir, jobs{j,3}, jobs{j,4});
    keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep);
    if jobs{j,5}, S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false); end
    for SP = 10:15
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(S, B, tVec, 13, SP, 0, 0);
        warning(ws);
        ok = seg(:,6) > 0;
        drop(end+1,:) = {jobs{j,1}, jobs{j,2}, SP, sum(~ok), size(seg,1)}; %#ok<SAGROW>
        seg = seg(ok,:);
        sgn = repmat({'Pos'}, size(seg,1), 1); sgn(seg(:,2) < 0) = {'Neg'};
        rows = [rows; [repmat(jobs(j,1:2), size(seg,1), 1), num2cell(repmat(SP, size(seg,1), 1)), sgn, ...
            num2cell(seg(:,1)), num2cell(abs(seg(:,5))), num2cell(seg(:,6)), num2cell(abs(seg(:,5)) ./ seg(:,6)), ...
            num2cell(seg(:,4) - seg(:,3))]]; %#ok<AGROW>
    end
    fprintf('%s %s done\n', jobs{j,1}, jobs{j,2});
end
writetable(cell2table(rows, 'VariableNames', {'dataset','measure','SP','sign','trace','absAmp','levelBefore','fraction','durationFrames'}), ...
    fullfile(outDir, 'fractions.csv'));
writetable(cell2table(drop, 'VariableNames', {'dataset','measure','SP','droppedNonPositiveLevel','bursts'}), ...
    fullfile(outDir, 'dropped.csv'));
disp('SHED_DONE');
