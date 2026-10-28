function reextract_steady(workRoot, thrWin)
% REEXTRACT_STEADY  Re-extraction with steady multithresh thresholds
% (collect_intensity opts.thrWin; RESULTS.md section 23). Same tracking, crop,
% offset 200 and multithresh classes as R2, but each track's four thresholds are
% the running median over thrWin frames (default 31 = ~5 min, fixed in advance)
% instead of being recomputed on every frame. In vivo only (no in vitro movies
% on the lab PC).
% Writes in_vivo_reextracted/after_R4_steadythr/<ds>_data.mat / _background.mat
% (spot mean, cell mean, as R2), <ds>_spotsum_data.mat / _spotsum_background.mat
% and <ds>_spotarea.mat (as after_R2_sums). The repo .mat files are untouched.
if nargin < 2, thrWin = 31; end
checksDir = fileparts(mfilename('fullpath')); repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir);
outDir = fullfile(repoDir, 'in_vivo_reextracted', 'after_R4_steadythr');
if ~exist(outDir, 'dir'), mkdir(outDir); end
for ds = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    [I, B, ~, S] = collect_intensity_tests(fullfile(workRoot, ds{1}, 'work'), struct('offset', 200, 'thrWin', thrWin));
    ref = load(fullfile(repoDir, 'in_vivo_reextracted', 'after_R2_offset', [ds{1} '_data.mat']));
    sameLen = isequal(cellfun(@numel, I), cellfun(@numel, ref.intensity));
    fprintf('%s: %d traces, same traces and lengths as R2: %d\n', ds{1}, numel(I) - 1, sameLen);
    background = B; intensity = I; %#ok<NASGU>
    save(fullfile(outDir, [ds{1} '_data.mat']), 'intensity'); save(fullfile(outDir, [ds{1} '_background.mat']), 'background');
    intensity = S.spotSum; %#ok<NASGU>
    save(fullfile(outDir, [ds{1} '_spotsum_data.mat']), 'intensity'); save(fullfile(outDir, [ds{1} '_spotsum_background.mat']), 'background');
    spotArea = S.spotArea; save(fullfile(outDir, [ds{1} '_spotarea.mat']), 'spotArea');
end
disp('STEADY_DONE');
end
