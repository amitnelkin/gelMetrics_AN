function reextract_fixed(workRoot, posWin)
% REEXTRACT_FIXED  Re-extraction with one fixed multithresh mask per granule
% (collect_intensity opts.fixedMask; RESULTS.md section 33). Same tracking, crop,
% offset 200 and multithresh classes as R2, but the mask is computed once from
% the granule's average crop (track position smoothed over posWin frames,
% default 15, fixed in advance) and used unchanged in every frame. In vivo only.
% Writes in_vivo_reextracted/after_R5_fixedmask/<ds>_data.mat / _background.mat
% (spot mean, cell mean, as R2), <ds>_spotsum_data.mat / _spotsum_background.mat,
% <ds>_spotarea.mat (as after_R2_sums) and <ds>_maskfrac.mat (share of the granule
% light inside the fixed mask, per frame). The repo .mat files are untouched.
if nargin < 2, posWin = 15; end
checksDir = fileparts(mfilename('fullpath')); repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir);
outDir = fullfile(repoDir, 'in_vivo_reextracted', 'after_R5_fixedmask');
if ~exist(outDir, 'dir'), mkdir(outDir); end
for ds = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    [I, B, ~, S] = collect_intensity_tests(fullfile(workRoot, ds{1}, 'work'), struct('offset', 200, 'fixedMask', true, 'posWin', posWin));
    ref = load(fullfile(repoDir, 'in_vivo_reextracted', 'after_R2_offset', [ds{1} '_data.mat']));
    sameLen = isequal(cellfun(@numel, I), cellfun(@numel, ref.intensity));
    fprintf('%s: %d traces, same traces and lengths as R2: %d\n', ds{1}, numel(I) - 1, sameLen);
    background = B; intensity = I; %#ok<NASGU>
    save(fullfile(outDir, [ds{1} '_data.mat']), 'intensity'); save(fullfile(outDir, [ds{1} '_background.mat']), 'background');
    intensity = S.spotSum; %#ok<NASGU>
    save(fullfile(outDir, [ds{1} '_spotsum_data.mat']), 'intensity'); save(fullfile(outDir, [ds{1} '_spotsum_background.mat']), 'background');
    spotArea = S.spotArea; save(fullfile(outDir, [ds{1} '_spotarea.mat']), 'spotArea');
    maskFrac = S.maskFrac; save(fullfile(outDir, [ds{1} '_maskfrac.mat']), 'maskFrac');
end
disp('FIXED_DONE');
end
