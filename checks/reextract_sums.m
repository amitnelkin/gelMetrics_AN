function reextract_sums(workRoot)
% REEXTRACT_SUMS  Integrated granule intensity (collect_intensity.m 4th output)
% for the PP7, PP7-FUS and Qbeta datasets, from the existing work folders
% (red-channel tracking, R2 settings: offset 200, multithresh spot). Used to
% express each burst as a fraction of the granule (checks/shed_fraction.m).
% Sanity check: the mean-based traces must equal in_vivo_reextracted/after_R2_offset.
% Writes in_vivo_reextracted/after_R2_sums/<ds>_spotsum_data.mat, _cropsum_data.mat
% (variable intensity) and matching _background.mat (the cell-class mean, as R2).
checksDir = fileparts(mfilename('fullpath')); repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir);
outDir = fullfile(repoDir, 'in_vivo_reextracted', 'after_R2_sums');
if ~exist(outDir, 'dir'), mkdir(outDir); end
for ds = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    [I, B, ~, S] = collect_intensity(fullfile(workRoot, ds{1}, 'work'), struct('offset', 200));
    ref = load(fullfile(repoDir, 'in_vivo_reextracted', 'after_R2_offset', [ds{1} '_data.mat']));
    same = isequal(numel(I), numel(ref.intensity)) && all(cellfun(@(a, b) isequal(a, b), I(2:end), ref.intensity(2:end)));
    fprintf('%s: %d traces, mean traces identical to R2: %d\n', ds{1}, numel(I) - 1, same);
    background = B; %#ok<NASGU>
    intensity = S.spotSum; save(fullfile(outDir, [ds{1} '_spotsum_data.mat']), 'intensity'); %#ok<NASGU>
    save(fullfile(outDir, [ds{1} '_spotsum_background.mat']), 'background');
    intensity = S.cropSum; save(fullfile(outDir, [ds{1} '_cropsum_data.mat']), 'intensity');
    save(fullfile(outDir, [ds{1} '_cropsum_background.mat']), 'background');
    spotArea = S.spotArea; save(fullfile(outDir, [ds{1} '_spotarea.mat']), 'spotArea');
end
disp('SUMS_DONE');
end
