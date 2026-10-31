% INSPECT_TIFF_CHANNELS  Read-only. ExampleTIF.tif is an ImageJ hyperstack
% with 2 channels x 361 frames, interleaved (page 1 = ch1 t1, page 2 = ch2
% t1, ...). generate_movArray.m reads every page as a frame, so this checks
% (a) the pixel floor per channel (camera dark offset), and (b) the lengths
% of the stored traces in the repo's .mat files, which tell whether the
% extraction ran on 360/361-frame or 722-page stacks.
clear; close all;
repoDir = '..';
outDir = fullfile(repoDir, 'analysis_outputs', 'baseline', 'tiff_check');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fn = fullfile(repoDir, 'ExampleTIF.tif');
nPages = numel(imfinfo(fn));
tIdx = unique(round(linspace(1, nPages/2, 25)));
for ch = 1:2
    mins = zeros(numel(tIdx),1); p01 = mins; p1 = mins; med = mins; low = [];
    for i = 1:numel(tIdx)
        I = double(imread(fn, 2*(tIdx(i)-1) + ch));
        mins(i) = min(I(:)); p01(i) = prctile(I(:), 0.1); p1(i) = prctile(I(:), 1); med(i) = median(I(:));
        low = [low; I(I <= prctile(I(:), 2))]; %#ok<AGROW>
    end
    fprintf('Channel %d (%d time points sampled):\n', ch, numel(tIdx));
    fprintf('  min        (min/med/max over frames): %s\n', mat2str(round(prctile(mins, [0 50 100]))));
    fprintf('  0.1 pct    : %s\n', mat2str(round(prctile(p01, [0 50 100]))));
    fprintf('  1 pct      : %s\n', mat2str(round(prctile(p1, [0 50 100]))));
    fprintf('  median     : %s\n', mat2str(round(prctile(med, [0 50 100]))));
    fprintf('  lowest-2%% pixels: mode %g, median %g, 5th pct %g\n', mode(low), median(low), prctile(low, 5));
    title(sprintf('Ch%d: lowest 2%% pixels', ch)); xlabel('pixel value');
end
for ch = 1:2
    I = double(imread(fn, ch));
end

% Stored trace lengths
fprintf('\nStored trace lengths (repo .mat files):\n');
fs = [dir(fullfile(repoDir, 'in_vitro_experiments-old', '*_data.mat')); dir(fullfile(repoDir, 'in_vivo_experiments-old', '*_data.mat'))];
for k = 1:numel(fs)
    S = struct2cell(load(fullfile(fs(k).folder, fs(k).name))); S = S{1}; S = S(2:end);
    L = cellfun(@numel, S);
    fprintf('  %-28s n=%4d  max=%d  n(>=360)=%d  n(>361)=%d  n(>=180)=%d  min=%d\n', ...
        fs(k).name, numel(L), max(L), sum(L >= 360), sum(L > 361), sum(L >= 180), min(L));
end
