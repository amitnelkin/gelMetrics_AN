% INSPECT_TIFF  Read-only look at ExampleTIF.tif to check the camera dark
% offset (claimed 200 for the Andor iXon Ultra) and the frame count.
% Prints metadata lines that mention offset/baseline/camera/gain and the
% low-end pixel statistics. (Figures dropped: MATLAB rendering hangs headless here.)
% Pixels that receive no light sit at the offset
% plus read noise, so the lowest pixel values bound the offset from above.
clear; close all;
repoDir = '..';
outDir = fullfile(repoDir, 'analysis_outputs', 'baseline', 'tiff_check');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fn = fullfile(repoDir, 'ExampleTIF.tif');
info = imfinfo(fn);
nF = numel(info);
fprintf('Frames (IFDs): %d, size %dx%d, BitDepth %d\n', nF, info(1).Height, info(1).Width, info(1).BitDepth);

% Metadata (OME-XML from bfconvert lives in ImageDescription of IFD 1)
desc = '';
if isfield(info, 'ImageDescription') && ~isempty(info(1).ImageDescription)
    desc = info(1).ImageDescription;
end
fid = fopen(fullfile(outDir, 'ImageDescription.txt'), 'w'); fwrite(fid, desc); fclose(fid);
fprintf('ImageDescription length: %d characters\n', numel(desc));
keys = {'offset','baseline','clamp','gain','camera','iXon','Andor','readout','binning','Exposure','TimeIncrement','DeltaT','PhysicalSize','SizeT'};
for k = 1:numel(keys)
    hits = regexpi(desc, ['[^<>"]{0,60}' keys{k} '[^<>]{0,80}'], 'match');
    for h = 1:min(numel(hits), 5)
        fprintf('  [%s] %s\n', keys{k}, strtrim(hits{h}));
    end
end

% Pixel statistics over up to 40 frames
idx = unique(round(linspace(1, nF, min(nF, 40))));
lowVals = []; frameMins = zeros(numel(idx),1); frameP01 = frameMins; frameMed = frameMins;
for i = 1:numel(idx)
    I = double(imread(fn, idx(i)));
    frameMins(i) = min(I(:)); frameP01(i) = prctile(I(:), 0.1); frameMed(i) = median(I(:));
    lowVals = [lowVals; I(I <= prctile(I(:), 5))]; %#ok<AGROW>
end
fprintf('\nPer-frame min (min/median/max): %s\n', mat2str(prctile(frameMins, [0 50 100])));
fprintf('Per-frame 0.1 pct: %s\n', mat2str(round(prctile(frameP01, [0 50 100]))));
fprintf('Per-frame median: %s\n', mat2str(round(prctile(frameMed, [0 50 100]))));
fprintf('Lowest-5%% pixels: mode %g, median %g, sd %.1f\n', mode(lowVals), median(lowVals), std(lowVals));

xline([100 200], '--', {'100','200'});
