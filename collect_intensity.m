function [intensity, background, info, sums] = collect_intensity(file_path, opts)
% COLLECT_INTENSITY  Body of Collect_intensity_data.m as a function, so the
% extraction can be rerun headless. For every <name>.tif in file_path it
% loads <name>_xy.mat (read_tables.m output) and returns, per track of >= 24
% points, the spot and surrounding-cell intensity on the full frame axis.
%
% Mapping (revision, analysis_outputs/zigzag_diagnosis):
%   - only the fluorescence channel is read (generate_movArray.m);
%   - Mosaic frames are 0-based time points: page = frame + 1 (was: page =
%     frame, first point dropped; HEAD~ added 1 only when some track started
%     at frame 0);
%   - Mosaic x/y are 0-based pixel coordinates: MATLAB column = x + 1,
%     row = y + 1 (measured: spot centroid - Mosaic = +0.87 col, +0.80 row).
%   - TIFF and _xy.mat are paired by name, not by directory order.
%
% opts.offset - camera dark offset subtracted from the spot and cell means
%               (default 200: Andor iXon EMCCD, fluorescence floor 193-200 in
%               ExampleTIF.tif). The pixel classes are computed on the raw
%               pixels, so only the levels change. 0 = no subtraction (R1).
%
% sums (optional 4th output; the other outputs are unchanged): integrated
% granule intensity per frame, for expressing bursts as a fraction of the
% granule (fraction shed per burst):
%   sums.spotSum  sum over the multithresh spot pixels (top two classes) of
%                 (pixel - mean of the cell class); the offset cancels
%   sums.cropSum  sum over the whole 15x15 crop of (pixel - median of the crop
%                 border); mask-free check
%   sums.spotArea number of spot pixels
% Diagnostics (opts.sumDiag = true; checks/crop_sum_diagnosis.m), all with the
% offset cancelling:
%   sums.cropSumW{k}  box sums at half-widths opts.cropHW(k) (default 7 10 14 20),
%                     each minus the median of its own border
%   sums.ringSum      sum within 6 px of the spot centre minus the median of the
%                     ring 8-11 px (granule-only, local background)
%   sums.nearCell     mean of the ring 4-8 px around the centre minus the ring
%                     8-11 px median (cell light just around the granule)
%
% opts.thrWin  steady thresholds (revision test, RESULTS.md section 23): 0 (default,
%              published) = multithresh recomputed on every frame. > 0 = the four
%              multithresh thresholds of each track are first computed per frame,
%              then replaced by their running median over opts.thrWin frames, so
%              the spot/cell masks no longer jump with each frame's noise but still
%              follow bleaching.

if nargin < 2, opts = struct(); end
if ~isfield(opts,'offset'), opts.offset = 200; end
if ~isfield(opts,'sumDiag'), opts.sumDiag = false; end
if ~isfield(opts,'cropHW'), opts.cropHW = [7 10 14 20]; end
if ~isfield(opts,'thrWin'), opts.thrWin = 0; end
file_path = [char(file_path), filesep];
tiff_files = dir([file_path,'*.tif']);

% Half-width of the subframe: the crop is (2*7+1) = 15 x 15 pixels
sub_frame_length = 7;

intensity = cell(1,1);
background = cell(1,1);
info = struct('movie',{},'channel',{},'nChannels',{},'nFrames',{},'nTracks',{},'nKept',{});
sums = struct('spotSum', {cell(1,1)}, 'cropSum', {cell(1,1)}, 'spotArea', {cell(1,1)});
if opts.sumDiag
    sums.cropSumW = repmat({cell(1,1)}, 1, numel(opts.cropHW));
    sums.ringSum = cell(1,1); sums.nearCell = cell(1,1);
end

for f = 1:numel(tiff_files)
    [~, base] = fileparts(tiff_files(f).name);
    xyName = [file_path, base, '_xy.mat'];
    if ~isfile(xyName)
        error('collect_intensity:noXY', 'No %s_xy.mat for %s', base, tiff_files(f).name);
    end
    [movArray, channel, nChannels] = generate_movArray(file_path,tiff_files(f).name);
    xydata = load(xyName);
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    nKept = 0;

    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)) + 1;   % Mosaic 0-based -> page index
        spot_x = cell2mat(xdata(t)) + 1;             % Mosaic 0-based -> MATLAB column
        spot_y = cell2mat(ydata(t)) + 1;             % Mosaic 0-based -> MATLAB row
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            frameIdx = spot_frames - spot_frames(1) + 1;  % position of each detection on the full frame axis
            spot_intensity = nan(frameIdx(end),1);       % NaN where the tracker skipped frames
            background_intensity = nan(frameIdx(end),1);
            spot_sum = nan(frameIdx(end),1); crop_sum = spot_sum; spot_area = spot_sum;
            if opts.sumDiag
                cw = nan(frameIdx(end), numel(opts.cropHW)); ring_sum = spot_sum; near_cell = spot_sum;
            end

            if opts.thrWin > 0   % steady thresholds: running median of the per-frame ones
                TH = nan(frameIdx(end), 4);
                for n = 1:length(spot_x)
                    TH(frameIdx(n),:) = double(multithresh(crop15(movArray(spot_frames(n)).data, spot_x(n), spot_y(n), sub_frame_length), 4));
                end
                TH = movmedian(TH, opts.thrWin, 1, 'omitnan');
            end

            for n = 1:length(spot_x)
               frame = movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordinates
               sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                   min(spot_y(n)+sub_frame_length,size(frame,1)),...
                                   max(spot_x(n)-sub_frame_length,1),...
                                   min(spot_x(n)+sub_frame_length,size(frame,2))];
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));

               % Separate to spot, and 2 levels of background
               th = multithresh(sub_frame,4);
               if opts.thrWin > 0, th = TH(frameIdx(n),:); end
               spot_frame = sub_frame>=th(3);
               cell_frame = sub_frame>=th(1) & sub_frame<th(3);

               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame)) - opts.offset;
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame)) - opts.offset;
               sf = double(sub_frame); cellMean = mean(sf(cell_frame));
               border = [sf(1,:), sf(end,:), sf(2:end-1,1)', sf(2:end-1,end)'];
               spot_sum(frameIdx(n),1) = sum(sf(spot_frame) - cellMean);
               crop_sum(frameIdx(n),1) = sum(sf(:) - median(border));
               spot_area(frameIdx(n),1) = nnz(spot_frame);
               if opts.sumDiag
                   for k = 1:numel(opts.cropHW)
                       h = opts.cropHW(k);
                       r = max(round(spot_y(n))-h,1):min(round(spot_y(n))+h,size(frame,1));
                       c = max(round(spot_x(n))-h,1):min(round(spot_x(n))+h,size(frame,2));
                       W = double(frame(r, c));
                       bd = [W(1,:), W(end,:), W(2:end-1,1)', W(2:end-1,end)'];
                       cw(frameIdx(n),k) = sum(W(:) - median(bd));
                   end
                   r = max(round(spot_y(n))-12,1):min(round(spot_y(n))+12,size(frame,1));
                   c = max(round(spot_x(n))-12,1):min(round(spot_x(n))+12,size(frame,2));
                   W = double(frame(r, c)); [cc, rr] = meshgrid(c, r);
                   d = hypot(cc - spot_x(n), rr - spot_y(n));
                   bgR = median(W(d >= 8 & d <= 11));
                   ring_sum(frameIdx(n),1) = sum(W(d <= 6) - bgR);
                   near_cell(frameIdx(n),1) = mean(W(d > 4 & d <= 8)) - bgR;
               end
            end

            spot_intensity = fillmissing(spot_intensity,'linear');           % fill skipped frames
            background_intensity = fillmissing(background_intensity,'linear');
            intensity = [intensity,{spot_intensity}]; %#ok<AGROW>
            background = [background,{background_intensity}]; %#ok<AGROW>
            sums.spotSum = [sums.spotSum, {fillmissing(spot_sum,'linear')}];
            sums.cropSum = [sums.cropSum, {fillmissing(crop_sum,'linear')}];
            sums.spotArea = [sums.spotArea, {fillmissing(spot_area,'linear')}];
            if opts.sumDiag
                for k = 1:numel(opts.cropHW)
                    sums.cropSumW{k} = [sums.cropSumW{k}, {fillmissing(cw(:,k),'linear')}];
                end
                sums.ringSum = [sums.ringSum, {fillmissing(ring_sum,'linear')}];
                sums.nearCell = [sums.nearCell, {fillmissing(near_cell,'linear')}];
            end
            nKept = nKept + 1;
        end
    end
    info(end+1) = struct('movie',base,'channel',channel,'nChannels',nChannels, ...
        'nFrames',numel(movArray),'nTracks',numel(xdata),'nKept',nKept); %#ok<AGROW>
end
end

function sub_frame = crop15(frame, x, y, hw)
% the same crop as the main loop
b = round([max(y-hw,1), min(y+hw,size(frame,1)), max(x-hw,1), min(x+hw,size(frame,2))]);
sub_frame = frame(b(1):b(2), b(3):b(4));
end
