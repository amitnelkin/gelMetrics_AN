function [intensity, background, info] = collect_intensity(file_path, opts)
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
% opts.aperture - radius (px) of a fixed circular aperture centred on the
%               Mosaic sub-pixel position; the spot intensity is the mean of
%               the pixels whose centres lie within it. Default 5 (PSF radius,
%               BatchParticleTracker.ijm). It replaces the mean of the top
%               two multithresh classes, whose pixel set changes with the
%               brightness, so it is not additive in molecule number.
%               0 = thresholded spot mean (R1/R2). The cell (background)
%               measure is unchanged.

if nargin < 2, opts = struct(); end
if ~isfield(opts,'offset'), opts.offset = 200; end
if ~isfield(opts,'aperture'), opts.aperture = 5; end
file_path = [char(file_path), filesep];
tiff_files = dir([file_path,'*.tif']);

% Half-width of the subframe: the crop is (2*7+1) = 15 x 15 pixels
sub_frame_length = 7;

intensity = cell(1,1);
background = cell(1,1);
info = struct('movie',{},'channel',{},'nChannels',{},'nFrames',{},'nTracks',{},'nKept',{});

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
               spot_frame = sub_frame>=th(3);
               cell_frame = sub_frame>=th(1) & sub_frame<th(3);

               if opts.aperture > 0
                   R = ceil(opts.aperture) + 1;
                   cols = max(floor(spot_x(n))-R,1):min(ceil(spot_x(n))+R,size(frame,2));
                   rows = max(floor(spot_y(n))-R,1):min(ceil(spot_y(n))+R,size(frame,1));
                   [cc, rr] = meshgrid(cols, rows);
                   inAp = (cc - spot_x(n)).^2 + (rr - spot_y(n)).^2 <= opts.aperture^2;
                   win = double(frame(rows, cols));
                   spot_intensity(frameIdx(n),1) = mean(win(inAp)) - opts.offset;
               else
                   spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame)) - opts.offset;
               end
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame)) - opts.offset;
            end

            spot_intensity = fillmissing(spot_intensity,'linear');           % fill skipped frames
            background_intensity = fillmissing(background_intensity,'linear');
            intensity = [intensity,{spot_intensity}]; %#ok<AGROW>
            background = [background,{background_intensity}]; %#ok<AGROW>
            nKept = nKept + 1;
        end
    end
    info(end+1) = struct('movie',base,'channel',channel,'nChannels',nChannels, ...
        'nFrames',numel(movArray),'nTracks',numel(xdata),'nKept',nKept); %#ok<AGROW>
end
end
