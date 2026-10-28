function [intensity, background, info, sums] = collect_intensity(file_path, offset)
% COLLECT_INTENSITY  The body of the original Collect_intensity_data.m, as a function.
% For every movie <name>.tif in file_path it loads <name>_xy.mat (the tracks from
% read_tables.m) and collects, for every track of at least 24 points, the intensity
% of the granule (spot) and of its surrounding cell, frame by frame.
%
% offset  camera dark level subtracted from the spot and cell means (default 200,
%         Andor iXon EMCCD). The pixel classes are found on the raw pixels, so only
%         the levels change. 0 = no subtraction.
%
% Outputs:
%   intensity   spot mean per frame (one cell per track; the first cell is empty, as before)
%   background  mean of the surrounding cell pixels per frame
%   info        per movie: name, channel read, number of channels, frames, tracks, tracks kept
%   sums        sums.spotSum: the integrated granule intensity per frame, i.e. the sum over
%               the spot pixels of (pixel - mean of the cell pixels); the offset cancels.
%               sums.spotArea: number of spot pixels per frame.
%
% Changes from the original script (RESULTS.md section 64), marked [fix n] below:
%   1. only the fluorescence channel of 2-channel movies is read (generate_movArray.m)
%   2. Mosaic frames start at 0: page = frame + 1, and the first point is kept
%   3. Mosaic x/y start at 0: column = x + 1, row = y + 1
%   4. every value is stored at its own frame; frames the tracker skipped are filled in
%   5. the camera offset is subtracted from the spot and cell means
%   6. each movie is paired with its track file by name, not by folder order
%   7. crop limits use the number of rows and columns (length() gives the larger one)
% The revision test options (steady thresholds, fixed mask, crop sums) are in
% checks/collect_intensity_tests.m.

if nargin < 2, offset = 200; end
file_path = [char(file_path), filesep];
tiff_files = dir([file_path,'*.tif']);

% Length of the subframe (the crop is 2*7+1 = 15 x 15 pixels)
sub_frame_length = 7;

% Initialize variables
intensity = cell(1,1);
background = cell(1,1);
sums.spotSum = cell(1,1);
sums.spotArea = cell(1,1);
info = struct('movie',{},'channel',{},'nChannels',{},'nFrames',{},'nTracks',{},'nKept',{});

% Iterate over files, convert to array of images
for f = 1:numel(tiff_files)
    [~, base] = fileparts(tiff_files(f).name);
    xyName = [file_path, base, '_xy.mat'];                          % [fix 6]
    if ~isfile(xyName)
        error('collect_intensity:noXY', 'No %s_xy.mat for %s', base, tiff_files(f).name);
    end
    [movArray, channel, nChannels] = generate_movArray(file_path,tiff_files(f).name);   % [fix 1]
    % load (x,y) coordinates of spots detected in imageJ
    xydata = load(xyName);
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    nKept = 0;

    % Start iterating through spots and gather intensity
    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)) + 1;   % [fix 2] Mosaic frame 0 is page 1
        spot_x = cell2mat(xdata(t)) + 1;             % [fix 3]
        spot_y = cell2mat(ydata(t)) + 1;             % [fix 3]
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            % [fix 4] position of each detection on the track's own frame axis;
            % frames the tracker skipped stay NaN and are filled in at the end
            frameIdx = spot_frames - spot_frames(1) + 1;
            spot_intensity = nan(frameIdx(end),1);
            background_intensity = nan(frameIdx(end),1);
            spot_sum = nan(frameIdx(end),1);
            spot_area = nan(frameIdx(end),1);

            % iterate over frames and get intensity data
            for n = 1:length(spot_x)
               frame = movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordinates
               sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                   min(spot_y(n)+sub_frame_length,size(frame,1)),...   % [fix 7]
                                   max(spot_x(n)-sub_frame_length,1),...
                                   min(spot_x(n)+sub_frame_length,size(frame,2))];     % [fix 7]
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
               % sub frame is ready from this point

               % Separate to spot, and 2 levels of background
               th = multithresh(sub_frame,4);
               spot_frame = sub_frame>=th(3);
               cell_frame = sub_frame>=th(1) & sub_frame<th(3);
               % frame partition to spot,background is done

               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame)) - offset;          % [fix 5]
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame)) - offset;    % [fix 5]
               % integrated granule intensity: spot pixels above the cell level
               sf = double(sub_frame);
               spot_sum(frameIdx(n),1) = sum(sf(spot_frame) - mean(sf(cell_frame)));
               spot_area(frameIdx(n),1) = nnz(spot_frame);
            end

            % [fix 4] fill the frames the tracker skipped
            intensity = [intensity,{fillmissing(spot_intensity,'linear')}]; %#ok<AGROW>
            background = [background,{fillmissing(background_intensity,'linear')}]; %#ok<AGROW>
            sums.spotSum = [sums.spotSum,{fillmissing(spot_sum,'linear')}];
            sums.spotArea = [sums.spotArea,{fillmissing(spot_area,'linear')}];
            nKept = nKept + 1;
        end
    end
    info(end+1) = struct('movie',base,'channel',channel,'nChannels',nChannels, ...
        'nFrames',numel(movArray),'nTracks',numel(xdata),'nKept',nKept); %#ok<AGROW>
end
% Intensity signal gathering is done.
end
