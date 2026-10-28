%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
%
% V2: the original Collect_intensity_data.m with the bug fixes of the
% revision, each marked [fix n] below (RESULTS.md section 64):
%   1. only the fluorescence channel of 2-channel movies is read (generate_movArray.m)
%   2. Mosaic frames start at 0: page = frame + 1, and the first point is kept
%   3. Mosaic x/y start at 0: column = x + 1, row = y + 1
%   4. every value is stored at its own frame; frames the tracker skipped are filled in
%   5. the camera offset is subtracted from the spot and cell means
%   6. each movie is paired with its track file (<name>_xy.mat) by name, not by folder order
%   7. crop limits use the number of rows and columns (length() gives the larger one)
% Also saved, for ExperimentalAnalysis.m:
%   spotsum.mat  the integrated granule intensity (spot pixels above the cell level)
%   movies.mat   the movie of every trace (movies are the replicates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
clc;

% Path to folder containing tif and mat files
file_path = [uigetdir("*", ''), filesep];

% load all .tiff files (each one is paired with its own <name>_xy.mat below)
tiff_files = dir([file_path,'*.tif']);

% Length of the subframe (the crop is 2*7+1 = 15 x 15 pixels)
sub_frame_length = 7;

% Camera dark level (Andor iXon EMCCD), subtracted from the spot and cell means   [fix 5]
offset = 200;

% Initialize variables
intensity = cell(1,1);
background = cell(1,1);
spotSum = cell(1,1);
movieNames = {};
movieId = [];

% Iterate over files, convert to array of images
for f = 1:size(tiff_files,1)
    [~, base] = fileparts(tiff_files(f).name);
    movArray = generate_movArray(file_path,tiff_files(f).name);   % [fix 1] fluorescence channel only
    % load (x,y) coordinates of spots detected in imageJ
    xydata = load([file_path, base, '_xy.mat']);                  % [fix 6] paired by name
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    movieNames{end+1} = base; %#ok<SAGROW>

    % Start iterating through spots and gather intensity
    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)) + 1;   % [fix 2] Mosaic frame 0 is page 1 (first point kept)
        spot_x = cell2mat(xdata(t)) + 1;             % [fix 3]
        spot_y = cell2mat(ydata(t)) + 1;             % [fix 3]
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            % [fix 4] position of each detection on the track's own frame axis;
            % frames the tracker skipped stay NaN and are filled in below
            frameIdx = spot_frames - spot_frames(1) + 1;
            spot_intensity = nan(frameIdx(end),1);
            background_intensity = nan(frameIdx(end),1);
            spot_sum = nan(frameIdx(end),1);

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
               spot_frame = zeros(size(sub_frame));
               cell_frame = zeros(size(sub_frame));
               spot_frame(sub_frame>=th(3)) = 1;
               cell_frame(sub_frame>=th(1) & sub_frame<th(3)) = 1;
               % frame partition to spot,background is done

               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame==1)) - offset;          % [fix 5]
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1)) - offset;    % [fix 5]
               % integrated granule intensity: the spot pixels above the cell level
               % (the offset cancels)
               pixels = double(sub_frame);
               spot_sum(frameIdx(n),1) = sum(pixels(spot_frame==1) - mean(pixels(cell_frame==1)));
            end

            % [fix 4] fill in the frames the tracker skipped
            intensity = [intensity,{fillmissing(spot_intensity,'linear')}]; %#ok<SAGROW>
            background = [background,{fillmissing(background_intensity,'linear')}]; %#ok<SAGROW>
            spotSum = [spotSum,{fillmissing(spot_sum,'linear')}]; %#ok<SAGROW>
            movieId(end+1,1) = f; %#ok<SAGROW>
        end
    end
end
% Intensity signal gathering is done.

save([file_path, 'data.mat'],'intensity');           % spot mean (camera offset removed)
save([file_path, 'background.mat'],'background');    % surrounding cell / dilute phase
save([file_path, 'spotsum.mat'],'spotSum');          % integrated granule intensity
save([file_path, 'movies.mat'],'movieId','movieNames');   % movie of every trace
