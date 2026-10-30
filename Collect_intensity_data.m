%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% offset  camera dark level subtracted from the spot and cell means (default 200,
%         Andor iXon EMCCD).
%
% Outputs:
%   signal.mat      spot mean per frame (intensity)
%   background.mat  mean of the surrounding cell pixels per frame (background)
%   spotsum.mat     the integrated granule intensity (spot pixels above the cell level)
%   movies.mat      the movie of every trace (movies are the replicates)
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

% Camera dark level; IN OUR SETUP (Andor iXon EMCCD) THIS IS EQUALS TO 200
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
    movArray = generate_movArray(file_path,tiff_files(f).name);
    
    % load (x,y) coordinates of spots detected in imageJ
    xydata = load([file_path, base, '_xy.mat']);           
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    movieNames{end+1} = base;

    % Start iterating through spots and gather intensity
    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)) + 1;   
        spot_x = cell2mat(xdata(t)) + 1;            
        spot_y = cell2mat(ydata(t)) + 1;            
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            frameIdx = spot_frames - spot_frames(1) + 1;
            spot_intensity = nan(frameIdx(end),1);
            background_intensity = nan(frameIdx(end),1);
            spot_sum = nan(frameIdx(end),1);

            % iterate over frames and get intensity data
            for n = 1:length(spot_x)
               frame = movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordinates
               sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                   min(spot_y(n)+sub_frame_length,size(frame,1)),... 
                                   max(spot_x(n)-sub_frame_length,1),...
                                   min(spot_x(n)+sub_frame_length,size(frame,2))];   
               sub_frame_bounds = round(sub_frame_bounds);
               sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
               % sub frame is ready from this point

               % Separate to spot, and 2 levels of background
               th = multithresh(sub_frame,4);
               spot_frame = zeros(size(sub_frame));
               cell_frame = zeros(size(sub_frame));
               back_frame = zeros(size(sub_frame));
               spot_frame(sub_frame>=th(3)) = 1;
               cell_frame(sub_frame>=th(1) & sub_frame<th(3)) = 1;
               back_frame(sub_frame<th(1)) = 1;
               % frame partition to spot,background is done

               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame==1)) - offset;    
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1)) - offset;  
               % integrated granule intensity: the spot pixels above the cell level
               pixels = double(sub_frame);
               spot_sum(frameIdx(n),1) = sum(pixels(spot_frame==1) - mean(pixels(cell_frame==1)));
            end

            intensity = [intensity,{fillmissing(spot_intensity,'linear')}]; 
            background = [background,{fillmissing(background_intensity,'linear')}];
            spotSum = [spotSum,{fillmissing(spot_sum,'linear')}]; 
            movieId(end+1,1) = f; 
        end
    end
end
% Intensity signal gathering is done.

save([file_path, 'data.mat'],'intensity');           % spot mean (camera offset removed)
save([file_path, 'background.mat'],'background');    % surrounding cell / dilute phase
save([file_path, 'spotsum.mat'],'spotSum');          % integrated granule intensity
save([file_path, 'movies.mat'],'movieId','movieNames');   % movie of every trace
