%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; 
close all;
clc;

% Path to folder containing tif and mat files
file_path = [uigetdir("*", ''), filesep];

% load all .tiff and all .mat files
tiff_files = dir([file_path,'*.tif']);
mat_files = dir([file_path,'*_xy.mat']); % only read_tables outputs (was *.mat: picked up data.mat/background.mat)

% Half-width of the subframe: the crop is (2*7+1) = 15 x 15 pixels
sub_frame_length = 7;

% Initialize variables
intensity = cell(1,1);
background = cell(1,1);

% Iterate over files, convert to array of images
for f = 1:size(tiff_files,1)
    movArray = generate_movArray(file_path,tiff_files(f).name);
    % load (x,y) coordinates of spots detected in imageJ
    xydata = load([file_path,mat_files(f).name]);
    xdata = xydata.x;
    ydata = xydata.y;
    framedata = xydata.frame;
    allFrames = cell2mat(framedata(:));
    frameOffset = double(min(allFrames) == 0); % Mosaic frame numbers start at 0; MATLAB indices at 1
    
    % Start iterating through spots and gather intensity
    for t = 1:length(xdata)
        spot_frames = cell2mat(framedata(t)); 
        spot_frames = spot_frames + frameOffset; % convert to 1-based (was: dropped first point, kept 0-based index)
        spot_x = cell2mat(xdata(t));
        spot_y = cell2mat(ydata(t));
        if length(spot_x)>=24 % continue if spot exists for >24 frames
            % iterate over frames and get intensity data
            frameIdx = spot_frames - spot_frames(1) + 1;  % position of each detection on the full frame axis
            spot_intensity = nan(frameIdx(end),1);       % NaN where the tracker skipped frames
            background_intensity = nan(frameIdx(end),1);
            
            for n = 1:length(spot_x)
               frame = movArray(spot_frames(n)).data; % complete current frame
               % Extract sub frame around the spot coordiantes
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
               
               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame==1));
               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1));
               
            end
                      
           spot_intensity = fillmissing(spot_intensity,'linear');           % fill skipped frames
           background_intensity = fillmissing(background_intensity,'linear');
           intensity = [intensity,{spot_intensity}];
           background = [background,{background_intensity}];
        end
    end
end
% Intensity signal gathering is done. 

save([file_path, 'data.mat'],'intensity');
save([file_path, 'background.mat'],'background');
