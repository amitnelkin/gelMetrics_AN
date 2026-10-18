%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
% The work is done by collect_intensity.m (fluorescence channel only,
% Mosaic frame/coordinate mapping, TIFF <-> _xy.mat paired by name).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
clc;

% Path to folder containing tif and mat files
file_path = [uigetdir("*", ''), filesep];

[intensity, background] = collect_intensity(file_path);
% Intensity signal gathering is done.

save([file_path, 'data.mat'],'intensity');
save([file_path, 'background.mat'],'background');
