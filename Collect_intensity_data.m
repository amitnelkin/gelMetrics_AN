%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads all .tiff files containing microscopy movies, and .mat
% files containing (x,y,frame) data about detected fluorescent granules.
% The data from each movie is loaded, converted to an array of images, and
% the fluorescence data (intensity of the granule, and
% intensity of its surrounding environment) from each granule is collected
% and saved for further analysis.
% The work is done by collect_intensity.m: the original loop, with the fixes
% listed at its top (fluorescence channel only, Mosaic frame/coordinate mapping,
% skipped frames, camera offset, TIFF <-> _xy.mat paired by name).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
clc;

% Path to folder containing tif and mat files
file_path = [uigetdir("*", ''), filesep];

[intensity, background, info, sums] = collect_intensity(file_path);
% Intensity signal gathering is done.

save([file_path, 'data.mat'],'intensity');           % spot mean (camera offset removed)
save([file_path, 'background.mat'],'background');    % surrounding cell / dilute phase
% For ExperimentalAnalysis.m (burstMeasures.m): the integrated granule intensity, and the movie of every trace
% (movies are the independent replicates of the confidence intervals).
spotSum = sums.spotSum;
save([file_path, 'spotsum.mat'],'spotSum');
movieNames = {info.movie};
movieId = repelem((1:numel(info))', [info.nKept]');
save([file_path, 'movies.mat'],'movieId','movieNames');
