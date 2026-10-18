% Function generate_movArray reads a tiff file containing a sequence of
% images into a matlab array prior to reading.
%
% Multi-channel TIFFs (bfconvert/ImageJ hyperstacks of NIS .nd2 movies with
% red fluorescence + DIC, pages interleaved c1,c2,c1,c2,...) are reduced to
% one channel, so movArray(t) is time point t. Reading every page as a time
% point made traces alternate between red and DIC pages (the PP7 "zigzag",
% analysis_outputs/zigzag_diagnosis).
%   channel - which channel to read; default [] picks the fluorescence
%             channel automatically (lowest pixel minimum on its first page:
%             EMCCD offset ~200 vs DIC >~1400).

function [movArray, channel, nChannels] = generate_movArray(filePath,fileName,channel)

if nargin < 3, channel = []; end
movArray = struct('data',[]);
fname = fullfile(filePath,fileName);
info = imfinfo(fname);
num_images = numel(info);

nChannels = 1;
if isfield(info,'ImageDescription') && ischar(info(1).ImageDescription)
    tok = regexp(info(1).ImageDescription,'channels=(\d+)','tokens','once');
    if ~isempty(tok), nChannels = str2double(tok{1}); end
end
if isempty(channel)
    channel = 1;
    if nChannels > 1
        mins = zeros(1,nChannels);
        for c = 1:nChannels
            I = imread(fname, c);
            mins(c) = min(I(:));
        end
        [~, channel] = min(mins);
    end
end

pages = channel:nChannels:num_images;
for k = 1:numel(pages)
    movArray(k).data = imread(fname, pages(k));
end


end
