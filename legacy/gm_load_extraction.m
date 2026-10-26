function D = gm_load_extraction(folder, name)
% GM_LOAD_EXTRACTION  Loads the files written by Collect_intensity_data.m in one
% folder (one condition) into the input of GelMetricsReport.m:
%   data.mat (intensity), background.mat (background), and, if present,
%   spotsum.mat (spotSum) and movies.mat (movieId, movieNames).
%   D = gm_load_extraction(folder)
%   D = gm_load_extraction(folder, name)
if nargin < 2, [~, name] = fileparts(char(folder)); end
f = @(x) fullfile(char(folder), x);
D.name = name;
D.spotMean = dropFirst(load(f('data.mat')));
D.background = dropFirst(load(f('background.mat')));
if exist(f('spotsum.mat'), 'file'), D.spotSum = dropFirst(load(f('spotsum.mat'))); end
if exist(f('movies.mat'), 'file')
    m = load(f('movies.mat')); D.movie = m.movieId(:); D.movieNames = m.movieNames;
else
    warning('gm_load_extraction:noMovies', ['No movies.mat in %s: every granule is treated as ' ...
        'its own replicate. Re-run Collect_intensity_data.m to save the movie of every trace.'], folder);
end
end


function c = dropFirst(s)
% the extraction stores a cell array whose first cell is an empty placeholder
v = struct2cell(s); c = v{1}; c = c(:);
if ~isempty(c) && isempty(c{1}), c = c(2:end); end
end
