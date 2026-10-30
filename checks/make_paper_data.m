% MAKE_PAPER_DATA  Copies the primary datasets of the paper into paper_data/<condition>/,
% in the same four files that Collect_intensity_data.m writes for a new experiment:
%   data.mat (intensity), background.mat (background), spotsum.mat (spotSum),
%   movies.mat (movieId, movieNames)
% so ExperimentalAnalysis.m and compareConditions.m load the paper data and new data the same way.
% Sources: in vivo  in_vivo_reextracted/after_R2_offset (spot mean, background) and
%                   after_R2_sums (spot sum); in vitro  in_vitro_reextracted/same_granules_R2.
% The traces are copied unchanged.
cd(fileparts(fileparts(mfilename('fullpath'))));
vivo = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'};
vitro = {'1_1', '10_1', '100_1'};
infoVivo = [readInfo('in_vivo_reextracted/after_R2_offset/extraction_info.csv'); ...
            readInfo('in_vivo_reextracted/after_R2_offset/extraction_info_5Qb_10Qb.csv')];
infoVitro = readInfo('in_vitro_reextracted/same_granules_R2/extraction_info_same.csv');
for d = 1:numel(vivo)
    src = 'in_vivo_reextracted/after_R2_offset/';
    writeCondition(vivo{d}, [src vivo{d} '_data.mat'], [src vivo{d} '_background.mat'], ...
        ['in_vivo_reextracted/after_R2_sums/' vivo{d} '_spotsum_data.mat'], infoVivo(infoVivo.dataset == vivo{d}, :));
end
for d = 1:numel(vitro)
    src = ['in_vitro_reextracted/same_granules_R2/same_' vitro{d}];
    writeCondition(['vitro_' vitro{d}], [src '_data.mat'], [src '_background.mat'], [src '_spotsum_data.mat'], ...
        infoVitro(infoVitro.dataset == ['same_' vitro{d}], :));
end

function writeCondition(name, dataFile, backgroundFile, sumFile, info)
out = fullfile('paper_data', name);
if ~exist(out, 'dir'), mkdir(out); end
s = struct2cell(load(dataFile)); intensity = s{1};
s = struct2cell(load(backgroundFile)); background = s{1};
s = struct2cell(load(sumFile)); spotSum = s{1};
movieId = repelem((1:height(info))', info.nKept);
movieNames = cellstr(info.movie)';
% the first cell is the empty placeholder of the extraction, as in Collect_intensity_data.m
if numel(intensity) ~= numel(movieId) + 1 || numel(background) ~= numel(intensity) || numel(spotSum) ~= numel(intensity)
    error('%s: trace counts do not match', name);
end
save(fullfile(out, 'data.mat'), 'intensity');
save(fullfile(out, 'background.mat'), 'background');
save(fullfile(out, 'spotsum.mat'), 'spotSum');
save(fullfile(out, 'movies.mat'), 'movieId', 'movieNames');
fprintf('%s: %d traces, %d movies\n', name, numel(movieId), height(info));
end

function info = readInfo(file)
opts = detectImportOptions(file, 'Delimiter', ',');
opts = setvartype(opts, {'movie', 'dataset'}, 'string');
info = readtable(file, opts);
end
