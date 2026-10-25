function D = gm_paper_dataset(name)
% GM_PAPER_DATASET  The datasets of the revised paper, ready for GelMetricsReport.m.
%   In vivo: '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%            (in_vivo_reextracted/: corrected extraction, camera offset removed,
%            spot mean + spot sum, movie of every trace)
%   In vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1' (protein:RNA)
%            (in_vitro_reextracted/same_granules_R2/: the original granules re-read
%            with the corrected extraction, movie of every trace)
% The original stored traces (in_vivo_experiments/, in_vitro_experiments/) are
% kept for reference only; they contain the errors described in the paper.
repo = fileparts(mfilename('fullpath'));
L = @(p) dropFirst(load(fullfile(repo, p)));
vitro = struct('vitro_1_1', 'same_1_1', 'vitro_10_1', 'same_10_1', 'vitro_100_1', 'same_100_1');
D.name = name;
if isfield(vitro, strrep(name, '-', '_'))
    f = vitro.(name); dir0 = 'in_vitro_reextracted/same_granules_R2';
    D.spotMean = L([dir0 '/' f '_data.mat']); D.background = L([dir0 '/' f '_background.mat']);
    info = readtable(fullfile(repo, dir0, 'extraction_info_same.csv'), 'Delimiter', ',', 'TextType', 'string');
    info = info(info.dataset == f, :);
else
    r2 = 'in_vivo_reextracted/after_R2_offset'; su = 'in_vivo_reextracted/after_R2_sums';
    D.spotMean = L([r2 '/' name '_data.mat']); D.background = L([r2 '/' name '_background.mat']);
    D.spotSum = L([su '/' name '_spotsum_data.mat']);
    info = [readtable(fullfile(repo, r2, 'extraction_info.csv'), 'TextType', 'string'); ...
            readtable(fullfile(repo, r2, 'extraction_info_5Qb_10Qb.csv'), 'TextType', 'string')];
    info = info(info.dataset == name, :);
end
if isempty(info), error('gm_paper_dataset:unknown', 'Unknown dataset "%s".', name); end
D.movie = repelem((1:height(info))', info.nKept);
assert(numel(D.movie) == numel(D.spotMean), 'movie list and traces do not match for %s', name);
end


function c = dropFirst(s)
v = struct2cell(s); c = v{1}; c = c(:);
if ~isempty(c) && isempty(c{1}), c = c(2:end); end
end
