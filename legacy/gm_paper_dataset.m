function D = gm_paper_dataset(name, measure)
% GM_PAPER_DATASET  The datasets of the revised paper, ready for GelMetricsReport.m.
%   In vivo: '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%            (in_vivo_reextracted/: corrected extraction, camera offset removed,
%            spot mean + spot sum, movie of every trace)
%   In vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1' (protein:RNA)
%            (in_vitro_reextracted/same_granules_R2/: the original granules re-read
%            with the corrected extraction, movie of every trace)
% measure (in vitro only): 'spotsum' (default since 2026-10-26, the user's decision, RESULTS.md
% section 53: D.spotSum = integrated granule light above the dilute phase, as in vivo; burst
% fractions, turnover, frame noise and net change use it) or 'spotmean' (the earlier in vitro
% version: fractions relative to spot mean - background; granules spread over time, so its
% 'net loss' reflects the growing area, not lost material).
% D.paperSeed: the random seed of the burst-free copies used for the paper's numbers
% (GelMetricsReport with opts.seed = D.paperSeed and opts.nCopies = 1 reproduces them).
% The original stored traces (in_vivo_experiments/, in_vitro_experiments/) are
% kept for reference only; they contain the errors described in the paper.
repo = fileparts(fileparts(mfilename('fullpath')));   % this file is in legacy/
L = @(p) dropFirst(load(fullfile(repo, p)));
vitro = struct('vitro_1_1', 'same_1_1', 'vitro_10_1', 'same_10_1', 'vitro_100_1', 'same_100_1');
D.name = name;
if isfield(vitro, strrep(name, '-', '_'))
    f = vitro.(name); dir0 = 'in_vitro_reextracted/same_granules_R2';
    D.spotMean = L([dir0 '/' f '_data.mat']); D.background = L([dir0 '/' f '_background.mat']);
    if nargin < 2 || ~strcmpi(measure, 'spotmean'), D.spotSum = L([dir0 '/' f '_spotsum_data.mat']); end
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
seeds = struct('x4PP7', 18, 'x4PP7_FUS', 19, 'x30PP7', 20, 'x30PP7_FUS', 21, 'x5Qb', 22, 'x10Qb', 23, ...
               'vitro_1_1', 18, 'vitro_10_1', 19, 'vitro_100_1', 20);
D.paperSeed = seeds.(matlab.lang.makeValidName(name));
assert(numel(D.movie) == numel(D.spotMean), 'movie list and traces do not match for %s', name);
end


function c = dropFirst(s)
v = struct2cell(s); c = v{1}; c = c(:);
if ~isempty(c) && isempty(c{1}), c = c(2:end); end
end
