% TREND_SHAPE_TEST2  Effect of the trend method (cubic vs trenddecomp) on the noise-controlled
% measures, all paper datasets, 20 copies, seed 1 (RESULTS.md section 60).
% Writes analysis_outputs/trend_shape/measures.csv.
cd(fileparts(fileparts(mfilename('fullpath')))); addpath('checks');
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'};
T = table();
for d = 1:numel(names)
    [S, B, G, mv] = loadData(names{d});
    for method = 1:2
        if method == 1
            M = burstMeasures_ref(S, B, G, mv, 20, 1); m = "poly3";
        else
            M = burstMeasures_td(S, B, G, mv, 20, 1); m = "trenddecomp";
        end
        fprintf('%-12s %-11s size beyond %.2f [%.2f, %.2f]  turnover %.1f [%.1f, %.1f]  copies %.2f  noise share %.2f\n', names{d}, m, ...
            M.sizeBeyond, M.sizeBeyond_lo, M.sizeBeyond_hi, M.turnover, M.turnover_lo, M.turnover_hi, M.sizeCopies, M.noiseShare);
        T = [T; table(string(names{d}), m, M.sizeGranules, M.sizeCopies, M.sizeBeyond, M.sizeBeyond_lo, M.sizeBeyond_hi, ...
            M.turnover, M.turnover_lo, M.turnover_hi, M.noiseShare, 'VariableNames', ...
            {'dataset', 'method', 'sizeGranules', 'sizeCopies', 'sizeBeyond', 'sizeBeyond_lo', 'sizeBeyond_hi', 'turnover', 'turnover_lo', 'turnover_hi', 'noiseShare'})]; %#ok<AGROW>
    end
end
writetable(T, fullfile('analysis_outputs', 'trend_shape', 'measures.csv'));

function [S, B, G, mv] = loadData(name)
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/'; file = strrep(name, 'vitro', 'same');
    info = readtable([folder 'extraction_info_same.csv'], 'Delimiter', ',', 'TextType', 'string'); info = info(info.dataset == file, :);
    S = loadTraces([folder file '_data.mat']); B = loadTraces([folder file '_background.mat']); G = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    info = [readtable([folder 'extraction_info.csv'], 'TextType', 'string'); readtable([folder 'extraction_info_5Qb_10Qb.csv'], 'TextType', 'string')];
    info = info(info.dataset == name, :);
    S = loadTraces([folder name '_data.mat']); B = loadTraces([folder name '_background.mat']);
    G = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
mv = repelem((1:height(info))', info.nKept);
end

function traces = loadTraces(file)
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end
