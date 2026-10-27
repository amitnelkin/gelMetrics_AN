% SEGAU_LEVEL_TEST  How many negative bursts in intensity units (segAU in burstMeasures.m)
% start at a level <= 0 (granule not above its background), and how much the mean and
% median burst size change if they are dropped, as segReal already does (RESULTS.md section 62).
cd(fileparts(fileparts(mfilename('fullpath'))));
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'};
tVec = (0:400) * 10 / 60;
T = table();
for d = 1:numel(names)
    [S, B] = loadData(names{d});
    keepS = {}; keepB = {};
    for i = 1:numel(S)
        if numel(S{i}) >= 180, keepS{end+1, 1} = double(S{i}(:)); keepB{end+1, 1} = double(B{i}(:)); end %#ok<AGROW>
    end
    oldWarnings = warning('off', 'all');
    [~, ~, ~, ~, ~, ~, seg] = analysis(keepS, keepB, tVec, 13, 12, 0, 0);
    warning(oldWarnings);
    neg = seg(seg(:, 2) < 0, :);
    low = neg(:, 6) <= 0;
    a = abs(neg(:, 5)); aKeep = a(~low);
    fprintf('%-12s %4d negative bursts, %3d with level <= 0 (%.1f%%); mean %.0f -> %.0f a.u., median %.0f -> %.0f a.u.\n', ...
        names{d}, numel(a), sum(low), 100 * mean(low), mean(a), mean(aKeep), median(a), median(aKeep));
    T = [T; table(string(names{d}), numel(a), sum(low), mean(a), mean(aKeep), median(a), median(aKeep), ...
        'VariableNames', {'dataset', 'bursts', 'levelAtOrBelow0', 'meanAll', 'meanKept', 'medianAll', 'medianKept'})]; %#ok<AGROW>
end
writetable(T, fullfile('analysis_outputs', 'trend_shape', 'segau_level.csv'));

function [S, B] = loadData(name)
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/'; file = strrep(name, 'vitro', 'same');
    S = loadTraces([folder file '_data.mat']); B = loadTraces([folder file '_background.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    S = loadTraces([folder name '_data.mat']); B = loadTraces([folder name '_background.mat']);
end
end

function traces = loadTraces(file)
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end
