% EXPERIMENTALANALYSIS  GelMetrics analysis of experimental data.
% For every condition: burst size, burst size beyond noise, turnover beyond noise
% and net change beyond photobleaching, each with a 95% confidence interval over
% movies (fields of view). Then every condition is tested against the first one.
% Figure 1: burst size distributions of the granules and of their burst-free copies
%           (as MS Figs. 6-9). Figure 2: mean burst size, turnover and net change of
%           all conditions (as MS Figs. 6 and 10).
% The manuscript figures themselves are made by PaperFigures_revised.m.
% The original lambda/k0 fit is in legacy/: lambda and k0 cannot be estimated
% separately; their product is the mean burst size reported here.
close all; clearvars; clc;

nCopies = 20;          % burst-free copies per granule (the paper used 1; 20 gives a stable noise floor)
saveFigures = false;   % true: save the figures as PNG in figures/

%% Load the conditions (the first one is the reference for the tests)
% Option 1: datasets of the paper
%   in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%   in vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'  (protein:RNA)
names = {'4PP7', '4PP7-FUS'};
data = cell(1, numel(names));
for c = 1:numel(names)
    data{c} = loadPaperDataset(names{c});
end
% Option 2: your own data, one folder per condition, with the files written by
% Collect_intensity_data.m (data.mat, background.mat, spotsum.mat, movies.mat)
% names = {'condition 1', 'condition 2'};
% data = {loadFolder('C:\data\condition1'), loadFolder('C:\data\condition2')};

%% Analysis
nCond = numel(names);
M = cell(1, nCond);
for c = 1:nCond
    d = data{c};
    M{c} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, 1);
    fprintf('\n%s: %d granules, %d movies, %d negative bursts\n', names{c}, M{c}.nGranules, M{c}.nMovies, M{c}.nBursts);
    fprintf('  brightness (median over movies)      %.0f\n', M{c}.brightness);
    fprintf('  frame noise (%% of granule)           %.1f\n', M{c}.frameNoise);
    fprintf('  mean burst size, lambda*k0 (a.u.)    %.0f [%.0f, %.0f]\n', M{c}.meanBurst, M{c}.meanBurst_lo, M{c}.meanBurst_hi);
    fprintf('  burst size, granules / copies (%%)    %.1f / %.1f\n', M{c}.sizeGranules, M{c}.sizeCopies);
    fprintf('  burst size beyond noise (%% points)   %.1f [%.1f, %.1f]\n', M{c}.sizeBeyond, M{c}.sizeBeyond_lo, M{c}.sizeBeyond_hi);
    fprintf('  turnover beyond noise (%%/h)          %.1f [%.1f, %.1f]\n', M{c}.turnover, M{c}.turnover_lo, M{c}.turnover_hi);
    fprintf('  net change beyond bleaching (%%/h)    %.1f [%.1f, %.1f]\n', M{c}.net, M{c}.net_lo, M{c}.net_hi);
    fprintf('  first vs second 30 min: turnover %.1f -> %.1f (p = %.3f), net %.1f -> %.1f (p = %.3f)\n', ...
        M{c}.turnover1, M{c}.turnover2, M{c}.pHalvesTurnover, M{c}.net1, M{c}.net2, M{c}.pHalvesNet);
end
for c = 2:nCond
    compareConditions(M{c}, M{1}, names{c}, names{1});
end

%% Figure 1: burst size distributions
f1 = figure('Units', 'inches', 'Position', [1 1 2.2 * min(nCond, 3) 2.3 * ceil(nCond / 3)], 'Color', 'w');
for c = 1:nCond
    ax = subplot(ceil(nCond / 3), min(nCond, 3), c);
    plotBurstHist(ax, M{c});
    title(ax, sprintf('%s  (n = %d)', names{c}, numel(vertcat(M{c}.fracReal{:}))), 'Interpreter', 'none');
    if c == 1, ylabel(ax, 'Frequency'); end
    styleAxes(ax, char('a' + c - 1));
    if c == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 6.5); end
end

%% Figure 2: summary of all conditions (95% CI over movies, then granules)
f2 = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
measures = {'meanBurst', 'turnover', 'net'};
labels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
colors = lines(nCond);
for j = 1:3
    ax = subplot(1, 3, j); hold(ax, 'on');
    y = zeros(1, nCond); lo = y; hi = y;
    for c = 1:nCond
        y(c) = M{c}.(measures{j}); lo(c) = M{c}.([measures{j} '_lo']); hi(c) = M{c}.([measures{j} '_hi']);
        bar(ax, c, y(c), 0.6, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    errorbar(ax, 1:nCond, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    yline(ax, 0, 'k', 'LineWidth', 0.5);
    set(ax, 'XTick', 1:nCond, 'XTickLabel', strrep(names, '_', ' ')); xlim(ax, [0.4, nCond + 0.6]);
    ylabel(ax, labels{j});
    styleAxes(ax, char('a' + j - 1));
end

if saveFigures
    if ~exist('figures', 'dir'), mkdir('figures'); end
    exportgraphics(f1, fullfile('figures', 'burst_sizes.png'), 'Resolution', 300);
    exportgraphics(f2, fullfile('figures', 'summary.png'), 'Resolution', 300);
end


%% ------------------------------------------------------------------ functions
function d = loadPaperDataset(name)
% One dataset of the paper: spot mean, background, spot sum and the movie of every trace.
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/';
    file = strrep(name, 'vitro', 'same');                 % e.g. vitro_1_1 -> same_1_1
    info = readtable([folder 'extraction_info_same.csv'], 'Delimiter', ',', 'TextType', 'string');
    info = info(info.dataset == file, :);
    d.spotMean = loadTraces([folder file '_data.mat']);
    d.background = loadTraces([folder file '_background.mat']);
    d.spotSum = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    info = [readtable([folder 'extraction_info.csv'], 'TextType', 'string'); ...
            readtable([folder 'extraction_info_5Qb_10Qb.csv'], 'TextType', 'string')];
    info = info(info.dataset == name, :);
    d.spotMean = loadTraces([folder name '_data.mat']);
    d.background = loadTraces([folder name '_background.mat']);
    d.spotSum = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
d.movie = repelem((1:height(info))', info.nKept);        % the traces are stored movie by movie
end


function d = loadFolder(folder)
% Your own data: the files written by Collect_intensity_data.m in one folder.
d.spotMean = loadTraces(fullfile(folder, 'data.mat'));
d.background = loadTraces(fullfile(folder, 'background.mat'));
d.spotSum = loadTraces(fullfile(folder, 'spotsum.mat'));
m = load(fullfile(folder, 'movies.mat'));
d.movie = m.movieId(:);
end


function traces = loadTraces(file)
% The extraction stores a cell array of traces whose first cell is empty.
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end


function plotBurstHist(ax, M)
% Burst sizes (% of the granule) of the granules (red bars) and of their copies (black outline).
fracReal = 100 * vertcat(M.fracReal{:}); fracCopy = 100 * vertcat(M.fracCopy{:});
edges = 0:2.5:60;
histogram(ax, min(fracReal, 59.99), edges, 'Normalization', 'pdf', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
hold(ax, 'on');
histogram(ax, min(fracCopy, 59.99), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
xlim(ax, [0 60]); xlabel(ax, 'Burst size (% of granule)');
text(ax, 0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
    M.sizeGranules, M.sizeCopies, M.turnover, M.turnover_lo, M.turnover_hi), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 6.5);
end


function styleAxes(ax, panel)
% Manuscript style: boxed axes, ticks inward, Arial 8, bold panel label "(a)".
set(ax, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.6, 'Layer', 'top', ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'TitleFontWeight', 'normal');
ax.Title.FontSize = 7.5;
if isempty(ax.Title.String), title(ax, ' '); end           % leaves room for the panel label
text(ax, -0.2, 1.04, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
    'VerticalAlignment', 'bottom', 'FontName', 'Arial');
end
