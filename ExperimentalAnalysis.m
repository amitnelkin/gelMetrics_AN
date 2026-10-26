% EXPERIMENTALANALYSIS  GelMetrics analysis of experimental data (revised method).
% For every condition: burst size, burst size beyond noise, turnover beyond noise
% and net change beyond photobleaching, each with a 95% confidence interval over
% movies (fields of view), and a movie-level test of every condition against the
% first one. Draws the same kind of figures as the manuscript:
%   Figure 1 - burst size distributions of the granules and of their burst-free
%              copies (as MS Figs. 6(a-c), 7, 8, 9(a, b))
%   Figure 2 - mean burst size, turnover and net change of all conditions
%              (as MS Figs. 6(d-f) and 10(b-d))
% To redraw the manuscript figures exactly, run PaperFigures_revised.m.
% The original lambda/k0 fit is in legacy/: lambda and k0 cannot be estimated
% separately; their product is the mean burst size reported here.
close all; clearvars; clc;

% Analysis settings (the values used in the paper)
opts.window = 13;          % moving-average window (time points)
opts.SP = 12;              % stringency parameter: minimal run of monotone change (time points)
opts.frameInterval = 10;   % seconds between frames - set to your acquisition interval
opts.minFrames = 180;      % analyze traces of at least 30 min
opts.nCopies = 20;         % burst-free copies per granule (the paper used 1; 20 is stable)

% Conditions to analyze: the first one is the reference for the tests.
% Option 1: datasets of the paper
%   in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%   in vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'  (protein:RNA)
conditions = {gm_paper_dataset('4PP7'), gm_paper_dataset('4PP7-FUS')};
% Option 2: your own data, one folder per condition, with the files written by
% Collect_intensity_data.m (data.mat, background.mat, spotsum.mat, movies.mat)
% conditions = {gm_load_extraction('C:\data\condition1', 'condition 1'), ...
%               gm_load_extraction('C:\data\condition2', 'condition 2')};

saveFigures = false;       % true: save Figure 1 and 2 as PNG and PDF in figures/

% Analysis
R = cellfun(@(D) GelMetricsReport(D, opts), conditions, 'UniformOutput', false);
for i = 2:numel(R)
    GelMetricsCompare(R{i}, R{1});
end

% Figure 1: burst size distributions
n = numel(R); nc = min(n, 3);
f1 = figure('Units', 'inches', 'Position', [1 1 2.2 * nc 2.3 * ceil(n / nc)], 'Color', 'w');
t = tiledlayout(f1, ceil(n / nc), nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:n
    ax = nexttile(t);
    nb = gm_plot_bursthist(ax, R{k}, 60, mod(k - 1, nc) == 0);
    title(ax, sprintf('%s  (n = %d)', R{k}.name, nb), 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 7.5);
    gm_style(ax, char('a' + k - 1));
    if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
end

% Figure 2: summary of all conditions (95% CI over movies, then granules)
f2 = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
t = tiledlayout(f2, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
fields = {'meanBurst', 'turnover', 'net'};
labels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
colors = lines(n);
names = cellfun(@(r) strrep(r.name, '_', ' '), R, 'UniformOutput', false);
for k = 1:3
    ax = nexttile(t);
    gm_plot_bars(ax, R, fields{k}, labels{k}, colors, 1:n, names);
    gm_style(ax, char('a' + k - 1));
end

if saveFigures
    outDir = fullfile(fileparts(mfilename('fullpath')), 'figures'); if ~exist(outDir, 'dir'), mkdir(outDir); end
    exportgraphics(f1, fullfile(outDir, 'burst_sizes.png'), 'Resolution', 300);
    exportgraphics(f2, fullfile(outDir, 'summary.png'), 'Resolution', 300);
end
