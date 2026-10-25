% EXPERIMENTALANALYSIS  GelMetrics analysis of experimental data (revised method).
% For every condition: burst size, burst size beyond noise, turnover beyond noise
% and net change beyond photobleaching, each with a 95% confidence interval over
% movies (fields of view); and a movie-level test between two conditions.
% Noise floor: every granule trace is compared with burst-free copies of itself
% (same photobleaching trend and noise level, no bursts), analyzed identically.
% The original lambda/k0 fit (modifiedPoissonFit.m) is in legacy/: lambda and k0
% cannot be estimated separately; their product is the mean burst size reported here.
close all; clearvars; clc;

% Analysis settings (the values used in the paper)
opts.window = 13;          % moving-average window (time points)
opts.SP = 12;              % stringency parameter: minimal run of monotone change (time points)
opts.frameInterval = 10;   % seconds between frames - set to your acquisition interval
opts.minFrames = 180;      % analyze traces of at least 30 min
opts.nCopies = 20;         % burst-free copies per granule (more copies = less random variation)

% Option 1: the datasets of the paper
%   in vivo:  '4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'
%   in vitro: 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'  (protein:RNA)
A = gm_paper_dataset('4PP7');
B = gm_paper_dataset('4PP7-FUS');

% Option 2: your own data, one folder per condition, with the files written by
% Collect_intensity_data.m (data.mat, background.mat, spotsum.mat, movies.mat)
% A = gm_load_extraction('C:\data\condition1', 'condition 1');
% B = gm_load_extraction('C:\data\condition2', 'condition 2');

RA = GelMetricsReport(A, opts);
RB = GelMetricsReport(B, opts);
T = GelMetricsCompare(RB, RA);

% Burst size distributions of the granules and of their burst-free copies
figure('Color', 'w');
Rs = {RA, RB};
for k = 1:2
    subplot(1, 2, k);
    r = 100 * vertcat(Rs{k}.items.r0); c = 100 * vertcat(Rs{k}.items.c0);
    edges = 0:2.5:60;
    histogram(min(r, 59.9), edges, 'Normalization', 'probability', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none'); hold on;
    histogram(min(c, 59.9), edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1);
    xlabel('Burst size (% of granule)'); ylabel('Frequency');
    title(sprintf('%s: median %.1f%% vs %.1f%%', Rs{k}.name, Rs{k}.sizeGranules, Rs{k}.sizeCopies), 'Interpreter', 'none');
    legend('granules', 'burst-free copies'); box on;
end
