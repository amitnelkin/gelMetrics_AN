% SIMULATIONRUNWITHANALYSIS  The revised GelMetrics analysis on simulated granules
% with known dynamics (generateSim.m), through the same code as the experimental data.
% Bursty granules: material leaves and enters in Poisson-distributed numbers of units
% (lambda units of k0 each, mean burst size lambda*k0). Constant granules: no bursts
% at all, so everything detected in them is noise.
close all; clearvars; clc;
rng(1);                    % fixed seed: reproducible simulations

% Simulation parameters
numTracks = 200;           % granules
trackLength = 361;         % time points (60 min at 10 s)
expConst = 0.0026;         % photobleaching rate (1/min), median of the measured backgrounds
noiseStd = 40;             % frame noise (a.u.)
lambda_sim = 2;            % mean number of units per burst
k0_sim = 60;               % intensity of one unit (a.u.)

% Analysis settings (as for the experimental data)
opts.window = 13; opts.SP = 12; opts.frameInterval = 10; opts.minFrames = 180;
opts.nCopies = 20;

[S, Bg] = generateSim(numTracks, trackLength, 3, expConst, noiseStd, lambda_sim, k0_sim);
D.name = sprintf('bursty (lambda %g, k0 %g)', lambda_sim, k0_sim);
D.spotMean = num2cell(S, 2); D.background = num2cell(Bg, 2);
R = GelMetricsReport(D, opts);
fprintf('  true mean burst size lambda*k0 = %g a.u.\n', lambda_sim * k0_sim);

[S, Bg] = generateSim(numTracks, trackLength, 1, expConst, noiseStd, lambda_sim, k0_sim);
D0.name = 'constant (no bursts)';
D0.spotMean = num2cell(S, 2); D0.background = num2cell(Bg, 2);
R0 = GelMetricsReport(D0, opts);
fprintf('  constant granules: burst size beyond noise and turnover beyond noise should be ~0\n');

% Simulations have no movies: every granule is its own replicate here.
