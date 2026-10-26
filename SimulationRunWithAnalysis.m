% SIMULATIONRUNWITHANALYSIS  The GelMetrics analysis on simulated granules with known
% dynamics (generateSim.m), through the same code as the experimental data.
% Bursty granules: material leaves and enters in Poisson-distributed numbers of units
% (lambda units of k0 each, so the mean burst size is lambda*k0).
% Constant granules: no bursts at all, so everything detected in them is noise.
close all; clearvars; clc;
rng(1);                    % fixed seed: reproducible simulations

% Simulation parameters
numTracks = 200;           % granules
trackLength = 361;         % time points (60 min at 10 s)
expConst = 0.0026;         % photobleaching rate (1/min), median of the measured backgrounds
noiseStd = 40;             % frame noise (a.u.)
lambda_sim = 2;            % mean number of units per burst
k0_sim = 60;               % intensity of one unit (a.u.)
nCopies = 20;              % burst-free copies per granule

for simType = [3 1]        % 3 = bursty granules, 1 = constant granules (no bursts)
    [signals, backgrounds] = generateSim(numTracks, trackLength, simType, expConst, noiseStd, lambda_sim, k0_sim);
    spotMean = cell(numTracks, 1); background = cell(numTracks, 1);
    for i = 1:numTracks
        spotMean{i} = signals(i, :)'; background{i} = backgrounds(i, :)';
    end
    % simulations have no movies (every granule is its own replicate) and no spot sum
    M = burstMeasures(spotMean, background, [], [], nCopies, 1);
    if simType == 3
        fprintf('\nBursty granules (lambda %g, k0 %g): %d bursts\n', lambda_sim, k0_sim, M.nBursts);
        fprintf('  mean burst size %.0f [%.0f, %.0f] a.u.  (true lambda*k0 = %g)\n', M.meanBurst, M.meanBurst_lo, M.meanBurst_hi, lambda_sim * k0_sim);
    else
        fprintf('\nConstant granules (no bursts): %d bursts detected (all noise)\n', M.nBursts);
    end
    fprintf('  burst size, granules / copies (%%)    %.1f / %.1f\n', M.sizeGranules, M.sizeCopies);
    fprintf('  burst size beyond noise (%% points)   %.1f [%.1f, %.1f]\n', M.sizeBeyond, M.sizeBeyond_lo, M.sizeBeyond_hi);
    fprintf('  turnover beyond noise (%%/h)          %.1f [%.1f, %.1f]\n', M.turnover, M.turnover_lo, M.turnover_hi);
    fprintf('  net change beyond bleaching (%%/h)    %.1f [%.1f, %.1f]\n', M.net, M.net_lo, M.net_hi);
end
fprintf('\nFor constant granules, the burst size beyond noise and the turnover beyond noise should be ~0.\n');
