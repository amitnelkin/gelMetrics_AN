% ESTIMATOR_COMPARE_SIM  Which estimator recovers lambda and k0 best on
% simulations with known truth? Reviewer points R1 2 (biased estimator: weak
% prior on lambda), R2 4a (WS weights are arbitrary), R2 4b (lambda recovery).
% Same code path as the experimental data: generateSim -> analysis.m per trace
% -> amplitudes (top 10% trim, 15th-percentile left truncation, as the runner).
% Estimators, all on the same amplitudes:
%   MLE        gelmetricsMLE, no prior
%   MAP_s1     + lognormal prior on lambda, median 1, sdLog 1   (weak)
%   MAP_s05    + lognormal prior on lambda, median 1, sdLog 0.5 (moderate)
%   WS_<set>   grid method (lambda 1:5, k0 1:500) with weights [MSE QQ KS]:
%              paper .5/.25/.25, equal, MSE only, QQ only, KS only,
%              QQ-heavy .25/.5/.25, KS-heavy .25/.25/.5
% Grid: lambda_true {0.5,1,1.5,2,3,4,5} x noiseStd {20,40,80} x seeds 1:2,
% 150 tracks, k0 = 60, SP 15, window 13. Both burst signs.
% Writes analysis_outputs/estimators/sim_estimates.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'estimators');
if ~exist(outDir, 'dir'), mkdir(outDir); end
W = [.5 .25 .25; 1/3 1/3 1/3; 1 0 0; 0 1 0; 0 0 1; .25 .5 .25; .25 .25 .5];
Wn = {'WS_paper','WS_equal','WS_MSE','WS_QQ','WS_KS','WS_QQheavy','WS_KSheavy'};
priors = {[], [1 1], [1 0.5]}; Pn = {'MLE','MAP_s1','MAP_s05'};
lams = [0.5 1 1.5 2 3 4 5]; noises = [20 40 80]; seeds = 1:2; k0 = 60;
tVec = (0:360) * 10/60; win = 13; SP = 15; trunc = 15;
rows = {};
for nz = noises
    for lam = lams
        for s = seeds
            rng(1000*s + round(10*lam) + nz);
            [S, B] = generateSim(150, 361, 3, 0.0026, nz, lam, k0);
            T = gm_amplitudes_by_trace(S, B, tVec, win, SP);
            T = T([T.len] >= 180);
            sg = gm_sigma_from_ampvar(vertcat(T.ampVar));
            for sgn = {'Neg', 'Pos'}
                A = abs(vertcat(T.(['amp' sgn{1}])));
                for p = 1:numel(priors)
                    R = gelmetricsMLE(A, sg, struct('truncPercent', trunc, 'topTrim', 10, ...
                        'profile', false, 'lamPrior', priors{p}));
                    rows(end+1,:) = {nz, lam, s, sgn{1}, Pn{p}, numel(A), R.lambda, R.k0}; %#ok<SAGROW>
                end
                a = rmoutliers(A, 'percentiles', [0 90]);
                [wl, wk] = gm_ws_weights(a, W, trunc);
                for k = 1:size(W,1)
                    rows(end+1,:) = {nz, lam, s, sgn{1}, Wn{k}, numel(A), wl(k), wk(k)}; %#ok<SAGROW>
                end
            end
            fprintf('noise %d lam %.1f seed %d done\n', nz, lam, s);
        end
    end
end
E = cell2table(rows, 'VariableNames', {'noise','lamTrue','seed','sign','estimator','nBursts','lamHat','k0Hat'});
writetable(E, fullfile(outDir, 'sim_estimates.csv'));
disp('SIM_EST_DONE');
