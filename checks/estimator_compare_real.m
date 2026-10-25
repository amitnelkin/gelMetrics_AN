% ESTIMATOR_COMPARE_REAL  Real-data counterpart of estimator_compare_sim.m.
% Reviewer points R1 2 (biased estimator), R2 4a (weights, trims), R2 5.
% Data: in vitro and Qbeta stored data; PP7 from the R2 re-extraction
% (red channel, multithresh spot, offset 200). Per-dataset paper SP and
% truncation; window 13; same code path as mle_real_data.m.
% A. Biased estimator: MAP with lognormal prior on lambda (median 1, sdLog 1
%    and 0.5), profile CIs + 100-resample granule bootstrap (negative bursts;
%    positive bursts point estimates). Plain MLE results: mle_real.csv.
% B. WS weight sets (grid method), both signs, point estimates.
% C. Data actually used: traces total / >= 180 frames, bursts detected, after
%    the top-10% trim, after the left truncation; and fits (MLE, WS paper
%    weights) with the trims on/off (top 10% x left truncation).
% Writes analysis_outputs/estimators/real_map.csv, real_weights.csv,
% real_trims.csv, real_counts.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'estimators');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset';
D = {'protein_RNA_1_1','in_vitro_experiments',12,15; 'protein_RNA_10_1','in_vitro_experiments',12,10;
     'protein_RNA_100_1','in_vitro_experiments',12,10; '4PP7',R2,10,15; '4PP7-FUS',R2,10,15;
     '30PP7',R2,12,10; '30PP7-FUS',R2,12,10; '5Qb','in_vivo_experiments',14,15; '10Qb','in_vivo_experiments',14,15};
W = [.5 .25 .25; 1/3 1/3 1/3; 1 0 0; 0 1 0; 0 0 1; .25 .5 .25; .25 .25 .5];
Wn = {'WS_paper','WS_equal','WS_MSE','WS_QQ','WS_KS','WS_QQheavy','WS_KSheavy'};
priors = {[1 1], [1 0.5]}; Pn = {'MAP_s1','MAP_s05'};
tVec = (0:360) * 10/60; win = 13;
mapRows = {}; wRows = {}; tRows = {}; cRows = {};
for d = 1:size(D,1)
    [S, B] = gm_loadDataset(repoDir, D{d,2}, D{d,1});
    T = gm_amplitudes_by_trace(S, B, tVec, win, D{d,3});
    T180 = T([T.len] >= 180);
    sg = gm_sigma_from_ampvar(vertcat(T180.ampVar));
    for sgn = {'Neg', 'Pos'}
        A = abs(vertcat(T180.(['amp' sgn{1}])));
        % C. counts
        a90 = A(A <= prctile(A, 90));
        cRows(end+1,:) = {D{d,1}, sgn{1}, numel(S), numel(T180), numel(A), numel(a90), ...
            sum(a90 >= prctile(a90, D{d,4}))}; %#ok<SAGROW>
        % A. MAP
        for p = 1:numel(priors)
            nb = 100 * strcmp(sgn{1}, 'Neg');
            o = gm_fit_dataset(T180, struct('sign', sgn{1}, 'truncPercent', D{d,4}, 'nBoot', nb, ...
                'seed', d, 'grid', false, 'lamPrior', priors{p}));
            m = o.mle;
            mapRows(end+1,:) = {D{d,1}, sgn{1}, Pn{p}, m.lambda, m.lamCI(1), m.lamCI(2), o.boot.lamCI(1), o.boot.lamCI(2), ...
                m.k0, m.k0CI(1), m.k0CI(2), o.boot.k0CI(1), o.boot.k0CI(2)}; %#ok<SAGROW>
        end
        % B. weights
        [wl, wk] = gm_ws_weights(rmoutliers(A, 'percentiles', [0 90]), W, D{d,4});
        for k = 1:size(W,1), wRows(end+1,:) = {D{d,1}, sgn{1}, Wn{k}, wl(k), wk(k)}; end %#ok<SAGROW>
        % C. trims on/off
        for top = [10 0]
            for tr = [D{d,4} 0]
                Rm = gelmetricsMLE(A, sg, struct('truncPercent', tr, 'topTrim', top, 'profile', false));
                a = A; if top > 0, a = rmoutliers(A, 'percentiles', [0 100-top]); end
                [gl, gk] = gm_ws_weights(a, W(1,:), tr);
                tRows(end+1,:) = {D{d,1}, sgn{1}, top, tr, Rm.n, Rm.lambda, Rm.k0, gl, gk}; %#ok<SAGROW>
            end
        end
    end
    fprintf('%s done\n', D{d,1});
end
writetable(cell2table(mapRows, 'VariableNames', {'dataset','sign','estimator','lambda','lamLoProf','lamHiProf', ...
    'lamLoBoot','lamHiBoot','k0','k0LoProf','k0HiProf','k0LoBoot','k0HiBoot'}), fullfile(outDir, 'real_map.csv'));
writetable(cell2table(wRows, 'VariableNames', {'dataset','sign','weights','lambda','k0'}), fullfile(outDir, 'real_weights.csv'));
writetable(cell2table(tRows, 'VariableNames', {'dataset','sign','topTrimPct','leftTruncPct','nFitMLE', ...
    'mleLambda','mleK0','gridLambda','gridK0'}), fullfile(outDir, 'real_trims.csv'));
writetable(cell2table(cRows, 'VariableNames', {'dataset','sign','tracesTotal','traces180','bursts', ...
    'afterTopTrim','afterLeftTrunc'}), fullfile(outDir, 'real_counts.csv'));
disp('REAL_EST_DONE');
