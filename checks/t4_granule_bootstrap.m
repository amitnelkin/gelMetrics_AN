% T4  Granule-level bootstrap at the repo's settings: resample whole
% traces (not bursts) with replacement, rerun detection and fitting, and
% report how often each lambda wins and the spread of k0.
% Bursts from one trace are not independent, so resampling bursts would
% understate the uncertainty.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');

datasets = { ...   % name, folder, SP, GOF
    '30PP7-FUS', 'in_vivo_experiments', 12, 10;
    '30PP7',     'in_vivo_experiments', 12, 10};
nBoot = 40;   % increase to 200+ for the paper
lamGrid = 1:5; K0Vector = 1:500; win = 13; tVec = linspace(0, 60, 360);

rng(0);
for d = 1:size(datasets, 1)
    [S, B] = gm_loadDataset(repoDir, datasets{d,2}, datasets{d,1}, 180);
    SP = datasets{d,3}; G = datasets{d,4};
    lamB = zeros(nBoot, 1); kB = zeros(nBoot, 1);
    for b = 1:nBoot
        idx = randi(numel(S), numel(S), 1);
        [~, ~, ap] = analysis(S(idx), B(idx), tVec, win, SP, 0, 0);
        ap = ap(2:end);
        a = rmoutliers(ap, "percentiles", [0, 90]);
        [lo, kb] = gm_mpfit(a, lamGrid, K0Vector, 10, G);
        lamB(b) = lo(1); kB(b) = kb(1);
    end
    fprintf('\n%s (%d traces, SP=%d, GOF=%d, top 10%% trimmed, %d resamples)\n', ...
        datasets{d,1}, numel(S), SP, G, nBoot);
    for l = lamGrid
        fprintf('  lambda=%d wins %d/%d\n', l, sum(lamB == l), nBoot);
    end
    fprintf('  k0 2.5 / 50 / 97.5 percentiles: %s\n', mat2str(round(prctile(kB, [2.5 50 97.5]))));
end
warning(ws);
