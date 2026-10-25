% T3  Sensitivity of the top-ranked lambda and its k0 to the stringency
% parameter (SP), the left-truncation percentile (GOF) and the top trim,
% on the experimental data in the repo. k0 grid 1:500 as in
% ExperimentalAnalysis.m. k0 = 1 or 500 means the fit hit a grid edge.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir); addpath(fullfile(repoDir, 'legacy'));
ws = warning('off', 'all');
fprintf('Scorer: %s (set env GM_SCORER=orig|patched; see gm_score.m)\n', gm_scorer());

datasets = { ...   % name, folder, SP and GOF listed in ExperimentalAnalysis.m
    '4PP7',            'in_vivo_experiments',  10, 15;
    '30PP7-FUS',       'in_vivo_experiments',  12, 10;
    'protein_RNA_1_1', 'in_vitro_experiments', 12, 15};
% Add more rows as needed, e.g. '30PP7', '4PP7-FUS', '5Qb', '10Qb',
% 'protein_RNA_10_1', 'protein_RNA_100_1'.

SPs = [10 12 14]; GOFs = [0 10 15 25]; topTrims = [0 10 20];   % % removed from top
lamGrid = 1:5; K0Vector = 1:500; win = 13; tVec = linspace(0, 60, 360);

rows = {};
for d = 1:size(datasets, 1)
    [S, B] = gm_loadDataset(repoDir, datasets{d,2}, datasets{d,1});
    nLong = sum(cellfun(@numel, S) >= 180);
    fprintf('\n== %s: %d traces, %d with >= 180 frames (repo settings SP=%d, GOF=%d)\n', ...
        datasets{d,1}, numel(S), nLong, datasets{d,3}, datasets{d,4});
    for SP = SPs
        [~, ~, ap] = analysis(S, B, tVec, win, SP, 0, 0);
        ap = ap(2:end);
        fprintf('  SP=%d  nPos=%d |', SP, numel(ap));
        for g = GOFs
            for tt = topTrims
                a = ap;
                if tt > 0, a = rmoutliers(ap, "percentiles", [0, 100 - tt]); end
                [lo, kb] = gm_score(a, lamGrid, K0Vector, 10, g);
                fprintf(' g%d/t%d:lam%d,k%d', g, tt, lo(1), kb(1));
                rows(end+1, :) = {datasets{d,1}, SP, numel(ap), g, tt, lo(1), kb(1)}; %#ok<SAGROW>
            end
        end
        fprintf('\n');
    end
end
warning(ws);

T = cell2table(rows, 'VariableNames', {'dataset','SP','nPos','GOFpct','topTrimPct','lamBest','k0Best'});
writetable(T, ['t3_parameter_sweep_' gm_scorer() '.csv']);
