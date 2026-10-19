% MLE_REAL_DATA  Likelihood fit (gelmetricsMLE) of every dataset, with profile
% CIs, granule-bootstrap CIs (200 resamples of traces) and the lambda-k0 NLL
% surface, next to the grid-score fit. Reviewer points R2 4a/4b, R2 5, R1 2.
% Data: PP7 from the R2 re-extraction (in_vivo_reextracted/after_R2_offset; multithresh spot, offset 200),
% in vitro and Qbeta from the repo .mat files (not re-extractable here).
% Per-dataset paper settings (SP, truncation GOF) as in run_runners.m; window
% 13; top trim 10%. Same code path as checks/mle_sim_validation.m.
% Writes analysis_outputs/mle_real/mle_real.csv, surfaces.mat, boot.mat.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'mle_real');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset';
D = {'protein_RNA_1_1',   'in_vitro_experiments', 12, 15;
     'protein_RNA_10_1',  'in_vitro_experiments', 12, 10;
     'protein_RNA_100_1', 'in_vitro_experiments', 12, 10;
     '4PP7',              R2,                     10, 15;
     '4PP7-FUS',          R2,                     10, 15;
     '30PP7',             R2,                     12, 10;
     '30PP7-FUS',         R2,                     12, 10;
     '5Qb',               'in_vivo_experiments',  14, 15;
     '10Qb',              'in_vivo_experiments',  14, 15};
tVec = (0:360) * 10/60; win = 13; nBoot = 200;
rows = {}; surfaces = struct(); boots = struct();
for d = 1:size(D,1)
    [S, B] = gm_loadDataset(repoDir, D{d,2}, D{d,1});
    T = gm_amplitudes_by_trace(S, B, tVec, win, D{d,3});
    for sg = {'Neg', 'Pos'}
        tic;
        o = gm_fit_dataset(T, struct('sign', sg{1}, 'truncPercent', D{d,4}, 'nBoot', nBoot, ...
            'seed', d, 'surface', true));
        m = o.mle;
        rows(end+1,:) = {D{d,1}, D{d,2}, sg{1}, D{d,3}, D{d,4}, o.nTraces, o.nBursts, m.n, o.sigma, ...
            m.lambda, m.lamCI(1), m.lamCI(2), m.lamCIopen(1), m.lamCIopen(2), o.boot.lamCI(1), o.boot.lamCI(2), ...
            m.k0, m.k0CI(1), m.k0CI(2), o.boot.k0CI(1), o.boot.k0CI(2), ...
            m.lamK0, o.boot.lamK0CI(1), o.boot.lamK0CI(2), o.grid.lambda, o.grid.k0, toc}; %#ok<SAGROW>
        key = matlab.lang.makeValidName([D{d,1} '_' sg{1}]);
        surfaces.(key) = m.surf;
        boots.(key) = [o.boot.lambda o.boot.k0];
        fprintf('%-18s %s: n=%d sigma=%.0f | MLE lam %.2f [prof %.2f-%.2f | boot %.2f-%.2f] k0 %.1f [boot %.1f-%.1f] | grid %d/%d (%.0f s)\n', ...
            D{d,1}, sg{1}, m.n, o.sigma, m.lambda, m.lamCI, o.boot.lamCI, m.k0, o.boot.k0CI, o.grid.lambda, o.grid.k0, toc);
    end
end
R = cell2table(rows, 'VariableNames', {'dataset','source','sign','SP','truncPct','nTraces180','nBursts','nFit','sigma', ...
    'lambda','lamLoProf','lamHiProf','lamLoOpen','lamHiOpen','lamLoBoot','lamHiBoot', ...
    'k0','k0LoProf','k0HiProf','k0LoBoot','k0HiBoot','lamK0','lamK0LoBoot','lamK0HiBoot','gridLam','gridK0','seconds'});
writetable(R, fullfile(outDir, 'mle_real.csv'));
save(fullfile(outDir, 'surfaces.mat'), '-struct', 'surfaces');
save(fullfile(outDir, 'boot.mat'), '-struct', 'boots');
disp('MLE_REAL_DONE');
