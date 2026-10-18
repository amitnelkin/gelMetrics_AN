% MLE_SIM_VALIDATION  Validates gelmetricsMLE.m on seeded simulations with
% known lambda and k0 (including non-integer lambda), through exactly the
% experimental code path: generateSim -> analysis.m (per trace) -> sigma from
% ampVar -> gelmetricsMLE (profile CIs) + granule bootstrap, side by side with
% the grid-score fit. Reviewer points R2 4b (integer lambda; 2 vs 4) and 4a.
% Settings = experimental runner: window 13, truncation 15%, top trim 10%,
% 10 s frames, 361-frame tracks; SP = 15 (simulation runner). Simulation:
% fixed generateSim (symmetric bursts, bleached baseline, expConst 0.0026),
% noiseStd 40, k0 = 60. Writes analysis_outputs/mle_validation/sim_validation.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'mle_validation');
if ~exist(outDir, 'dir'), mkdir(outDir); end

lams = [0.5 1 1.5 2 3 4 5]; k0 = 60; seeds = 1:3; nTracks = 150; nBoot = 100;  % bootstrap on negative bursts only (the manuscript's fits)
win = 13; SP = 15; tVec = (0:360) * 10/60;
fitSigma = strcmp(getenv('GM_FITSIGMA'), '1');   % 1: sigma free (no bootstrap: 20x slower)
sfx = ''; if fitSigma, sfx = '_fitsigma'; nBoot = 0; end
rows = {};
for lam = lams
    for s = seeds
        rng(100*s + round(10*lam));
        [S, B] = generateSim(nTracks, 361, 3, 0.0026, 40, lam, k0);
        T = gm_amplitudes_by_trace(S, B, tVec, win, SP);
        for sg = {'Neg', 'Pos'}
            o = gm_fit_dataset(T, struct('sign', sg{1}, 'nBoot', nBoot*strcmp(sg{1},'Neg'), 'seed', s, 'fitSigma', fitSigma));
            m = o.mle;
            rows(end+1,:) = {lam, s, sg{1}, o.nBursts, m.sigma, m.lambda, m.lamCI(1), m.lamCI(2), ...
                o.boot.lamCI(1), o.boot.lamCI(2), m.k0, m.k0CI(1), m.k0CI(2), o.boot.k0CI(1), o.boot.k0CI(2), ...
                m.lamK0, o.grid.lambda, o.grid.k0}; %#ok<SAGROW>
            fprintf('lam %.1f seed %d %s: n=%d sigma=%.0f | MLE lam %.2f [%.2f %.2f] k0 %.1f [%.1f %.1f] | grid %d/%d\n', ...
                lam, s, sg{1}, o.nBursts, m.sigma, m.lambda, m.lamCI, m.k0, m.k0CI, o.grid.lambda, o.grid.k0);
        end
    end
end
V = cell2table(rows, 'VariableNames', {'lamTrue','seed','sign','nBursts','sigma', ...
    'lamHat','lamLoProf','lamHiProf','lamLoBoot','lamHiBoot','k0Hat','k0LoProf','k0HiProf', ...
    'k0LoBoot','k0HiBoot','lamK0Hat','gridLam','gridK0'});
V.lamCovProf = V.lamLoProf <= V.lamTrue & V.lamTrue <= V.lamHiProf;
V.lamCovBoot = V.lamLoBoot <= V.lamTrue & V.lamTrue <= V.lamHiBoot;
V.k0CovProf = V.k0LoProf <= k0 & k0 <= V.k0HiProf;
V.k0CovBoot = V.k0LoBoot <= k0 & k0 <= V.k0HiBoot;
writetable(V, fullfile(outDir, ['sim_validation' sfx '.csv']));
G = groupsummary(V, {'sign','lamTrue'}, {'median','mean'}, ...
    {'nBursts','lamHat','k0Hat','lamK0Hat','gridLam','gridK0','lamCovProf','lamCovBoot','k0CovProf','k0CovBoot'});
writetable(G, fullfile(outDir, ['sim_validation_summary' sfx '.csv']));
disp(G(:, {'sign','lamTrue','median_nBursts','median_lamHat','median_k0Hat','median_lamK0Hat', ...
    'median_gridLam','median_gridK0','mean_lamCovProf','mean_lamCovBoot','mean_k0CovProf','mean_k0CovBoot'}));
