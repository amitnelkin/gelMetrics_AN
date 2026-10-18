% MLE_UNIT_TEST  gelmetricsMLE.m recovers lambda and k0 from amplitudes drawn
% directly from its own model (R2 4b: non-integer lambda). Seeded. This tests
% the fitter only; the full-pipeline validation is checks/mle_sim_validation.m.
% Writes analysis_outputs/mle_validation/unit_test.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoDir, 'analysis'));
outDir = fullfile(repoDir, 'analysis_outputs', 'mle_validation');
if ~exist(outDir, 'dir'), mkdir(outDir); end
lams = [0.5 1 1.5 2 3 4 5]; k0 = 60; sigmas = [15 30]; nA = 300; seeds = 1:5;
rows = {};
for sg = sigmas
    for lam = lams
        for s = seeds
            rng(1000*s + round(10*lam) + sg);
            N = zeros(0,1);
            while numel(N) < nA
                x = poissrnd(lam, 4*nA, 1); N = [N; x(x >= 1)]; %#ok<AGROW>
            end
            A = N(1:nA)*k0 + sg*randn(nA,1);
            R = gelmetricsMLE(A, sg, struct('truncPercent', 15, 'topTrim', 10));
            rows(end+1,:) = {sg, lam, s, R.n, R.lambda, R.lamCI(1), R.lamCI(2), ...
                R.k0, R.k0CI(1), R.k0CI(2), R.lamCI(1) <= lam && lam <= R.lamCI(2), ...
                R.k0CI(1) <= k0 && k0 <= R.k0CI(2)}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'sigma','lamTrue','seed','n','lamHat','lamLo','lamHi', ...
    'k0Hat','k0Lo','k0Hi','lamCovered','k0Covered'});
writetable(T, fullfile(outDir, 'unit_test.csv'));
G = groupsummary(T, {'sigma','lamTrue'}, {'median','mean'}, {'lamHat','k0Hat','lamCovered','k0Covered'});
disp(G(:, {'sigma','lamTrue','median_lamHat','median_k0Hat','mean_lamCovered','mean_k0Covered'}));
