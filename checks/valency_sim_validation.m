% VALENCY_SIM_VALIDATION  Does linking k0 across RNA designs by valency
% (analysis/gelmetricsJoint.m) recover lambda where separate fits cannot?
% Reviewer points R2 4b, R1 2. Seeded; same code path as the real data
% (generateSim -> analysis.m per trace -> amplitudes; negative bursts; top 10%
% trim; 15th-percentile truncation; SP 15; window 13; 150 tracks per dataset).
% Cases: pairs (low, high valency) with known lambda_low, lambda_high, k0_low,
% true brightness ratio, and the ratio ASSUMED from valency:
%   ratio 2   (Qb-like 5x vs 10x):  lambda (1,1) (2,2) (4,4) (1,3)
%   ratio 7.5 (PP7-like 4x vs 30x): lambda (1,1) (2,2)
%   partial occupancy: assumed 2, true 1.5, lambda (2,2)
% noise SD 40 and 80; seeds 1:2. Fits: separate MLE per dataset; joint
% 'fixed' (ratio = assumed); joint 'bounded' (1 <= ratio <= assumed).
% Writes analysis_outputs/valency/sim_validation.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'valency');
if ~exist(outDir, 'dir'), mkdir(outDir); end
% lamLo lamHi k0Lo trueRatio assumedRatio
C = [1 1 60 2 2; 2 2 60 2 2; 4 4 60 2 2; 1 3 60 2 2; 1 1 30 7.5 7.5; 2 2 30 7.5 7.5; 2 2 60 1.5 2];
tVec = (0:360) * 10/60; rows = {};
for nz = [40 80]
    for c = 1:size(C,1)
        for s = 1:2
            A = cell(1,2); sg = zeros(1,2);
            lamT = C(c,1:2); k0T = C(c,3) * [1 C(c,4)];
            for i = 1:2
                rng(10000*s + 100*c + 10*i + nz);
                [S, B] = generateSim(150, 361, 3, 0.0026, nz, lamT(i), k0T(i));
                T = gm_amplitudes_by_trace(S, B, tVec, 13, 15);
                T = T([T.len] >= 180);
                A{i} = abs(vertcat(T.ampNeg)); sg(i) = gm_sigma_from_ampvar(vertcat(T.ampVar));
            end
            for i = 1:2
                Rs = gelmetricsMLE(A{i}, sg(i), struct('profile', false));
                rows(end+1,:) = {nz, c, s, 'separate', i, lamT(i), k0T(i), Rs.lambda, Rs.k0}; %#ok<SAGROW>
            end
            for md = {'fixed', 'bounded'}
                Rj = gelmetricsJoint(A, sg, [1 C(c,5)], md{1});
                for i = 1:2
                    rows(end+1,:) = {nz, c, s, md{1}, i, lamT(i), k0T(i), Rj.lambda(i), Rj.k0(i)}; %#ok<SAGROW>
                end
            end
            fprintf('noise %d case %d seed %d done\n', nz, c, s);
        end
    end
end
V = cell2table(rows, 'VariableNames', {'noise','case','seed','fit','dataset','lamTrue','k0True','lamHat','k0Hat'});
writetable(V, fullfile(outDir, 'sim_validation.csv'));
disp('VALENCY_SIM_DONE');
