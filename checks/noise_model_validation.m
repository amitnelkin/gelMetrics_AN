function noise_model_validation(k0list, tag, only)
% NOISE_MODEL_VALIDATION  Which noise / false-burst model recovers lambda?
% R2 4b (lambda recovery), R2 6 (noise floor). Held-out seeds (calibration used
% 101-102). Same code path as the data: generateSim -> analysis.m per trace ->
% amplitudes (top 10% trim). lambda {0.5,1,1.5,2,3,4,5} x noise SD {20,40,80} x
% k0 (k0list) x seed 1, both burst signs, SP 15.
% Estimators (gelmetricsMLE options):
%   fixed          sigma from ampVar, 15th-percentile left cut (the current method)
%   free / growN / prop   noise models with one fitted parameter
%   calib+floor    sigma_n^2 = (a sigma)^2 + (b n k0)^2 (a, b calibrated) + noise floor
%   mix, mix+free, mix+prop   + false-burst component (share and width fitted)
%   floor, floor+free         left cut at floorK*sigma (calibrated noise floor)
%   mixH+free, mixH+prop      false-burst component with its width fixed at
%                             h0K*sigma (calibrated, noise_calibration_h0.m)
% only (optional): cell of model names to run (e.g. the mixH models).
% Writes analysis_outputs/noise_models/validation_<tag>.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_models');
C = readtable(fullfile(outDir, 'calibration.csv'));
SP = 15; fk = C.floorK(C.SP == SP); ab = [C.a(1) C.b(1)]; hk = C.h0K(C.SP == SP);
E = {'fixed',       struct('noise','fixed');
     'free',        struct('noise','free');
     'growN',       struct('noise','growN');
     'prop',        struct('noise','prop');
     'calib+floor', struct('noise','calib','calib',ab,'floorK',fk);
     'mix',         struct('noise','fixed','falseBursts','mix');
     'mix+free',    struct('noise','free','falseBursts','mix');
     'mix+prop',    struct('noise','prop','falseBursts','mix');
     'floor',       struct('noise','fixed','floorK',fk);
     'floor+free',  struct('noise','free','floorK',fk);
     'mixH+free',   struct('noise','free','falseBursts','mix','h0K',hk);
     'mixH+prop',   struct('noise','prop','falseBursts','mix','h0K',hk)};
if nargin >= 3, E = E(ismember(E(:,1), only), :); end
tVec = (0:360) * 10/60; rows = {};
for k0 = k0list
    for nz = [20 40 80]
        for lam = [0.5 1 1.5 2 3 4 5]
            rng(1000 + round(10*lam) + k0 + nz);
            [S, B] = generateSim(150, 361, 3, 0.0026, nz, lam, k0);
            T = gm_amplitudes_by_trace(S, B, tVec, 13, SP);
            T = T([T.len] >= 180);
            sg = gm_sigma_from_ampvar(vertcat(T.ampVar));
            for sgn = {'Neg', 'Pos'}
                A = abs(vertcat(T.(['amp' sgn{1}])));
                for e = 1:size(E,1)
                    o = E{e,2}; o.profile = false;
                    R = gelmetricsMLE(A, sg, o);
                    rows(end+1,:) = {k0, nz, lam, sgn{1}, E{e,1}, R.n, R.lambda, R.k0}; %#ok<AGROW>
                end
            end
            fprintf('k0 %d noise %d lam %.1f done\n', k0, nz, lam);
        end
    end
end
writetable(cell2table(rows, 'VariableNames', {'k0True','noise','lamTrue','sign','model','nFit','lamHat','k0Hat'}), ...
    fullfile(outDir, ['validation_' tag '.csv']));
disp('NOISE_VAL_DONE');
end
