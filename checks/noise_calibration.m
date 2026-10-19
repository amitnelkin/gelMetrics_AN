% NOISE_CALIBRATION  Calibrations for the noise models of gelmetricsMLE.m
% (R2 6 noise floor; R2 4b lambda recovery). Seeds 101-102 only, so the
% validation (checks/noise_model_validation.m, seeds 1-2) is held out.
% 1. Noise floor: burst-free simulations, noise SD 20/40/80, 150 tracks; every
%    detected "burst" is false. floorK(SP) = 95th percentile of |false
%    amplitude| / sigma(ampVar), per SP 10..15, from type 1 (constant) signals.
%    Type 2 (linear drift over the whole movie) is reported separately: a slow
%    drift is detected as one huge 'burst', which no amplitude floor can remove.
% 2. 'calib' noise model sigma_n^2 = (a*sigma)^2 + (b*n*k0)^2: burst
%    simulations (lambda 0.5..5, k0 60/120, noise 20/40/80, SP 15) with the noise
%    floor applied; (a, b) maximise the pooled likelihood with lambda fixed at
%    the truth and k0 free per simulation.
% Writes analysis_outputs/noise_models/calibration.csv (floorK per SP; a; b).
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_models');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tVec = (0:360) * 10/60;
% 1. noise floor
rat = cell(1, 15); ratDrift = cell(1, 15);
for typ = [1 2]
    for nz = [20 40 80]
        for s = 101:102
            rng(s * 1000 + nz + typ);
            [S, B] = generateSim(150, 361, typ, 0.0026, nz, 0);
            for SP = 10:15
                T = gm_amplitudes_by_trace(S, B, tVec, 13, SP);
                sg = gm_sigma_from_ampvar(vertcat(T.ampVar));
                a = abs([vertcat(T.ampNeg); vertcat(T.ampPos)]);
                if typ == 1, rat{SP} = [rat{SP}; a / sg]; else, ratDrift{SP} = [ratDrift{SP}; a / sg]; end
            end
        end
    end
end
floorK = cellfun(@(x) prctile(x, 95), rat(10:15));
nFalse = cellfun(@numel, rat(10:15));
fprintf('floorK (SP 10..15): %s  (false bursts pooled: %s)\n', mat2str(floorK, 3), mat2str(nFalse));
% 2. calib model
sims = {};
for lam = [0.5 1 2 3 5]
    for k0 = [60 120]
        for nz = [20 40 80]
            rng(101 * 1000 + round(10*lam) + k0 + nz);
            [S, B] = generateSim(150, 361, 3, 0.0026, nz, lam, k0);
            T = gm_amplitudes_by_trace(S, B, tVec, 13, 15);
            T = T([T.len] >= 180);
            sg = gm_sigma_from_ampvar(vertcat(T.ampVar));
            A = abs(vertcat(T.ampNeg)); Amax = prctile(A, 90); A = A(A <= Amax);
            Amin = floorK(15 - 9) * sg; A = A(A >= Amin);
            sims(end+1,:) = {lam, A, sg, Amin, Amax, k0}; %#ok<SAGROW>
        end
    end
end
pooled = @(p) sum(cellfun(@(lam, A, sg, Amin, Amax, k0) innerK0(lam, A, exp(p(1)) * sg, Amin, Amax, exp(p(2)), k0), ...
    sims(:,1), sims(:,2), sims(:,3), sims(:,4), sims(:,5), sims(:,6)));
best = struct('v', Inf);
for a0 = [1 1.5 2]
    for b0 = [0.05 0.15 0.3]
        [p, v] = fminsearch(pooled, log([a0 b0]), optimset('TolX', 1e-3, 'TolFun', 1e-3, 'MaxFunEvals', 200));
        if v < best.v, best = struct('v', v, 'p', p); end
    end
end
ab = exp(best.p);
fprintf('calib: a = %.3f, b = %.3f\n', ab);
C = table((10:15)', floorK', nFalse', 'VariableNames', {'SP', 'floorK', 'nFalseBursts'});
C.a = repmat(ab(1), height(C), 1); C.b = repmat(ab(2), height(C), 1);
writetable(C, fullfile(outDir, 'calibration.csv'));
disp('CALIB_DONE');

function v = innerK0(lam, A, sg, Amin, Amax, b, k0guess)
g = @(lk) gm_mle_nll(lam, exp(lk), A, sg, Amin, Amax, struct('f', b));
[~, v] = fminbnd(g, log(k0guess / 3), log(k0guess * 3));
end
