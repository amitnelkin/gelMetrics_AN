% NOISE_CALIBRATION_H0  Width of the false-burst component, in units of the
% quiescent noise sigma (ampVar), from burst-free simulations (generateSim
% type 1, noise SD 20/40/80, 150 tracks, seeds 101-102 as in
% noise_calibration.m). Half-normal MLE: h0K = sqrt(mean((|A|/sigma)^2)).
% Used by gelmetricsMLE opts.h0K (false-burst width fixed at h0K*sigma), so the
% false-burst component can only absorb noise-sized events (R2 6).
% Adds column h0K to analysis_outputs/noise_models/calibration.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_models');
tVec = (0:360) * 10/60; rat = cell(1, 15);
for nz = [20 40 80]
    for s = 101:102
        rng(s * 1000 + nz + 1);
        [S, B] = generateSim(150, 361, 1, 0.0026, nz, 0);
        for SP = 10:15
            T = gm_amplitudes_by_trace(S, B, tVec, 13, SP);
            sg = gm_sigma_from_ampvar(vertcat(T.ampVar));
            rat{SP} = [rat{SP}; abs([vertcat(T.ampNeg); vertcat(T.ampPos)]) / sg];
        end
    end
end
C = readtable(fullfile(outDir, 'calibration.csv'));
C.h0K = cellfun(@(x) sqrt(mean(x.^2)), rat(10:15))';
writetable(C, fullfile(outDir, 'calibration.csv'));
disp(C); disp('H0_DONE');
