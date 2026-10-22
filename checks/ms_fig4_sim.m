% MS_FIG4_SIM  Data for the revised manuscript Fig. 4 and Sec. III numbers
% (Reviewer 2 minor point: colour-blind-safe Fig. 4; Sec. III counts redone
% with the corrected simulation). Same code path as SimulationRunWithAnalysis.m:
% generateSim (types 3 = bursts, 1 = constant, 2 = gradually increasing),
% 1000 traces x 360 points, noise SD 40, expConst 0.0026 /min, lambda 2,
% analysis.m with window 13 and SP 15. Also counts false events per flat
% trace at SP 10, 12, 15 (the Sec. III calibration). Seeded.
% Output: analysis_outputs/ms_revision/fig4_*.csv
clear; close all;
repoDir = '..'; addpath(repoDir);
outDir = fullfile(repoDir, 'analysis_outputs', 'ms_revision'); if ~exist(outDir, 'dir'), mkdir(outDir); end
nTr = 1000; T = 360; noiseSD = 40; expC = 0.0026; lam = 2; win = 13;
tVec = linspace(0, 60, T);
types = [3 1 2]; names = {'bursty', 'constant', 'increasing'};
A = []; C = []; X = [];
for k = 1:3
    rng(10 + k);
    [S, B] = generateSim(nTr, T, types(k), expC, noiseSD, lam);
    [~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, 15, 0, 0);
    ap = ap(2:end); an = an(2:end);
    A = [A; repmat(k, numel(ap), 1), ap(:); repmat(k, numel(an), 1), an(:)]; %#ok<AGROW>
    C = [C; k, numel(ap), numel(an), sum(ev) / nTr]; %#ok<AGROW>
    % example trace: base, + bleaching, + noise (same seed, same draws)
    rng(100 + k); s0 = generateSim(1, T, types(k), 0, 0, lam);
    rng(100 + k); s1 = generateSim(1, T, types(k), expC, 0, lam);
    rng(100 + k); s2 = generateSim(1, T, types(k), expC, noiseSD, lam);
    X = [X; repmat(k, T, 1), tVec(:), s0(:), s1(:), s2(:)]; %#ok<AGROW>
    fprintf('%s: pos %d, neg %d, events/trace %.2f\n', names{k}, numel(ap), numel(an), sum(ev) / nTr);
end
% calibration on flat signals
rng(21); [S, B] = generateSim(nTr, T, 1, expC, noiseSD, lam);
cal = [];
for SP = [10 12 15]
    [~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, SP, 0, 0);
    cal = [cal; SP, sum(ev) / nTr, numel(ap) - 1, numel(an) - 1]; %#ok<AGROW>
    fprintf('flat SP %d: %.2f events/trace\n', SP, sum(ev) / nTr);
end
writetable(array2table(A, 'VariableNames', {'type', 'amp'}), fullfile(outDir, 'fig4_amplitudes.csv'));
writetable(array2table(C, 'VariableNames', {'type', 'nPos', 'nNeg', 'eventsPerTrace'}), fullfile(outDir, 'fig4_counts.csv'));
writetable(array2table(X, 'VariableNames', {'type', 't', 'base', 'bleach', 'noisy'}), fullfile(outDir, 'fig4_traces.csv'));
writetable(array2table(cal, 'VariableNames', {'SP', 'eventsPerTrace', 'nPos', 'nNeg'}), fullfile(outDir, 'fig4_calibration.csv'));
