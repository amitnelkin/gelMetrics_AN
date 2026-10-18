% T1  Lambda/k0 recovery on simulated data, original vs corrected simulation.
% Uses the repo's own analysis.m; scoring via gm_mpfit (same math as
% modifiedPoissonFit.m, no plots). k0 grid built as in
% SimulationRunWithAnalysis.m (khat - 8 SE : khat + 18 SE).
% True k0 = 60. Watch the lam*k0 column: if it is stable while lam and k0
% swing, the two parameters are not separately identified.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');
fprintf('Scorer: %s (set env GM_SCORER=orig|patched; see gm_score.m)\n', gm_scorer());

nTracks = 300; seeds = 0:2; lamSims = 1:5; lamGrid = 1:5;
win = 13; SP = 15; truncPct = 15; tVec = linspace(0, 60, 360);

rows = {};
for fixBleach = [false true]
    negFrac = 0.6; if fixBleach, negFrac = 0.5; end
    for lamTrue = lamSims
        for s = seeds
            rng(s);
            [S, B] = gm_generateSim(nTracks, 3, 0.0005, 40, lamTrue, fixBleach, negFrac);
            [~, ~, ap, an] = analysis(S, B, tVec, win, SP, 0, 0);
            ap = ap(2:end); an = abs(an(2:end));
            [lp, kp] = gm_score(ap, lamGrid, k0grid(ap), 10, truncPct);
            [ln, kn] = gm_score(an, lamGrid, k0grid(an), 10, truncPct);
            rows(end+1, :) = {fixBleach, lamTrue, s, numel(ap), lp(1), kp(1), lp(1)*kp(1), ...
                              numel(an), ln(1), kn(1), ln(1)*kn(1)}; %#ok<SAGROW>
            fprintf('fix=%d lamTrue=%d seed=%d | pos: lam=%d k0=%d | neg: lam=%d k0=%d\n', ...
                fixBleach, lamTrue, s, lp(1), kp(1), ln(1), kn(1));
        end
    end
end
warning(ws);

T = cell2table(rows, 'VariableNames', {'fixedSim','lamTrue','seed', ...
    'nPos','lamPos','k0Pos','lamK0Pos','nNeg','lamNeg','k0Neg','lamK0Neg'});
disp(T);
writetable(T, ['t1_lambda_recovery_' gm_scorer() '.csv']);

function K = k0grid(a)
    a = abs(a(:));
    kh = var(a) / mean(a);
    se = kh * sqrt(2 / numel(a));
    K = floor(kh - 8*se):ceil(kh + 18*se);
    K = K(K > 0);
end
