% FUS_TEST_4X  Is the 4x PP7 vs PP7-FUS difference in burst amplitudes real?
% Manuscript claim: FUS shifts lambda up and k0 down for the 4x RNA (Figs. 7, 9).
% Reviewer points R2 4a/5 (error bars, replicates); answers why separate CIs on
% lambda and k0 overlap while the histograms look different.
% Negative bursts, same code path as checks/mle_real_data.m (window 13, SP 10).
% For each data version (R2 re-extraction; stored MS data):
%  1. descriptives (traces, bursts, quantiles, mean);
%  2. difference in mean and median |amplitude|: granule-level permutation
%     (whole traces are shuffled between the two conditions; 2000 perms);
%  3. two-sample KS on bursts (bursts within a trace are not independent, so
%     the permutation p is the one to trust);
%  4. joint likelihood-ratio test of H0: same (lambda, k0) vs separate, with a
%     COMMON fit window (pooled 90th-percentile top trim, pooled 15th-percentile
%     left truncation) and each condition's own sigma (ampVar). p from the
%     granule permutation (500 perms) and, for reference, chi2(2).
% Writes analysis_outputs/fus_test_4x/fus_test_4x.csv and amplitudes_4x.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'fus_test_4x');
if ~exist(outDir, 'dir'), mkdir(outDir); end
versions = {'R2_reextracted', 'in_vivo_reextracted/after_R2_offset';
            'stored_MS',      'in_vivo_experiments-old'};
tVec = (0:360) * 10/60; win = 13; SP = 10; nPermMean = 2000; nPermLR = 500;
rows = {}; ampRows = {};
for v = 1:size(versions,1)
    rng(7);
    T = cell(1,2); names = {'4PP7','4PP7-FUS'};
    for g = 1:2
        [S, B] = gm_loadDataset(repoDir, versions{v,2}, names{g});
        Tg = gm_amplitudes_by_trace(S, B, tVec, win, SP);
        T{g} = Tg([Tg.len] >= 180);
        for i = 1:numel(T{g})
            a = abs(T{g}(i).ampNeg);
            for k = 1:numel(a), ampRows(end+1,:) = {versions{v,1}, names{g}, i, a(k)}; end %#ok<SAGROW>
        end
    end
    A = cellfun(@(t) abs(vertcat(t.ampNeg)), T, 'UniformOutput', false);
    sg = cellfun(@(t) gm_sigma_from_ampvar(vertcat(t.ampVar)), T);
    % 2. permutation of traces for mean / median difference
    Tall = [T{1}; T{2}]; n1 = numel(T{1});
    stat = @(Ta, Tb) [mean(abs(vertcat(Ta.ampNeg))) - mean(abs(vertcat(Tb.ampNeg))), ...
                      median(abs(vertcat(Ta.ampNeg))) - median(abs(vertcat(Tb.ampNeg)))];
    obs = stat(T{1}, T{2});
    perm = zeros(nPermMean, 2);
    for p = 1:nPermMean
        ix = randperm(numel(Tall));
        perm(p,:) = stat(Tall(ix(1:n1)), Tall(ix(n1+1:end)));
    end
    pMean = mean(abs(perm(:,1)) >= abs(obs(1))); pMed = mean(abs(perm(:,2)) >= abs(obs(2)));
    % 3. KS (burst level)
    [~, pKS, ksD] = kstest2(A{1}, A{2});
    % 4. joint LRT with a common window
    [LR, est] = lrStat(A, sg);
    LRp = zeros(nPermLR, 1);
    for p = 1:nPermLR
        ix = randperm(numel(Tall)); Ta = Tall(ix(1:n1)); Tb = Tall(ix(n1+1:end));
        Ap = {abs(vertcat(Ta.ampNeg)), abs(vertcat(Tb.ampNeg))};
        sp = [gm_sigma_from_ampvar(vertcat(Ta.ampVar)), gm_sigma_from_ampvar(vertcat(Tb.ampVar))];
        LRp(p) = lrStat(Ap, sp, true);
    end
    pLRperm = mean(LRp >= LR); pLRchi2 = 1 - chi2cdf(LR, 2);
    q = @(x) prctile(x, [25 50 75]);
    rows(end+1,:) = {versions{v,1}, numel(T{1}), numel(T{2}), numel(A{1}), numel(A{2}), ...
        mean(A{1}), mean(A{2}), q(A{1}), q(A{2}), sg(1), sg(2), obs(1), pMean, obs(2), pMed, ksD, pKS, ...
        est.Amin, est.Amax, est.lamSep(1), est.k0Sep(1), est.lamSep(2), est.k0Sep(2), est.lam0, est.k00, ...
        LR, pLRperm, pLRchi2}; %#ok<SAGROW>
    fprintf(['%s: traces %d/%d bursts %d/%d | mean %.0f vs %.0f (perm p %.4f) | median %.0f vs %.0f (perm p %.4f) | KS p %.2g\n' ...
             '   common window [%.0f %.0f]: sep 4PP7 %.2f/%.1f, FUS %.2f/%.1f; shared %.2f/%.1f; LR %.1f perm p %.4f chi2 p %.2g\n'], ...
        versions{v,1}, numel(T{1}), numel(T{2}), numel(A{1}), numel(A{2}), mean(A{1}), mean(A{2}), pMean, ...
        median(A{1}), median(A{2}), pMed, pKS, est.Amin, est.Amax, est.lamSep(1), est.k0Sep(1), ...
        est.lamSep(2), est.k0Sep(2), est.lam0, est.k00, LR, pLRperm, pLRchi2);
end
R = cell2table(rows, 'VariableNames', {'data','nTraces4PP7','nTracesFUS','nBursts4PP7','nBurstsFUS', ...
    'mean4PP7','meanFUS','q4PP7','qFUS','sigma4PP7','sigmaFUS','dMean','pMeanPerm','dMedian','pMedianPerm', ...
    'KS_D','KS_p','Amin','Amax','lam4PP7','k04PP7','lamFUS','k0FUS','lamShared','k0Shared','LR','pLRperm','pLRchi2'});
writetable(R, fullfile(outDir, 'fus_test_4x.csv'));
writetable(cell2table(ampRows, 'VariableNames', {'data','dataset','trace','amp'}), fullfile(outDir, 'amplitudes_4x.csv'));
disp('FUS_TEST_DONE');

function [LR, est] = lrStat(A, sg, fast)
% 2*(NLL shared - NLL separate) on a common window
if nargin < 3, fast = false; end
pool = [A{1}; A{2}];
Amax = prctile(pool, 90); pool = pool(pool <= Amax);
Amin = prctile(pool, 15);
W = cellfun(@(a) a(a >= Amin & a <= Amax), A, 'UniformOutput', false);
lb = [log(0.02), log(Amax/400)]; ub = [log(50), log(2*Amax)];
f1 = @(th) gm_mle_nll(exp(th(1)), exp(th(2)), W{1}, sg(1), Amin, Amax);
f2 = @(th) gm_mle_nll(exp(th(1)), exp(th(2)), W{2}, sg(2), Amin, Amax);
f0 = @(th) f1(th) + f2(th);
ns = 30; if fast, ns = 15; end
[th1, n1] = fitLK(f1, lb, ub, ns); [th2, n2] = fitLK(f2, lb, ub, ns);
[th0, n0] = fitLK(f0, lb, ub, ns, [th1; th2]);
LR = max(2*(n0 - n1 - n2), 0);
est = struct('Amin', Amin, 'Amax', Amax, 'lamSep', exp([th1(1) th2(1)]), 'k0Sep', exp([th1(2) th2(2)]), ...
    'lam0', exp(th0(1)), 'k00', exp(th0(2)));
end

function [thBest, nBest] = fitLK(f, lb, ub, ns, extra)
if nargin < 5, extra = zeros(0,2); end
gl = linspace(lb(1), ub(1), ns); gk = linspace(lb(2), ub(2), ns);
G = inf(ns);
for i = 1:ns, for j = 1:ns, G(i,j) = f([gl(i) gk(j)]); end, end
[~, ord] = sort(G(:)); starts = extra;
for q = ord(1:4)', [i, j] = ind2sub(size(G), q); starts(end+1,:) = [gl(i) gk(j)]; end %#ok<AGROW>
starts = min(max(starts, lb + 1e-6), ub - 1e-6);
fo = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
nBest = inf; thBest = starts(1,:);
for s = 1:size(starts,1)
    [th, nll] = fmincon(f, starts(s,:), [], [], [], [], lb, ub, [], fo);
    if nll < nBest, nBest = nll; thBest = th; end
end
end
