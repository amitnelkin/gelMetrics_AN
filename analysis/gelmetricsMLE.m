function R = gelmetricsMLE(A, sigma, opts)
% GELMETRICSMLE  Maximum-likelihood fit of continuous lambda and k0 to burst
% amplitudes. Reviewer points R2 4a/4b (error bars; lambda restricted to
% integers; simulated vs estimated lambda) and R1 point 2. Runs alongside the
% grid-score fit (modifiedPoissonFit.m), which is kept for comparison.
%
% Model: a burst moves N >= 1 units, N ~ zero-truncated Poisson(lambda), and
% its measured amplitude is A ~ Normal(N*k0, sigma^2):
%     p(A) = sum_{n>=1} P(n | lambda, n>=1) * phi((A - n*k0)/sigma) / sigma
% sigma is fixed, from the quiescent-segment amplitudes (ampVar). The fitted
% amplitudes are the ones the grid method fits: |A| after the top trim, left-
% truncated at the truncation percentile. Both cuts are modelled: the
% likelihood is divided by P(Amin < A < Amax | lambda, k0).
%
% Inputs
%   A      burst amplitudes (one sign; absolute values are used)
%   sigma  amplitude noise SD (gm_sigma_from_ampvar(ampVar))
%   opts   .truncPercent  left-truncation percentile (default 15)
%          .topTrim       top percent removed before fitting (default 10)
%          .profile       profile likelihoods + 95% CIs (default true)
%          .surface       lambda-k0 NLL surface (default false)
%          .nGrid         surface grid size (default 60)
% Output R: lambda, k0, lamK0, nll, n, Amin, Amax, sigma, exitflag, lamAtBound;
%   if profile: lamCI, k0CI (95%: 2*dNLL <= 3.84), lamCIopen / k0CIopen
%               ([low high] true = the interval reaches the search bound, so
%               the data do not bound it), profLam / profK0 (tables);
%   if surface: surf.lam, surf.k0, surf.nll (NLL - min).
% Optimisation is in log(lambda), log(k0): a 40 x 40 grid scan seeds fmincon
% from the best grid cells plus fixed starts (the likelihood is multimodal
% when sigma << k0: k0/m sub-harmonics).

if nargin < 3, opts = struct(); end
opts = setdef(opts, 'truncPercent', 15);
opts = setdef(opts, 'topTrim', 10);
opts = setdef(opts, 'profile', true);
opts = setdef(opts, 'surface', false);
opts = setdef(opts, 'nGrid', 60);
opts = setdef(opts, 'scanN', 40);       % seed-grid size per axis
opts = setdef(opts, 'start', []);       % extra start [lambda k0] (bootstrap warm start)

A = abs(A(:));
A = A(isfinite(A));
if opts.topTrim > 0
    Amax = prctile(A, 100 - opts.topTrim);
    A = A(A <= Amax);
else
    Amax = Inf;
end
Amin = prctile(A, opts.truncPercent);
A = A(A >= Amin);

R = struct('n', numel(A), 'Amin', Amin, 'Amax', Amax, 'sigma', sigma, ...
    'lambda', NaN, 'k0', NaN, 'lamK0', NaN, 'nll', NaN, 'exitflag', -99, 'lamAtBound', false);
if numel(A) < 5, return, end

lb = [log(0.02), log(max(A)/400)];
ub = [log(50),   log(2*max(A))];
f = @(th) gm_mle_nll(exp(th(1)), exp(th(2)), A, sigma, Amin, Amax);

% grid scan
ns = opts.scanN;
gl = linspace(lb(1), ub(1), ns); gk = linspace(lb(2), ub(2), ns);
G = inf(ns);
for i = 1:ns, for j = 1:ns, G(i,j) = f([gl(i) gk(j)]); end, end
[~, ord] = sort(G(:));
starts = zeros(0,2);
for q = ord(1:min(6, ceil(ns/6)))'
    [i, j] = ind2sub(size(G), q); starts(end+1,:) = [gl(i) gk(j)]; %#ok<AGROW>
end
if isempty(opts.start)
    for lam0 = [0.5 1 2 4 8]
        starts(end+1,:) = [log(lam0), log(mean(A)*(1-exp(-lam0))/lam0)]; %#ok<AGROW>
    end
else
    starts(end+1,:) = log(opts.start(:)');
end
starts = min(max(starts, lb + 1e-6), ub - 1e-6);

fo = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'OptimalityTolerance', 1e-9, 'StepTolerance', 1e-10);
best = struct('nll', Inf, 'th', starts(1,:), 'ef', -99);
for s = 1:size(starts,1)
    [th, nll, ef] = fmincon(f, starts(s,:), [], [], [], [], lb, ub, [], fo);
    if nll < best.nll, best = struct('nll', nll, 'th', th, 'ef', ef); end
end
R.lambda = exp(best.th(1)); R.k0 = exp(best.th(2)); R.lamK0 = R.lambda * R.k0;
R.nll = best.nll; R.exitflag = best.ef;
R.lamAtBound = any(abs(best.th(1) - [lb(1) ub(1)]) < 1e-3);

crit = chi2inv(0.95, 1) / 2;      % 1.92
if opts.profile
    % profile over lambda (k0 free) and over k0 (lambda free)
    [R.lamCI, R.lamCIopen, PL, nll1] = profileCI(f, best, 1, lb, ub, crit);
    [R.k0CI,  R.k0CIopen,  PK, nll2] = profileCI(f, best, 2, lb, ub, crit);
    if min(nll1, nll2) < R.nll - 1e-6    % profile found a better optimum
        R.note = 'profile found a lower NLL than the multistart optimum';
    end
    R.profLam = array2table([exp(PL(:,1)) PL(:,2)], 'VariableNames', {'lambda','dNLL'});
    R.profK0  = array2table([exp(PK(:,1)) PK(:,2)], 'VariableNames', {'k0','dNLL'});
end
if opts.surface
    lamS = exp(linspace(lb(1), ub(1), opts.nGrid));
    k0S = exp(linspace(log(R.lamK0/60), log(min(R.lamK0*4, 2*max(A))), opts.nGrid));
    S = zeros(numel(lamS), numel(k0S));
    for i = 1:numel(lamS)
        for j = 1:numel(k0S)
            S(i,j) = gm_mle_nll(lamS(i), k0S(j), A, sigma, Amin, Amax);
        end
    end
    R.surf = struct('lam', lamS, 'k0', k0S, 'nll', S - R.nll);
end
end

function [ci, open, P, nllBest] = profileCI(f, best, p, lb, ub, crit)
% March outward from the MLE in parameter p (log scale), minimising over the
% other parameter q at each step (warm start from the previous step plus a
% coarse scan), until dNLL > crit or the bound; then root-find the crossing.
q = 3 - p;
step = 0.04;
P = [best.th(p), 0]; nllBest = best.nll;
ci = exp([best.th(p) best.th(p)]); open = [false false];
for dirn = [-1 1]
    x = best.th(p); qPrev = best.th(q); dPrev = 0;
    side = (dirn > 0) + 1;
    while true
        xn = x + dirn*step;
        if xn < lb(p) || xn > ub(p)
            open(side) = true; ci(side) = exp(min(max(xn, lb(p)), ub(p))); break
        end
        [v, qv] = innerMin(f, p, xn, qPrev, lb(q), ub(q));
        nllBest = min(nllBest, v);
        d = v - best.nll;
        P(end+1,:) = [xn d]; %#ok<AGROW>
        if d > crit
            g = @(xx) innerMin(f, p, xx, qv, lb(q), ub(q)) - best.nll - crit;
            try
                xr = fzero(g, sort([x xn]));
            catch
                xr = x + dirn*step*(crit - dPrev)/max(d - dPrev, eps);
            end
            ci(side) = exp(xr); break
        end
        x = xn; qPrev = qv; dPrev = d;
        if abs(x - best.th(p)) > 0.5, step = 0.1; end
    end
end
P = sortrows(P, 1);
end

function [v, qv] = innerMin(f, p, xp, q0, lo, hi)
th = zeros(1,2); th(p) = xp;
    function y = g(xq), th(3-p) = xq; y = f(th); end
% local search around the warm start, plus a coarse scan for mode jumps
a = max(q0 - 0.3, lo); b = min(q0 + 0.3, hi);
[qv, v] = fminbnd(@g, a, b, optimset('TolX', 1e-7));
xs = linspace(lo, hi, 30); gs = arrayfun(@g, xs);
[gmin, i] = min(gs);
if gmin < v
    a = xs(max(i-1,1)); b = xs(min(i+1,numel(xs)));
    [q2, v2] = fminbnd(@g, a, b, optimset('TolX', 1e-7));
    if v2 < v, v = v2; qv = q2; end
end
end

function s = setdef(s, f, v)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = v; end
end
