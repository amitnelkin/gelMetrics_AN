function R = gelmetricsJoint(A, sigma, ratioMax, mode, opts)
% GELMETRICSJOINT  Fits several datasets of the SAME protein together, with the
% unit brightness k0 linked across RNA designs by their valency (number of
% protein binding sites per RNA). Reviewer points R2 4b / R1 2 (lambda and k0
% are not separately identified per dataset; valency is independent
% information that ties the k0s together).
% Model per dataset i: the likelihood of gelmetricsMLE (zero-truncated Poisson
% mixture of Gaussians, own sigma, own fit window), own lambda_i, and
%   k0_i = k0_1 * r_i,   r_1 = 1.
%   mode 'fixed'   : r_i = ratioMax(i)              (full occupancy)
%   mode 'bounded' : r_{i-1} <= r_i <= ratioMax(i), i.e. brightness never
%                    decreases with valency and never exceeds full occupancy
%                    (datasets ordered by valency / protein amount).
% Inputs: A, sigma: cell / vector per dataset (amplitudes, noise SD);
% ratioMax(i) = maximum brightness of dataset i relative to dataset 1
% (valency ratio, ratioMax(1) = 1); opts.truncPercent (per dataset or scalar,
% default 15), opts.topTrim (10), opts.start (struct from a previous fit).
% Output R: lambda (1 x m), k0 (1 x m), ratio (1 x m), nll, n (1 x m).
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'truncPercent'), opts.truncPercent = 15; end
if ~isfield(opts, 'topTrim'), opts.topTrim = 10; end
m = numel(A);
tp = opts.truncPercent; if isscalar(tp), tp = repmat(tp, 1, m); end
W = cell(1, m); Amin = zeros(1, m); Amax = zeros(1, m); n = zeros(1, m);
for i = 1:m
    a = abs(A{i}(:)); a = a(isfinite(a));
    Amax(i) = prctile(a, 100 - opts.topTrim); a = a(a <= Amax(i));
    Amin(i) = prctile(a, tp(i)); W{i} = a(a >= Amin(i)); n(i) = numel(W{i});
end
fixed = strcmp(mode, 'fixed');
% parameters (log): lambda_1..m, k0_1, then (bounded) r_2..m
amax = max(Amax);
lb = [log(0.02) * ones(1, m), log(amax / 400)];
ub = [log(50) * ones(1, m), log(2 * amax)];
if ~fixed
    lb = [lb, zeros(1, m - 1)];                    % log r >= 0
    ub = [ub, log(ratioMax(2:m))];                 % log r <= log ratioMax
end
% bounded: r non-decreasing along the dataset order -> linear inequalities
Aineq = []; bineq = [];
if ~fixed && m > 2
    for i = 3:m
        row = zeros(1, numel(lb)); row(m + 1 + i - 2) = 1; row(m + 1 + i - 1) = -1;   % log r_{i-1} - log r_i <= 0
        Aineq = [Aineq; row]; bineq = [bineq; 0]; %#ok<AGROW>
    end
end
f = @(th) jointNLL(th, W, sigma, Amin, Amax, ratioMax, fixed, m);
% starts: grid over k0_1 and a common lambda, plus opts.start
starts = zeros(0, numel(lb));
for lam0 = [0.3 1 2 4]
    for k0f = [0.25 0.5 1 2]
        k0s = k0f * median(W{1}) * (1 - exp(-lam0)) / lam0;
        th = [log(lam0) * ones(1, m), log(k0s)];
        if ~fixed, th = [th, log(max(1, ratioMax(2:m) / 2))]; end %#ok<AGROW>
        starts(end + 1, :) = th; %#ok<AGROW>
    end
end
if isfield(opts, 'start') && ~isempty(opts.start)
    s0 = [log(opts.start.lambda), log(opts.start.k0(1))];
    if ~fixed, s0 = [s0, log(opts.start.ratio(2:m))]; end
    starts = [s0; starts];
end
starts = min(max(starts, lb + 1e-6), ub - 1e-6);
G = arrayfun(@(i) f(starts(i, :)), 1:size(starts, 1));
[~, ord] = sort(G); starts = starts(ord(1:min(6, numel(ord))), :);
fo = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
best = struct('nll', Inf, 'th', starts(1, :));
for s = 1:size(starts, 1)
    [th, v] = fmincon(f, starts(s, :), Aineq, bineq, [], [], lb, ub, [], fo);
    if v < best.nll, best = struct('nll', v, 'th', th); end
end
th = best.th;
R.lambda = exp(th(1:m)); r = ratioMax(:)';
if ~fixed, r = [1, exp(th(m + 2:end))]; end
R.ratio = r; R.k0 = exp(th(m + 1)) * r; R.nll = best.nll; R.n = n; R.mode = mode;
end

function v = jointNLL(th, W, sigma, Amin, Amax, ratioMax, fixed, m)
lam = exp(th(1:m)); k1 = exp(th(m + 1));
if fixed, r = ratioMax(:)'; else, r = [1, exp(th(m + 2:end))]; end
v = 0;
for i = 1:m
    v = v + gm_mle_nll(lam(i), k1 * r(i), W{i}, sigma(i), Amin(i), Amax(i));
end
end
