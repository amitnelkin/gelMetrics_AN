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
% default 15), opts.topTrim (10), opts.start (struct from a previous fit),
% opts.mixFree (false): per-dataset fitted burst scatter sigma_i and false-burst
%   component (share w_i, half-normal width h_i), the model chosen in
%   checks/noise_model_validation.m ('mix+free'); opts.sep (struct array of
%   separate gelmetricsMLE mix+free fits, used as starting values);
% opts.h0K (0): > 0 fixes each false-burst width at h0K*sigma_i (calibrated),
%   so only sigma_i and the share w_i (<= 0.5) are fitted per dataset.
% Output R: lambda (1 x m), k0 (1 x m), ratio (1 x m), nll, n (1 x m).
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'truncPercent'), opts.truncPercent = 15; end
if ~isfield(opts, 'topTrim'), opts.topTrim = 10; end
if ~isfield(opts, 'mixFree'), opts.mixFree = false; end
if ~isfield(opts, 'h0K'), opts.h0K = 0; end
nex = 3 - (opts.h0K > 0);   % extras per dataset
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
nb = numel(lb);                                   % index of the last core parameter
if opts.mixFree   % per dataset: log sigma_i, logit w_i, log h_i
    lgt = @(p) log(p/(1-p));
    for i = 1:m
        if opts.h0K > 0
            lb = [lb, log(sigma(i)/5), lgt(0.001)]; ub = [ub, log(5*sigma(i)), lgt(0.5)]; %#ok<AGROW>
        else
            lb = [lb, log(sigma(i)/5), lgt(0.001), log(sigma(i)/3)]; %#ok<AGROW>
            ub = [ub, log(5*sigma(i)), lgt(0.6),   log(5*sigma(i))]; %#ok<AGROW>
        end
    end
end
% bounded: r non-decreasing along the dataset order -> linear inequalities
Aineq = []; bineq = [];
if ~fixed && m > 2
    for i = 3:m
        row = zeros(1, numel(lb)); row(m + 1 + i - 2) = 1; row(m + 1 + i - 1) = -1;   % log r_{i-1} - log r_i <= 0
        Aineq = [Aineq; row]; bineq = [bineq; 0]; %#ok<AGROW>
    end
end
f = @(th) jointNLL(th, W, sigma, Amin, Amax, ratioMax, fixed, m, opts.mixFree, nb, opts.h0K, nex);
% starts: grid over k0_1 and a common lambda, plus opts.start
starts = zeros(0, numel(lb));
for lam0 = [0.3 1 2 4]
    for k0f = [0.25 0.5 1 2]
        k0s = k0f * median(W{1}) * (1 - exp(-lam0)) / lam0;
        th = [log(lam0) * ones(1, m), log(k0s)];
        if ~fixed, th = [th, log(max(1, ratioMax(2:m) / 2))]; end %#ok<AGROW>
        if opts.mixFree, th = [th, extraStart(opts, sigma, m, nex)]; end %#ok<AGROW>
        starts(end + 1, :) = th; %#ok<AGROW>
    end
end
if isfield(opts, 'start') && ~isempty(opts.start)
    s0 = [log(opts.start.lambda), log(opts.start.k0(1))];
    if ~fixed, s0 = [s0, log(opts.start.ratio(2:m))]; end
    if opts.mixFree, s0 = [s0, opts.start.extra]; end
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
if ~fixed, r = [1, exp(th(m + 2:nb))]; end   % ratios only (extras follow)
R.ratio = r; R.k0 = exp(th(m + 1)) * r; R.nll = best.nll; R.n = n; R.mode = mode;
if opts.mixFree
    R.extra = th(nb + 1:end); X = reshape(R.extra, nex, m);
    R.sigmaFit = exp(X(1,:)); R.w0 = 1 ./ (1 + exp(-X(2,:)));
    if nex == 3, R.h0 = exp(X(3,:)); else, R.h0 = opts.h0K * sigma(:)'; end
end
end

function v = jointNLL(th, W, sigma, Amin, Amax, ratioMax, fixed, m, mixFree, nb, h0K, nex)
lam = exp(th(1:m)); k1 = exp(th(m + 1));
if fixed, r = ratioMax(:)'; else, r = [1, exp(th(m + 2:nb))]; end
v = 0;
for i = 1:m
    if mixFree
        x = th(nb + nex*(i-1) + (1:nex));
        if nex == 3, h = exp(x(3)); else, h = h0K * sigma(i); end
        v = v + gm_mle_nll(lam(i), k1 * r(i), W{i}, exp(x(1)), Amin(i), Amax(i), ...
            struct('w0', 1/(1+exp(-x(2))), 'h0', h));
    else
        v = v + gm_mle_nll(lam(i), k1 * r(i), W{i}, sigma(i), Amin(i), Amax(i));
    end
end
end

function e = extraStart(opts, sigma, m, nex)
% starting extras: from separate mix+free fits when given, else defaults
e = zeros(1, nex*m);
for i = 1:m
    if isfield(opts, 'sep') && ~isempty(opts.sep)
        s = opts.sep(i); v = [log(s.sigmaFit), log(s.w0/(1-s.w0)), log(s.h0)];
    else
        v = [log(sigma(i)), log(0.15/0.85), log(1.5*sigma(i))];
    end
    e(nex*(i-1) + (1:nex)) = v(1:nex);
end
end
