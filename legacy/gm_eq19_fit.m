function [lamWS, k0WS, lamNLL, k0NLL] = gm_eq19_fit(a)
% GM_EQ19_FIT  Fits the manuscript's Eq. (1.9) to burst amplitudes a, as in
% checks/ms_eq19_fit.py (MS Fig. 5(d); RESULTS.md section 30):
%   p(I) = (1/k0) lambda^(I/k0) exp(-lambda) / Gamma(I/k0 + 1), normalised over I >= 0.
% Grid: lambda 0.05-20 (log-spaced) plus 1..5; k0 2-600 (log-spaced).
% WS  = the manuscript's weighted rank score 0.5 rank(MSE*) + 0.25 rank(QQ) + 0.25 rank(KS)
%       (MSE* on the peak-normalised 10-bin histogram, as in the original code);
% NLL = maximum likelihood (continuous lambda).
% Returns the best lambda and k0 for each score.
LAMS = unique([exp(linspace(log(0.05), log(20), 80)), 1:5]);
K0S = exp(linspace(log(2), log(600), 140));
a = sort(double(a(:)))'; n = numel(a); pq = ((1:n) - 0.5) / n;
edges = linspace(a(1), a(end), 11); cnt = histcounts(a, edges);
ctr = (edges(1:end-1) + edges(2:end)) / 2;
ehi = (1:n) / n; elo = (0:n-1) / n; K = K0S(:);
nL = numel(LAMS); nK = numel(K0S);
NLL = zeros(nL, nK); KS = NLL; QQ = NLL; MSEpk = NLL;
for i = 1:nL
    lam = LAMS(i);
    xg = linspace(0, lam + 15 * sqrt(lam) + 30, 6000);
    pdf = exp(xg * log(lam) - lam - gammaln(xg + 1)); Z = trapz(xg, pdf);
    cdf = [0, cumsum((pdf(2:end) + pdf(1:end-1)) / 2 .* diff(xg))] / Z; pdf = pdf / Z;
    X = a ./ K;                                                     % k0 x n
    F = interp1(xg, cdf, X, 'linear', 1);
    f = interp1(xg, pdf, X, 'linear', 0) ./ K;
    NLL(i, :) = -sum(log(f + 1e-300), 2)';
    KS(i, :) = max(max(abs(F - ehi), [], 2), max(abs(F - elo), [], 2))';
    [cu, iu] = unique(cdf); xu = xg(iu);
    qx = interp1(cu, xu, min(max(pq, cu(1)), cu(end)));
    QQ(i, :) = sqrt(mean((a - K .* qx) .^ 2, 2))';
    m = interp1(xg, pdf, ctr ./ K, 'linear', 0) ./ K;
    mp = m ./ max(max(m, [], 2), 1e-300) * max(cnt);
    MSEpk(i, :) = mean((mp - cnt) .^ 2, 2)';
end
WS = 0.5 * rk(MSEpk) + 0.25 * rk(QQ) + 0.25 * rk(KS);
[i, j] = best(WS); lamWS = LAMS(i); k0WS = K0S(j);
[i, j] = best(NLL); lamNLL = LAMS(i); k0NLL = K0S(j);
end


function R = rk(M)
% rank of every grid cell among all cells, scaled to (0, 1]
[~, o] = sort(M(:)); r = zeros(numel(M), 1); r(o) = 1:numel(M);
R = reshape(r, size(M)) / numel(M);
end


function [i, j] = best(M)
% first minimum in row-major order (lambda first), as numpy argmin
Mt = M'; [~, q] = min(Mt(:)); [j, i] = ind2sub(size(Mt), q);
end
