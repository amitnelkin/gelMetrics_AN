
function [lamWS, k0WS, WS, scores] = PoissonFit2WS(data, lambda, K0Vector, numBins)
%   p(I) = (1/k0) lambda^(I/k0) exp(-lambda) / Gamma(I/k0 + 1)
%   on a grid of lambda and k0, scored by the weighted score (WS)
%   WS = 0.5 rank(MSE) + 0.25 rank(QQ) + 0.25 rank(KS)
%   MSE on a `numBins`-bin histogram, QQ and KS on the amplitudes. 
%   Returns the best lambda and k0.
%   Optional: WS, the weighted score of every (lambda, k0) pair (rows lambda, columns k0),
%   and scores, the single scores on the same grid (MSE, QQ, KS, and CvM, the Cramer-von Mises
%   score, used only in AppendixFigures.m).

data = sort(data(:))'; 
n = numel(data);
binEdges = linspace(data(1), data(end), numBins+1); 
binCounts = histcounts(data, binEdges);         % 10-bin histogram
centers = (binEdges(1:end-1) + binEdges(2:end)) / 2;

% Initialize single scores matrices
KS = zeros(numel(lambda), numel(K0Vector)); 
QQ = zeros(numel(lambda), numel(K0Vector));
MSE = zeros(numel(lambda), numel(K0Vector));
CvM = zeros(numel(lambda), numel(K0Vector));

for i = 1:numel(lambda)
    lam = lambda(i);
    % the model (in units of k0), normalised, and its cumulative distribution
    x = linspace(0, lam + 15 * sqrt(lam) + 30, 6000); % An x grid
    pdf = exp(x * log(lam) - lam - gammaln(x + 1));
    pdf(x < 0.5) = 0;                                 % Zero-truncated: a burst of 0 units cannot be detected
    pdf = pdf / trapz(x, pdf);                        % Rescaling to 100% (proper distribution)
    cdf = [0, cumsum((pdf(2:end) + pdf(1:end-1)) / 2 .* diff(x))];
    [cdfUnique, iu] = unique(cdf, 'last');

    for j = 1:numel(K0Vector)
        k0 = K0Vector(j);
        F = interp1(x, cdf, data / k0, 'linear', 1);                         % Model CDF at the data
        KS(i, j) = max(max(abs(F - (1:n) / n)), max(abs(F - (0:n-1) / n)));  % Kolmogorov–Smirnov test score
        quantiles = interp1(cdfUnique, x(iu), min(max(((1:n) - 0.5) / n, cdfUnique(1)), cdfUnique(end)));  
        QQ(i, j) = sqrt(mean((data - k0 * quantiles) .^ 2));                 % Quantile-quanitle RMSE
        m = interp1(x, pdf, centers / k0, 'linear', 0);                      % Model at the histogram bins,
        m = m / max(max(m), 1e-300) * max(binCounts);                        % Scaled to the histogram peak
        MSE(i, j) = mean((m - binCounts) .^ 2);                              % MSE
        CvM(i, j) = 1 / (12 * n) + sum((F - ((1:n) - 0.5) / n) .^ 2);       % Cramer-von Mises score
    end
end

% Ranking Normalization
M = numel(lambda)*numel(K0Vector);
rankKS = tiedrank(KS(:)); rankKS = (rankKS-1)/(M-1); rankKS = reshape(rankKS, size(KS));
rankMSE = tiedrank(MSE(:)); rankMSE = (rankMSE-1)/(M-1); rankMSE = reshape(rankMSE, size(MSE));
rankQQ = tiedrank(QQ(:)); rankQQ = (rankQQ-1)/(M-1); rankQQ = reshape(rankQQ, size(QQ));
WS = (0.25 * rankQQ) + (0.25 * rankKS) + (0.5 * rankMSE); % Weighted scores
[~, best] = min(WS(:));
[lam_id, k0_id] = ind2sub(size(WS), best); % Ties go to the smallest k0 and then smallest lambda
lamWS = lambda(lam_id);
k0WS = K0Vector(k0_id);
scores = struct('MSE', MSE, 'QQ', QQ, 'KS', KS, 'CvM', CvM);
end