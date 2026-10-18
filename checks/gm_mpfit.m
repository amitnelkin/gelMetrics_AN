function [lamOrder, kBest, scoreSorted, w_score] = gm_mpfit(data, lambda, K0Vector, numBins, percent)
% GM_MPFIT  Same scoring math as modifiedPoissonFit.m, with three changes
% so it can be used inside parameter sweeps:
%   1. the truncation percentile is an argument (it is hardcoded to 15 in
%      the original function),
%   2. nothing is plotted,
%   3. the lambdas are returned ranked by their best weighted score.
%
% Outputs
%   lamOrder    - lambda values sorted from best to worst score
%   kBest       - best k0 for each lambda, in the same order
%   scoreSorted - best weighted score for each lambda, in the same order
%   w_score     - full weighted-score matrix (lambda x k0)

if nargin < 4 || isempty(numBins), numBins = 10; end
if nargin < 5 || isempty(percent), percent = 15; end
lambda = lambda(:)';
K0Vector = K0Vector(:)';
data = data(:)';

[binCounts, binEdges] = histcounts(data, numBins, 'Normalization', 'pdf');
if binEdges(1) < 0
    binEdges = fliplr(binEdges);
    binCounts = fliplr(binCounts);
end
nn_full = binCounts;
hist_x_axis = abs(mean([binEdges(1:end-1); binEdges(2:end)]));

% As in the original: percentile of the BIN CENTRES, not of the data.
% percent = 0 means no truncation (the original comment says to set both
% GOF_threshold and percent to 0 for that).
if percent > 0
    GOF_threshold = prctile(hist_x_axis, percent);
else
    GOF_threshold = 0;
end

MSE      = zeros(numel(lambda), numel(K0Vector));
QQ_error = zeros(numel(lambda), numel(K0Vector));
KS_dist  = zeros(numel(lambda), numel(K0Vector));

w_MSE = 0.5; w_QQ = 0.25; w_KS = 0.25;

for j = 1:numel(lambda)
    for i = 1:numel(K0Vector)
        normalizer = K0Vector(i);
        xx = 0:max(hist_x_axis)/normalizer+1;
        yy = poisspdf(xx, lambda(j));
        yq = interp1(xx, yy./max(yy), linspace(0, max(hist_x_axis)/normalizer, length(binCounts)));

        sel = hist_x_axis > GOF_threshold;
        yq_thres = yq(sel);
        binCounts_thres = binCounts(sel);
        nn_full_thres = nn_full(sel);

        if all(~yq_thres) || all(~nn_full_thres) || all(~binCounts_thres)
            yq_thres = yq;
            binCounts_thres = binCounts;
            nn_full_thres = nn_full;
        end

        % goodnessOfFit(x, xref, 'MSE') = mean((xref - x).^2); written
        % out here so the System Identification Toolbox is not needed
        ok = ~isnan(yq_thres);
        xref = binCounts_thres(ok) ./ max(binCounts_thres);
        MSE(j,i) = mean((xref - yq_thres(ok)).^2);

        n = numel(nn_full_thres);
        p = ((1:n) - 0.5) / n;
        q1 = quantile(nn_full_thres, p);
        q2 = quantile(yq_thres*max(binCounts_thres)./max(yq_thres), p);
        QQ_error(j,i) = sqrt(mean((q2 - q1).^2));
        [~, ~, KS_dist(j,i)] = kstest2(yq_thres*max(binCounts_thres)./max(yq_thres), nn_full_thres, 0.01);
    end
end

M = numel(MSE);
rankKS  = reshape((tiedrank(KS_dist(:))  - 1) / (M - 1), size(KS_dist));
rankMSE = reshape((tiedrank(MSE(:))      - 1) / (M - 1), size(MSE));
rankQQ  = reshape((tiedrank(QQ_error(:)) - 1) / (M - 1), size(QQ_error));
w_score = w_QQ*rankQQ + w_KS*rankKS + w_MSE*rankMSE;

[bestPerLam, kIdx] = min(w_score, [], 2);
[scoreSorted, ord] = sort(bestPerLam);
lamOrder = lambda(ord);
kBest = K0Vector(kIdx(ord));
end
