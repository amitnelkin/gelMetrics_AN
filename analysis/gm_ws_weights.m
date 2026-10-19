function [lam, k0, comp] = gm_ws_weights(data, W, truncPercent, comp)
% GM_WS_WEIGHTS  Grid/WS fit (modifiedPoissonFit.m, current scoring) for
% several weight sets at once. Reviewer point R2 4a (arbitrary choices): does
% the 0.5 MSE / 0.25 QQ / 0.25 KS weighting matter?
% data: amplitudes after the top trim (as the runner passes them); W: K x 3
% weights [wMSE wQQ wKS]; truncPercent: left-truncation percentile.
% The three component matrices are computed once (lambda = 1:5, k0 = 1:500,
% 10 bins) and re-ranked per weight set exactly as modifiedPoissonFit.m does
% (tiedrank, scaled to [0,1], weighted sum). comp can be passed back in to
% skip the fit. Returns lam(k), k0(k) of the minimum WS for each weight set.
persistent noplotDir
lambda = 1:5; K0 = 1:500;
if nargin < 4 || isempty(comp)
    if isempty(noplotDir) || ~exist(noplotDir, 'dir'), noplotDir = gm_noplot(); end
    MSE = []; QQ = []; KS = [];
    evalc('[~, MSE, QQ, ~, KS] = modifiedPoissonFit(data(:)'', lambda, K0, 10, ''w'', truncPercent);');
    comp = struct('MSE', MSE, 'QQ', QQ, 'KS', KS);
end
M = numel(comp.MSE);
rk = @(X) reshape((tiedrank(X(:)) - 1) / (M - 1), size(X));
rM = rk(comp.MSE); rQ = rk(comp.QQ); rK = rk(comp.KS);
lam = zeros(size(W,1),1); k0 = lam;
for k = 1:size(W,1)
    S = W(k,1)*rM + W(k,2)*rQ + W(k,3)*rK;
    [~, i] = min(S(:));
    [li, ki] = ind2sub(size(S), i);
    lam(k) = lambda(li); k0(k) = K0(ki);
end
end
