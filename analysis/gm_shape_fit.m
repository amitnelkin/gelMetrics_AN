function [lam, k0, lo, hi] = gm_shape_fit(amps, trace, nBoot)
% GM_SHAPE_FIT  Shape-only fit of lambda and k0 to all burst amplitudes, as in
% checks/lambda_k0_shape.py --notrim (MS Fig. 5(e); RESULTS.md sections 26-28).
% Model: the original modifiedPoissonFit curve, i.e. the Poisson pmf at whole
% numbers of units (P(0) = 0) joined by straight lines and evaluated at A/k0,
% normalised over [0, 1.001 * max(A)]; maximum likelihood (Nelder-Mead).
%   [lam, k0] = gm_shape_fit(amps)
%   [lam, k0, lo, hi] = gm_shape_fit(amps, trace, nBoot)   95% CIs of
%       [lambda, k0, lambda*k0] from resampling traces (granules), nBoot times.
amps = double(amps(:));
[lam, k0] = fit1(amps);
lo = []; hi = [];
if nargin >= 3 && nBoot > 0
    [~, ~, g] = unique(trace(:)); G = accumarray(g, (1:numel(amps))', [], @(v) {v});
    v = zeros(nBoot, 3);
    for b = 1:nBoot
        idx = vertcat(G{randi(numel(G), numel(G), 1)});
        [l, k] = fit1(amps(idx)); v(b, :) = [l, k, l * k];
    end
    q = prctile(v, [2.5 97.5]); lo = q(1, :); hi = q(2, :);
end
end


function [lam, k0] = fit1(a)
hiA = max(a) * 1.001; loK = hiA / 300; hiK = 3 * hiA;
NN = 0:399;
f = @(th) nll(th, a, hiA, loK, hiK, NN);
o = optimset('TolX', 1e-4, 'TolFun', 1e-4, 'MaxIter', 3000, 'MaxFunEvals', 6000, 'Display', 'off');
bestF = inf; bx = [0 0];
for l0 = [0.3 1 3 10]
    for fct = [0.7 1.5]
        [x, fv] = fminsearch(f, [log(l0), log(fct * mean(a) / max(l0, 1))], o);
        if fv < bestF, bestF = fv; bx = x; end
    end
end
lam = exp(bx(1)); k0 = exp(bx(2));
end


function v = nll(th, a, hiA, loK, hiK, NN)
lam = exp(th(1)); k0 = exp(th(2));
if lam < 0.05 || lam > 60 || k0 < loK || k0 > hiK, v = 1e12; return; end
pm = exp(NN * log(lam) - lam - gammaln(NN + 1)); pm(1) = 0;
shape = @(x) interp1(NN, pm, x, 'linear', 0);
g = linspace(0, hiA / k0, 2000);
z = trapz(g, shape(g)) * k0;
v = -sum(log(shape(a / k0) + 1e-300) - log(z));
end
