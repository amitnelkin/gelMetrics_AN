function [a_grid, F_model] = modelCDF(lambda_val, k0, threshold, a_max)
% MODELCDF  CDF of the fitted model in intensity units, on (threshold, a_max].
% Same piecewise-linear interpolated Poisson shape as the curve used in
% modifiedPoissonFit.m.
a_grid = linspace(threshold, a_max, 2000);
xx = 0:ceil(a_max/k0)+1;
yy = poisspdf(xx, lambda_val);
% poisspdf underflows to 0 far in the tail (small k0), which made F_model 0/0.
% Same shape computed in log space, rescaled; the CDF is normalised below anyway.
lp = xx.*log(lambda_val) - lambda_val - gammaln(xx+1);
inRange = xx >= floor(threshold/k0) & xx <= ceil(a_max/k0)+1;
yy = exp(lp - max(lp(inRange)));
dens = interp1(xx, yy, a_grid./k0);
F_model = cumtrapz(a_grid, dens);
F_model = F_model ./ F_model(end);
[F_model, iu] = unique(F_model);          % interp1 needs strictly increasing values
a_grid = a_grid(iu);
end
