function [lamOrder, kBest, scoreSorted, w_score] = gm_mpfit_patched(data, lambda, K0Vector, numBins, percent)
% GM_MPFIT_PATCHED  Same outputs as gm_mpfit.m, but scores with the repo's
% current modifiedPoissonFit.m (after the bug fixes: sample-based QQ/KS,
% model on the data's bin centres, percentile of the amplitudes).
% Figures and console output from modifiedPoissonFit are suppressed.
if nargin < 4 || isempty(numBins), numBins = 10; end
if nargin < 5 || isempty(percent), percent = 15; end
assert(nargin(@modifiedPoissonFit) >= 6, ...
    'modifiedPoissonFit.m has no percent argument: the bug-fix commit is not applied');
lambda = lambda(:)';
K0Vector = K0Vector(:)';
before = findall(0, 'Type', 'figure');
dv = get(0, 'DefaultFigureVisible'); set(0, 'DefaultFigureVisible', 'off');
w_score = [];
evalc('[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(data(:)'', lambda, K0Vector, numBins, ''check'', percent);');
set(0, 'DefaultFigureVisible', dv);
close(setdiff(findall(0, 'Type', 'figure'), before));
[bestPerLam, kIdx] = min(w_score, [], 2);
[scoreSorted, ord] = sort(bestPerLam);
lamOrder = lambda(ord);
kBest = K0Vector(kIdx(ord));
end
