function [lamOrder, kBest, scoreSorted, w_score] = gm_score(data, lambda, K0Vector, numBins, percent)
% GM_SCORE  Scorer switch for the checks. Set the environment variable
% GM_SCORER before running a check script (the scripts call clear, so a
% workspace variable would not survive):
%   GM_SCORER = 'orig'    (default) gm_mpfit.m, the ORIGINAL scoring
%   GM_SCORER = 'patched' gm_mpfit_patched.m, the repo's patched
%                         modifiedPoissonFit.m
switch lower(gm_scorer())
    case 'patched'
        [lamOrder, kBest, scoreSorted, w_score] = gm_mpfit_patched(data, lambda, K0Vector, numBins, percent);
    otherwise
        [lamOrder, kBest, scoreSorted, w_score] = gm_mpfit(data, lambda, K0Vector, numBins, percent);
end
end
