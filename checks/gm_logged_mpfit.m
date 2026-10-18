function varargout = gm_logged_mpfit(data, lambda, K0Vector, numBins, title_str, varargin)
% GM_LOGGED_MPFIT  Calls the repo's modifiedPoissonFit.m unchanged and
% records the ranked (lambda, k0) result in the global GM_FITLOG, so the
% runners' results can be tabulated (SimulationRunWithAnalysis.m only shows
% them as text on figures). Used by run_experimental.m / run_simulation.m.
global GM_FITLOG
out = cell(1, 6);
[out{:}] = modifiedPoissonFit(data, lambda, K0Vector, numBins, title_str, varargin{:});
w_score = out{6};
[bestPerLam, kIdx] = min(w_score, [], 2);
[scoreSorted, ord] = sort(bestPerLam);
entry.title    = title_str;
entry.nFit     = numel(data);
entry.lamOrder = lambda(ord);
entry.k0Order  = K0Vector(kIdx(ord));
entry.score    = scoreSorted';
entry.K0min    = K0Vector(1);
entry.K0max    = K0Vector(end);
entry.data     = data(:)';
entry.lambda   = lambda;
entry.K0Vector = K0Vector;
entry.w_score  = w_score;
entry.extraArgs = varargin;
GM_FITLOG{end+1} = entry;
varargout = out(1:max(nargout, 0));
end
