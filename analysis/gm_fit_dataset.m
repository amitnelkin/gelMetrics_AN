function out = gm_fit_dataset(T, opts)
% GM_FIT_DATASET  Fits one dataset (struct array from gm_amplitudes_by_trace)
% with the likelihood fitter and, for comparison, the grid-score method.
% Reviewer points R2 4a/4b (CIs, continuous lambda) and R2 5 (granule-level
% bootstrap: traces are resampled, never bursts).
% opts: .sign 'Neg'|'Pos' (default 'Neg'), .truncPercent (15), .topTrim (10),
%       .nBoot (200), .seed (1), .surface (false), .grid (true: also run the
%       grid-score fit, gm_mpfit_patched = current modifiedPoissonFit.m).
% Output: point estimate + profile CIs (out.mle), bootstrap percentile CIs
% (out.boot), grid result (out.grid), counts.
if nargin < 2, opts = struct(); end
d = struct('sign','Neg','truncPercent',15,'topTrim',10,'nBoot',200,'seed',1,'surface',false,'grid',true,'fitSigma',false,'lamPrior',[],'noise','fixed','calib',[1 0],'falseBursts','none','floorK',0,'h0K',0);
fn = fieldnames(d);
for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = d.(fn{i}); end, end
amp = @(TT) abs(vertcat(TT.(['amp' opts.sign])));
keep = [T.len] >= 180;
T = T(keep);
A = amp(T); sigma = gm_sigma_from_ampvar(vertcat(T.ampVar));
mo = struct('truncPercent', opts.truncPercent, 'topTrim', opts.topTrim, 'profile', true, ...
    'surface', opts.surface, 'fitSigma', opts.fitSigma, 'lamPrior', opts.lamPrior, ...
    'noise', opts.noise, 'calib', opts.calib, 'falseBursts', opts.falseBursts, 'floorK', opts.floorK, 'h0K', opts.h0K);
out.mle = gelmetricsMLE(A, sigma, mo);
out.nTraces = numel(T); out.nBursts = numel(A); out.sigma = sigma;

% granule bootstrap: resample traces with replacement; sigma re-estimated
nb = opts.nBoot; bl = nan(nb,1); bk = nan(nb,1);
if nb > 0
    rng(opts.seed);
    idx = randi(numel(T), numel(T), nb);
    mb = mo; mb.profile = false; mb.surface = false;
    mb.scanN = 12; mb.start = [out.mle.lambda out.mle.k0];   % warm start: speed only
    for b = 1:nb
        Tb = T(idx(:,b));
        Rb = gelmetricsMLE(amp(Tb), gm_sigma_from_ampvar(vertcat(Tb.ampVar)), mb);
        bl(b) = Rb.lambda; bk(b) = Rb.k0;
    end
end
out.boot = struct('lambda', bl, 'k0', bk, 'lamCI', prctile(bl, [2.5 97.5]), ...
    'k0CI', prctile(bk, [2.5 97.5]), 'lamK0CI', prctile(bl.*bk, [2.5 97.5]));

out.grid = struct('lambda', NaN, 'k0', NaN);
if opts.grid
    a = rmoutliers(A, 'percentiles', [0 100 - opts.topTrim]);
    [lo, kb] = gm_mpfit_patched(a, 1:5, 1:500, 10, opts.truncPercent);
    out.grid = struct('lambda', lo(1), 'k0', kb(1));
end
end
