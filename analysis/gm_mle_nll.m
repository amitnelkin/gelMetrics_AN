function nll = gm_mle_nll(lambda, k0, A, sigma, Amin, Amax, nz)
% GM_MLE_NLL  Negative log-likelihood of amplitudes A under the zero-truncated
% Poisson mixture of Gaussians, truncated to [Amin, Amax] (gelmetricsMLE.m).
% Reviewer points R2 4a/4b. Computed in log space (logsumexp over n).
% Noise of a burst of n units (optional nz, default: constant sigma):
%   sigma_n^2 = sigma^2 + n*nz.s1^2 + (nz.f*n*k0)^2
%   s1: extra scatter per unit (growN model); f: scatter proportional to the
%   burst size (prop / calib models).
% False-burst component (optional nz.w0, nz.h0): a fraction w0 of the
% detected bursts are noise ramps or step fragments with |A| ~ half-normal of
% scale h0, not n*k0 events:  p(A) = w0*h(A) + (1-w0)*mixture(A).
if nargin < 7 || isempty(nz), nz = struct('s1', 0, 'f', 0); end
if ~isfield(nz, 's1'), nz.s1 = 0; end
if ~isfield(nz, 'f'), nz.f = 0; end
if ~isfield(nz, 'w0'), nz.w0 = 0; end
if ~isfield(nz, 'h0'), nz.h0 = 1; end
A = A(:);
if ~(lambda > 0 && k0 > 0 && sigma > 0 && nz.s1 >= 0 && nz.f >= 0), nll = Inf; return, end
amax = max(A); if isfinite(Amax), amax = max(amax, Amax); end
N = max(ceil((amax + 8*sigma) / k0), ceil(lambda + 10*sqrt(lambda) + 10));
N = min(N, 5000);
n = 1:N;
% zero-truncated Poisson weights, log space
logw = n*log(lambda) - lambda - gammaln(n+1) - log(-expm1(-lambda));
mu = n * k0;
sg = sqrt(sigma^2 + n*nz.s1^2 + (nz.f*n*k0).^2);          % 1 x N
% log density per amplitude: logsumexp_n [logw + log N(A; mu, sg)]
L = logw - 0.5*((A - mu)./sg).^2 - log(sg) - 0.5*log(2*pi);   % numel(A) x N
m = max(L, [], 2);
logp = m + log(sum(exp(L - m), 2));
% probability mass inside the truncation window
% Phi(zb) - Phi(za) without cancellation: when the window lies above the
% component mean (za > 0) use upper tails, Q(za) - Q(zb).
za = (Amin - mu)./sg; zb = (Amax - mu)./sg;
up = za > 0;
inWin = zeros(size(mu));
inWin(up)  = 0.5*(erfc(za(up)/sqrt(2))  - erfc(zb(up)/sqrt(2)));
inWin(~up) = 0.5*(erfc(-zb(~up)/sqrt(2)) - erfc(-za(~up)/sqrt(2)));
Z = sum(exp(logw) .* inWin);
if nz.w0 > 0
    if ~(nz.w0 < 1 && nz.h0 > 0), nll = Inf; return, end
    logh = log(2) - 0.5*(A/nz.h0).^2 - log(nz.h0) - 0.5*log(2*pi);        % half-normal density
    a1 = log(nz.w0) + logh; a2 = log1p(-nz.w0) + logp;
    mm = max(a1, a2); logp = mm + log(exp(a1 - mm) + exp(a2 - mm));
    hWin = erfc(Amin/(nz.h0*sqrt(2))) - erfc(min(Amax, 1e300)/(nz.h0*sqrt(2)));   % mass of the half-normal in the window
    Z = nz.w0 * hWin + (1 - nz.w0) * Z;
end
if ~(Z > 0), nll = Inf; return, end
nll = -sum(logp) + numel(A)*log(Z);
if ~isfinite(nll), nll = Inf; end
end
