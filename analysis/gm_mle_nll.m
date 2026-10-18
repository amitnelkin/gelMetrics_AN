function nll = gm_mle_nll(lambda, k0, A, sigma, Amin, Amax)
% GM_MLE_NLL  Negative log-likelihood of amplitudes A under the zero-truncated
% Poisson mixture of Gaussians, truncated to [Amin, Amax] (gelmetricsMLE.m).
% Reviewer points R2 4a/4b. Computed in log space (logsumexp over n).
A = A(:);
if ~(lambda > 0 && k0 > 0 && sigma > 0), nll = Inf; return, end
amax = max(A); if isfinite(Amax), amax = max(amax, Amax); end
N = max(ceil((amax + 8*sigma) / k0), ceil(lambda + 10*sqrt(lambda) + 10));
N = min(N, 5000);
n = 1:N;
% zero-truncated Poisson weights, log space
logw = n*log(lambda) - lambda - gammaln(n+1) - log(-expm1(-lambda));
mu = n * k0;
% log density per amplitude: logsumexp_n [logw + log N(A; mu, sigma)]
L = logw - 0.5*((A - mu)/sigma).^2 - log(sigma) - 0.5*log(2*pi);   % numel(A) x N
m = max(L, [], 2);
logp = m + log(sum(exp(L - m), 2));
% probability mass inside the truncation window
% Phi(zb) - Phi(za) without cancellation: when the window lies above the
% component mean (za > 0) use upper tails, Q(za) - Q(zb).
za = (Amin - mu)/sigma; zb = (Amax - mu)/sigma;
up = za > 0;
inWin = zeros(size(mu));
inWin(up)  = 0.5*(erfc(za(up)/sqrt(2))  - erfc(zb(up)/sqrt(2)));
inWin(~up) = 0.5*(erfc(-zb(~up)/sqrt(2)) - erfc(-za(~up)/sqrt(2)));
Z = sum(exp(logw) .* inWin);
if ~(Z > 0), nll = Inf; return, end
nll = -sum(logp) + numel(A)*log(Z);
if ~isfinite(nll), nll = Inf; end
end
