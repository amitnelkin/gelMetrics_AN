"""LAMBDA_K0_DEGENERACY  Why lambda and k0 are not reported: one burst-amplitude
distribution is explained equally well by small units exchanged in large
numbers, large units in small numbers, or anything in between. RESULTS.md
section 25.
Model: the likelihood model of analysis/gelmetricsMLE.m in its basic form:
amplitude A = n*k0 + Gaussian noise (sigma), n >= 1 zero-truncated Poisson(lambda),
fitted on the amplitudes between the 15th percentile (detection truncation, as
the manuscript) and the 90th percentile (top-10% trim), normalised to that
window. Data: 4PP7-FUS negative-burst amplitudes, R2, SP 10, window 13
(analysis_outputs/fus_test_4x/amplitudes_4x.csv); sigma = 71 (section 9, from
the quiescent segments).
Panels: (a) data + best fits at fixed lambda 0.3, 1, 2, 4 (k0 refitted);
(b) log-likelihood surface: the ridge; (c) data simulated from three different
(lambda, k0) truths on the ridge: hard to tell apart; (d) 200 datasets
simulated from a known truth (lambda 1), each fitted: the estimates scatter
along the ridge. Seeded.
Usage: python checks/lambda_k0_degeneracy.py
Writes analysis_outputs/lambda_k0_degeneracy/{degeneracy.png/.pdf, summary.csv}."""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm, poisson, ks_2samp, kstest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'lambda_k0_degeneracy'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(26)
A0 = pd.read_csv('analysis_outputs/fus_test_4x/amplitudes_4x.csv')
A0 = A0[(A0.data == 'R2_reextracted') & (A0.dataset == '4PP7-FUS')].amp.abs().values
LO, HI = np.percentile(A0, 15), np.percentile(A0, 90)
A = A0[(A0 >= LO) & (A0 <= HI)]
SIG = 71.0; NMAX = 40; ns = np.arange(1, NMAX + 1)
LAMS = np.exp(np.linspace(np.log(0.05), np.log(8), 60)); K0S = np.linspace(5, 300, 100)


def weights(lam):
    w = poisson.pmf(ns, lam); return w / w.sum()


def loglik_grid(a):
    W = np.array([weights(l) for l in LAMS])                       # lambda x n
    mu = K0S[:, None] * ns[None, :]                                # k0 x n
    phi = np.exp(-0.5 * ((a[None, None, :] - mu[:, :, None]) / SIG) ** 2) / (SIG * np.sqrt(2 * np.pi))   # k0 x n x data
    dens = np.einsum('ln,knd->lkd', W, phi)                        # lambda x k0 x data
    Zc = norm.cdf((HI - mu) / SIG) - norm.cdf((LO - mu) / SIG)     # k0 x n
    Z = W @ Zc.T                                                   # lambda x k0
    return np.log(dens + 1e-300).sum(2) - len(a) * np.log(Z)


def pdf(x, lam, k0):
    w = weights(lam); mu = k0 * ns
    d = (w[:, None] * norm.pdf((x[None, :] - mu[:, None]) / SIG) / SIG).sum(0)
    Z = (w * (norm.cdf((HI - mu) / SIG) - norm.cdf((LO - mu) / SIG))).sum()
    return d / Z


def simulate(lam, k0, n):
    s = []
    while len(s) < n:
        m = rng.poisson(lam, 4 * n); m = m[m >= 1]
        x = m * k0 + SIG * rng.standard_normal(len(m)); s.extend(x[(x >= LO) & (x <= HI)])
    return np.array(s[:n])


L = loglik_grid(A); best = np.unravel_index(L.argmax(), L.shape)
prof = {lam: K0S[L[np.argmin(abs(LAMS - lam))].argmax()] for lam in [0.3, 1, 2, 4]}
xx = np.linspace(LO, HI, 400)
cdf = lambda lam, k0: (lambda q: np.interp(q, xx, np.cumsum(pdf(xx, lam, k0)) / np.sum(pdf(xx, lam, k0))))
rows = []
for lam, k0 in prof.items():
    i = np.argmin(abs(LAMS - lam)); j = np.argmin(abs(K0S - k0))
    rows.append(dict(lam=lam, k0=k0, mean_units_x_k0=lam / (1 - np.exp(-lam)) * k0, dlogL=L[i, j] - L.max(), KS_p=kstest(A, cdf(lam, k0)).pvalue))

fig, axs = plt.subplots(1, 4, figsize=(22, 5.2))
ax = axs[0]
ax.hist(A, bins=14, density=True, color='tab:green', alpha=0.45, edgecolor='k', label=f'4PP7-FUS bursts (n = {len(A)})')
for r, c in zip(rows, ['tab:purple', 'tab:orange', 'tab:blue', 'tab:red']):
    ax.plot(xx, pdf(xx, r['lam'], r['k0']), color=c, lw=2,
            label=f"λ = {r['lam']}, k0 = {r['k0']:.0f}   (ΔlogL {r['dlogL']:.1f}, KS p {r['KS_p']:.2f})")
ax.set_xlabel('burst amplitude (a.u.)'); ax.set_ylabel('probability density'); ax.legend(fontsize=7.5)
ax.set_title('(a) one histogram, four very different explanations', fontsize=10, loc='left')
ax = axs[1]
D = L - L.max()
cs = ax.contourf(K0S, LAMS, np.clip(D, -20, 0), levels=np.linspace(-20, 0, 21), cmap='viridis')
ax.contour(K0S, LAMS, D, levels=[-3.0], colors='w', linewidths=2)
for m in [100, 150, 200]:
    lamc = LAMS; k0c = m * (1 - np.exp(-lamc)) / lamc
    ax.plot(k0c, lamc, 'w:', lw=1)
ax.set_yscale('log'); ax.set_xlabel('k0 (intensity of one unit)'); ax.set_ylabel('λ (units per burst)')
ax.set_title('(b) likelihood: a long ridge (white: 95% region)\ndotted: constant mean burst size', fontsize=10, loc='left'); fig.colorbar(cs, ax=ax, label='log-likelihood vs best')
ax = axs[2]
truths = [(0.3, prof[0.3]), (1, prof[1]), (4, prof[4])]; sims = [simulate(l, k, len(A)) for l, k in truths]
bins = np.linspace(LO, HI, 15)
for (l, k), s, c in zip(truths, sims, ['tab:purple', 'tab:orange', 'tab:red']):
    ax.hist(s, bins=bins, histtype='step', lw=2, color=c, density=True, label=f'simulated truth λ = {l}, k0 = {k:.0f}')
p01 = ks_2samp(sims[0], sims[1]).pvalue; p02 = ks_2samp(sims[0], sims[2]).pvalue; p12 = ks_2samp(sims[1], sims[2]).pvalue
ax.text(0.98, 0.70, f'pairwise KS p:\n{p01:.2f}, {p02:.2f}, {p12:.2f}', transform=ax.transAxes, ha='right', fontsize=8)
ax.set_xlabel('burst amplitude (a.u.)'); ax.legend(fontsize=7.5)
ax.set_title(f'(c) three different truths, same size data (n = {len(A)}):\nhard to tell apart', fontsize=10, loc='left')
ax = axs[3]
lt, kt = 1, prof[1]; est = []
for it in range(200):
    Ls = loglik_grid(simulate(lt, kt, len(A))); i, j = np.unravel_index(Ls.argmax(), Ls.shape); est.append((LAMS[i], K0S[j]))
    if it % 50 == 0: print('repeat', it, flush=True)
est = np.array(est)
ax.plot(est[:, 1], est[:, 0], '.', color='0.4', alpha=0.6, label='fit to one simulated dataset')
ax.plot(kt, lt, 'r*', ms=18, label=f'truth (λ = {lt}, k0 = {kt:.0f})')
ax.set_yscale('log'); ax.set_xlabel('fitted k0'); ax.set_ylabel('fitted λ'); ax.legend(fontsize=8)
within = np.mean((est[:, 0] > 0.5 * lt) & (est[:, 0] < 2 * lt))
ax.set_title(f'(d) 200 repeats of a known truth: fitted λ within 2x of truth\nin {100 * within:.0f}% (λ range {np.percentile(est[:, 0], 5):.2f}-{np.percentile(est[:, 0], 95):.1f}, 90%)', fontsize=10, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(out, 'degeneracy.png'), dpi=130); fig.savefig(os.path.join(out, 'degeneracy.pdf'))
S = pd.DataFrame(rows); S['sim_KS_pairwise'] = f'{p01:.3f},{p02:.3f},{p12:.3f}'; S['recovery_within2x'] = within
S['lam90lo'], S['lam90hi'] = np.percentile(est[:, 0], [5, 95]); S['k090lo'], S['k090hi'] = np.percentile(est[:, 1], [5, 95])
S.to_csv(os.path.join(out, 'summary.csv'), index=False); print(S.round(3).to_string(index=False))
print('best overall: lambda %.2f, k0 %.0f' % (LAMS[best[0]], K0S[best[1]]))
