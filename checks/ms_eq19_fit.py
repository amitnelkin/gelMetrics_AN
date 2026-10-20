"""MS_EQ19_FIT  Fit the manuscript's Eq. (1.9) to simulated burst amplitudes with
several scores, and check recovery of lambda, k0 and lambda*k0 (RESULTS.md
section 30). Data: checks/ms_sim_test.m (MS section III simulations, k0 = 60,
negative bursts, nothing removed) and, for comparison, the true jump sizes of
the noise-free twin (what Eq. 1.9 is meant to describe).
Eq. (1.9) exactly as written: p(I) = (1/k0) lambda^(I/k0) exp(-lambda) / (I/k0)!,
with (I/k0)! = Gamma(I/k0 + 1); normalised numerically over I >= 0 where a
normalised density is needed. No measurement-noise term (the baseline model).
Scores, each minimised over the same grid (lambda 0.05-20 log-spaced plus
1..5 exactly; k0 2-600 log-spaced):
  MSE        10-bin histogram density (data range, as the manuscript) vs model density
  NLL        negative log-likelihood (maximum likelihood)
  KS         max |ECDF - model CDF|
  CvM        Cramer-von Mises: sum (F(x_i) - (2i-1)/2n)^2 (whole-CDF distance)
  QQ         RMSE between sorted data and model quantiles
  WS         the manuscript's weighted rank score 0.5 rank(MSE*) + 0.25 rank(QQ)
             + 0.25 rank(KS), MSE* on the peak-normalised curve as in the code;
             on the full grid, and restricted to integer lambda 1..5 (WS-int, as the MS)
Goodness of fit of the best-NLL curve: KS test p value.
Usage: python checks/ms_eq19_fit.py
Writes analysis_outputs/ms_sim_test/{fits.csv, recovery.png/.pdf, fits_examples.png/.pdf, traces.png/.pdf}."""
import os
import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import kstest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'ms_sim_test')
A = pd.read_csv(os.path.join(d, 'amplitudes.csv')); J = pd.read_csv(os.path.join(d, 'true_jumps.csv'))
LAMS = np.unique(np.concatenate([np.exp(np.linspace(np.log(0.05), np.log(20), 80)), [1, 2, 3, 4, 5]]))
K0S = np.exp(np.linspace(np.log(2), np.log(600), 140))
SC = ['MSE', 'NLL', 'KS', 'CvM', 'QQ', 'WS', 'WS-int']
COL = dict(zip(SC, ['tab:orange', 'tab:blue', 'tab:green', 'tab:purple', 'tab:brown', 'tab:red', 'k']))
XG = {}   # per lambda: x grid, normalised pdf, cdf (x = I/k0)
for lam in LAMS:
    xg = np.linspace(0, lam + 15 * np.sqrt(lam) + 30, 6000)
    pdf = np.exp(xg * np.log(lam) - lam - gammaln(xg + 1)); Z = np.trapezoid(pdf, xg)
    cdf = np.concatenate([[0], np.cumsum((pdf[1:] + pdf[:-1]) / 2 * np.diff(xg))]) / Z
    XG[lam] = (xg, pdf / Z, cdf, pdf.max() / Z)


def score_grid(a):
    a = np.sort(a); n = len(a); pq = (np.arange(1, n + 1) - 0.5) / n
    cnt, edges = np.histogram(a, bins=10); ctr = (edges[:-1] + edges[1:]) / 2; dens = cnt / (n * np.diff(edges))
    S = {s: np.zeros((len(LAMS), len(K0S))) for s in ['MSE', 'MSEpk', 'NLL', 'KS', 'CvM', 'QQ']}
    ecdf_hi = np.arange(1, n + 1) / n; ecdf_lo = np.arange(0, n) / n
    for i, lam in enumerate(LAMS):
        xg, pdf, cdf, pk = XG[lam]
        X = a[None, :] / K0S[:, None]                                   # k0 x n
        F = np.interp(X, xg, cdf, right=1.0)
        f = np.interp(X, xg, pdf, right=0.0) / K0S[:, None]
        S['NLL'][i] = -np.log(f + 1e-300).sum(1)
        S['KS'][i] = np.maximum(np.abs(F - ecdf_hi).max(1), np.abs(F - ecdf_lo).max(1))
        S['CvM'][i] = 1 / (12 * n) + ((F - (2 * np.arange(1, n + 1) - 1) / (2 * n)) ** 2).sum(1)
        qx = np.interp(pq, cdf, xg); S['QQ'][i] = np.sqrt(((a[None, :] - K0S[:, None] * qx[None, :]) ** 2).mean(1))
        m = np.interp(ctr[None, :] / K0S[:, None], xg, pdf, right=0.0) / K0S[:, None]
        S['MSE'][i] = ((m - dens[None, :]) ** 2).mean(1)
        mp = m / m.max(1, keepdims=True).clip(1e-300) * cnt.max(); S['MSEpk'][i] = ((mp - cnt[None, :]) ** 2).mean(1)
    rk = lambda M: (np.argsort(np.argsort(M.ravel())).reshape(M.shape) + 1) / M.size
    S['WS'] = 0.5 * rk(S['MSEpk']) + 0.25 * rk(S['QQ']) + 0.25 * rk(S['KS'])
    integ = np.isin(LAMS, [1, 2, 3, 4, 5])
    Wi = np.full_like(S['WS'], np.inf); sub = S['MSEpk'][integ], S['QQ'][integ], S['KS'][integ]
    Wi[integ] = 0.5 * rk(sub[0]) + 0.25 * rk(sub[1]) + 0.25 * rk(sub[2]); S['WS-int'] = Wi
    return S


rows, keep = [], {}
for src, D, col in [('detected', A, 'amp'), ('true jumps', J, 'jump')]:
    for (i_, lam, nz, rep), g in D.groupby(['id', 'lambda', 'noise', 'rep']):
        a = g[col].values.astype(float)
        if src == 'true jumps' and len(a) > 3000:
            a = np.random.default_rng(int(i_)).choice(a, 3000, replace=False)
        S = score_grid(a)
        for s in SC:
            i, j = np.unravel_index(np.argmin(S[s]), S[s].shape)
            r = dict(source=src, lambdaTrue=lam, noise=nz, rep=rep, n=len(a), score=s, lam=LAMS[i], k0=K0S[j], lamk0=LAMS[i] * K0S[j],
                     at_edge=(i in (0, len(LAMS) - 1)) or (j in (0, len(K0S) - 1)))
            if s == 'NLL':
                xg, pdf, cdf, _ = XG[LAMS[i]]; r['KS_p'] = kstest(a, lambda q: np.interp(q / K0S[j], xg, cdf, right=1.0)).pvalue
            rows.append(r)
        if rep == 1:
            keep[(src, lam, nz)] = (a, {s: np.unravel_index(np.argmin(S[s]), S[s].shape) for s in SC})
    print(src, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(d, 'fits.csv'), index=False)
pd.set_option('display.width', 250)
summ = R.groupby(['source', 'noise', 'lambdaTrue', 'score']).agg(lam=('lam', 'median'), k0=('k0', 'median'), lamk0=('lamk0', 'median'),
                                                                   edge=('at_edge', 'mean'), KS_p=('KS_p', 'median'), n=('n', 'median')).reset_index()
summ.to_csv(os.path.join(d, 'fits_summary.csv'), index=False)
print(summ.round(2).to_string(index=False))

# ---------- figure 1: example traces
E = pd.read_csv(os.path.join(d, 'example_traces.csv')); EB = pd.read_csv(os.path.join(d, 'example_bursts.csv'))
fig, axs = plt.subplots(3, 2, figsize=(14, 9), sharex=True)
for k, tr in enumerate(sorted(E.trace.unique())):
    e = E[E.trace == tr]; b = EB[EB.trace == tr]; tm = (e.frame - 1) / 6
    axs[k, 0].plot(tm, e.raw, color='0.6', lw=0.6, label='simulated signal (noise 40)'); axs[k, 0].plot(tm, e.noiseFree, color='k', lw=1.2, label='same signal without noise')
    axs[k, 1].plot(tm, e.processed, color='0.2', lw=1)
    for _, s in b.iterrows():
        axs[k, 1].plot([(s.firstFrame - 1) / 6, (s.lastFrame - 1) / 6], [s.levelStart, s.levelEnd], color='tab:red' if s.sign < 0 else 'tab:blue', lw=3)
    axs[k, 1].text(0.01, 0.92, f"{(b.sign < 0).sum()} negative (red), {(b.sign > 0).sum()} positive (blue) bursts detected", transform=axs[k, 1].transAxes, fontsize=8)
axs[0, 0].legend(fontsize=8); axs[0, 0].set_title('simulated signal, λ = 2, k0 = 60 (MS section III)', loc='left', fontsize=10)
axs[0, 1].set_title('what the detector sees (smoothed, bleaching-corrected) and the bursts it calls', loc='left', fontsize=10)
for ax in axs[-1]:
    ax.set_xlabel('time (min)')
fig.tight_layout(); fig.savefig(os.path.join(d, 'traces.png'), dpi=130); fig.savefig(os.path.join(d, 'traces.pdf'))

# ---------- figure 2: fits for lambda 1, 2, 4 (noise 40) and the true jumps
fig, axs = plt.subplots(2, 3, figsize=(18, 9))
for c, lam in enumerate([1, 2, 4]):
    for r_, src in enumerate(['detected', 'true jumps']):
        ax = axs[r_, c]; a, best = keep[(src, lam, 40)]
        hi = np.percentile(a, 99.5); bins = np.linspace(0, hi, 40)
        ax.hist(a, bins=bins, density=True, color='0.8', edgecolor='0.5', label=f'{src} (n = {len(a)})')
        xx = np.linspace(0.5, hi, 600)
        for s in SC:
            i, j = best[s]; xg, pdf, _, _ = XG[LAMS[i]]
            ax.plot(xx, np.interp(xx / K0S[j], xg, pdf) / K0S[j], color=COL[s], lw=1.8, ls='--' if s == 'WS-int' else '-',
                    label=f'{s}: λ {LAMS[i]:.2g}, k0 {K0S[j]:.0f}')
        ax.set_title(f"{'detected bursts' if src == 'detected' else 'true jumps (no noise, no detection)'}: true λ = {lam}, k0 = 60", fontsize=10, loc='left')
        ax.set_xlabel('burst amplitude (a.u.)'); ax.legend(fontsize=7)
axs[0, 0].set_ylabel('probability density'); axs[1, 0].set_ylabel('probability density')
fig.tight_layout(); fig.savefig(os.path.join(d, 'fits_examples.png'), dpi=130); fig.savefig(os.path.join(d, 'fits_examples.pdf'))

# ---------- figure 3: recovery
fig, axs = plt.subplots(3, 3, figsize=(18, 13))
for c, (src, nz, lab) in enumerate([('detected', 40, 'detected bursts, noise 40 (MS)'), ('detected', 4, 'detected bursts, noise 4'), ('true jumps', 40, 'true jumps (no noise, no detection)')]):
    g = R[(R.source == src) & (R.noise == nz)]
    for r_, (q, truth, ylab) in enumerate([('lam', lambda l: l, 'fitted λ'), ('k0', lambda l: 60 + 0 * l, 'fitted k0'), ('lamk0', lambda l: 60 * l, 'fitted λ·k0')]):
        ax = axs[r_, c]
        for k, s in enumerate(SC):
            m = g[g.score == s].groupby('lambdaTrue')[q].agg(['median', 'min', 'max'])
            xs = m.index.values * (1 + (k - 3) * 0.02)
            ax.errorbar(xs, m['median'], yerr=[m['median'] - m['min'], m['max'] - m['median']], fmt='o-', color=COL[s], ms=4, capsize=2, lw=1, label=s,
                        ls='--' if s == 'WS-int' else '-')
        lt = np.array([0.5, 1, 2, 3, 4, 5]); ax.plot(lt, truth(lt), 'k:', lw=2, label='truth')
        if q != 'lamk0':
            ax.set_yscale('log')
        ax.set_xlabel('true λ'); ax.set_ylabel(ylab)
        if r_ == 0:
            ax.set_title(lab, fontsize=11, loc='left')
axs[0, 0].legend(fontsize=8, ncol=2)
fig.suptitle('Eq. (1.9) fitted to simulated data (k0 = 60) with six scores: median and range over 5 independent datasets', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(d, 'recovery.png'), dpi=130); fig.savefig(os.path.join(d, 'recovery.pdf'))
