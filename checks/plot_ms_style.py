"""PLOT_MS_STYLE  Why the fits looked better in the manuscript (R2 4a/4b).
Same amplitudes as checks/mle_real_data.m (negative bursts), drawn the way
modifiedPoissonFit.m draws the manuscript figures (Figs. 6-9, panels c-f):
  - 10-bin 'pdf' histogram of |A| after the top-10% trim, bin centres as points;
  - red dashed = grid-score fit: zero-truncated Poisson pmf at n = 0,1,2,..,
    linearly interpolated at (bin centre / k0), divided by its max and
    multiplied by the tallest bin (only the SHAPE is compared, at 10 points);
  - blue dashed = the MLE (lambda, k0) drawn in exactly the same style;
  - orange solid = the MLE's actual fitted density: Gaussian mixture with sigma
    from ampVar, renormalised on the fit window and scaled by the fraction of
    amplitudes in the window (true normalisation, no rescaling to the bars).
Usage (repo root): python checks/plot_ms_style.py
Writes analysis_outputs/mle_real/ms_style_Neg.png/.pdf."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import gammaln

d = os.path.join('analysis_outputs', 'mle_real')
A = pd.read_csv(os.path.join(d, 'amplitudes.csv'))
R = pd.read_csv(os.path.join(d, 'mle_real.csv'))
BAR, PTS, GRID, MLEMS, MLED, INK, MUTED = '#f3c2b0', '#2a78d6', '#d62a2a', '#2a78d6', '#eb6834', '#1f1e1c', '#898781'


def ms_curve(c, lam, k0, peak):
    """modifiedPoissonFit.m plotting: ZTP pmf interpolated at c/k0, scaled to the tallest bin"""
    xx = np.arange(0, np.floor(c.max() / k0 + 1) + 1)
    lp = xx * np.log(lam) - lam - gammaln(xx + 1)
    yy = np.exp(lp - lp.max()); yy[0] = 0
    yq = np.interp(c / k0, xx, yy / yy.max(), right=np.nan)
    return yq * peak / np.nanmax(yq)


def mle_density(x, lam, k0, sig, lo, hi):
    N = int(min(max(np.ceil((hi + 8 * sig) / k0), lam + 10 * np.sqrt(lam) + 10), 5000))
    n = np.arange(1, N + 1)
    w = np.exp(n * np.log(lam) - lam - gammaln(n + 1) - np.log(-np.expm1(-lam)))
    dens = (w[None, :] * stats.norm.pdf(x[:, None], n * k0, sig)).sum(1)
    Z = np.sum(w * (stats.norm.cdf(hi, n * k0, sig) - stats.norm.cdf(lo, n * k0, sig)))
    return dens / Z


fig, axs = plt.subplots(3, 3, figsize=(12, 10))
for ax, (_, r) in zip(axs.ravel(), R[R.sign == 'Neg'].iterrows()):
    a = A[(A.dataset == r.dataset) & (A.sign == 'Neg')].amp.values
    amax = np.percentile(a, 90, method='hazen')
    a = a[a <= amax]                                     # rmoutliers(...,'percentiles',[0 90])
    cnt, edges = np.histogram(a, bins=10, density=True)
    c = 0.5 * (edges[:-1] + edges[1:])
    ax.bar(c, cnt, width=np.diff(edges), color=BAR, edgecolor='white')
    ax.plot(c, cnt, 'o', mfc='none', mec=PTS, ms=6, label='data (bin centres)')
    ax.plot(c, ms_curve(c, r.gridLam, r.gridK0, cnt.max()), '--', color=GRID, lw=1.6,
            label='grid fit, MS style')
    ax.plot(c, ms_curve(c, r['lambda'], r.k0, cnt.max()), '--', color=MLEMS, lw=1.2,
            label='MLE params, MS style')
    amin = np.percentile(a, r.truncPct, method='hazen')
    xs = np.linspace(amin, amax, 500)
    frac = np.mean(a >= amin)
    ax.plot(xs, frac * mle_density(xs, r['lambda'], r.k0, r.sigma, amin, amax), color=MLED, lw=1.2,
            label='MLE fitted density')
    ax.set_title(f'{r.dataset} (n={len(a)})\ngrid {int(r.gridLam)}/{r.gridK0:.0f}   MLE {r["lambda"]:.2f}/{r.k0:.0f}  '
                 f'sigma {r.sigma:.0f}', fontsize=8.5, loc='left', color=INK)
    ax.tick_params(labelsize=7, colors=MUTED)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
for ax in axs[-1]:
    ax.set_xlabel('Intensity (a.u.)', fontsize=8)
for ax in axs[:, 0]:
    ax.set_ylabel('Probability', fontsize=8)
axs[0, 0].legend(fontsize=7, frameon=False)
fig.suptitle('Negative bursts drawn as in the manuscript (10 bins, curve rescaled to the tallest bin) vs the MLE density',
             fontsize=10, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'ms_style_Neg.{ext}'), dpi=150)
print('saved')
