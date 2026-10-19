"""PLOT_MLE_HISTOGRAMS  Burst-amplitude histograms with the fitted likelihood
model overlaid (R2 4a/4b: does the model describe the data?).
Inputs: analysis_outputs/mle_real/amplitudes.csv (checks/export_amplitudes.m)
and mle_real.csv (checks/mle_real_data.m: lambda, k0, sigma per dataset).
The fitter's cuts are redone exactly (top 10% removed, then left-truncated at
the dataset's percentile; MATLAB prctile = numpy 'hazen'); the fitted n is
checked against mle_real.csv. Model curve = zero-truncated Poisson mixture of
Gaussians N(n*k0, sigma^2), renormalised to the fit window, scaled to counts;
thin dashed curves = the n = 1, 2, 3, ... components. Grey bars = amplitudes
outside the fit window. KS test of the windowed data against the truncated
model CDF (parameters estimated from the same data, so p is optimistic).
Usage (repo root): python checks/plot_mle_histograms.py
Writes analysis_outputs/mle_real/mle_hist_<Neg|Pos>.png/.pdf and mle_hist_gof.csv."""
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
BAR, OUT, MODEL, COMP, INK, MUTED = '#9fc2ec', '#d9d7d2', '#eb6834', '#eb6834', '#1f1e1c', '#898781'


def mixture(x, lam, k0, sig):
    """density of the untruncated model and its components (n = 1..N)"""
    N = int(max(np.ceil((x.max() + 8 * sig) / k0), lam + 10 * np.sqrt(lam) + 10))
    N = min(N, 5000)
    n = np.arange(1, N + 1)
    logw = n * np.log(lam) - lam - gammaln(n + 1) - np.log(-np.expm1(-lam))
    w = np.exp(logw)
    comp = w[None, :] * stats.norm.pdf(x[:, None], n[None, :] * k0, sig)
    return comp.sum(1), comp, w, n


def window_mass(lam, k0, sig, lo, hi):
    _, _, w, n = mixture(np.array([hi]), lam, k0, sig)
    return np.sum(w * (stats.norm.cdf(hi, n * k0, sig) - stats.norm.cdf(lo, n * k0, sig)))


rows = []
for sg in ['Neg', 'Pos']:
    fig, axs = plt.subplots(3, 3, figsize=(12, 10))
    for ax, (_, r) in zip(axs.ravel(), R[R.sign == sg].iterrows()):
        a = np.sort(A[(A.dataset == r.dataset) & (A.sign == sg)].amp.values)
        amax = np.percentile(a, 90, method='hazen')
        a1 = a[a <= amax]
        amin = np.percentile(a1, r.truncPct, method='hazen')
        fit = a1[a1 >= amin]
        assert len(fit) == r.nFit, (r.dataset, sg, len(fit), r.nFit)
        lam, k0, sig = r['lambda'], r.k0, r.sigma
        Z = window_mass(lam, k0, sig, amin, amax)
        # histogram: fit window (blue) and excluded amplitudes (grey), same bins
        top = min(np.percentile(a, 99, method='hazen'), 1.5 * amax)   # long tails would squash the window
        nb = 20
        w = (amax - amin) / nb
        edges = np.arange(amin - w * np.ceil(amin / w), top + w, w)
        ax.hist(a[(a < amin) | (a > amax)], bins=edges, color=OUT, label='outside fit window')
        ax.hist(fit, bins=edges, color=BAR, label='fitted amplitudes')
        xs = np.linspace(max(edges[0], 1e-6), edges[-1], 600)
        dens, comp, _, n = mixture(xs, lam, k0, sig)
        inwin = (xs >= amin) & (xs <= amax)
        scale = len(fit) * w / Z
        ax.plot(xs[inwin], dens[inwin] * scale, color=MODEL, lw=2, label='fitted model (window)')
        ax.plot(xs, dens * scale, color=MODEL, lw=1, ls=':', alpha=0.7)
        for j in range(min(4, comp.shape[1])):
            ax.plot(xs, comp[:, j] * scale, color=COMP, lw=0.7, ls='--', alpha=0.6)
        ax.axvline(amin, color=MUTED, lw=0.8); ax.axvline(amax, color=MUTED, lw=0.8)
        # KS against the window-truncated model CDF
        cdf_grid = np.linspace(amin, amax, 2000)
        dg, _, _, _ = mixture(cdf_grid, lam, k0, sig)
        cdf = np.cumsum(dg); cdf /= cdf[-1]
        ks = stats.kstest(fit, lambda x: np.interp(x, cdf_grid, cdf))
        rows.append(dict(dataset=r.dataset, sign=sg, n=len(fit), Amin=amin, Amax=amax, lam=lam, k0=k0,
                         sigma=sig, KS_D=ks.statistic, KS_p=ks.pvalue))
        ax.set_title(f'{r.dataset}  n={len(fit)}\nlambda={lam:.2f}  k0={k0:.0f}  sigma={sig:.0f}  KS p={ks.pvalue:.2f}',
                     fontsize=8.5, loc='left', color=INK)
        ax.set_xlim(edges[0], edges[-1])
        ax.tick_params(labelsize=7, colors=MUTED)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    for ax in axs[-1]:
        ax.set_xlabel('|burst amplitude| (a.u.)', fontsize=8)
    for ax in axs[:, 0]:
        ax.set_ylabel('bursts', fontsize=8)
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle(f'{"Negative" if sg == "Neg" else "Positive"} bursts: amplitude histogram and fitted likelihood model '
                 '(dashed: n = 1-4 unit components; dotted: model outside the fit window)', fontsize=10, color=INK)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(d, f'mle_hist_{sg}.{ext}'), dpi=150)
G = pd.DataFrame(rows)
G.to_csv(os.path.join(d, 'mle_hist_gof.csv'), index=False)
print(G.round(3).to_string(index=False))
