"""PLOT_RUNNER_OUTPUTS  Figures for the runner outputs of one stage.

MATLAB figure rendering hangs in headless mode on this machine, so the
runners save their data (checks/run_experimental.m, run_simulation.m) and
this script draws the figures from it.

Usage (repo root):  py -3 checks/plot_runner_outputs.py <stage>
Reads  analysis_outputs/<stage>/*_run.mat
Writes analysis_outputs/<stage>/figs/<dataset>.png and .pdf, one per dataset:
  row 1: segment amplitudes (positive / negative / quiescent)
  row 2: per sign, the amplitudes that were fitted (normalised histogram,
         10 bins as in modifiedPoissonFit) with the rank-1 model's unit
         positions n*k0 (n >= 1) weighted by Poisson(n; lambda)
  row 3: per sign, the weighted-score surface (lambda x k0; lower = better)
"""
import glob
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.stats import poisson

stage = sys.argv[1]
base = os.path.join('analysis_outputs', stage)
fig_dir = os.path.join(base, 'figs')
os.makedirs(fig_dir, exist_ok=True)

for f in sorted(glob.glob(os.path.join(base, '*_run.mat'))):
    R = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
    name = str(R['name'])
    fits = np.atleast_1d(R['fits'])
    fits = {str(e.title): e for e in fits}
    fig, ax = plt.subplots(3, 2, figsize=(11, 11))
    gs = ax[0, 0].get_gridspec()
    for a in ax[0, :]:
        a.remove()
    a0 = fig.add_subplot(gs[0, :])
    amps = [np.atleast_1d(R['ampPos']), np.atleast_1d(R['ampNeg']), np.atleast_1d(R['ampVar'])]
    amps = [x[np.isfinite(x)] for x in amps]
    allv = np.concatenate(amps)
    lim = np.percentile(np.abs(allv), 99) if allv.size else 1
    bins = np.linspace(-lim, lim, 80)
    for x, c, lab in zip(amps, ['tab:green', 'tab:red', 'tab:blue'], ['positive', 'negative', 'quiescent']):
        a0.hist(x, bins=bins, color=c, alpha=0.55, label=f'{lab} (n={x.size})')
    a0.set_xlabel('segment amplitude [A.U.]'); a0.set_ylabel('count'); a0.legend()
    a0.set_title(f'{name}  ({stage})')
    for j, t in enumerate(['Positive amplitudes', 'Negative amplitudes']):
        if t not in fits:
            continue
        e = fits[t]
        d = np.atleast_1d(e.data).astype(float)
        lam1 = float(np.atleast_1d(e.lamOrder)[0]); k1 = float(np.atleast_1d(e.k0Order)[0])
        a = ax[1, j]
        cnt, edges = np.histogram(d, bins=10, density=True)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        a.bar(ctr, cnt / cnt.max(), width=np.diff(edges), color='0.75', edgecolor='k', label='fitted amplitudes')
        n = np.arange(1, int(d.max() / k1) + 2)
        pm = poisson.pmf(n, lam1)
        a.stem(n * k1, pm / pm.max(), linefmt='r-', markerfmt='ro', basefmt=' ',
               label=f'rank 1: λ={lam1:g}, k0={k1:g}')
        a.set_xlim(0, d.max() * 1.05); a.set_xlabel('|amplitude| [A.U.]'); a.set_ylabel('normalised density')
        a.set_title(f'{t[:3]}: n={d.size}'); a.legend(fontsize=8)
        W = np.atleast_2d(e.w_score); K = np.atleast_1d(e.K0Vector).astype(float); L = np.atleast_1d(getattr(e, 'lambda'))
        a = ax[2, j]
        im = a.imshow(W, aspect='auto', origin='lower', cmap='viridis',
                      extent=[K[0], K[-1], L[0] - 0.5, L[-1] + 0.5], vmin=0, vmax=1)
        for li, lv in enumerate(L):
            ki = int(np.argmin(W[li])); a.plot(K[ki], lv, 'w*', ms=9)
        a.set_xlabel('k0'); a.set_ylabel('λ'); a.set_title('weighted score (★ = best k0 per λ)')
        fig.colorbar(im, ax=a)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(fig_dir, f'{name}.{ext}'), dpi=110)
    plt.close(fig)
    print('wrote', name)
