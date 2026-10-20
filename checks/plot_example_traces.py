"""PLOT_EXAMPLE_TRACES  Real trace vs its burst-free surrogate, raw and as the
detector sees it, with detected bursts (checks/example_traces.m output).
Usage: python checks/plot_example_traces.py
Writes analysis_outputs/noise_background/example_traces.png/.pdf."""
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'noise_background')
T = pd.read_csv(os.path.join(d, 'example_traces.csv')); B = pd.read_csv(os.path.join(d, 'example_bursts.csv'))
NAMES = {'4PP7-FUS': '4PP7-FUS (in vivo, spot sum)', 'protein_RNA_10_1': 'in vitro 10:1 (spot mean)'}
fig, axs = plt.subplots(4, 2, figsize=(12, 11), sharex=True)
for r, ds in enumerate(NAMES):
    for c, kind in enumerate(['real', 'surrogate']):
        t = T[(T.dataset == ds) & (T.kind == kind)]; b = B[(B.dataset == ds) & (B.kind == kind)]
        tm = (t.frame - 1) * 10 / 60
        ax = axs[2 * r, c]
        ax.plot(tm, t.raw, color='0.35', lw=0.7)
        ttl = 'real trace' if kind == 'real' else 'simulated, no bursts: typical bleaching of the dataset + this trace\'s noise'
        ax.set_title(f'{NAMES[ds]} - {ttl}', fontsize=8.5, loc='left')
        ax.set_ylabel('raw intensity')
        ax = axs[2 * r + 1, c]
        ax.plot(tm, t.processed, color='0.2', lw=1)
        for _, s in b.iterrows():
            col = 'tab:red' if s.sign < 0 else 'tab:blue'
            ax.plot([(s.firstFrame - 1) / 6, (s.lastFrame - 1) / 6], [s.levelStart, s.levelEnd], color=col, lw=3, alpha=0.8)
        nn, npos = (b.sign < 0).sum(), (b.sign > 0).sum()
        ax.text(0.01, 0.95, f'{nn} negative (red), {npos} positive (blue) bursts detected', transform=ax.transAxes,
                va='top', fontsize=8, bbox=dict(fc='w', ec='none', alpha=0.8))
        ax.set_ylabel('smoothed, bleaching-\ncorrected (detector input)')
    for k in (2 * r, 2 * r + 1):
        lo = min(axs[k, 0].get_ylim()[0], axs[k, 1].get_ylim()[0]); hi = max(axs[k, 0].get_ylim()[1], axs[k, 1].get_ylim()[1])
        axs[k, 0].set_ylim(lo, hi); axs[k, 1].set_ylim(lo, hi)
for ax in axs[-1]:
    ax.set_xlabel('time (min)')
fig.suptitle('Same detector (window 13, SP 12) on a real trace and on a burst-free simulated copy of it', fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(d, 'example_traces.png'), dpi=140); fig.savefig(os.path.join(d, 'example_traces.pdf'))
