"""PLOT_FUS_TEST_4X  Figure for checks/fus_test_4x.m: 4x PP7 vs PP7-FUS negative
burst amplitudes, R2 re-extraction vs the stored (manuscript) data. Histograms
(same bins per row) and empirical CDFs, with the permutation / LRT p-values.
Usage (repo root): python checks/plot_fus_test_4x.py
Writes analysis_outputs/fus_test_4x/fus_test_4x.png/.pdf."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'fus_test_4x')
A = pd.read_csv(os.path.join(d, 'amplitudes_4x.csv'))
R = pd.read_csv(os.path.join(d, 'fus_test_4x.csv'))
COL = {'4PP7': '#2a78d6', '4PP7-FUS': '#eb6834'}
INK, MUTED = '#1f1e1c', '#898781'
fig, axs = plt.subplots(2, 2, figsize=(11, 8))
for i, (data, label) in enumerate([('R2_reextracted', 'R2 re-extraction (multithresh, red channel, offset)'),
                                   ('stored_MS', 'stored data used in the manuscript')]):
    r = R[R.data == data].iloc[0]
    sub = A[A.data == data]
    top = np.percentile(sub.amp, 95)
    bins = np.linspace(0, top, 21)
    for ds in ['4PP7', '4PP7-FUS']:
        a = sub[sub.dataset == ds].amp.values
        axs[i, 0].hist(np.clip(a, 0, top), bins=bins, density=True, histtype='stepfilled', alpha=0.35,
                       color=COL[ds], label=f'{ds} (n={len(a)}, median {np.median(a):.0f})')
        axs[i, 0].hist(np.clip(a, 0, top), bins=bins, density=True, histtype='step', color=COL[ds], lw=1.5)
        s = np.sort(a)
        axs[i, 1].step(s, np.arange(1, len(s) + 1) / len(s), where='post', color=COL[ds], lw=1.5, label=ds)
    axs[i, 0].set_title(label, fontsize=10, loc='left', color=INK)
    axs[i, 0].legend(fontsize=8, frameon=False)
    axs[i, 1].set_xlim(0, top)
    axs[i, 1].set_title(f'median: perm p {r.pMedianPerm:.4f}   mean: perm p {r.pMeanPerm:.3f}   KS p {r.KS_p:.1e}\n'
                        f'joint (lambda,k0) LRT: perm p {r.pLRperm:.3f}  |  4PP7 {r.lam4PP7:.2f}/{r.k04PP7:.0f}, '
                        f'FUS {r.lamFUS:.2f}/{r.k0FUS:.0f}', fontsize=8.5, loc='left', color=INK)
    axs[i, 1].set_ylabel('cumulative fraction', fontsize=8)
    for ax in axs[i]:
        ax.set_xlabel('|negative burst amplitude| (a.u.; last bin = everything above)', fontsize=8)
        ax.tick_params(labelsize=7, colors=MUTED)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
fig.suptitle('4x RNA: PP7 vs PP7-FUS negative burst amplitudes (granule-level permutation tests)', fontsize=10, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, 'fus_test_4x.' + ext), dpi=150)
print('saved')
