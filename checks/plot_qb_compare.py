"""PLOT_QB_COMPARE  Figure for checks/qb_compare.m: stored MS Qbeta data vs the
re-extraction of the Naor .nd2 movies (5Qb\\250ms, 10Qb). No fitting in the
histograms. Row 1: negative-burst |amplitude| histograms at the paper SP (14),
stored vs re-extracted on one common axis (frequency). Row 2: same for
positive bursts. Row 3: median |negative amplitude| and bursts per trace vs SP.
Usage (repo root): python checks/plot_qb_compare.py
Writes analysis_outputs/qb_compare/qb_compare.png/.pdf and prints the tables."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'qb_compare')
A = pd.read_csv(os.path.join(d, 'amplitudes.csv'))
Bu = pd.read_csv(os.path.join(d, 'bursts.csv'))
Tr = pd.read_csv(os.path.join(d, 'traces.csv'))
Fi = pd.read_csv(os.path.join(d, 'fits.csv'))
COL = {'stored_MS': '#d62a2a', 'R2_reextracted': '#2a78d6'}
LAB = {'stored_MS': 'stored (MS)', 'R2_reextracted': 're-extracted'}
INK, MUTED = '#1f1e1c', '#898781'
fig, axs = plt.subplots(3, 2, figsize=(12, 12))
for j, ds in enumerate(['5Qb', '10Qb']):
    for i, sg in enumerate(['Neg', 'Pos']):
        ax = axs[i, j]
        sub = A[(A.dataset == ds) & (A.SP == 14) & (A.sign == sg)]
        top = np.percentile(sub.amp, 97.5)
        edges = np.linspace(0, top, 21)
        for v in ['stored_MS', 'R2_reextracted']:
            a = sub[sub.data == v].amp.values
            if len(a) == 0:
                continue
            ax.hist(np.clip(a, 0, top * 0.999), bins=edges, weights=np.ones(len(a)) / len(a), color=COL[v], alpha=0.45,
                    label=f'{LAB[v]}: n={len(a)}, median {np.median(a):.0f}')
            ax.axvline(np.median(a), color=COL[v], ls='--', lw=1.2)
        ax.set_title(f'{ds}, {"negative" if sg == "Neg" else "positive"} bursts, SP 14', fontsize=10, loc='left', color=INK)
        ax.set_xlabel('|amplitude| (a.u.; last bin = above)', fontsize=8); ax.set_ylabel('Frequency', fontsize=8)
        ax.legend(fontsize=8, frameon=False)
    ax = axs[2, j]
    for v in ['stored_MS', 'R2_reextracted']:
        b = Bu[(Bu.dataset == ds) & (Bu.data == v)]
        ax.plot(b.SP, b.medianAbsNeg, 'o-', color=COL[v], label=f'{LAB[v]} median |neg amp|')
        ax2 = ax.twinx() if v == 'stored_MS' else ax2
        ax2.plot(b.SP, b.burstsPerTrace, 's:', color=COL[v], alpha=0.7, label=f'{LAB[v]} bursts/trace')
    ax.set_xlabel('SP', fontsize=8); ax.set_ylabel('median |negative amplitude|', fontsize=8)
    ax2.set_ylabel('bursts per trace (dotted)', fontsize=8)
    ax.legend(fontsize=7, frameon=False, loc='upper left'); ax2.legend(fontsize=7, frameon=False, loc='upper right')
    ax.set_title(f'{ds}: vs SP', fontsize=10, loc='left', color=INK)
for ax in axs.ravel():
    ax.tick_params(labelsize=7, colors=MUTED)
    for s in ('top',):
        ax.spines[s].set_visible(False)
fig.suptitle('Qbeta: stored MS data vs re-extraction of the Naor .nd2 movies (red channel, multithresh spot, offset 200); no fitting',
             fontsize=11, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'qb_compare.{ext}'), dpi=140)
pd.set_option('display.width', 250)
print(Tr.round(0).to_string(index=False)); print()
print(Bu[Bu.SP == 14].round(1).to_string(index=False)); print()
print(Fi.round(2).to_string(index=False))
