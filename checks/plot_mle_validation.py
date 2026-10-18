"""PLOT_MLE_VALIDATION  Figure for checks/mle_sim_validation.m (R2 4b):
estimated vs true lambda and k0 (MLE with profile and granule-bootstrap 95% CIs,
and the grid-score method), negative and positive bursts.
Usage (repo root): python checks/plot_mle_validation.py
Writes analysis_outputs/mle_validation/sim_validation.png and .pdf."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'mle_validation')
V = pd.read_csv(os.path.join(d, 'sim_validation.csv'))
MLE, BOOT, GRID, INK, MUTED = '#2a78d6', '#1baf7a', '#eb6834', '#1f1e1c', '#898781'
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
for c, sg in enumerate(['Neg', 'Pos']):
    v = V[V.sign == sg]
    for r, (par, true) in enumerate([('lam', None), ('k0', 60)]):
        ax = axs[r, c]
        for j, (_, row) in enumerate(v.iterrows()):
            x = row.lamTrue + (row.seed - 2) * 0.08
            hat = row[par + 'Hat']
            ax.errorbar(x - 0.03, hat, yerr=[[hat - row[par + 'LoProf']], [row[par + 'HiProf'] - hat]],
                        fmt='o', ms=4, color=MLE, lw=1, label='MLE, profile CI' if j == 0 else None)
            if np.isfinite(row[par + 'LoBoot']):
                ax.errorbar(x + 0.03, hat, yerr=[[max(hat - row[par + 'LoBoot'], 0)], [max(row[par + 'HiBoot'] - hat, 0)]],
                            fmt='none', color=BOOT, lw=1, label='granule bootstrap CI' if j == 0 else None)
            g = row.gridLam if par == 'lam' else row.gridK0
            ax.plot(x, g, 'x', color=GRID, ms=6, label='grid-score method' if j == 0 else None)
        lt = np.array([0.5, 1, 1.5, 2, 3, 4, 5])
        ax.plot(lt, lt if par == 'lam' else np.full_like(lt, true), color=MUTED, lw=1, ls='--', label='truth')
        ax.set_xlabel('true lambda', color=INK, fontsize=9)
        ax.set_ylabel('estimated ' + ('lambda' if par == 'lam' else 'k0'), color=INK, fontsize=9)
        ax.set_title(f"{'negative' if sg == 'Neg' else 'positive'} bursts: {'lambda' if par == 'lam' else 'k0 (true 60)'}",
                     fontsize=10, color=INK, loc='left')
        if par == 'k0':
            ax.set_ylim(0, 300)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=8, colors=MUTED)
        ax.grid(color='#e6e4df', lw=0.5)
axs[0, 0].legend(fontsize=8, frameon=False)
fig.suptitle('Full-pipeline simulation (generateSim -> analysis.m -> fit), 150 tracks, noise 40, 3 seeds per lambda',
             fontsize=10, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, 'sim_validation.' + ext), dpi=150)
print('saved')
