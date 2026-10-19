"""PLOT_MLE_REAL  Figures for checks/mle_real_data.m (R2 4a/4b, R2 5, R1 2).
1) mle_real_estimates.png/.pdf: MLE lambda and k0 per dataset with granule-
   bootstrap 95% CIs (thick) and profile CIs (thin), grid-score estimate as x.
2) mle_real_surfaces.png/.pdf: lambda-k0 NLL surface per dataset (negative
   bursts), contours at 2*dNLL = 2.30 / 5.99 (joint 68% / 95%), MLE marked.
Usage (repo root): python checks/plot_mle_real.py"""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'mle_real')
R = pd.read_csv(os.path.join(d, 'mle_real.csv'))
NEG, POS, GRID, INK, MUTED = '#2a78d6', '#eb6834', '#1f1e1c', '#1f1e1c', '#898781'
names = list(dict.fromkeys(R.dataset))
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
for ax, par, lo, hi, plo, phi, g, lab in [
        (axs[0], 'lambda', 'lamLoBoot', 'lamHiBoot', 'lamLoProf', 'lamHiProf', 'gridLam', 'lambda'),
        (axs[1], 'k0', 'k0LoBoot', 'k0HiBoot', 'k0LoProf', 'k0HiProf', 'gridK0', 'k0 (a.u.)')]:
    for i, nm in enumerate(names):
        for j, (sg, col) in enumerate([('Neg', NEG), ('Pos', POS)]):
            r = R[(R.dataset == nm) & (R.sign == sg)].iloc[0]
            y = i + (j - 0.5) * 0.3
            if par == 'lambda':
                ax.plot([r[plo], r[phi]], [y, y], color=col, lw=0.8, alpha=0.6)
            ax.plot([r[lo], r[hi]], [y, y], color=col, lw=3, solid_capstyle='round',
                    label=(('negative' if sg == 'Neg' else 'positive') + ' bursts, MLE + bootstrap CI') if i == 0 else None)
            ax.plot(r[par], y, 'o', color='white', ms=5, mec=col, mew=1.5)
            ax.plot(r[g], y, 'x', color=GRID, ms=5, label='grid-score method' if (i == 0 and j == 0) else None)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlabel(lab + ' (log scale)', fontsize=9, color=INK)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.grid(axis='x', color='#e6e4df', lw=0.5)
    ax.tick_params(labelsize=8, colors=MUTED)
axs[0].legend(fontsize=7, frameon=False, loc='lower right')
axs[0].set_title('lambda (thin line: profile CI)', fontsize=10, loc='left', color=INK)
axs[1].set_title('k0', fontsize=10, loc='left', color=INK)
fig.suptitle('Likelihood fit per dataset; PP7 = R2 re-extracted, in vitro / Qb = stored data', fontsize=10, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, 'mle_real_estimates.' + ext), dpi=150)

S = sio.loadmat(os.path.join(d, 'surfaces.mat'), squeeze_me=True, struct_as_record=False)
fig, axs = plt.subplots(3, 3, figsize=(11, 10))
for ax, nm in zip(axs.ravel(), names):
    want = nm.replace('-', '_')
    if want[0].isdigit():
        want = 'x' + want          # matlab.lang.makeValidName prefix
    key = [k for k in S if k == want + '_Neg']
    s = S[key[0]]
    L, K = np.meshgrid(s.lam, s.k0, indexing='ij')
    z = 2 * np.asarray(s.nll)
    ax.contourf(L, K, np.minimum(z, 20), levels=20, cmap='Blues_r')
    ax.contour(L, K, z, levels=[2.30, 5.99], colors=[INK, INK], linewidths=[1.2, 0.8], linestyles=['-', '--'])
    r = R[(R.dataset == nm) & (R.sign == 'Neg')].iloc[0]
    ax.plot(r['lambda'], r['k0'], 'o', color='#eb6834', ms=5)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(nm, fontsize=9, loc='left', color=INK)
    ax.tick_params(labelsize=7, colors=MUTED)
for ax in axs[-1]:
    ax.set_xlabel('lambda', fontsize=8)
for ax in axs[:, 0]:
    ax.set_ylabel('k0', fontsize=8)
fig.suptitle('Negative bursts: 2*dNLL surface (solid 68%, dashed 95% joint region); the ridge is lambda*k0 = const',
             fontsize=10, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, 'mle_real_surfaces.' + ext), dpi=150)
print('saved')
