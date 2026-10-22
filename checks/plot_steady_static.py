"""PLOT_STEADY_STATIC  The per-frame (R2) vs steady-threshold (R4) comparison of
section 23, with the static-granule control drawn next to every real value
(RESULTS.md section 23, addendum). The static control is a synthetic granule
that never changes - each real granule's average image kept on its real track,
bleached, with fresh camera + photon noise (checks/static_control.m) - measured
with the same extraction and the same downstream steps as the real data.
Panels: (1) frame noise, sqrt(D1/2) of the relative granule intensity (median
over granules; D1 = mean squared change over one frame); (2) false / real negative
bursts (window 13, SP 12, burst-free copies as in section 18); (3) turnover beyond
noise (%/h) with 95% CI (nested movie -> granule bootstrap, 1000x).
Usage: python checks/plot_steady_static.py
Writes analysis_outputs/steady/steady_with_static.png/.pdf and steady_with_static.csv."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'steady'); rng = np.random.default_rng(32); H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
VER = {'v0': ('per-frame (R2)', 'analysis_outputs/noise_methods/real_sf.csv', '0.55'), 'v31': ('steady (R4)', 'analysis_outputs/steady/real_sf.csv', 'tab:green')}
D = pd.read_csv(os.path.join(d, 'detection.csv')); D = D[(D.window == 13) & (D.SP == 12)]
SC = pd.concat([pd.read_csv(os.path.join(d, 'static_control.csv')), pd.read_csv(os.path.join(d, 'static_control_more.csv'))])


def excess_ci(t):
    f = lambda x: 100 * (x.fracReal13.sum() - x.fracSurr13.sum()) / (x.lengthFrames.sum() * H)
    movs = {m: g.index.values for m, g in t.groupby('movie')}; v = []
    for _ in range(1000):
        pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        v.append(f(t.loc[pick]))
    return f(t), np.percentile(v, [2.5, 97.5])


rows = []
for v, (lab, sf, _) in VER.items():
    R = pd.read_csv(sf); R = R[R.measure == 'spotsum']
    for ds in VIVO:
        Dr = np.array([float(x.strip('[]').split()[0]) for x in R[R.dataset == ds].Dreal])
        s = SC[(SC.dataset == ds) & (SC.version == v)].reset_index(drop=True)
        dr = D[(D.version == lab) & (D.dataset == ds)].iloc[0]
        se, sci = excess_ci(s)
        rows.append(dict(version=lab, dataset=ds, noise_real=100 * np.sqrt(np.median(Dr) / 2), noise_static=100 * np.sqrt(s.D1.dropna().median() / 2),
                         FDR_real=dr.FDR, FDR_static=s.nSurr13.sum() / max(s.nReal13.sum(), 1),
                         excess_real=dr.excess, excess_real_lo=dr.excess - 1.96 * dr.excess_se, excess_real_hi=dr.excess + 1.96 * dr.excess_se,
                         excess_static=se, excess_static_lo=sci[0], excess_static_hi=sci[1]))
T = pd.DataFrame(rows); T.to_csv(os.path.join(d, 'steady_with_static.csv'), index=False)
pd.set_option('display.width', 250); print(T.round(2).to_string(index=False))

fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
x = np.arange(len(VIVO)); wd = 0.2
for k, (v, (lab, _, col)) in enumerate(VER.items()):
    t = T[T.version == lab].set_index('dataset').reindex(VIVO)
    for q, (key, hatch, alpha, lab2) in enumerate([('real', None, 1.0, 'real granules'), ('static', '///', 0.35, 'unchanging granule (control)')]):
        off = (2 * k + q - 1.5) * wd
        axs[0].bar(x + off, t[f'noise_{key}'], wd, color=col, alpha=alpha, hatch=hatch, edgecolor='k', lw=0.4, label=f'{lab}: {lab2}')
        axs[1].bar(x + off, t[f'FDR_{key}'], wd, color=col, alpha=alpha, hatch=hatch, edgecolor='k', lw=0.4, label=f'{lab}: {lab2}')
        m = t[f'excess_{key}']
        axs[2].errorbar(x + off, m, yerr=[m - t[f'excess_{key}_lo'], t[f'excess_{key}_hi'] - m], fmt='o', color=col, mfc=col if key == 'real' else 'w',
                        ms=7, capsize=3, label=f'{lab}: {lab2}')
axs[0].set_ylabel('frame-to-frame noise (% of granule)')
axs[1].set_ylabel('false / real negative bursts (window 13, SP 12)'); axs[1].axhline(1, color='k', lw=0.8, ls=':')
axs[2].set_ylabel('turnover beyond noise (%/h, 95% CI)'); axs[2].axhline(0, color='k', lw=0.8)
for ax, ttl in zip(axs, ['(1) noise: the unchanging granule matches the real noise', '(2) share of bursts that noise also produces',
                         '(3) turnover beyond noise: compare filled (real) with open (control)']):
    ax.set_xticks(x); ax.set_xticklabels(VIVO, rotation=30, fontsize=9); ax.set_title(ttl, fontsize=10, loc='left')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, fontsize=9, frameon=False)
fig.tight_layout(rect=[0, 0.07, 1, 1]); fig.savefig(os.path.join(d, 'steady_with_static.png'), dpi=140); fig.savefig(os.path.join(d, 'steady_with_static.pdf'))
