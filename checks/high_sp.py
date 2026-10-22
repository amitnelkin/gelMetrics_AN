"""HIGH_SP  Real traces vs their burst-free copies at increasing stringency SP
(RESULTS.md section 34). Question: if SP is raised until noise alone produces
(almost) no bursts, is anything real left?
Inputs: checks/example_traces_sp.m (example traces at SP 12-25) and
checks/detection_grid.m run as detection_grid('in_vivo_reextracted/after_R2_sums',
'detection_grid_highSP', [10 12 14 16 18 20 22 25], 13, false): in vivo R2 spot sum,
window 13, traces >= 180 frames, negative bursts, burst-free copies = the
dataset's typical bleaching shape + white noise at each trace's own level.
Per dataset and SP: negative bursts per granule-hour in the real traces and in
the copies, share of turnover reproduced by the copies, turnover beyond noise
(%/h) with a nested movie -> granule bootstrap 95% CI (1000x).
Usage: python checks/high_sp.py
Writes analysis_outputs/detection_grid_highSP/{summary.csv, high_sp.png/.pdf}
and analysis_outputs/noise_background/example_traces_sp_<dataset>.png."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.default_rng(34); H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
COL = {'4PP7': 'tab:blue', '4PP7-FUS': 'tab:green', '30PP7': 'navy', '30PP7-FUS': 'darkgreen', '5Qb': 'tab:red', '10Qb': 'darkred'}
d = os.path.join('analysis_outputs', 'detection_grid_highSP'); nb = os.path.join('analysis_outputs', 'noise_background')

# ---- examples
E = pd.read_csv(os.path.join(nb, 'example_traces_sp.csv')); EB = pd.read_csv(os.path.join(nb, 'example_bursts_sp.csv'))
for ds in E.dataset.unique():
    sps = sorted(EB[EB.dataset == ds].SP.unique()) or [12]
    sps = sorted(set(sps) | {12, 15, 18, 21, 25})
    fig, axs = plt.subplots(2, len(sps), figsize=(4.2 * len(sps), 6.2), sharex=True, sharey=True)
    for r, kind in enumerate(['real', 'surrogate']):
        e = E[(E.dataset == ds) & (E.kind == kind)]; tm = (e.frame - 1) / 6
        for c, sp in enumerate(sps):
            ax = axs[r, c]; ax.plot(tm, e.processed, color='0.25', lw=0.9)
            b = EB[(EB.dataset == ds) & (EB.kind == kind) & (EB.SP == sp)]
            for _, s in b.iterrows():
                ax.plot([(s.firstFrame - 1) / 6, (s.lastFrame - 1) / 6], [s.levelStart, s.levelEnd], color='tab:red' if s.sign < 0 else 'tab:blue', lw=3)
            ax.text(0.02, 0.96, f"{(b.sign < 0).sum()} neg, {(b.sign > 0).sum()} pos", transform=ax.transAxes, va='top', fontsize=9,
                    bbox=dict(fc='w', ec='none', alpha=0.8))
            if r == 0:
                ax.set_title(f'SP = {sp}', fontsize=11)
            if c == 0:
                ax.set_ylabel(('real trace' if kind == 'real' else 'burst-free copy') + '\n(detector input)')
            if r == 1:
                ax.set_xlabel('time (min)')
    fig.suptitle(f'{ds}: same real trace and same burst-free copy, detector at increasing SP (window 13). Red = loss, blue = gain', fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(nb, f'example_traces_sp_{ds}.png'), dpi=120)

# ---- all granules
B = pd.read_csv(os.path.join(d, 'bursts.csv')); B = B[(B.levelBefore > 0)].copy(); B['f'] = B.amp.abs() / B.levelBefore
T = pd.read_csv(os.path.join(d, 'traces.csv'))
rows = []
for ds in VIVO:
    t = T[T.dataset == ds].reset_index(drop=True); L = t.lengthFrames.values.astype(float); mov = t.movie.values
    movs = {m: np.where(mov == m)[0] for m in np.unique(mov)}
    boots = [np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))]) for _ in range(1000)]
    for sp in sorted(B.SP.unique()):
        b = B[(B.dataset == ds) & (B.SP == sp)]
        fr = t.trace.map(b[b.kind == 'real'].groupby('trace').f.sum()).fillna(0).values
        fs = t.trace.map(b[b.kind == 'surr'].groupby('trace').f.sum()).fillna(0).values
        nr = t.trace.map(b[b.kind == 'real'].groupby('trace').size()).fillna(0).values
        ns = t.trace.map(b[b.kind == 'surr'].groupby('trace').size()).fillna(0).values
        ex = lambda i: 100 * (fr[i].sum() - fs[i].sum()) / (L[i].sum() * H)
        bs = np.array([ex(i) for i in boots]); hours = L.sum() * H
        rows.append(dict(dataset=ds, SP=sp, granules=len(t), realBursts=int(nr.sum()), copyBursts=int(ns.sum()),
                         real_per_granule_h=nr.sum() / hours, copy_per_granule_h=ns.sum() / hours,
                         noise_share=fs.sum() / fr.sum() if fr.sum() else np.nan, turnover_beyond=ex(slice(None)),
                         lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
S = pd.DataFrame(rows); S.to_csv(os.path.join(d, 'summary.csv'), index=False)
pd.set_option('display.width', 220)
print(S.pivot(index='SP', columns='dataset', values='realBursts')[VIVO].to_string(), '\n')
print(S.pivot(index='SP', columns='dataset', values='copyBursts')[VIVO].to_string(), '\n')
print(S.pivot(index='SP', columns='dataset', values='noise_share')[VIVO].round(2).to_string(), '\n')
S['v'] = S.turnover_beyond.round(1).astype(str) + ' [' + S.lo.round(1).astype(str) + ', ' + S.hi.round(1).astype(str) + ']'
print(S.pivot(index='SP', columns='dataset', values='v')[VIVO].to_string())

fig, axs = plt.subplots(1, 3, figsize=(18, 5.2))
for ds in VIVO:
    s = S[S.dataset == ds]
    axs[0].plot(s.SP, s.real_per_granule_h, 'o-', color=COL[ds], label=f'{ds} real')
    axs[0].plot(s.SP, s.copy_per_granule_h, 'x--', color=COL[ds], alpha=0.7)
    axs[1].plot(s.SP, s.noise_share, 'o-', color=COL[ds], label=ds)
    axs[2].errorbar(s.SP + (VIVO.index(ds) - 2.5) * 0.12, s.turnover_beyond, yerr=[s.turnover_beyond - s.lo, s.hi - s.turnover_beyond], fmt='o-', color=COL[ds], capsize=2, label=ds)
axs[0].set_yscale('log'); axs[0].set_ylabel('negative bursts per granule-hour\n(solid: real traces, dashed: burst-free copies)')
axs[1].axhline(1, color='k', ls=':', lw=0.8); axs[1].set_ylabel('share of turnover that noise also produces'); axs[1].set_ylim(0, 1.6)
axs[2].axhline(0, color='k', lw=0.8); axs[2].set_ylabel('turnover beyond noise (%/h, 95% CI)')
for ax, ttl in zip(axs, ['(a) bursts detected', '(b) share explained by noise', '(c) what remains beyond noise']):
    ax.set_xlabel('stringency SP'); ax.set_title(ttl, fontsize=10, loc='left')
axs[1].legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(d, 'high_sp.png'), dpi=130); fig.savefig(os.path.join(d, 'high_sp.pdf'))
