"""WIN_SP  Smoothing window x higher SP (RESULTS.md section 35).
In vivo: checks/detection_grid.m run as detection_grid('in_vivo_reextracted/
after_R2_sums', 'detection_grid_win_sp', [12 15 18 21], [5 9 13 21 31], false):
real traces vs burst-free copies; per window, SP and dataset: negative bursts,
share of turnover the copies reproduce, turnover beyond noise (%/h, nested movie
-> granule bootstrap 95% CI, 1000x); 4x and 30x FUS vs PP7 z on turnover beyond noise.
Simulations: checks/sim_high_sp.m run as sim_high_sp([5 9 13 21 31], [12 15 18 21],
'sim_win_sp.csv'): share of true losses found and false / detected bursts.
Usage: python checks/win_sp.py
Writes analysis_outputs/detection_grid_win_sp/{summary.csv, fus_z.csv, win_sp.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'detection_grid_win_sp'); rng = np.random.default_rng(37); H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
B = pd.read_csv(os.path.join(d, 'bursts.csv')); B = B[B.levelBefore > 0].copy(); B['f'] = B.amp.abs() / B.levelBefore
T = pd.read_csv(os.path.join(d, 'traces.csv'))
rows = []
for ds in VIVO:
    t = T[T.dataset == ds].reset_index(drop=True); L = t.lengthFrames.values.astype(float); mov = t.movie.values
    movs = {m: np.where(mov == m)[0] for m in np.unique(mov)}
    boots = [np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))]) for _ in range(1000)]
    for (w, sp), g in B[B.dataset == ds].groupby(['window', 'SP']):
        fr = t.trace.map(g[g.kind == 'real'].groupby('trace').f.sum()).fillna(0).values
        fs = t.trace.map(g[g.kind == 'surr'].groupby('trace').f.sum()).fillna(0).values
        ex = lambda i: 100 * (fr[i].sum() - fs[i].sum()) / (L[i].sum() * H)
        bs = np.array([ex(i) for i in boots])
        rows.append(dict(dataset=ds, window=w, SP=sp, real=(g.kind == 'real').sum(), copies=(g.kind == 'surr').sum(),
                         noise_share=fs.sum() / fr.sum() if fr.sum() else np.nan, beyond=ex(slice(None)), se=bs.std(),
                         lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
S = pd.DataFrame(rows); S.to_csv(os.path.join(d, 'summary.csv'), index=False)
pd.set_option('display.width', 250)
S['rc'] = S.real.astype(str) + '/' + S.copies.astype(str)
for ds in VIVO:
    print(ds, '(real / copy bursts)'); print(S[S.dataset == ds].pivot(index='SP', columns='window', values='rc').to_string())
Z = []
for (w, sp), g in S.groupby(['window', 'SP']):
    g = g.set_index('dataset')
    for a, b in [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7')]:
        if a not in g.index or b not in g.index:
            continue
        Z.append(dict(window=w, SP=sp, compare=f'{a} vs {b}', z=(g.loc[a, 'beyond'] - g.loc[b, 'beyond']) / np.hypot(g.loc[a, 'se'], g.loc[b, 'se'])))
Z = pd.DataFrame(Z); Z.to_csv(os.path.join(d, 'fus_z.csv'), index=False)
print('\nFUS vs PP7 z (turnover beyond noise):'); print(Z.pivot_table(index=['compare', 'SP'], columns='window', values='z').round(1).to_string())
M = pd.read_csv(os.path.join('analysis_outputs', 'detection_grid_highSP', 'sim_win_sp.csv'))
print('\nsimulations:'); print(M.pivot_table(index=['setting', 'SP'], columns='window', values=['recall', 'false_over_detected']).round(2).to_string())

pp7 = S[S.dataset.isin(VIVO[:4])].groupby(['SP', 'window']).noise_share.median().unstack()
qb = S[S.dataset.isin(VIVO[4:])].groupby(['SP', 'window']).noise_share.median().unstack()
sim = M[M.setting.str.startswith('PP7')]
panels = [(pp7, 'in vivo PP7 (median of 4): share of turnover noise also produces', 'viridis_r', '{:.2f}'),
          (qb, 'in vivo Qβ (median of 2): share of turnover noise also produces', 'viridis_r', '{:.2f}'),
          (sim.pivot(index='SP', columns='window', values='recall') * 100, 'PP7-like simulation: % of true bursts found', 'viridis', '{:.0f}'),
          (sim.pivot(index='SP', columns='window', values='false_over_detected'), 'PP7-like simulation: false / detected bursts', 'viridis_r', '{:.2f}')]
fig, axs = plt.subplots(1, 4, figsize=(24, 4.8))
for ax, (V, ttl, cm, fmt) in zip(axs, panels):
    im = ax.imshow(V.values, cmap=cm, aspect='auto', origin='lower'); fig.colorbar(im, ax=ax)
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            v = V.values[i, j]; ax.text(j, i, fmt.format(v) if v == v else '-', ha='center', va='center', fontsize=9, color='w')
    ax.set_xticks(range(V.shape[1])); ax.set_xticklabels(V.columns); ax.set_yticks(range(V.shape[0])); ax.set_yticklabels(V.index)
    ax.set_xlabel('smoothing window (frames)'); ax.set_ylabel('SP'); ax.set_title(ttl, fontsize=10, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(d, 'win_sp.png'), dpi=130); fig.savefig(os.path.join(d, 'win_sp.pdf'))
