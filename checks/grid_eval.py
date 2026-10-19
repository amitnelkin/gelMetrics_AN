"""GRID_EVAL  Evaluates checks/detection_grid.m (SP x smoothing window, real and
burst-free surrogate). RESULTS.md section 20; reviewer points R2 6, R2 2b, R2 4a.
Per setting and dataset: real and false negative-burst counts, false-discovery
share (false / real), uncorrected turnover and excess turnover (%/h, real minus
its paired surrogate), nested (movie -> granule) bootstrap SE, 500x.
Setting choice rule, fixed before looking at any comparison: the setting with
the lowest median false-discovery share over the six in vivo datasets among
settings that keep >= 50 real negative bursts in every in vivo dataset. The
comparisons (z of excess turnover) are then read off at that setting; the
whole grid is shown so that the choice can be checked.
Usage: python checks/grid_eval.py
Writes analysis_outputs/detection_grid/{grid.csv, comparisons.csv, chosen.json, grid.png}."""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'detection_grid')
rng = np.random.default_rng(5)
H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
ORDER = VIVO + ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7'),
         ('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1')]

B = pd.read_csv(os.path.join(d, 'bursts.csv')); B = B[B.levelBefore > 0].copy(); B['frac'] = B.amp.abs() / B.levelBefore
T = pd.read_csv(os.path.join(d, 'traces.csv'))
TR = {ds: g for ds, g in T.groupby('dataset')}


def stats_for(ds, w, sp, nb=500):
    t = TR[ds]; L = dict(zip(t.trace, t.lengthFrames)); idx = t.trace.values; mov = t.movie.values
    b = B[(B.dataset == ds) & (B.window == w) & (B.SP == sp)]
    R = {k: g.frac.values for k, g in b[b.kind == 'real'].groupby('trace')}
    S = {k: g.frac.values for k, g in b[b.kind == 'surr'].groupby('trace')}
    def exc(p):
        return 100 * (sum(R.get(k, np.zeros(0)).sum() for k in p) - sum(S.get(k, np.zeros(0)).sum() for k in p)) / (sum(L[k] for k in p) * H)
    movs = {m: idx[mov == m] for m in np.unique(mov)}
    bs = []
    for _ in range(nb):
        if len(movs) > 1:
            pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        else:
            pick = rng.choice(idx, len(idx))
        bs.append(exc(pick))
    nr, ns = (b.kind == 'real').sum(), (b.kind == 'surr').sum()
    hours = t.lengthFrames.sum() * H
    return dict(dataset=ds, window=w, SP=sp, nReal=nr, nFalse=ns, FDR=ns / max(nr, 1),
                turnover=100 * b[b.kind == 'real'].frac.sum() / hours, excess=exc(idx), excess_se=np.std(bs))


rows = [stats_for(ds, w, sp) for w in sorted(B.window.unique()) for sp in sorted(B.SP.unique()) for ds in ORDER]
G = pd.DataFrame(rows); G.to_csv(os.path.join(d, 'grid.csv'), index=False)
vivo = G[G.dataset.isin(VIVO)]
summ = vivo.groupby(['window', 'SP']).agg(medianFDR=('FDR', 'median'), worstFDR=('FDR', 'max'), minReal=('nReal', 'min')).reset_index()
ok = summ[summ.minReal >= 50]
best = ok.sort_values('medianFDR').iloc[0]
json.dump(dict(window=int(best.window), SP=int(best.SP), medianFDR=float(best.medianFDR), worstFDR=float(best.worstFDR),
               rule='lowest median in vivo false-discovery share with >= 50 real negative bursts per in vivo dataset'),
          open(os.path.join(d, 'chosen.json'), 'w'), indent=1)
pd.set_option('display.width', 250)
print('median in vivo false-discovery share (rows SP, cols window); [min real bursts]:')
print(summ.assign(v=summ.medianFDR.round(2).astype(str) + ' [' + summ.minReal.astype(int).astype(str) + ']')
      .pivot(index='SP', columns='window', values='v').to_string())
print(f"\nchosen: window {int(best.window)}, SP {int(best.SP)} (median FDR {best.medianFDR:.2f}, worst {best.worstFDR:.2f})")
C = []
for (w, sp), g in G.groupby(['window', 'SP']):
    g = g.set_index('dataset')
    for a, b in COMPS:
        C.append(dict(window=w, SP=sp, compare=f'{a} vs {b}', excessA=g.loc[a, 'excess'], excessB=g.loc[b, 'excess'],
                      z=(g.loc[a, 'excess'] - g.loc[b, 'excess']) / np.hypot(g.loc[a, 'excess_se'], g.loc[b, 'excess_se'])))
C = pd.DataFrame(C); C.to_csv(os.path.join(d, 'comparisons.csv'), index=False)
sel = C[(C.window == best.window) & (C.SP == best.SP)]
print('\nat the chosen setting (excess turnover %/h, nested z):'); print(sel.round(2).to_string(index=False))
print('\n4PP7-FUS vs 4PP7 z across the grid (rows SP, cols window):')
print(C[C['compare'] == '4PP7-FUS vs 4PP7'].pivot(index='SP', columns='window', values='z').round(1).to_string())
print('\nexcess turnover at the chosen setting per dataset:')
print(G[(G.window == best.window) & (G.SP == best.SP)].set_index('dataset').reindex(ORDER)[['nReal', 'nFalse', 'FDR', 'turnover', 'excess', 'excess_se']].round(2).to_string())
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
for ax, (val, lab, cmap) in zip(axs, [(summ.pivot(index='SP', columns='window', values='medianFDR'), 'median in vivo false-discovery share', 'viridis'),
                                       (C[C['compare'] == '4PP7-FUS vs 4PP7'].pivot(index='SP', columns='window', values='z'), '4PP7-FUS vs 4PP7, z (excess turnover)', 'coolwarm')]):
    im = ax.imshow(val.values, cmap=cmap, aspect='auto', origin='lower')
    ax.set_xticks(range(val.shape[1])); ax.set_xticklabels(val.columns); ax.set_yticks(range(val.shape[0])); ax.set_yticklabels(val.index)
    for i in range(val.shape[0]):
        for j in range(val.shape[1]):
            ax.text(j, i, f'{val.values[i, j]:.2f}' if 'share' in lab else f'{val.values[i, j]:.1f}', ha='center', va='center', fontsize=7, color='w')
    ax.set_xlabel('smoothing window (frames)'); ax.set_ylabel('SP'); ax.set_title(lab, fontsize=10, loc='left'); fig.colorbar(im, ax=ax)
    ax.plot(list(val.columns).index(best.window), list(val.index).index(best.SP), 's', ms=18, mfc='none', mec='k', mew=2)
fig.tight_layout(); fig.savefig(os.path.join(d, 'grid.png'), dpi=140); fig.savefig(os.path.join(d, 'grid.pdf'))
print('\nfalse-discovery share, 4PP7 / 4PP7-FUS (rows SP, cols window):')
for ds in ['4PP7', '4PP7-FUS']:
    print(ds); print(G[G.dataset == ds].pivot(index='SP', columns='window', values='FDR').round(2).to_string())
print('\nall comparisons, z at every window for SP 12 and SP 5:')
for sp in [12, 5]:
    print('SP', sp); print(C[C.SP == sp].pivot(index='compare', columns='window', values='z').round(1).to_string())
