"""PLOT_CROP_DIAGNOSIS  Summary and figure for checks/crop_sum_diagnosis.m.
1. Radial profiles (normalised to the centre) per dataset.
2. Median fraction shed per negative burst for each intensity measure, and
   granule-level permutation tests (2000 shuffles) Qb vs PP7 (pooled) and
   5Qb vs 4PP7 (same RNA).
3. Where the shed light goes: for spot-sum negative bursts, the change of the
   cell ring 4-8 px around the granule (positive = the cell gets brighter).
Usage (repo root): python checks/plot_crop_diagnosis.py
Writes analysis_outputs/crop_diagnosis/summary.csv, tests.csv, crop_diagnosis.png/.pdf."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'crop_diagnosis')
R = pd.read_csv(os.path.join(d, 'radial.csv'))
F = pd.read_csv(os.path.join(d, 'fractions.csv'))
N = pd.read_csv(os.path.join(d, 'nearcell.csv'))
rng = np.random.default_rng(1)
COL = {'4PP7': '#d62a2a', '30PP7': '#f08a8a', '4PP7-FUS': '#2a5bd6', '30PP7-FUS': '#8aa8f0', '5Qb': '#1b8a4a', '10Qb': '#7fcf9f'}
INK, MUTED = '#1f1e1c', '#898781'
measures = ['spotsum', 'ringsum', 'box7', 'box10', 'box14', 'box20']


def groups(sub, keys):
    return {(ds, t): g.fraction.values for ds in keys for t, g in sub[sub.dataset == ds].groupby('trace')}


def perm(ga, gb, nb=2000):
    a = list(ga.values()); b = list(gb.values()); allg = a + b; na = len(a)
    obs = np.median(np.concatenate(a)) - np.median(np.concatenate(b)); cnt = 0
    for _ in range(nb):
        ix = rng.permutation(len(allg))
        cnt += abs(np.median(np.concatenate([allg[i] for i in ix[:na]])) - np.median(np.concatenate([allg[i] for i in ix[na:]]))) >= abs(obs)
    return obs, (cnt + 1) / (nb + 1)


S = F.groupby(['dataset', 'measure']).fraction.median().unstack('measure')[measures]
L = F.groupby(['dataset', 'measure']).levelBefore.median().unstack('measure')[measures]
S.to_csv(os.path.join(d, 'summary.csv'))
print('median fraction shed (negative bursts, SP 12):'); print(S.round(3).to_string()); print()
print('median granule level before burst:'); print(L.round(0).to_string()); print()
T = []
for m in measures:
    sub = F[F.measure == m]
    for A, B in [(('5Qb', '10Qb'), ('4PP7', '30PP7')), (('5Qb',), ('4PP7',)), (('4PP7-FUS',), ('4PP7',))]:
        o, p = perm(groups(sub, A), groups(sub, B))
        T.append(dict(measure=m, compare='+'.join(A) + ' vs ' + '+'.join(B), diff=o, p=p))
T = pd.DataFrame(T); T.to_csv(os.path.join(d, 'tests.csv'), index=False)
print(T.pivot(index='measure', columns='compare', values='p').reindex(measures).round(4).to_string()); print()
N['cellUp'] = N.dNearCellMean > 0
print('where the shed light goes (spot-sum negative bursts):')
print(N.groupby('dataset').agg(bursts=('trace', 'size'), medianSpotDrop=('dSpotSum', 'median'),
                               medianNearCellChange=('dNearCellMean', 'median'), shareCellUp=('cellUp', 'mean'),
                               nearCellLevel=('nearCellBefore', 'median')).round(2).to_string())

fig, axs = plt.subplots(1, 3, figsize=(17, 5))
ax = axs[0]
for ds, g in R.groupby('dataset'):
    ax.plot(g.radius, g.medianIntensity / g.medianIntensity.iloc[0], 'o-', ms=3, color=COL[ds], label=ds)
for x, lab in [(7, 'box 15x15'), (10, 'box 21'), (14, 'box 29')]:
    ax.axvline(x, color=MUTED, lw=0.8, ls=':'); ax.text(x + 0.2, 0.9, lab, fontsize=7, color=MUTED)
ax.set_xlabel('distance from granule centre (px)', fontsize=8); ax.set_ylabel('intensity / centre (bg = 18-20 px ring)', fontsize=8)
ax.set_title('radial profile', fontsize=10, loc='left', color=INK); ax.legend(fontsize=7, frameon=False)
ax = axs[1]
for k, ds in enumerate(COL):
    ax.plot(range(len(measures)), S.loc[ds, measures].values, 'o-', color=COL[ds], label=ds)
ax.set_xticks(range(len(measures))); ax.set_xticklabels(measures, fontsize=8)
ax.set_ylabel('median fraction shed per negative burst (SP 12)', fontsize=8)
ax.set_title('fraction shed by intensity measure', fontsize=10, loc='left', color=INK)
ax = axs[2]
for ds in COL:
    n = N[N.dataset == ds]
    ax.scatter(-n.dSpotSum / n.dSpotSum.abs().median(), n.dNearCellMean / n.nearCellBefore.abs().median(), s=6, alpha=0.4, color=COL[ds], label=ds)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlabel('granule drop (spot sum, relative to its median drop)', fontsize=8)
ax.set_ylabel('change of the 4-8 px cell ring (relative to its level)', fontsize=8)
ax.set_title('does the surrounding cell brighten when the granule sheds?', fontsize=10, loc='left', color=INK)
for a in axs:
    a.tick_params(labelsize=7, colors=MUTED)
    for s in ('top', 'right'):
        a.spines[s].set_visible(False)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'crop_diagnosis.{ext}'), dpi=140)
