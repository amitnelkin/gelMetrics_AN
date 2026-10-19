"""PLOT_SHED_FRACTION  Statistics and figures for checks/shed_fraction.m.
Fraction shed per negative burst = |amplitude| / granule level just before it.
- medians.csv: median fraction per dataset x measure x SP, granule-bootstrap
  95% CI (traces resampled, 1000x).
- tests.csv: granule-level permutation tests (2000 shuffles of whole traces)
  of the difference in median fraction: 5Qb vs 4PP7 (same RNA), 10Qb vs 30PP7,
  Qb vs PP7 (pooled), 4PP7 vs 4PP7-FUS, 30PP7 vs 30PP7-FUS.
- shed_hist.png/pdf: histograms at SP 12 for each measure.
- shed_vs_sp.png/pdf: median fraction vs SP per measure.
Usage (repo root): python checks/plot_shed_fraction.py"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'shed_fraction')
F = pd.read_csv(os.path.join(d, 'fractions.csv'))
F = F[F.sign == 'Neg']
rng = np.random.default_rng(1)
COL = {'4PP7': '#d62a2a', '30PP7': '#f08a8a', '4PP7-FUS': '#2a5bd6', '30PP7-FUS': '#8aa8f0',
       '5Qb': '#1b8a4a', '10Qb': '#7fcf9f', 'protein_RNA_1_1': '#7a4fd6', 'protein_RNA_10_1': '#a98be6',
       'protein_RNA_100_1': '#d3c3f5'}
INK, MUTED = '#1f1e1c', '#898781'


def groups(sub):
    return {t: g.fraction.values for t, g in sub.groupby('trace')}


def boot_median(g, nb=1000):
    keys = np.array(list(g))
    m = [np.median(np.concatenate([g[k] for k in rng.choice(keys, len(keys))])) for _ in range(nb)]
    return np.percentile(m, [2.5, 97.5])


def perm(ga, gb, nb=2000):
    a = list(ga.values()); b = list(gb.values()); allg = a + b; na = len(a)
    obs = np.median(np.concatenate(a)) - np.median(np.concatenate(b))
    cnt = 0
    for _ in range(nb):
        ix = rng.permutation(len(allg))
        x = np.concatenate([allg[i] for i in ix[:na]]); y = np.concatenate([allg[i] for i in ix[na:]])
        cnt += abs(np.median(x) - np.median(y)) >= abs(obs)
    return obs, (cnt + 1) / (nb + 1)


med = []
for (ds, ms, sp), sub in F.groupby(['dataset', 'measure', 'SP']):
    g = groups(sub); lo, hi = boot_median(g)
    med.append(dict(dataset=ds, measure=ms, SP=sp, n=len(sub), nTraces=len(g), median=np.median(sub.fraction), lo=lo, hi=hi))
M = pd.DataFrame(med); M.to_csv(os.path.join(d, 'medians.csv'), index=False)

pairs = [('5Qb', '4PP7'), ('10Qb', '30PP7'), (('5Qb', '10Qb'), ('4PP7', '30PP7')), ('4PP7', '4PP7-FUS'), ('30PP7', '30PP7-FUS')]
tests = []
for ms in ['mean', 'spotsum', 'cropsum']:
    for sp in sorted(F.SP.unique()):
        sub = F[(F.measure == ms) & (F.SP == sp)]
        for a, b in pairs:
            A = a if isinstance(a, tuple) else (a,); B = b if isinstance(b, tuple) else (b,)
            ga = {(ds, t): v for ds in A for t, v in groups(sub[sub.dataset == ds]).items()}
            gb = {(ds, t): v for ds in B for t, v in groups(sub[sub.dataset == ds]).items()}
            if not ga or not gb:
                continue
            obs, p = perm(ga, gb)
            tests.append(dict(measure=ms, SP=sp, compare=' vs '.join(['+'.join(A), '+'.join(B)]),
                              medianA=np.median(np.concatenate(list(ga.values()))),
                              medianB=np.median(np.concatenate(list(gb.values()))), diff=obs, p=p))
T = pd.DataFrame(tests); T.to_csv(os.path.join(d, 'tests.csv'), index=False)
pd.set_option('display.width', 220)
print(M[M.SP == 12].round(3).to_string(index=False)); print()
print(T[T.SP == 12].round(4).to_string(index=False))

# histograms at SP 12
fig, axs = plt.subplots(3, 3, figsize=(15, 11))
bins = np.linspace(0, 1.2, 25)
for i, ms in enumerate(['mean', 'spotsum', 'cropsum']):
    for j, grp in enumerate([('4PP7', '4PP7-FUS', '5Qb'), ('30PP7', '30PP7-FUS', '10Qb'),
                             ('protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1')]):
        ax = axs[i, j]
        for ds in grp:
            x = F[(F.dataset == ds) & (F.measure == ms) & (F.SP == 12)].fraction.values
            if len(x) == 0:
                continue
            ax.hist(np.clip(x, 0, 1.199), bins=bins, weights=np.ones(len(x)) / len(x), histtype='step', lw=1.8,
                    color=COL[ds], label=f'{ds}: n={len(x)}, median {np.median(x):.2f}')
        ax.set_title(f'{ms}, SP 12', fontsize=9, loc='left', color=INK)
        ax.set_xlabel('fraction of granule intensity shed per negative burst (last bin = above)', fontsize=7.5)
        ax.set_ylabel('Frequency', fontsize=8)
        ax.legend(fontsize=7, frameon=False)
        ax.tick_params(labelsize=7, colors=MUTED)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
fig.suptitle('Fraction of the granule shed per negative burst (no model, no k0): PP7 red, PP7-FUS blue, Qbeta green, in vitro purple',
             fontsize=11, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'shed_hist.{ext}'), dpi=140)

fig, axs = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
for ax, ms in zip(axs, ['mean', 'spotsum', 'cropsum']):
    for k, ds in enumerate(COL):
        m = M[(M.dataset == ds) & (M.measure == ms)].sort_values('SP')
        if m.empty:
            continue
        off = (k - 4) * 0.05
        ax.errorbar(m.SP + off, m['median'], yerr=[m['median'] - m.lo, m.hi - m['median']], fmt='o-', ms=3.5,
                    lw=1.2, capsize=2, color=COL[ds], label=ds)
    ax.set_title(ms, fontsize=10, loc='left', color=INK)
    ax.set_xlabel('SP', fontsize=8); ax.set_ylabel('median fraction shed (granule bootstrap 95% CI)', fontsize=8)
    ax.tick_params(labelsize=7, colors=MUTED)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
axs[0].legend(fontsize=7, frameon=False)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'shed_vs_sp.{ext}'), dpi=140)
