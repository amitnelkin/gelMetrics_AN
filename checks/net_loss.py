"""NET_LOSS  Net loss or gain of granule material, measured without burst
detection (RESULTS.md section 22). Per granule (traces >= 180 frames):
  cell fade rate     slope of ln(cell background mean) vs time
  granule fade rate  slope of ln(granule) vs time
  net change rate    slope of ln(granule / cell) = granule rate - cell rate
The cell background is the same fluorescent protein under the same light, so it
bleaches at the same rate. A granule fading faster than its surroundings is
losing material (net < 0); slower means it is gaining. Slopes: Theil-Sen (robust
to outlier frames and to bursts), reported in %/h (100 * slope * 360 frames/h).
Granule measures:
  spotsum   integrated granule intensity above the cell (R2 sums; in vivo only)
  excess    spot mean - cell mean (R2 offset data in vivo; stored in vitro data
            with the 200 offset subtracted - assumed, in vitro camera not confirmed)
Frames with granule or cell <= 0 are dropped; traces with > 10% dropped frames
are skipped. Dataset value = mean over granules; CI and z from a nested (movie
-> granule) bootstrap, 2000x (in vitro: granule bootstrap, no movie ids).
Usage: python checks/net_loss.py
Writes analysis_outputs/net_loss/{per_trace.csv, summary.csv, comparisons.csv, net_loss.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import theilslopes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'net_loss'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(22)
R2, SU, IV = 'in_vivo_reextracted/after_R2_offset', 'in_vivo_reextracted/after_R2_sums', 'in_vitro_experiments'
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
VITRO = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7'), ('30PP7', '4PP7'),
         ('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1')]
PER_H = 100 * 360


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def rate(y):
    t = np.arange(len(y))
    ok = y > 0
    if ok.mean() < 0.9:
        return np.nan
    return theilslopes(np.log(y[ok]), t[ok])[0] * PER_H


info = pd.concat([pd.read_csv(f'{R2}/extraction_info.csv'), pd.read_csv(f'{R2}/extraction_info_5Qb_10Qb.csv')])
rows = []
for ds in VIVO + VITRO:
    if ds in VIVO:
        S, B = load(f'{R2}/{ds}_data.mat'), load(f'{R2}/{ds}_background.mat')
        SS = load(f'{SU}/{ds}_spotsum_data.mat'); off = 0
        I = info[info.dataset == ds]; mov = np.repeat(np.arange(1, len(I) + 1), I.nKept.values)
    else:
        S, B = load(f'{IV}/{ds}_data.mat'), load(f'{IV}/{ds}_background.mat'); SS = [None] * len(S); off = 200
        mov = np.zeros(len(S), int)
    for i, (s, b, ss) in enumerate(zip(S, B, SS)):
        s = np.ravel(s).astype(float) - off; b = np.ravel(b).astype(float) - off
        if len(s) < 180:
            continue
        r = dict(dataset=ds, trace=i + 1, movie=mov[i], lengthFrames=len(s), cell=rate(b), excessGranule=rate(s - b), excessNet=rate((s - b) / b))
        if ss is not None:
            ss = np.ravel(ss).astype(float); r.update(sumGranule=rate(ss), sumNet=rate(ss / b))
        rows.append(r)
P = pd.DataFrame(rows); P.to_csv(os.path.join(out, 'per_trace.csv'), index=False)


def boot(t, col, nb=2000):
    t = t.dropna(subset=[col]); v = t[col].values; m = t.movie.values; movs = np.unique(m)
    idx = {k: np.where(m == k)[0] for k in movs}; res = []
    for _ in range(nb):
        if len(movs) > 1:
            pick = np.concatenate([rng.choice(idx[k], len(idx[k])) for k in rng.choice(movs, len(movs))])
        else:
            pick = rng.choice(len(v), len(v))
        res.append(v[pick].mean())
    return v.mean(), np.std(res), np.percentile(res, [2.5, 97.5]), len(v)


MEAS = ['cell', 'sumGranule', 'sumNet', 'excessGranule', 'excessNet']
summ, est = [], {}
for ds in VIVO + VITRO:
    t = P[P.dataset == ds]
    for col in MEAS:
        if t[col].notna().sum() < 5:
            continue
        m, se, ci, n = boot(t, col); est[(ds, col)] = (m, se)
        summ.append(dict(dataset=ds, measure=col, mean_pct_per_h=m, ci_lo=ci[0], ci_hi=ci[1], se=se, granules=n,
                         movies=t.movie.nunique(), share_negative=(t[col].dropna() < 0).mean()))
Sm = pd.DataFrame(summ); Sm.to_csv(os.path.join(out, 'summary.csv'), index=False)
comp = []
for a, b in COMPS:
    for col in ['sumNet', 'excessNet']:
        if (a, col) in est and (b, col) in est:
            (ma, sa), (mb, sb) = est[(a, col)], est[(b, col)]
            comp.append(dict(compare=f'{a} vs {b}', measure=col, A=ma, B=mb, diff=ma - mb, z=(ma - mb) / np.hypot(sa, sb)))
Cp = pd.DataFrame(comp); Cp.to_csv(os.path.join(out, 'comparisons.csv'), index=False)
pd.set_option('display.width', 250)
Sm['v'] = Sm.mean_pct_per_h.round(1).astype(str) + ' [' + Sm.ci_lo.round(1).astype(str) + ', ' + Sm.ci_hi.round(1).astype(str) + ']'
print('mean rate %/h [95% CI, nested bootstrap]:')
print(Sm.pivot(index='dataset', columns='measure', values='v').reindex(VIVO + VITRO)[MEAS].to_string())
print('\nshare of granules with net < 0 (sumNet / excessNet):')
print(Sm[Sm.measure.isin(['sumNet', 'excessNet'])].pivot(index='dataset', columns='measure', values='share_negative').reindex(VIVO + VITRO).round(2).to_string())
print('\ncomparisons (z):'); print(Cp.round(2).to_string(index=False))

fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
for ax, col, lab in [(axs[0], 'sumNet', 'net change, granule sum vs cell (in vivo)'), (axs[1], 'excessNet', 'net change, granule mean above cell vs cell')]:
    dss = [d for d in VIVO + VITRO if (d, col) in est]
    for k, ds in enumerate(dss):
        v = P[P.dataset == ds][col].dropna()
        ax.plot(k + rng.uniform(-0.25, 0.25, len(v)), v.clip(-300, 300), '.', ms=3, alpha=0.35, color='0.5')
        r = Sm[(Sm.dataset == ds) & (Sm.measure == col)].iloc[0]
        ax.errorbar(k, r.mean_pct_per_h, yerr=[[r.mean_pct_per_h - r.ci_lo], [r.ci_hi - r.mean_pct_per_h]], fmt='o', color='tab:red', capsize=4)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(range(len(dss))); ax.set_xticklabels([d.replace('protein_RNA_', 'in vitro ').replace('_', ':') for d in dss], rotation=30, fontsize=8)
    ax.set_ylabel('%/h  (< 0: granule loses material faster than it bleaches)'); ax.set_title(lab, fontsize=10, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(out, 'net_loss.png'), dpi=140); fig.savefig(os.path.join(out, 'net_loss.pdf'))
