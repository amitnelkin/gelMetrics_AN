"""SLOW_CHANGES  Large slow climbs and drops, counted without the run-length rule,
in real traces vs their burst-free copies (RESULTS.md section 35).
Signal as the detector sees it (analysis.m): x = movmean(S, 13) / Bfit, with S =
spot sum + cell background (R2, as detection_grid.m) and Bfit the cubic fit of the
cell background. An event = x changes by >= thr (20, 30, 40%) relative to its
starting value within <= 90 frames (15 min). Events are counted without overlap:
scanning forward, the first frame that reaches the threshold ends the event and
the scan resumes from there. Rises and drops are counted separately.
Burst-free copies: the dataset's typical bleaching shape (median of the per-trace
cubic trends of S, normalised) scaled to each trace's mean, plus white noise at
the trace's robust frame-to-frame level (as detection_grid.m); 5 copies per trace,
seeded; same Bfit as the real trace.
Rates per granule-hour; real minus copies with a nested movie -> granule bootstrap
95% CI (1000x). Traces >= 180 frames.
Usage: python checks/slow_changes.py [--absolute]
  --absolute: the same amount of light for every granule instead of a percentage of
  its own level. Signal y = movmean(S, 13) / Bfit * Bfit(1) (bleaching-corrected
  intensity); event = |change| >= thr x the median PP7 granule level (pooled over the
  four PP7 datasets); output folder slow_changes_absolute.
Writes analysis_outputs/slow_changes/{summary.csv, slow_changes.png/.pdf}."""
import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ABS = '--absolute' in sys.argv
out = os.path.join('analysis_outputs', 'slow_changes' + ('_absolute' if ABS else '')); os.makedirs(out, exist_ok=True)
CABS = None   # set below in absolute mode
rng = np.random.default_rng(35); H = 10 / 3600; NCOPY = 5; MAXK = 90; THR = [0.2, 0.3, 0.4]
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
R2 = 'in_vivo_reextracted/after_R2_offset'
info = pd.concat([pd.read_csv(f'{R2}/extraction_info.csv'), pd.read_csv(f'{R2}/extraction_info_5Qb_10Qb.csv')])


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def movmean(v, w=13):
    return uniform_filter1d(v, w, mode='nearest')


def events(x, thr):
    up = dn = 0; t = 0; n = len(x)
    while t < n - 1:
        seg = (x[t + 1:min(n, t + 1 + MAXK)] - x[t]) / CABS if ABS else x[t + 1:min(n, t + 1 + MAXK)] / x[t] - 1
        iu = np.argmax(seg >= thr) if np.any(seg >= thr) else None
        idn = np.argmax(seg <= -thr) if np.any(seg <= -thr) else None
        if iu is None and idn is None:
            t += 1; continue
        if idn is None or (iu is not None and iu < idn):
            up += 1; t += iu + 1
        else:
            dn += 1; t += idn + 1
    return up, dn


if ABS:   # pooled median bleaching-corrected level of PP7 granules
    lv = []
    for ds in VIVO[:4]:
        S0 = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_data.mat'); B0 = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_background.mat')
        for s0, b0 in zip(S0, B0):
            s0 = np.ravel(s0).astype(float); b0 = np.ravel(b0).astype(float)
            if len(s0) >= 180:
                t0 = np.arange(len(s0)); bf0 = np.polyval(np.polyfit(t0, b0, 3), t0); lv.append(np.median((s0 + b0) / bf0 * bf0[0]))
    CABS = np.median(lv); print('absolute unit (median PP7 granule level):', round(CABS))
rows, per = [], []
for ds in VIVO:
    S = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_data.mat'); B = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_background.mat')
    I = info[info.dataset == ds]; mov = np.repeat(np.arange(1, len(I) + 1), I.nKept.values)
    keep = [i for i, s in enumerate(S) if np.size(s) >= 180]
    Ss = [np.ravel(S[i]).astype(float) + np.ravel(B[i]).astype(float) for i in keep]; Bs = [np.ravel(B[i]).astype(float) for i in keep]; mk = mov[keep]
    Lm = min(len(s) for s in Ss)
    shp = np.array([np.polyval(np.polyfit(np.arange(len(s)), s, 3), np.arange(len(s)))[:Lm] for s in Ss]); shp /= shp.mean(1, keepdims=True)
    shape = np.median(shp, 0)
    for i, (s, b) in enumerate(zip(Ss, Bs)):
        n = len(s); t = np.arange(n); bf = np.polyval(np.polyfit(t, b, 3), t)
        sh = np.concatenate([shape, np.repeat(shape[-1], max(n - Lm, 0))])[:n]
        sg = 1.4826 * np.median(np.abs(np.diff(s) - np.median(np.diff(s)))) / np.sqrt(2)
        rec = dict(dataset=ds, trace=i + 1, movie=mk[i], hours=n * H)
        for thr in THR:
            u, d = events(movmean(s) / bf * (bf[0] if ABS else 1), thr); rec[f'real_up_{thr}'] = u; rec[f'real_dn_{thr}'] = d
            cu = cd = 0
            for c in range(NCOPY):
                sc = s.mean() / sh.mean() * sh + sg * rng.standard_normal(n)
                a, b_ = events(movmean(sc) / bf * (bf[0] if ABS else 1), thr); cu += a; cd += b_
            rec[f'copy_up_{thr}'] = cu / NCOPY; rec[f'copy_dn_{thr}'] = cd / NCOPY
        per.append(rec)
    print(ds, 'done', flush=True)
P = pd.DataFrame(per); P.to_csv(os.path.join(out, 'per_trace.csv'), index=False)

for ds in VIVO:
    t = P[P.dataset == ds].reset_index(drop=True); movs = {m: np.where(t.movie.values == m)[0] for m in t.movie.unique()}
    boots = [np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))]) for _ in range(1000)]
    for thr in THR:
        for dirn in ['up', 'dn']:
            r, c = t[f'real_{dirn}_{thr}'].values, t[f'copy_{dirn}_{thr}'].values; h = t.hours.values
            f = lambda i: (r[i].sum() - c[i].sum()) / h[i].sum()
            bs = np.array([f(i) for i in boots])
            rows.append(dict(dataset=ds, threshold=int(thr * 100), direction='rise' if dirn == 'up' else 'drop', real_per_granule_h=r.sum() / h.sum(),
                             copy_per_granule_h=c.sum() / h.sum(), excess=f(np.arange(len(t))), lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5),
                             granules_with_event_real=(r > 0).mean(), granules_with_event_copy=(c > 0).mean()))
S_ = pd.DataFrame(rows); S_.to_csv(os.path.join(out, 'summary.csv'), index=False)
pd.set_option('display.width', 250)
S_['v'] = S_.real_per_granule_h.round(2).astype(str) + ' vs ' + S_.copy_per_granule_h.round(2).astype(str) + '  [' + S_.lo.round(2).astype(str) + ', ' + S_.hi.round(2).astype(str) + ']'
print(S_.pivot_table(index=['threshold', 'direction'], columns='dataset', values='v', aggfunc='first')[VIVO].to_string())

fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=False)
COL = {'4PP7': 'tab:blue', '4PP7-FUS': 'tab:green', '30PP7': 'navy', '30PP7-FUS': 'darkgreen', '5Qb': 'tab:red', '10Qb': 'darkred'}
x = np.arange(len(VIVO))
for k, thr in enumerate([20, 30, 40]):
    ax = axs[k]
    for j, dirn in enumerate(['rise', 'drop']):
        s = S_[(S_.threshold == thr) & (S_.direction == dirn)].set_index('dataset').reindex(VIVO)
        off = (j - 0.5) * 0.3
        ax.bar(x + off, s.real_per_granule_h, 0.28, color=[COL[d] for d in VIVO], alpha=0.85 if dirn == 'rise' else 0.5, edgecolor='k', lw=0.4,
               label=f'{dirn}s: real (bar)')
        ax.plot(x + off, s.copy_per_granule_h, 'k_', ms=18, mew=2.5, label=f'{dirn}s: noise-only copies (line)' if k == 0 else None)
    ax.set_xticks(x); ax.set_xticklabels(VIVO, rotation=30, fontsize=9); ax.set_title(f'changes of ≥ {thr}% within 15 min (left bar rises, right bar drops)', fontsize=10, loc='left')
    ax.set_ylabel('events per granule-hour')
axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out, 'slow_changes.png'), dpi=130); fig.savefig(os.path.join(out, 'slow_changes.pdf'))
