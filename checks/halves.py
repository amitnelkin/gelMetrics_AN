"""HALVES  First vs second 30 minutes of each movie: are the reported parameters
stable over the hour (granule ageing, drifting measurement)? (RESULTS.md section 36;
reviewer 2 point 5, time-uniformity.)
Bursts: checks/detection_grid.m output (R2 spot sum in vivo, mean in vitro;
window 13, SP 12; negative bursts on the real traces and on their burst-free
copies), assigned to the half in which they start (frame <= 180 = first half).
Per half: median burst size (fraction of the granule shed) in the real traces and
in the copies, burst size beyond noise, turnover beyond noise (%/h, per granule-hour
of that half).
Net change beyond bleaching per half: robust (Theil-Sen) slope of
ln(granule / cell background) fitted separately to frames 1-180 and 181-end
(in vivo: spot sum / cell mean, R2 sums; in vitro: (spot - cell) / cell, stored data,
offset 200 assumed); second halves shorter than 60 frames are skipped.
Second minus first half, with a nested movie -> granule bootstrap 95% CI and
two-sided p (2000x); granules keep both their halves (paired). Traces >= 180 frames.
Usage: python checks/halves.py
Writes analysis_outputs/halves/{summary.csv, halves.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import theilslopes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'halves'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(38); H = 10 / 3600; NB = 2000; CUT = 180
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']; VITRO = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
LAB = {'protein_RNA_1_1': 'in vitro 1:1', 'protein_RNA_10_1': 'in vitro 10:1', 'protein_RNA_100_1': 'in vitro 100:1'}
G = pd.read_csv('analysis_outputs/detection_grid/bursts.csv')
G = G[(G.window == 13) & (G.SP == 12) & (G.levelBefore > 0)].copy(); G['f'] = G.amp.abs() / G.levelBefore
T = pd.read_csv('analysis_outputs/detection_grid/traces.csv')


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def slope(y):
    y = np.asarray(y, float); ok = y > 0
    if ok.sum() < 30:
        return np.nan
    t = np.arange(len(y)); return theilslopes(np.log(y[ok]), t[ok])[0] * 360 * 100


def net_halves(ds):
    if ds in VIVO:
        S = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_data.mat'); B = load(f'in_vivo_reextracted/after_R2_sums/{ds}_spotsum_background.mat')
        sig = [np.ravel(s).astype(float) / np.ravel(b).astype(float) for s, b in zip(S, B) if np.size(s) >= 180]
    else:
        S = load(f'in_vitro_experiments-old/{ds}_data.mat'); B = load(f'in_vitro_experiments-old/{ds}_background.mat')
        sig = [(np.ravel(s).astype(float) - np.ravel(b).astype(float)) / (np.ravel(b).astype(float) - 200) for s, b in zip(S, B) if np.size(s) >= 180]
    return [(slope(x[:CUT]), slope(x[CUT:]) if len(x) - CUT >= 60 else np.nan) for x in sig]


rows, items = [], {}
for ds in VIVO + VITRO:
    t = T[T.dataset == ds].reset_index(drop=True); b = G[G.dataset == ds]; nets = net_halves(ds)
    assert len(nets) == len(t), (ds, len(nets), len(t))
    gr = []
    for k, (tr, L, m) in enumerate(zip(t.trace, t.lengthFrames, t.movie)):
        bb = b[b.trace == tr]; e = {}
        for h, sel in [(1, bb.startFrame <= CUT), (2, bb.startFrame > CUT)]:
            e[f'r{h}'] = bb[sel & (bb.kind == 'real')].f.values; e[f'c{h}'] = bb[sel & (bb.kind == 'surr')].f.values
        e['h1'] = min(L, CUT) * H; e['h2'] = max(L - CUT, 0) * H; e['n1'], e['n2'] = nets[k]; e['movie'] = m
        gr.append(e)
    items[ds] = gr


def stats(it):
    o = {}
    for h in (1, 2):
        r = np.concatenate([x[f'r{h}'] for x in it]); c = np.concatenate([x[f'c{h}'] for x in it]); hrs = sum(x[f'h{h}'] for x in it)
        o[f'size{h}'] = 100 * np.median(r) if len(r) else np.nan
        o[f'beyond{h}'] = 100 * (np.median(r) - np.median(c)) if len(r) and len(c) else np.nan
        o[f'turn{h}'] = 100 * (r.sum() - c.sum()) / hrs if hrs else np.nan
        n = np.array([x[f'n{h}'] for x in it], float); o[f'net{h}'] = np.nanmean(n)
    return o


KEYS = [('size', 'median burst size, granules (% of granule)'), ('beyond', 'burst size beyond noise (% points)'),
        ('turn', 'turnover beyond noise (%/h)'), ('net', 'net change beyond bleaching (%/h)')]
for ds, it in items.items():
    obs = stats(it); movs = {}
    for x in it:
        movs.setdefault(x['movie'], []).append(x)
    keys = list(movs); bs = []
    for _ in range(NB):
        res = []
        for kk in rng.choice(len(keys), len(keys)):
            g = movs[keys[kk]]; res.extend(g[i] for i in rng.integers(0, len(g), len(g)))
        bs.append(stats(res))
    bs = pd.DataFrame(bs)
    for key, lab in KEYS:
        d = bs[f'{key}2'] - bs[f'{key}1']; d = d.dropna()
        diff = obs[f'{key}2'] - obs[f'{key}1']
        p = 2 * min((d <= 0).mean(), (d >= 0).mean()) if len(d) else np.nan
        rows.append(dict(dataset=ds, measure=key, label=lab, first=obs[f'{key}1'], second=obs[f'{key}2'], diff=diff,
                         lo=np.percentile(d, 2.5) if len(d) else np.nan, hi=np.percentile(d, 97.5) if len(d) else np.nan, p=max(p, 1 / NB)))
    print(ds, 'done', flush=True)
R = pd.DataFrame(rows); R.to_csv(os.path.join(out, 'summary.csv'), index=False)
pd.set_option('display.width', 250)
R['v'] = R['first'].round(1).astype(str) + ' -> ' + R['second'].round(1).astype(str) + '  (' + R['diff'].round(1).astype(str) + ' [' + R.lo.round(1).astype(str) + ', ' + R.hi.round(1).astype(str) + '], p ' + R.p.round(3).astype(str) + ')'
print(R.pivot(index='dataset', columns='measure', values='v').reindex(VIVO + VITRO)[[k for k, _ in KEYS]].to_string())

fig, axs = plt.subplots(1, 4, figsize=(22, 5))
DS = VIVO + VITRO; x = np.arange(len(DS))
for ax, (key, lab) in zip(axs, KEYS):
    r = R[R.measure == key].set_index('dataset').reindex(DS)
    ax.bar(x - 0.2, r['first'], 0.38, color='0.7', edgecolor='k', lw=0.4, label='first 30 min')
    ax.bar(x + 0.2, r['second'], 0.38, color='tab:purple', alpha=0.75, edgecolor='k', lw=0.4, label='second 30 min')
    for i, (_, q) in enumerate(r.iterrows()):
        if q.p < 0.05:
            ax.text(i, max(q['first'], q['second']), '*', ha='center', va='bottom', fontsize=14)
    ax.axhline(0, color='k', lw=0.7); ax.set_xticks(x); ax.set_xticklabels([LAB.get(d, d) for d in DS], rotation=35, ha='right', fontsize=8)
    ax.set_title(lab, fontsize=10, loc='left')
axs[0].legend(fontsize=8)
fig.suptitle('First vs second 30 minutes (* = difference p < 0.05, movie-level bootstrap)', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, 'halves.png'), dpi=130); fig.savefig(os.path.join(out, 'halves.pdf'))
