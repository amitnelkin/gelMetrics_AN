"""CONTROL_TABLE  Minimal per-dataset table of publication-worthy parameters,
each with its noise control (RESULTS.md section 27).
Columns (95% CI: nested movie -> granule bootstrap in vivo, granule bootstrap
in vitro, 2000x, seeded):
  granules (movies)
  granule brightness       median spot mean per granule (offset 200 removed)
  frame noise              robust frame-to-frame noise of the granule signal,
                           % of its level (spot sum in vivo, spot mean in vitro)
  mean burst size          mean negative-burst amplitude, all bursts, no trim
                           (spot mean units, window 13, SP 12), model-free;
                           also the median (robust to the long tail)
  burst size, granules / copies  median fraction of the granule shed per burst,
                           granules and their burst-free copies (spot sum), and
                           the difference = burst size beyond noise, with CI
  lambda*k0, shape fit     product of lambda and k0 from the shape-only fit of the
                           manuscript's model curve on all bursts (checks/
                           lambda_k0_shape.py --notrim, granule bootstrap CI)
  turnover beyond noise    %/h (sections 18-24), and the share of the detected
                           turnover that burst-free copies reproduce (noise share)
  net change               %/h, granule vs cell-background fade (section 22)
Usage: python checks/control_table.py
Writes analysis_outputs/control_table/{control_table.csv, control_table.md, control_table.png}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'control_table'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(28); NB = 2000; H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
VITRO = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
LAB = {'4PP7': 'PP7 + QCP-5x/PCP-4x', '4PP7-FUS': 'PP7-FUS + QCP-5x/PCP-4x', '30PP7': 'PP7 + PCP-30x', '30PP7-FUS': 'PP7-FUS + PCP-30x',
       '5Qb': 'Qβ + QCP-5x/PCP-4x', '10Qb': 'Qβ + QCP-10x', 'protein_RNA_1_1': 'in vitro, protein:RNA 1:1',
       'protein_RNA_10_1': 'in vitro, protein:RNA 10:1', 'protein_RNA_100_1': 'in vitro, protein:RNA 100:1'}


def nested(units, stat):
    """units: {movie: list of granule items}; percentile CI of stat(flat list)"""
    keys = list(units); v = []
    for _ in range(NB):
        it = []
        for k in rng.choice(len(keys), len(keys)):
            g = units[keys[k]]; it.extend(g[i] for i in rng.integers(0, len(g), len(g)))
        v.append(stat(it))
    return np.percentile(v, [2.5, 97.5])


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def rs(x):
    d = np.diff(x); return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.median(x)


AMP = pd.read_csv('analysis_outputs/lambda_k0_shape/amplitudes.csv')
BR = pd.read_csv('analysis_outputs/noise_background/traces.csv'); BR = BR[BR.lengthFrames >= 180]
GB = pd.read_csv('analysis_outputs/detection_grid/bursts.csv')
GB = GB[(GB.window == 13) & (GB.SP == 12) & (GB.levelBefore > 0)].copy(); GB['f'] = GB.amp.abs() / GB.levelBefore
GT = pd.read_csv('analysis_outputs/detection_grid/traces.csv')
NL = pd.read_csv('analysis_outputs/net_loss/per_trace.csv')
SH = pd.read_csv('analysis_outputs/lambda_k0_shape_notrim/fits.csv')
rows = []
for ds in VIVO + VITRO:
    vivo = ds in VIVO
    r = dict(dataset=ds, label=LAB[ds])
    # granules / movies, brightness, noise
    if vivo:
        b = BR[BR.dataset == ds]
        r['brightness'] = b.spotMean.median(); r['noise'] = 100 * (b.sdSum_meas / b.spotSum).median()
        U = {m: list(g.spotMean.values) for m, g in b.groupby('movie')}
        r['brightness_lo'], r['brightness_hi'] = nested(U, np.median)
    else:
        S = load(f'in_vitro_experiments-old/{ds}_data.mat')
        S = [np.ravel(s).astype(float) - 200 for s in S if np.size(s) >= 180]
        lv = np.array([np.median(s) for s in S]); r['brightness'] = np.median(lv); r['noise'] = 100 * np.median([rs(s) for s in S])
        r['brightness_lo'], r['brightness_hi'] = nested({0: list(lv)}, np.median)
    t = GT[GT.dataset == ds]; r['granules'] = len(t); r['movies'] = t.movie.nunique() if vivo else np.nan
    # mean burst size (all bursts, no trim)
    a = AMP[AMP.dataset == ds]
    U = {m: [g2.amp.values for _, g2 in g.groupby('trace')] for m, g in a.groupby('movie')}
    r['bursts'] = len(a); r['meanBurst'] = a.amp.mean(); r['medianBurst'] = a.amp.median()
    r['meanBurst_lo'], r['meanBurst_hi'] = nested(U, lambda it: np.concatenate(it).mean() if len(it) else np.nan)
    r['medianBurst_lo'], r['medianBurst_hi'] = nested(U, lambda it: np.median(np.concatenate(it)) if len(it) else np.nan)
    # turnover beyond noise + noise share
    bb = GB[GB.dataset == ds]
    R = bb[bb.kind == 'real'].groupby('trace').f.sum(); Su = bb[bb.kind == 'surr'].groupby('trace').f.sum()
    items = {m: [(R.get(k, 0.0), Su.get(k, 0.0), L) for k, L in zip(g.trace, g.lengthFrames)] for m, g in t.groupby('movie')}
    turn = lambda it: 100 * (sum(x[0] for x in it) - sum(x[1] for x in it)) / (sum(x[2] for x in it) * H)
    flat = [x for g in items.values() for x in g]
    r['turnover'] = turn(flat); r['turnover_lo'], r['turnover_hi'] = nested(items, turn)
    r['noiseShare'] = sum(x[1] for x in flat) / sum(x[0] for x in flat)
    # burst size of the granules vs their burst-free copies (median % of granule shed), and the part beyond noise
    Rb = {k: g.f.values for k, g in bb[bb.kind == 'real'].groupby('trace')}; Sb = {k: g.f.values for k, g in bb[bb.kind == 'surr'].groupby('trace')}
    e0 = np.zeros(0)
    sz = {m: [(Rb.get(k, e0), Sb.get(k, e0)) for k in g.trace] for m, g in t.groupby('movie')}
    med = lambda it, j: 100 * np.median(np.concatenate([x[j] for x in it]))
    fl = [x for g in sz.values() for x in g]
    r['sizeReal'], r['sizeCopy'] = med(fl, 0), med(fl, 1); r['sizeBeyond'] = r['sizeReal'] - r['sizeCopy']
    r['sizeBeyond_lo'], r['sizeBeyond_hi'] = nested(sz, lambda it: med(it, 0) - med(it, 1))
    # lambda*k0 from the shape fit on all bursts (checks/lambda_k0_shape.py --notrim)
    f = SH[SH.dataset == ds].iloc[0]
    r['lamk0'], r['lamk0_lo'], r['lamk0_hi'] = f['mean'], f['mean_lo'], f['mean_hi']
    # net change
    col = 'sumNet' if vivo else 'excessNet'
    n = NL[(NL.dataset == ds)].dropna(subset=[col])
    U = {m: list(g[col].values) for m, g in n.groupby('movie')}
    r['net'] = n[col].mean(); r['net_lo'], r['net_hi'] = nested(U, np.mean)
    rows.append(r); print(ds, 'done', flush=True)
T = pd.DataFrame(rows); T.to_csv(os.path.join(out, 'control_table.csv'), index=False)

ci = lambda r, c, f='.1f': f"{r[c]:{f}} [{r[c + '_lo']:{f}}, {r[c + '_hi']:{f}}]"
hdr = ['dataset', 'granules (movies)', 'granule brightness (a.u.)', 'frame noise (% of granule)',
       'burst size, granules / noise copies (% of granule)', 'burst size beyond noise (% points)', 'median burst size (a.u.)', 'mean burst size (a.u.)', 'λ·k0, shape fit (a.u.)',
       'turnover beyond noise (%/h)', 'net change beyond bleaching (%/h)']
cells = []
for _, r in T.iterrows():
    cells.append([r.label, f"{r.granules}" + (f" ({int(r.movies)})" if r.dataset in VIVO else ''), ci(r, 'brightness', '.0f'), f"{r.noise:.1f}",
                  f"{r.sizeReal:.1f} / {r.sizeCopy:.1f}", ci(r, 'sizeBeyond'), ci(r, 'medianBurst', '.0f'), ci(r, 'meanBurst', '.0f'), ci(r, 'lamk0', '.0f'), ci(r, 'turnover'), ci(r, 'net')])
with open(os.path.join(out, 'control_table.md'), 'w', encoding='utf-8', newline=chr(13) + chr(10)) as fh:
    fh.write('| ' + ' | '.join(hdr) + ' |' + chr(10) + '|' + '---|' * len(hdr) + chr(10))
    for c in cells:
        fh.write('| ' + ' | '.join(c) + ' |' + chr(10))
print(open(os.path.join(out, 'control_table.md'), encoding='utf-8').read())
fig, ax = plt.subplots(figsize=(26, 4.8)); ax.axis('off')
tb = ax.table(cellText=cells, colLabels=[h.replace(', ', ',\n').replace(' (', '\n(').replace('turnover explained', 'turnover\nexplained') for h in hdr], loc='center', cellLoc='center')
tb.auto_set_font_size(False); tb.set_fontsize(8.5); tb.scale(1, 2.1)
for (i, j), c in tb.get_celld().items():
    if i == 0:
        c.set_facecolor('#e8e8e8'); c.set_height(0.2); c.set_text_props(weight='bold')
    elif j == 0:
        c.set_text_props(ha='left')
    if i == len(VIVO) + 1 and i > 0:
        c.visible_edges = 'BRTL'
fig.text(0.01, -0.03, 'In vivo: R2 extraction (published multithresh, offset 200), traces ≥ 30 min; turnover and noise from the granule spot sum; brightness and mean burst size in spot-mean units. '
         'In vitro: stored data (offset 200 assumed). Window 13, SP 12, all bursts (no trim). 95% CI: movie → granule bootstrap (in vitro: granules).', fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(out, 'control_table.png'), dpi=150, bbox_inches='tight')
