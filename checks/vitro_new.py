"""VITRO_NEW  In vitro parameters from the re-extracted titration movies
(RESULTS.md section 38), with movies as replicates for the first time.
Mapping, from trace matching (checks/vitro_compare.py): protein_RNA_1_1 <-
1000ng_ul_RNA_0.1mg_ml_protein, protein_RNA_10_1 <- 1000ng_RNA_10ug_protein,
protein_RNA_100_1 <- 1000ng_ul_RNA_1mg_ml_protein; extra_100ngul_1mgml <-
100ng_ul_RNA_1mg_ml_protein (not in the stored data).
Data: in_vitro_reextracted/after_R2_offset (fluorescence channel, offset 200 removed,
published multithresh spot), traces >= 180 frames. Bursts: detection_grid.m run as
detection_grid('', 'detection_grid_vitro_new', 12, 13, true,
'in_vitro_reextracted/after_R2_offset', {...}, false): spot mean, window 13, SP 12,
negative bursts on the traces and on burst-free copies.
Per dataset: granules (movies), brightness (median spot mean), frame noise (robust,
% of level), median burst amplitude, median burst size (fraction of granule shed)
real / copies and the part beyond noise, turnover beyond noise (%/h), net change
beyond bleaching (%/h; Theil-Sen slope of ln((spot - cell) / cell)), and turnover /
net change per 30-min half. 95% CI: nested movie -> granule bootstrap (2000x).
Comparisons: movie-level permutation p (10000x) for 10:1 vs 1:1, 100:1 vs 1:1 and
extra vs 1:1. The stored-trace values (sections 22, 27-29, 36) are listed alongside.
Usage: python checks/vitro_new.py
Writes analysis_outputs/vitro_new/{summary.csv, tests.csv, summary.md, vitro_new.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import theilslopes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'vitro_new'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(39); H = 10 / 3600; NB, NP, CUT = 2000, 10000, 180
DS = {'protein_RNA_1_1': 'vitro_1000ngul_0p1mgml', 'protein_RNA_10_1': 'vitro_1000ng_10ug',
      'protein_RNA_100_1': 'vitro_1000ngul_1mgml', 'extra_100ngul_1mgml': 'vitro_100ngul_1mgml'}
LAB = {'protein_RNA_1_1': '1:1', 'protein_RNA_10_1': '10:1', 'protein_RNA_100_1': '100:1', 'extra_100ngul_1mgml': 'extra (100 ng/µl RNA, 1 mg/ml protein)'}
R2 = 'in_vitro_reextracted/after_R2_offset'
G = pd.read_csv('analysis_outputs/detection_grid_vitro_new/bursts.csv')
G = G[(G.levelBefore > 0)].copy(); G['f'] = G.amp.abs() / G.levelBefore
T = pd.read_csv('analysis_outputs/detection_grid_vitro_new/traces.csv')


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return [np.ravel(x).astype(float) for x in m[[k for k in m if not k.startswith('__')][0]][1:]]


def slope(y):
    y = np.asarray(y, float); ok = np.isfinite(y) & (y > 0)
    if ok.sum() < 30:
        return np.nan
    t = np.arange(len(y)); return theilslopes(np.log(y[ok]), t[ok])[0] * 360 * 100


def rs(x):
    d = np.diff(x); return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.median(x)


items = {}
for lab, stem in DS.items():
    S = [s for s in load(f'{R2}/{stem}_data.mat')]; B = load(f'{R2}/{stem}_background.mat')
    keep = [i for i, s in enumerate(S) if len(s) >= 180]
    t = T[T.dataset == lab].reset_index(drop=True); assert len(t) == len(keep), (lab, len(t), len(keep))
    b = G[G.dataset == lab]; it = []
    for k, i in enumerate(keep):
        s, c = S[i], B[i]; tr = t.trace[k]; L = len(s); bb = b[b.trace == tr]
        e = dict(movie=t.movie[k], L=L, level=np.median(s), noise=rs(s), amps=bb[bb.kind == 'real'].amp.abs().values)
        for h, sel in [(0, slice(None)), (1, bb.startFrame <= CUT), (2, bb.startFrame > CUT)]:
            sb = bb if h == 0 else bb[sel]
            e[f'r{h}'] = sb[sb.kind == 'real'].f.values; e[f'c{h}'] = sb[sb.kind == 'surr'].f.values
        e['h0'], e['h1'], e['h2'] = L * H, min(L, CUT) * H, max(L - CUT, 0) * H
        y = (s - c) / c
        e['n0'] = slope(y); e['n1'] = slope(y[:CUT]); e['n2'] = slope(y[CUT:]) if L - CUT >= 60 else np.nan
        it.append(e)
    items[lab] = it


def stats(it):
    o = {}
    cat = lambda key: np.concatenate([x[key] for x in it]) if it else np.zeros(0)
    o['brightness'] = np.median([x['level'] for x in it]); o['noise'] = 100 * np.median([x['noise'] for x in it])
    a = cat('amps'); o['median_burst'] = np.median(a) if len(a) else np.nan; o['mean_burst'] = a.mean() if len(a) else np.nan
    r, c = cat('r0'), cat('c0')
    o['size_real'] = 100 * np.median(r) if len(r) else np.nan; o['size_copy'] = 100 * np.median(c) if len(c) else np.nan
    o['beyond'] = o['size_real'] - o['size_copy']
    for h in (0, 1, 2):
        hrs = sum(x[f'h{h}'] for x in it); rr, cc = cat(f'r{h}'), cat(f'c{h}')
        o[f'turn{h}'] = 100 * (rr.sum() - cc.sum()) / hrs if hrs else np.nan
        o[f'net{h}'] = np.nanmean([x[f'n{h}'] for x in it])
    o['noise_share'] = c.sum() / r.sum() if r.sum() else np.nan
    return o


def by_movie(it):
    g = {}
    for x in it:
        g.setdefault(x['movie'], []).append(x)
    return g


KEYS = ['brightness', 'noise', 'median_burst', 'mean_burst', 'size_real', 'size_copy', 'beyond', 'turn0', 'net0', 'turn1', 'turn2', 'net1', 'net2', 'noise_share']
rows, boots = [], {}
for lab, it in items.items():
    obs = stats(it); g = by_movie(it); keys = list(g); bs = []
    for _ in range(NB):
        res = []
        for k in rng.choice(len(keys), len(keys)):
            m = g[keys[k]]; res.extend(m[i] for i in rng.integers(0, len(m), len(m)))
        bs.append(stats(res))
    bs = pd.DataFrame(bs); boots[lab] = bs
    r = dict(dataset=lab, label=LAB[lab], granules=len(it), movies=len(keys))
    for k in KEYS:
        r[k] = obs[k]; r[k + '_lo'], r[k + '_hi'] = np.nanpercentile(bs[k], [2.5, 97.5])
    r['halves_turn_p'] = 2 * min((bs.turn2 - bs.turn1 <= 0).mean(), (bs.turn2 - bs.turn1 >= 0).mean())
    r['halves_net_p'] = 2 * min((bs.net2 - bs.net1 <= 0).mean(), (bs.net2 - bs.net1 >= 0).mean())
    rows.append(r); print(lab, 'done', flush=True)
S = pd.DataFrame(rows); S.to_csv(os.path.join(out, 'summary.csv'), index=False)

tests = []
for a, b in [('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1'), ('extra_100ngul_1mgml', 'protein_RNA_1_1'), ('extra_100ngul_1mgml', 'protein_RNA_100_1')]:
    ga, gb = by_movie(items[a]), by_movie(items[b]); units = list(ga.values()) + list(gb.values()); na = len(ga)
    oa, ob = stats(items[a]), stats(items[b])
    for k in ['beyond', 'turn0', 'net0', 'median_burst']:
        obs = oa[k] - ob[k]; null = np.empty(NP)
        for i in range(NP):
            p = rng.permutation(len(units))
            null[i] = stats([x for j in p[:na] for x in units[j]])[k] - stats([x for j in p[na:] for x in units[j]])[k]
        tests.append(dict(compare=f'{LAB[a]} vs {LAB[b]}', measure=k, A=oa[k], B=ob[k], diff=obs,
                          p_movie_perm=(np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (NP + 1)))
    print(a, 'vs', b, 'done', flush=True)
Tt = pd.DataFrame(tests); Tt.to_csv(os.path.join(out, 'tests.csv'), index=False)

ci = lambda r, k, f='.1f': f"{r[k]:{f}} [{r[k + '_lo']:{f}}, {r[k + '_hi']:{f}}]"
with open(os.path.join(out, 'summary.md'), 'w', encoding='utf-8', newline=chr(13) + chr(10)) as fh:
    fh.write('| condition | granules (movies) | brightness | frame noise (%) | median burst (a.u.) | burst size real / copies (%) | beyond noise (pts) | turnover beyond noise (%/h) | net change (%/h) | share of turnover from noise |' + chr(10) + '|---|---|---|---|---|---|---|---|---|---|' + chr(10))
    for _, r in S.iterrows():
        fh.write(f"| {r.label} | {r.granules} ({r.movies}) | {r.brightness:.0f} | {r.noise:.1f} | {r.median_burst:.0f} | {r.size_real:.1f} / {r.size_copy:.1f} | {ci(r, 'beyond')} | {ci(r, 'turn0')} | {ci(r, 'net0')} | {r.noise_share:.2f} |" + chr(10))
    fh.write(chr(10) + '| condition | turnover, first → second half (%/h) | p | net change, first → second half (%/h) | p |' + chr(10) + '|---|---|---|---|---|' + chr(10))
    for _, r in S.iterrows():
        fh.write(f"| {r.label} | {r.turn1:.1f} → {r.turn2:.1f} | {r.halves_turn_p:.3f} | {r.net1:.1f} → {r.net2:.1f} | {r.halves_net_p:.3f} |" + chr(10))
    fh.write(chr(10) + '| comparison | measure | A | B | p (movies shuffled) |' + chr(10) + '|---|---|---|---|---|' + chr(10))
    for _, r in Tt.iterrows():
        fh.write(f"| {r['compare']} | {r.measure} | {r.A:.1f} | {r.B:.1f} | {r.p_movie_perm:.4f} |" + chr(10))
print(open(os.path.join(out, 'summary.md'), encoding='utf-8').read())

fig, axs = plt.subplots(1, 4, figsize=(22, 5))
x = np.arange(len(S)); cols = ['0.75', '0.5', '0.25', 'tab:purple']
for ax, (k, t_) in zip(axs, [('beyond', 'burst size beyond noise (% points)'), ('turn0', 'turnover beyond noise (%/h)'), ('net0', 'net change beyond bleaching (%/h)'), ('median_burst', 'median burst amplitude (a.u.)')]):
    ax.bar(x, S[k], color=cols, edgecolor='k', lw=0.5)
    ax.errorbar(x, S[k], yerr=[S[k] - S[k + '_lo'], S[k + '_hi'] - S[k]], fmt='none', ecolor='k', capsize=4)
    ax.axhline(0, color='k', lw=0.7); ax.set_xticks(x); ax.set_xticklabels([l.split(' (')[0] for l in S.label], fontsize=9); ax.set_title(t_, fontsize=10, loc='left')
fig.suptitle('In vitro, re-extracted from the titration movies (95% CI over movies, then granules)', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, 'vitro_new.png'), dpi=130); fig.savefig(os.path.join(out, 'vitro_new.pdf'))
