"""FUS_BRIGHTNESS  Why 4x PP7 and PP7-FUS look consistently different while 30x
do not: how much of the burst-size gap is set by granule brightness (noise
level) and how much is beyond noise. RESULTS.md section 25.
Data: in vivo R2 (published extraction), spot sum, window 13, SP 12, negative
bursts on the traces and on their burst-free copies (detection_grid.m);
per-granule brightness (spot mean, offset removed) and frame noise
(noise_from_background.m); fluctuation G30 (noise_methods_real.m); net change
(net_loss.py). Traces >= 180 frames.
p values: movie-level permutation (whole movies reassigned between the two
conditions, 10000x, seeded): movies are the replicate unit. For contrast, the
burst-level Mann-Whitney p (every burst treated as independent, as in a naive
histogram comparison) is also given for the raw burst sizes.
CIs: nested movie -> granule bootstrap, 2000x.
Usage: python checks/fus_brightness.py
Writes analysis_outputs/fus_brightness/{pvalues.csv, pvalues.md, fus_brightness.png/.pdf}."""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'fus_brightness'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(25)
H = 10 / 3600; NPERM = 10000; NB = 2000
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7')]
LAB = {'4PP7': 'PP7, 4x RNA', '4PP7-FUS': 'PP7-FUS, 4x RNA', '30PP7': 'PP7, 30x RNA', '30PP7-FUS': 'PP7-FUS, 30x RNA', '5Qb': 'Qβ, 5x', '10Qb': 'Qβ, 10x'}
COL = {'4PP7': 'tab:blue', '30PP7': 'tab:blue', '4PP7-FUS': 'tab:green', '30PP7-FUS': 'tab:green', '5Qb': 'tab:red', '10Qb': 'tab:red'}

# ---- per-granule tables, each with (dataset, movie)
B = pd.read_csv('analysis_outputs/detection_grid/bursts.csv')
B = B[(B.window == 13) & (B.SP == 12) & (B.levelBefore > 0)].copy(); B['f'] = B.amp.abs() / B.levelBefore
T = pd.read_csv('analysis_outputs/detection_grid/traces.csv')
B = B.merge(T[['dataset', 'trace', 'movie']], on=['dataset', 'trace'])
BR = pd.read_csv('analysis_outputs/noise_background/traces.csv'); BR = BR[BR.lengthFrames >= 180].copy()
BR['noise'] = BR.sdSum_meas / BR.spotSum
SF = pd.read_csv('analysis_outputs/noise_methods/real_sf.csv'); SF = SF[SF.measure == 'spotsum'].copy()
SF['g'] = [float(x.strip('[]').split()[6]) - float(x.strip('[]').split()[0]) for x in SF.Dreal]
NL = pd.read_csv('analysis_outputs/net_loss/per_trace.csv').dropna(subset=['sumNet'])


def items(key, ds):
    """{movie: list of granule items}; an item is a scalar, or for bursts
    (real burst sizes, copy burst sizes, trace length)"""
    if key in ('brightness', 'noise', 'G30', 'net'):
        df, col = {'brightness': (BR, 'spotMean'), 'noise': (BR, 'noise'), 'G30': (SF, 'g'), 'net': (NL, 'sumNet')}[key]
        d = df[df.dataset == ds]
        return {m: list(g[col].values) for m, g in d.groupby('movie')}
    t = T[T.dataset == ds]; b = B[B.dataset == ds]
    R = {k: g.f.values for k, g in b[b.kind == 'real'].groupby('trace')}; S_ = {k: g.f.values for k, g in b[b.kind == 'surr'].groupby('trace')}
    e = np.zeros(0)
    return {m: [(R.get(k, e), S_.get(k, e), L) for k, L in zip(g.trace, g.lengthFrames)] for m, g in t.groupby('movie')}


def stat(key, it):
    """it: flat list of granule items"""
    if key == 'brightness': return np.median(it)
    if key == 'noise': return 100 * np.median(it)
    if key == 'G30': return 100 * np.sqrt(max(np.mean(it), 0))
    if key == 'net': return np.mean(it)
    r = np.concatenate([x[0] for x in it]); c = np.concatenate([x[1] for x in it])
    if key == 'size_real': return 100 * np.median(r)
    if key == 'size_noise': return 100 * np.median(c)
    if key == 'size_beyond': return 100 * (np.median(r) - np.median(c))
    if key == 'turnover': return 100 * (r.sum() - c.sum()) / (sum(x[2] for x in it) * H)


def boot(key, G):
    keys = list(G); v = []
    for _ in range(NB):
        it = []
        for k in rng.choice(len(keys), len(keys)):
            g = G[keys[k]]; it.extend(g[i] for i in rng.integers(0, len(g), len(g)))
        v.append(stat(key, it))
    return np.percentile(v, [2.5, 97.5])


MEAS = [('brightness', 'granule brightness (spot mean)'), ('noise', 'frame noise (% of granule)'),
        ('size_real', 'median burst size, granules (%)'), ('size_noise', 'median burst size, burst-free copies (%)'),
        ('size_beyond', 'burst size beyond noise (granules - copies, % points)'), ('turnover', 'turnover beyond noise (%/h)'),
        ('G30', 'fluctuation over 5 min (%)'), ('net', 'net change beyond bleaching (%/h)')]
DS = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
rows, cis = [], {}
for key, lab in MEAS:
    G = {ds: items(key, ds) for ds in DS}
    for ds in DS:
        flat = [x for g in G[ds].values() for x in g]
        cis[(key, ds)] = (stat(key, flat), boot(key, G[ds]))
    for a, b in COMPS:
        movies = list(G[a].values()) + list(G[b].values()); na = len(G[a])
        obs = cis[(key, a)][0] - cis[(key, b)][0]; cnt = 0
        for _ in range(NPERM):
            p = rng.permutation(len(movies))
            d = stat(key, [x for i in p[:na] for x in movies[i]]) - stat(key, [x for i in p[na:] for x in movies[i]])
            cnt += abs(d) >= abs(obs) - 1e-12
        r = dict(measure=key, label=lab, compare=f'{a} vs {b}', A=cis[(key, a)][0], B=cis[(key, b)][0], diff=obs,
                 p_movie_perm=(cnt + 1) / (NPERM + 1), moviesA=na, moviesB=len(G[b]))
        if key == 'size_real':
            r['p_burst_MW_naive'] = mannwhitneyu(B[(B.dataset == a) & (B.kind == 'real')].f, B[(B.dataset == b) & (B.kind == 'real')].f).pvalue
        rows.append(r)
    print(key, 'done', flush=True)
P = pd.DataFrame(rows); P.to_csv(os.path.join(out, 'pvalues.csv'), index=False)
pd.set_option('display.width', 250)
print(P.round(4).to_string(index=False))
LF, CRLF = chr(10), chr(13) + chr(10)
with open(os.path.join(out, 'pvalues.md'), 'w', encoding='utf-8', newline=CRLF) as fh:
    fh.write('| measure | ' + ' | '.join(f'{a} vs {b}' for a, b in COMPS) + ' |' + LF + '|---|' + '---|' * len(COMPS) + LF)
    for key, lab in MEAS:
        cells = []
        for a, b in COMPS:
            r = P[(P.measure == key) & (P['compare'] == f'{a} vs {b}')].iloc[0]
            c = f"{r.A:.3g} vs {r.B:.3g}, p = {r.p_movie_perm:.2g}"
            if key == 'size_real':
                c += f" (naive burst-level p = {r.p_burst_MW_naive:.1g})"
            cells.append(c)
        fh.write(f'| {lab} | ' + ' | '.join(cells) + ' |' + LF)

# ---- figure
pv = lambda key, a, b: P[(P.measure == key) & (P['compare'] == f'{a} vs {b}')].p_movie_perm.iloc[0]
ptxt = lambda p: f'p = {p:.2g}' if p >= 0.001 else 'p < 0.001'
fig, axs = plt.subplots(2, 4, figsize=(19, 9.5))
for r_, (a, b) in enumerate(COMPS[:2]):
    pair = [b, a]
    # (1) brightness
    ax = axs[r_, 0]
    for k, ds in enumerate(pair):
        v = BR[BR.dataset == ds].spotMean.values
        ax.plot(k + rng.uniform(-0.2, 0.2, len(v)), v, '.', color=COL[ds], alpha=0.35, ms=4)
        m, c = cis[('brightness', ds)]; ax.errorbar(k, m, yerr=[[m - c[0]], [c[1] - m]], fmt='o', color='k', capsize=5, ms=7)
    ax.set_xticks([0, 1]); ax.set_xticklabels([LAB[d] for d in pair]); ax.set_ylabel('granule brightness (spot mean)')
    ax.set_title(f"(a) brightness: {ptxt(pv('brightness', a, b))}\nnoise per frame {cis[('noise', b)][0]:.0f}% vs {cis[('noise', a)][0]:.0f}% ({ptxt(pv('noise', a, b))})", fontsize=9, loc='left')
    ax.set_yscale('log')
    # (2) burst size real vs copies
    ax = axs[r_, 1]
    for k, ds in enumerate(pair):
        for j, (key, mk, lab) in enumerate([('size_real', 'o', 'granules'), ('size_noise', 's', 'burst-free copies (noise only)')]):
            m, c = cis[(key, ds)]
            ax.errorbar(k + (j - 0.5) * 0.25, m, yerr=[[m - c[0]], [c[1] - m]], fmt=mk, color=COL[ds], mfc=COL[ds] if j == 0 else 'w', ms=9, capsize=4,
                        label=lab if k == 0 else None)
    ax.set_xticks([0, 1]); ax.set_xticklabels([LAB[d] for d in pair]); ax.set_ylabel('median burst size (% of granule shed)')
    ax.set_xlim(-0.6, 1.6); ax.legend(fontsize=8, loc='upper left')
    ax.set_title(f"(b) burst size: granules {ptxt(pv('size_real', a, b))}\n(naive burst-level test: p = {P[(P.measure == 'size_real') & (P['compare'] == f'{a} vs {b}')].p_burst_MW_naive.iloc[0]:.1g})", fontsize=9, loc='left')
    # (3) decomposition of the gap
    ax = axs[r_, 2]
    tot = cis[('size_real', a)][0] - cis[('size_real', b)][0]; noi = cis[('size_noise', a)][0] - cis[('size_noise', b)][0]; bey = tot - noi
    ax.bar([0], [noi], color='0.7', edgecolor='k', label=f"set by brightness / noise ({ptxt(pv('size_noise', a, b))})")
    ax.bar([0], [bey], bottom=[noi] if np.sign(bey) == np.sign(noi) else [0], color=COL[a], edgecolor='k', label=f"beyond noise ({ptxt(pv('size_beyond', a, b))})")
    ax.axhline(0, color='k', lw=0.8); ax.axhline(tot, color='k', ls=':', lw=1)
    ax.text(0.33, tot, f'total gap {tot:+.1f} points', va='bottom', fontsize=8)
    ax.set_xlim(-0.8, 1.4); ax.set_xticks([]); ax.set_ylabel('FUS minus PP7, median burst size (% points)')
    ax.set_title('(c) where the gap comes from', fontsize=9, loc='left'); ax.legend(fontsize=8, loc='lower right')
    # (4) dynamics beyond noise
    ax = axs[r_, 3]
    for k, ds in enumerate(pair):
        m, c = cis[('turnover', ds)]; ax.errorbar(k, m, yerr=[[m - c[0]], [c[1] - m]], fmt='o', color=COL[ds], ms=9, capsize=5)
    ax.axhline(0, color='k', lw=0.8); ax.set_xticks([0, 1]); ax.set_xticklabels([LAB[d] for d in pair]); ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel('turnover beyond noise (%/h, 95% CI)')
    ax.set_title(f"(d) dynamics beyond noise: turnover {ptxt(pv('turnover', a, b))}\nfluctuation {ptxt(pv('G30', a, b))}, net change {ptxt(pv('net', a, b))}", fontsize=9, loc='left')
fig.suptitle('Rows: 4x RNA (top), 30x RNA (bottom). p: movie-level permutation (movies are the replicates)', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fus_brightness.png'), dpi=130); fig.savefig(os.path.join(out, 'fus_brightness.pdf'))
