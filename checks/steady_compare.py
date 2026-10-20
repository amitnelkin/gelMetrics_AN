"""STEADY_COMPARE  Per-frame (R2, published) vs steady multithresh thresholds
(R4, checks/reextract_steady.m) on the in vivo data. RESULTS.md section 23.
1. Frame-to-frame jitter (robust, relative) of spot mean, spot sum, cell mean,
   spot area, per dataset.
2. Burst detection on real traces and burst-free surrogates (checks/detection_grid.m
   run on each version; SP 8/10/12 x window 13/17): false-discovery share and
   excess turnover (%/h, real minus surrogate fraction shed) with nested
   (movie -> granule) bootstrap SE, 1000x; comparison z.
3. Noise-free fluctuation G30 (M2b, checks/noise_methods_real.m on each version),
   nested bootstrap z and the data multiple needed for 80% power.
Usage: python checks/steady_compare.py
Writes analysis_outputs/steady/{jitter.csv, detection.csv, comparisons.csv, m2b.csv, steady.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'steady'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(23)
H = 10 / 3600
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7')]
VER = {'per-frame (R2)': dict(r2='in_vivo_reextracted/after_R2_offset', su='in_vivo_reextracted/after_R2_sums',
                              grid='analysis_outputs/detection_grid', sf='analysis_outputs/noise_methods/real_sf.csv'),
       'steady (R4)': dict(r2='in_vivo_reextracted/after_R4_steadythr', su='in_vivo_reextracted/after_R4_steadythr',
                           grid='analysis_outputs/steady/grid', sf='analysis_outputs/steady/real_sf.csv')}


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def rs(x):
    d = np.diff(x)
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.median(x)


def nested(stat, idx, mov, nb=1000):
    movs = {m: idx[mov == m] for m in np.unique(mov)}; v = []
    for _ in range(nb):
        pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        v.append(stat(pick))
    return np.std(v)


# 1. jitter
J = []
for ver, f in VER.items():
    for ds in VIVO:
        S, B = load(f"{f['r2']}/{ds}_data.mat"), load(f"{f['r2']}/{ds}_background.mat")
        SS, A = load(f"{f['su']}/{ds}_spotsum_data.mat"), load(f"{f['su']}/{ds}_spotarea.mat")
        r = [(rs(s), rs(ss), rs(b), rs(a)) for s, b, ss, a in
             ((np.ravel(s).astype(float), np.ravel(b).astype(float), np.ravel(ss).astype(float), np.ravel(a).astype(float)) for s, b, ss, a in zip(S, B, SS, A))
             if len(s) >= 180 and np.all(b > 0)]
        r = np.median(np.array(r), axis=0)
        J.append(dict(version=ver, dataset=ds, spotMean=r[0], spotSum=r[1], cellMean=r[2], spotArea=r[3]))
J = pd.DataFrame(J); J.to_csv(os.path.join(out, 'jitter.csv'), index=False)
pd.set_option('display.width', 250)
print('relative frame-to-frame jitter (median over traces):')
print(J.pivot(index='dataset', columns='version', values=['spotSum', 'spotMean', 'cellMean', 'spotArea']).reindex(VIVO).round(3).to_string())

# 2. detection
D, C = [], []
for ver, f in VER.items():
    B = pd.read_csv(os.path.join(f['grid'], 'bursts.csv')); B = B[(B.levelBefore > 0) & B.dataset.isin(VIVO)].copy()
    B['frac'] = B.amp.abs() / B.levelBefore
    T = pd.read_csv(os.path.join(f['grid'], 'traces.csv'))
    for (w, sp), g in B.groupby(['window', 'SP']):
        if w not in (13, 17) or sp not in (8, 10, 12):
            continue
        est = {}
        for ds in VIVO:
            t = T[T.dataset == ds]; L = dict(zip(t.trace, t.lengthFrames)); b = g[g.dataset == ds]
            R = b[b.kind == 'real'].groupby('trace').frac.sum().to_dict(); Su = b[b.kind == 'surr'].groupby('trace').frac.sum().to_dict()
            exc = lambda p: 100 * (sum(R.get(k, 0) for k in p) - sum(Su.get(k, 0) for k in p)) / (sum(L[k] for k in p) * H)
            idx = t.trace.values; m, se = exc(idx), nested(exc, idx, t.movie.values)
            nr, ns = (b.kind == 'real').sum(), (b.kind == 'surr').sum(); est[ds] = (m, se)
            D.append(dict(version=ver, window=w, SP=sp, dataset=ds, nReal=nr, nFalse=ns, FDR=ns / max(nr, 1), excess=m, excess_se=se,
                          medianFracReal=b[b.kind == 'real'].frac.median()))
        for a, bb in COMPS:
            C.append(dict(version=ver, window=w, SP=sp, compare=f'{a} vs {bb}', A=est[a][0], B=est[bb][0],
                          z=(est[a][0] - est[bb][0]) / np.hypot(est[a][1], est[bb][1])))
D = pd.DataFrame(D); D.to_csv(os.path.join(out, 'detection.csv'), index=False)
C = pd.DataFrame(C); C.to_csv(os.path.join(out, 'comparisons.csv'), index=False)
print('\nfalse-discovery share (false / real negative bursts):')
print(D.pivot_table(index=['dataset'], columns=['version', 'window', 'SP'], values='FDR').reindex(VIVO).round(2).to_string())
print('\nreal negative bursts:')
print(D.pivot_table(index=['dataset'], columns=['version', 'window', 'SP'], values='nReal').reindex(VIVO).astype(int).to_string())
print('\nexcess turnover %/h:')
print(D.pivot_table(index=['dataset'], columns=['version', 'window', 'SP'], values='excess').reindex(VIVO).round(1).to_string())
print('\ncomparison z (excess turnover, nested):')
print(C.pivot_table(index='compare', columns=['version', 'window', 'SP'], values='z').round(1).to_string())

# 3. M2b
M = []
for ver, f in VER.items():
    SF = pd.read_csv(f['sf']); SF = SF[(SF.measure == 'spotsum') & SF.dataset.isin(VIVO)]
    est = {}
    for ds, g in SF.groupby('dataset'):
        Dm = np.array([[float(v) for v in x.strip('[]').split()] for x in g.Dreal]); mov = g.movie.values; idx = np.arange(len(Dm))
        st = lambda p, Dm=Dm: 100 * np.sqrt(max(np.mean(Dm[p, 6] - Dm[p, 0]), 0))
        est[ds] = (st(idx), nested(st, idx, mov)); M.append(dict(version=ver, dataset=ds, G30=est[ds][0], se=est[ds][1]))
    for a, b in COMPS:
        z = (est[a][0] - est[b][0]) / np.hypot(est[a][1], est[b][1])
        M.append(dict(version=ver, dataset=f'{a} vs {b}', G30=np.nan, se=np.nan, z=z, k80=(2.8 / z) ** 2 if z else np.inf))
M = pd.DataFrame(M); M.to_csv(os.path.join(out, 'm2b.csv'), index=False)
print('\nnoise-free fluctuation G30 (%), and comparison z / data multiple for 80% power:')
print(M.round(2).to_string(index=False))

fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
x = np.arange(len(VIVO)); wd = 0.38
for k, (ver, col) in enumerate(zip(VER, ['0.6', 'tab:green'])):
    j = J[J.version == ver].set_index('dataset').reindex(VIVO)
    axs[0].bar(x + (k - 0.5) * wd, 100 * j.spotSum, wd, color=col, label=ver)
    d = D[(D.version == ver) & (D.window == 13) & (D.SP == 12)].set_index('dataset').reindex(VIVO)
    axs[1].bar(x + (k - 0.5) * wd, d.FDR, wd, color=col, label=ver)
    axs[2].errorbar(x + (k - 0.5) * 0.3, d.excess, yerr=1.96 * d.excess_se, fmt='o', color=col, capsize=3, label=ver)
axs[0].set_ylabel('spot-sum jitter per frame (% of level)'); axs[1].set_ylabel('false / real negative bursts (window 13, SP 12)')
axs[1].axhline(1, color='k', lw=0.8, ls=':'); axs[2].axhline(0, color='k', lw=0.8); axs[2].set_ylabel('excess turnover beyond noise (%/h, 95% CI)')
for ax in axs:
    ax.set_xticks(x); ax.set_xticklabels(VIVO, rotation=30, fontsize=8); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out, 'steady.png'), dpi=140); fig.savefig(os.path.join(out, 'steady.pdf'))
