"""SCORE_COMPARE  Candidate summary scores for distinguishing systems, computed
from the burst tables of checks/score_export.m (negative bursts unless stated).
Reviewer points R2 4a/4b/7 (interpretable, robust parameters with error bars).
Scores (frac = |amplitude| / granule level just before the burst):
  f_med      median frac
  f_mean     mean frac, top 10% of bursts removed (as the original method)
  A_mean     mean |amplitude|, top 10% removed (~lambda*k0; absolute units)
  shape      gamma shape of frac, truncated-gamma MLE on the 15-90 percentile
             window (low = mostly small bursts + long tail; high = bursts
             cluster around a typical size)
  lam_mom    moment 'number of units': mean^2 / (var - sigma^2) of |amplitude|
             (top 10% removed; sigma = quiescent noise); NaN if var <= sigma^2
  rate       negative bursts per granule per hour (traces >= 180 frames, 10 s)
  turnover   rate * f_mean = % of the granule shed per hour
  net        positive rate * positive f_mean - turnover (% per hour; > 0 growing)
Uncertainty: granule bootstrap (traces resampled within a dataset), 300x, SP 12.
Stability: CV of each score across SP 10..15; agreement between the mean and
spot-sum measures (in vivo). Discrimination: z = difference / sqrt(seA^2 + seB^2)
for the comparisons of interest. Simulations: Spearman correlation of each score
with the true lambda within each k0 x noise setting.
Usage (repo root): python checks/score_compare.py
Writes analysis_outputs/scores/{scores.csv, stability.csv, discrimination.csv,
sim_tracking.csv, scores.png, scores.pdf}."""
import os
import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'scores')
rng = np.random.default_rng(1)
FPS_H = 10 / 3600  # hours per frame


def trimmed_mean(x):
    x = np.asarray(x, float)
    return np.mean(x[x <= np.percentile(x, 90)]) if len(x) else np.nan


def gamma_shape(x):
    x = np.asarray(x, float); x = x[x > 0]
    if len(x) < 10:
        return np.nan
    lo, hi = np.percentile(x, [15, 90]); w = x[(x >= lo) & (x <= hi)]

    def nll(p):
        k, th = np.exp(p)
        z = stats.gamma.cdf(hi, k, scale=th) - stats.gamma.cdf(lo, k, scale=th)
        if not z > 0:
            return 1e12
        return -np.sum(stats.gamma.logpdf(w, k, scale=th)) + len(w) * np.log(z)
    k0, _, th0 = stats.gamma.fit(x, floc=0)
    best = None
    for kk in [k0, 1, 3, 8]:
        r = optimize.minimize(nll, np.log([max(kk, 0.05), np.mean(x) / max(kk, 0.05)]), method='Nelder-Mead',
                              options=dict(xatol=1e-4, fatol=1e-6, maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    return float(np.exp(best.x[0]))


def scores(neg, pos, hours, sigma):
    """neg/pos: DataFrames of bursts with columns frac, absamp."""
    out = {}
    out['f_med'] = np.median(neg.frac) if len(neg) else np.nan
    out['f_mean'] = trimmed_mean(neg.frac)
    out['A_mean'] = trimmed_mean(neg.absamp)
    out['shape'] = gamma_shape(neg.frac)
    a = neg.absamp.values; a = a[a <= np.percentile(a, 90)] if len(a) else a
    v = np.var(a) - sigma ** 2 if len(a) else np.nan
    out['lam_mom'] = np.mean(a) ** 2 / v if len(a) and v > 0 else np.nan
    out['rate'] = len(neg) / hours
    out['turnover'] = 100 * out['rate'] * out['f_mean']
    out['net'] = 100 * (len(pos) / hours) * trimmed_mean(pos.frac) - out['turnover']
    return out


NAMES = ['f_med', 'f_mean', 'A_mean', 'shape', 'lam_mom', 'rate', 'turnover', 'net']


def prep(B):
    B = B[B.levelBefore > 0].copy()
    B['frac'] = B.amp.abs() / B.levelBefore
    B['absamp'] = B.amp.abs()
    return B


def dataset_scores(B, T, sigma, boot=0):
    hours = T.lengthFrames.sum() * FPS_H
    neg, pos = B[B.sign < 0], B[B.sign > 0]
    point = scores(neg, pos, hours, sigma)
    if not boot:
        return point, None
    tr = T.trace.values; L = dict(zip(T.trace, T.lengthFrames))
    gN = {t: g for t, g in neg.groupby('trace')}; gP = {t: g for t, g in pos.groupby('trace')}
    empty = neg.iloc[:0]
    bs = []
    for _ in range(boot):
        pick = rng.choice(tr, len(tr))
        n = pd.concat([gN.get(t, empty) for t in pick]); p = pd.concat([gP.get(t, empty) for t in pick])
        bs.append(scores(n, p, sum(L[t] for t in pick) * FPS_H, sigma))
    return point, pd.DataFrame(bs)


if __name__ == '__main__':
    Bd = prep(pd.read_csv(os.path.join(d, 'data_bursts.csv')))
    Td = pd.read_csv(os.path.join(d, 'data_traces.csv'))
    Sd = pd.read_csv(os.path.join(d, 'data_sigma.csv'))
    rows, bootsd = [], {}
    for (ds, ms), Tg in Td.groupby(['dataset', 'measure']):
        for sp in sorted(Bd.SP.unique()):
            sig = Sd[(Sd.dataset == ds) & (Sd.measure == ms) & (Sd.SP == sp)].sigma.iloc[0]
            b = Bd[(Bd.dataset == ds) & (Bd.measure == ms) & (Bd.SP == sp)]
            point, bs = dataset_scores(b, Tg, sig, boot=300 if sp == 12 else 0)
            r = dict(dataset=ds, measure=ms, SP=sp, **point)
            if bs is not None:
                for k in NAMES:
                    r[k + '_lo'], r[k + '_hi'] = np.nanpercentile(bs[k], [2.5, 97.5])
                    r[k + '_se'] = np.nanstd(bs[k])
                bootsd[(ds, ms)] = bs
            rows.append(r)
        print(ds, ms, 'done', flush=True)
    R = pd.DataFrame(rows); R.to_csv(os.path.join(d, 'scores.csv'), index=False)
    pd.set_option('display.width', 250)
    R12 = R[R.SP == 12]
    for ms in ['spotsum', 'mean']:
        t = R12[R12.measure == ms].set_index('dataset')
        print(f'\n== SP 12, measure {ms}: score [95% CI]')
        print(pd.DataFrame({k: [f"{t.loc[i, k]:.3g} [{t.loc[i, k + '_lo']:.3g}-{t.loc[i, k + '_hi']:.3g}]" for i in t.index]
                            for k in NAMES}, index=t.index).to_string())
    # stability across SP and agreement between measures
    st = []
    for k in NAMES:
        cv = R.groupby(['dataset', 'measure'])[k].agg(lambda x: np.nanstd(x) / abs(np.nanmean(x)) if np.nanmean(x) else np.nan)
        a = R12[R12.measure == 'mean'].set_index('dataset')[k]; b = R12[R12.measure == 'spotsum'].set_index('dataset')[k]
        common = b.index.intersection(a.index)
        rho = stats.spearmanr(a[common], b[common]).statistic if len(common) > 2 else np.nan
        relwid = ((R12[k + '_hi'] - R12[k + '_lo']) / R12[k].abs()).median()
        st.append(dict(score=k, medianCV_acrossSP=cv.median(), rankAgreement_mean_vs_spotsum=rho, medianRelCIwidth=relwid))
    S = pd.DataFrame(st); S.to_csv(os.path.join(d, 'stability.csv'), index=False)
    print('\n== stability / precision'); print(S.round(3).to_string(index=False))
    # discrimination
    comps = [('4PP7-FUS', '4PP7', 'spotsum'), ('30PP7-FUS', '30PP7', 'spotsum'), ('5Qb', '4PP7', 'spotsum'),
             ('10Qb', '30PP7', 'spotsum'), ('4PP7-FUS', '4PP7', 'mean'), ('5Qb', '4PP7', 'mean'),
             ('protein_RNA_10_1', 'protein_RNA_1_1', 'mean'), ('protein_RNA_100_1', 'protein_RNA_10_1', 'mean'),
             ('protein_RNA_100_1', 'protein_RNA_1_1', 'mean')]
    dc = []
    for a, b, ms in comps:
        ra = R12[(R12.dataset == a) & (R12.measure == ms)].iloc[0]; rb = R12[(R12.dataset == b) & (R12.measure == ms)].iloc[0]
        row = dict(compare=f'{a} vs {b}', measure=ms)
        for k in NAMES:
            row[k] = (ra[k] - rb[k]) / np.sqrt(ra[k + '_se'] ** 2 + rb[k + '_se'] ** 2)
        dc.append(row)
    D = pd.DataFrame(dc); D.to_csv(os.path.join(d, 'discrimination.csv'), index=False)
    print('\n== discrimination (z = difference / combined bootstrap SE; |z| > 2 ~ p < 0.05)')
    print(D.round(1).to_string(index=False))
    # simulations
    if os.path.exists(os.path.join(d, 'sim_bursts.csv')):
        Bs = prep(pd.read_csv(os.path.join(d, 'sim_bursts.csv'))); Ts = pd.read_csv(os.path.join(d, 'sim_traces.csv'))
        Ss = pd.read_csv(os.path.join(d, 'sim_sigma.csv'))
        sr = []
        for ds, Tg in Ts.groupby('dataset'):
            p = dict(x.split('l')[-1:] and (x[0], x[1:]) for x in ds.split('_')[1:])
            sig = Ss[Ss.dataset == ds].sigma.iloc[0]
            pt, _ = dataset_scores(Bs[Bs.dataset == ds], Tg, sig)
            sr.append(dict(lamTrue=float(p['l']), k0=int(p['k']), noise=int(p['n']), seed=int(p['s']), **pt))
        Q = pd.DataFrame(sr)
        tk = []
        for (k0, nz), g in Q.groupby(['k0', 'noise']):
            tk.append(dict(k0=k0, noise=nz, **{k: stats.spearmanr(g.lamTrue, g[k]).statistic for k in NAMES}))
        TK = pd.DataFrame(tk); TK.to_csv(os.path.join(d, 'sim_tracking.csv'), index=False)
        Q.to_csv(os.path.join(d, 'sim_scores.csv'), index=False)
        print('\n== simulations: Spearman correlation of each score with the true lambda')
        print(TK.round(2).to_string(index=False))
        print('\n== simulations: lam_mom and shape vs true lambda (median over seeds)')
        print(Q.groupby(['k0', 'noise', 'lamTrue'])[['lam_mom', 'shape', 'f_mean', 'A_mean']].median().round(2).to_string())
    # figure: every score per dataset at SP 12 (spot sum in vivo, mean in vitro)
    order = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
    COL = ['#d62a2a', '#2a5bd6', '#f08a8a', '#8aa8f0', '#1b8a4a', '#7fcf9f', '#7a4fd6', '#a98be6', '#d3c3f5']
    fig, axs = plt.subplots(2, 4, figsize=(18, 8))
    for ax, k in zip(axs.ravel(), NAMES):
        for i, ds in enumerate(order):
            ms = 'mean' if ds.startswith('protein') else 'spotsum'
            r = R12[(R12.dataset == ds) & (R12.measure == ms)]
            if r.empty or not np.isfinite(r[k].iloc[0]):
                continue
            r = r.iloc[0]
            ax.errorbar(i, r[k], yerr=[[max(r[k] - r[k + '_lo'], 0)], [max(r[k + '_hi'] - r[k], 0)]], fmt='o', color=COL[i], capsize=3)
            m2 = R12[(R12.dataset == ds) & (R12.measure == 'mean')]
            if ms == 'spotsum' and not m2.empty:
                ax.plot(i + 0.25, m2[k].iloc[0], 'x', color=COL[i], ms=6)
        ax.set_xticks(range(len(order))); ax.set_xticklabels([o.replace('protein_RNA_', 'iv ') for o in order], rotation=45, fontsize=7)
        ax.set_title(k, fontsize=10, loc='left')
        ax.tick_params(labelsize=7)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    fig.suptitle('Candidate scores per dataset, SP 12 (dots + 95% granule-bootstrap CI: spot sum in vivo, mean in vitro; x = mean measure in vivo)',
                 fontsize=10)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(d, f'scores.{ext}'), dpi=140)
