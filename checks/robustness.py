"""ROBUSTNESS  Robustness analyses of the chosen scores (checks/score_export2.m).
Scores (negative bursts; frac = |amplitude| / granule level before the burst):
  f_mean   mean fraction shed per burst, top 10% of bursts removed  (MAIN)
  f_med    median fraction shed per burst
  A_mean   mean absolute burst size, top 10% removed
  turnover % of the granule shed per hour = rate x f_mean
Steps (python checks/robustness.py <step>):
  noise     false-burst rate from the burst-free surrogates / real rate = false-
            discovery share per dataset x measure x SP; chooses one SP (lowest
            SP with share <= 10% for every dataset with the main measure).
            Reviewer point R2 6.
  sens      sensitivity table at the chosen SP: measure, top trim, SP +-1,
            minimum trace length, smoothing window; ranking agreement and key
            comparisons. R2 4a.
  nested    granule-only vs movie-then-granule bootstrap CIs, counts, first vs
            second 30 min. R2 5.
  window    scores vs smoothing window + synthetic duration / spacing tests. R2 2b.
Outputs in analysis_outputs/robustness/."""
import json
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'robustness')
rng = np.random.default_rng(1)
FPS_H = 10 / 3600
MAIN = {'in vivo': 'spotsum', 'in vitro': 'mean'}
ORDER = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7'),
         ('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1')]
SCORES = ['f_mean', 'f_med', 'A_mean', 'turnover']


def main_measure(ds):
    return MAIN['in vitro'] if ds.startswith('protein') else MAIN['in vivo']


def load(tag):
    B = pd.read_csv(os.path.join(d, f'{tag}_bursts.csv'))
    B = B[B.levelBefore > 0].copy(); B['frac'] = B.amp.abs() / B.levelBefore; B['absamp'] = B.amp.abs()
    T = pd.read_csv(os.path.join(d, f'{tag}_traces.csv'))
    return B, T


def tmean(x, top=10):
    x = np.asarray(x, float)
    if not len(x):
        return np.nan
    return np.mean(x[x <= np.percentile(x, 100 - top)]) if top > 0 else np.mean(x)


def score(neg, hours, top=10):
    fm = tmean(neg.frac, top)
    return dict(f_mean=fm, f_med=np.median(neg.frac) if len(neg) else np.nan, A_mean=tmean(neg.absamp, top),
                rate=len(neg) / hours, turnover=100 * fm * len(neg) / hours, n=len(neg))


def scores_with_boot(B, T, nb=200, top=10, nested=False, half=None):
    neg = B[B.sign < 0]
    if half == 1:
        neg = neg[neg.startFrame <= 180]
    elif half == 2:
        neg = neg[neg.startFrame > 180]
    L = dict(zip(T.trace, T.lengthFrames))
    hours_of = lambda tr: sum(min(L[t], 180) if half == 1 else (max(L[t] - 180, 0) if half == 2 else L[t]) for t in tr) * FPS_H
    pt = score(neg, hours_of(T.trace), top)
    g = {t: x for t, x in neg.groupby('trace')}; empty = neg.iloc[:0]
    mv = T.groupby('movie').trace.apply(list).to_dict()
    bs = []
    for _ in range(nb):
        if nested and len(mv) > 1:
            pick = [t for m in rng.choice(list(mv), len(mv)) for t in rng.choice(mv[m], len(mv[m]))]
        else:
            pick = list(rng.choice(T.trace.values, len(T)))
        bs.append(score(pd.concat([g.get(t, empty) for t in pick]), hours_of(pick), top))
    bs = pd.DataFrame(bs)
    for k in SCORES:
        pt[k + '_lo'], pt[k + '_hi'] = np.nanpercentile(bs[k], [2.5, 97.5]); pt[k + '_se'] = np.nanstd(bs[k])
    return pt


def comps_table(R, label):
    rows = []
    for a, b in COMPS:
        ra, rb = R.loc[a], R.loc[b]
        row = dict(variant=label, compare=f'{a} vs {b}')
        for k in SCORES:
            row[k + '_z'] = (ra[k] - rb[k]) / np.sqrt(ra[k + '_se'] ** 2 + rb[k + '_se'] ** 2)
        rows.append(row)
    return rows


def step_noise():
    Br, Tr = load('base'); Bs, Ts = load('surrpoly')
    rows = []
    for (ds, ms), T in Tr.groupby(['dataset', 'measure']):
        hours = T.lengthFrames.sum() * FPS_H
        for sp in sorted(Br.SP.unique()):
            r = Br[(Br.dataset == ds) & (Br.measure == ms) & (Br.SP == sp) & (Br.sign < 0)]
            s = Bs[(Bs.dataset == ds) & (Bs.measure == ms) & (Bs.SP == sp) & (Bs.sign < 0)]
            rows.append(dict(dataset=ds, measure=ms, SP=sp, realRate=len(r) / hours, falseRate=len(s) / hours,
                             FDR=len(s) / max(len(r), 1), realFmed=np.median(r.frac) if len(r) else np.nan,
                             falseFmed=np.median(s.frac) if len(s) else np.nan))
    N = pd.DataFrame(rows); N.to_csv(os.path.join(d, 'noise_floor.csv'), index=False)
    main = N[N.apply(lambda r: r.measure == main_measure(r.dataset), axis=1)]
    worst = main.groupby('SP').FDR.max()
    ok = worst[worst <= 0.10]
    sp = int(ok.index.min()) if len(ok) else int(worst.idxmin())
    json.dump(dict(chosenSP=sp, rule='lowest SP with false-discovery share <= 10% in every dataset (main measure)',
                   worstFDRbySP={int(k): float(v) for k, v in worst.items()}), open(os.path.join(d, 'chosen_sp.json'), 'w'), indent=1)
    pd.set_option('display.width', 220)
    print(main.pivot(index='dataset', columns='SP', values='FDR').reindex(ORDER).round(3).to_string())
    print('\nworst FDR by SP:', worst.round(3).to_dict(), '-> chosen SP', sp)
    print('\nmedian fraction of false bursts vs real (SP chosen):')
    print(main[main.SP == sp][['dataset', 'realFmed', 'falseFmed', 'realRate', 'falseRate']].set_index('dataset').reindex(ORDER).round(3).to_string())
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for ds, g in main.groupby('dataset'):
        ax.plot(g.SP, 100 * g.FDR, 'o-', label=ds)
    ax.axhline(10, color='grey', ls='--', lw=1); ax.axvline(sp, color='k', lw=0.8)
    ax.set_xlabel('SP'); ax.set_ylabel('false bursts as % of detected bursts'); ax.legend(fontsize=7, frameon=False)
    ax.set_title(f'noise floor (burst-free surrogates); chosen SP = {sp}', fontsize=10, loc='left')
    fig.tight_layout(); fig.savefig(os.path.join(d, 'noise_floor.png'), dpi=140); fig.savefig(os.path.join(d, 'noise_floor.pdf'))


def per_dataset(B, T, sp, top=10, measure_override=None, nb=200):
    out = {}
    for ds in ORDER:
        ms = measure_override if (measure_override and not ds.startswith('protein')) else main_measure(ds)
        b = B[(B.dataset == ds) & (B.measure == ms) & (B.SP == sp)]; t = T[(T.dataset == ds) & (T.measure == ms)]
        if t.empty:
            continue
        out[ds] = scores_with_boot(b, t, nb=nb, top=top)
    return pd.DataFrame(out).T


def step_sens():
    sp = json.load(open(os.path.join(d, 'chosen_sp.json')))['chosenSP']
    B, T = load('base')
    variants = [('baseline', B, T, sp, 10, None), ('mean measure (in vivo)', B, T, sp, 10, 'mean'),
                ('top trim 0%', B, T, sp, 0, None), ('top trim 5%', B, T, sp, 5, None), ('top trim 20%', B, T, sp, 20, None),
                (f'SP {sp - 1}', B, T, sp - 1, 10, None), (f'SP {sp + 1}', B, T, sp + 1, 10, None)]
    for tag, lab in [('minlen120', 'min trace 120 frames'), ('minlen240', 'min trace 240 frames'),
                     ('win9', 'window 9'), ('win17', 'window 17')]:
        if os.path.exists(os.path.join(d, f'{tag}_bursts.csv')):
            Bv, Tv = load(tag); variants.append((lab, Bv, Tv, sp, 10, None))
    base = None; rows = []; zrows = []
    for lab, Bv, Tv, spv, top, mo in variants:
        R = per_dataset(Bv, Tv, spv, top, mo)
        if base is None:
            base = R
        r = dict(variant=lab)
        for k in SCORES:
            r[k + '_rho'] = stats.spearmanr(base[k], R.loc[base.index, k]).statistic
            r[k + '_medRelChange'] = np.median(np.abs(R.loc[base.index, k] / base[k] - 1))
        rows.append(r); zrows += comps_table(R, lab)
        R.assign(variant=lab).to_csv(os.path.join(d, f'sens_{lab.replace(" ", "_").replace("%", "pct")}.csv'))
    S = pd.DataFrame(rows); Z = pd.DataFrame(zrows)
    S.to_csv(os.path.join(d, 'sensitivity_summary.csv'), index=False); Z.to_csv(os.path.join(d, 'sensitivity_comparisons.csv'), index=False)
    pd.set_option('display.width', 250)
    print(S.round(2).to_string(index=False)); print()
    print(Z.pivot(index='variant', columns='compare', values='f_mean_z').reindex([v[0] for v in variants]).round(1).to_string())


def step_nested():
    sp = json.load(open(os.path.join(d, 'chosen_sp.json')))['chosenSP']
    B, T = load('base'); rows = []; hrows = []
    for ds in ORDER:
        ms = main_measure(ds)
        b = B[(B.dataset == ds) & (B.measure == ms) & (B.SP == sp)]; t = T[(T.dataset == ds) & (T.measure == ms)]
        flat = scores_with_boot(b, t, nb=300); nest = scores_with_boot(b, t, nb=300, nested=True)
        nmov = t.movie.nunique() if t.movie.max() > 0 else np.nan
        rows.append(dict(dataset=ds, movies=nmov, traces=len(t), negBursts=flat['n'],
                         **{f'{k}': flat[k] for k in SCORES},
                         **{f'{k}_CI_granule': f"{flat[k + '_lo']:.3g}-{flat[k + '_hi']:.3g}" for k in SCORES},
                         **{f'{k}_CI_nested': f"{nest[k + '_lo']:.3g}-{nest[k + '_hi']:.3g}" for k in SCORES},
                         f_mean_se_ratio=nest['f_mean_se'] / flat['f_mean_se']))
        h1 = scores_with_boot(b, t, nb=300, half=1); h2 = scores_with_boot(b, t, nb=300, half=2)
        hrows.append(dict(dataset=ds, f_mean_first=h1['f_mean'], f_mean_second=h2['f_mean'],
                          z_f_mean=(h2['f_mean'] - h1['f_mean']) / np.hypot(h1['f_mean_se'], h2['f_mean_se']),
                          turnover_first=h1['turnover'], turnover_second=h2['turnover'],
                          z_turnover=(h2['turnover'] - h1['turnover']) / np.hypot(h1['turnover_se'], h2['turnover_se']),
                          bursts_first=h1['n'], bursts_second=h2['n']))
    N = pd.DataFrame(rows); H = pd.DataFrame(hrows)
    N.to_csv(os.path.join(d, 'nested_bootstrap.csv'), index=False); H.to_csv(os.path.join(d, 'time_halves.csv'), index=False)
    pd.set_option('display.width', 250)
    print(N[['dataset', 'movies', 'traces', 'negBursts', 'f_mean', 'f_mean_CI_granule', 'f_mean_CI_nested', 'f_mean_se_ratio',
             'turnover', 'turnover_CI_granule', 'turnover_CI_nested']].round(3).to_string(index=False)); print()
    print(H.round(3).to_string(index=False))


def step_window():
    sp = json.load(open(os.path.join(d, 'chosen_sp.json')))['chosenSP']
    rows = []
    for w in range(7, 23, 2):
        tag = 'base' if w == 13 else f'win{w}'
        if not os.path.exists(os.path.join(d, f'{tag}_bursts.csv')):
            continue
        Bv, Tv = load(tag); R = per_dataset(Bv, Tv, sp, nb=100)
        for ds in R.index:
            rows.append(dict(window=w, dataset=ds, **{k: R.loc[ds, k] for k in SCORES + ['rate', 'n']}))
    W = pd.DataFrame(rows); W.to_csv(os.path.join(d, 'window_scores.csv'), index=False)
    pd.set_option('display.width', 250)
    for k in ['f_mean', 'turnover', 'rate']:
        print(f'\n{k} vs window:'); print(W.pivot(index='dataset', columns='window', values=k).reindex(ORDER).round(3).to_string())
    if os.path.exists(os.path.join(d, 'window_duration.csv')):
        A = pd.read_csv(os.path.join(d, 'window_duration.csv')); Sp = pd.read_csv(os.path.join(d, 'window_spacing.csv'))
        print('\nP(detected) by event duration (rows) and window (cols), amplitude 150 / 300:')
        for amp in sorted(A.amp.unique()):
            print(f'amp {amp}:'); print(A[A.amp == amp].pivot(index='durationFrames', columns='window', values='pDetected').round(2).to_string())
        print('\nP(two drops resolved) by gap (rows) and window (cols):')
        print(Sp.pivot(index='gapFrames', columns='window', values='pResolved').round(2).to_string())
        fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
        for amp, ax in zip(sorted(A.amp.unique()), axs[:2]):
            for w, g in A[A.amp == amp].groupby('window'):
                ax.plot(g.durationFrames, g.pDetected, 'o-', label=f'window {w}')
            ax.set_xlabel('event duration (frames, 10 s each)'); ax.set_ylabel('P(detected)')
            ax.set_title(f'single drop of {amp} (noise SD 40), SP {sp}', fontsize=9, loc='left')
        axs[0].legend(fontsize=7, frameon=False)
        for w, g in Sp.groupby('window'):
            axs[2].plot(g.gapFrames, g.pResolved, 'o-', label=f'window {w}')
        axs[2].set_xlabel('gap between two drops (frames)'); axs[2].set_ylabel('P(resolved as 2 bursts)')
        axs[2].set_title('closely spaced bursts', fontsize=9, loc='left')
        fig.tight_layout(); fig.savefig(os.path.join(d, 'window_tests.png'), dpi=140); fig.savefig(os.path.join(d, 'window_tests.pdf'))


if __name__ == '__main__':
    {'noise': step_noise, 'sens': step_sens, 'nested': step_nested, 'window': step_window}[sys.argv[1]]()
