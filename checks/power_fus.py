"""POWER_FUS  How many granules / movies would make the 4x PP7 vs PP7-FUS
difference (and the other in vivo comparisons) significant, assuming the same
effect and the same spread? RESULTS.md section 20; reviewer point R2 5.
Scores (spot sum, in vivo):
  M2b  noise-free fluctuation G(tau) = 100*sqrt(mean_traces(D(tau) - D(1 frame)))
       (checks/noise_methods_real.m), tau = 30 frames main; also 6..60 frames
  M1   excess turnover (%/h) = 100*(sum frac_real - sum frac_surr)/hours, shared-
       bleaching surrogate (robustness/surrbleach), SP 12, window 13
Standard errors: granule bootstrap (traces resampled) and nested bootstrap
(movies resampled, then traces within each movie), 1000x. If the nested SE is
larger, movies add their own variability and more MOVIES (not only more
granules per movie) are needed.
Sample size: SE scales as 1/sqrt(k) when every group is k times larger, so
z_k = z * sqrt(k). k needed for z = 1.96 (p = 0.05, 50% power) and z = 2.80
(80% power at p = 0.05).
Usage: python checks/power_fus.py   Writes analysis_outputs/power/power.csv, lags.csv."""
import os
import numpy as np
import pandas as pd
from scipy import stats

out = os.path.join('analysis_outputs', 'power'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(3)
TAUS = [1, 2, 3, 6, 12, 18, 30, 45, 60]
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7')]
H = 10 / 3600


def parse(s):
    return np.array([float(v) for v in s.strip('[]').split()])


SF = pd.read_csv(os.path.join('analysis_outputs', 'noise_methods', 'real_sf.csv'))
SF = SF[SF.measure == 'spotsum']
D = {ds: dict(D=np.array([parse(x) for x in g.Dreal]), mov=g.movie.values) for ds, g in SF.groupby('dataset')}
rob = os.path.join('analysis_outputs', 'robustness')
Bur = pd.read_csv(os.path.join(rob, 'base_bursts.csv')); Bsu = pd.read_csv(os.path.join(rob, 'surrbleach_bursts.csv'))
Tr = pd.read_csv(os.path.join(rob, 'base_traces.csv'))


def frac_by_trace(B, ds):
    b = B[(B.dataset == ds) & (B.measure == 'spotsum') & (B.SP == 12) & (B.sign < 0) & (B.levelBefore > 0)]
    return {t: (g.amp.abs() / g.levelBefore).values for t, g in b.groupby('trace')}


M1 = {}
for ds in D:
    t = Tr[(Tr.dataset == ds) & (Tr.measure == 'spotsum')]
    M1[ds] = dict(traces=t.trace.values, mov=t.movie.values, L=dict(zip(t.trace, t.lengthFrames)),
                  R=frac_by_trace(Bur, ds), S=frac_by_trace(Bsu, ds))


def m2b(Dm, k):
    return 100 * np.sqrt(max(np.mean(Dm[:, k] - Dm[:, 0]), 0))


def m1(x, pick):
    real = sum(x['R'].get(t, np.zeros(0)).sum() for t in pick); surr = sum(x['S'].get(t, np.zeros(0)).sum() for t in pick)
    return 100 * (real - surr) / (sum(x['L'][t] for t in pick) * H)


def boot(stat, idx, mov, nested, nb=1000):
    vals = []
    movs = {m: idx[mov == m] for m in np.unique(mov)}
    for _ in range(nb):
        if nested and len(movs) > 1:
            pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        else:
            pick = rng.choice(idx, len(idx))
        vals.append(stat(pick))
    return np.std(vals)


rows = []
for score in ['M2b_G30', 'M1_excess']:
    est, seg, sen = {}, {}, {}
    for ds in D:
        if score == 'M2b_G30':
            Dm = D[ds]['D']; idx = np.arange(len(Dm)); mov = D[ds]['mov']
            f = lambda p, Dm=Dm: m2b(Dm[p], 6)
            est[ds] = m2b(Dm, 6)
        else:
            x = M1[ds]; idx = x['traces']; mov = x['mov']
            f = lambda p, x=x: m1(x, p)
            est[ds] = m1(x, idx)
        seg[ds] = boot(f, idx, mov, False); sen[ds] = boot(f, idx, mov, True)
    for a, b in COMPS:
        diff = est[a] - est[b]
        for kind, se in [('granule', seg), ('nested (movies)', sen)]:
            z = diff / np.hypot(se[a], se[b])
            k50 = (1.96 / z) ** 2 if z else np.inf; k80 = (2.80 / z) ** 2 if z else np.inf
            na, nb_ = (len(D[a]['D']), len(D[b]['D'])) if score == 'M2b_G30' else (len(M1[a]['traces']), len(M1[b]['traces']))
            ma, mb = len(np.unique(D[a]['mov'])), len(np.unique(D[b]['mov']))
            rows.append(dict(score=score, compare=f'{a} vs {b}', se=kind, valueA=est[a], valueB=est[b], z=z,
                             power_now=stats.norm.cdf(abs(z) - 1.96), k_for_p05=k50, k_for_80power=k80,
                             granules_now=f'{na}/{nb_}', movies_now=f'{ma}/{mb}',
                             granules_80power=f'{int(np.ceil(na * k80))}/{int(np.ceil(nb_ * k80))}',
                             movies_80power=f'{int(np.ceil(ma * k80))}/{int(np.ceil(mb * k80))}',
                             nestedSE_over_granuleSE=np.hypot(sen[a], sen[b]) / np.hypot(seg[a], seg[b])))
P = pd.DataFrame(rows); P.to_csv(os.path.join(out, 'power.csv'), index=False)
pd.set_option('display.width', 260)
print(P.round(2).to_string(index=False))

lag = []
for a, b in COMPS:
    for k, tau in enumerate(TAUS):
        if tau == 1:
            continue
        ga, gb = m2b(D[a]['D'], k), m2b(D[b]['D'], k)
        sa = boot(lambda p: m2b(D[a]['D'][p], k), np.arange(len(D[a]['D'])), D[a]['mov'], True, 500)
        sb = boot(lambda p: m2b(D[b]['D'][p], k), np.arange(len(D[b]['D'])), D[b]['mov'], True, 500)
        lag.append(dict(compare=f'{a} vs {b}', lagMin=tau * 10 / 60, valueA=ga, valueB=gb, z_nested=(ga - gb) / np.hypot(sa, sb)))
Lg = pd.DataFrame(lag); Lg.to_csv(os.path.join(out, 'lags.csv'), index=False)
print('\nnoise-free fluctuation G by lag (nested bootstrap z):')
print(Lg.assign(v=Lg.valueA.round(1).astype(str) + ' vs ' + Lg.valueB.round(1).astype(str) + ' z ' + Lg.z_nested.round(1).astype(str))
      .pivot(index='lagMin', columns='compare', values='v').to_string())
