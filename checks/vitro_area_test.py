"""VITRO_AREA_TEST  Does the growth of the in vitro granule area differ between the
protein:RNA ratios? (RESULTS.md section 55.) Area N = spot sum / (spot mean - dilute phase),
the number of granule pixels; per granule (>= 180 frames, same_granules_R2) the Theil-Sen
slope of ln(N) in %/h, as checks/vitro_area_trend.py. Movies are the replicates:
dataset value = mean over granules; 95% CI from resampling movies, then granules (5000x);
p from shuffling whole movies between two ratios (10000x), and a trend test against
log10(protein:RNA) (slope of the movie-level values, movies shuffled across all ratios).
Usage: python checks/vitro_area_test.py -> analysis_outputs/vitro_spotsum/area_test.csv"""
import numpy as np, scipy.io as sio, pandas as pd
from scipy.stats import theilslopes

D = 'in_vitro_reextracted/same_granules_R2/'
rng = np.random.default_rng(55)


def L(p):
    m = sio.loadmat(p, squeeze_me=True); k = [k for k in m if not k.startswith('__')][0]
    return [np.ravel(x).astype(float) for x in m[k][1:]]


def sl(y):
    ok = np.isfinite(y) & (y > 0)
    if ok.mean() < 0.9: return np.nan
    t = np.arange(len(y)); return theilslopes(np.log(y[ok]), t[ok])[0] * 360 * 100


info = pd.read_csv(D + 'extraction_info_same.csv')
data = {}
for st, lab, lg in [('same_1_1', '1:1', 0), ('same_10_1', '10:1', 1), ('same_100_1', '100:1', 2)]:
    S, B, G = L(D + st + '_data.mat'), L(D + st + '_background.mat'), L(D + st + '_spotsum_data.mat')
    I = info[info.dataset == st]; mov = np.repeat(np.arange(len(I)), I.nKept.values)
    v, m = [], []
    for s, b, g, mv in zip(S, B, G, mov):
        if len(s) < 180: continue
        a = sl(g / (s - b))
        if np.isfinite(a): v.append(a); m.append(mv)
    data[lab] = (np.array(v), np.array(m), lg)


def boot(v, m, nb=5000):
    movies = np.unique(m); out = np.empty(nb)
    for b in range(nb):
        pick = []
        for k in rng.choice(movies, len(movies)):
            idx = np.where(m == k)[0]; pick.extend(rng.choice(idx, len(idx)))
        out[b] = v[pick].mean()
    return np.percentile(out, [2.5, 97.5])


rows = []
for lab, (v, m, lg) in data.items():
    lo, hi = boot(v, m)
    rows.append(dict(test=lab, value=v.mean(), lo=lo, hi=hi, granules=len(v), movies=len(np.unique(m))))


def perm(a, b, n=10000):
    (va, ma, _), (vb, mb, _) = data[a], data[b]
    units = [va[ma == k] for k in np.unique(ma)] + [vb[mb == k] for k in np.unique(mb)]; na = len(np.unique(ma))
    obs = va.mean() - vb.mean(); c = 0
    for _ in range(n):
        p = rng.permutation(len(units)); x = np.concatenate([units[i] for i in p[:na]]); y = np.concatenate([units[i] for i in p[na:]])
        c += abs(x.mean() - y.mean()) >= abs(obs) - 1e-12
    return obs, (c + 1) / (n + 1)


for a, b in [('100:1', '1:1'), ('10:1', '1:1'), ('100:1', '10:1')]:
    d, p = perm(a, b); rows.append(dict(test=f'{a} vs {b}', value=d, p=p))
# trend with log10(protein:RNA): slope over granules; null = whole movies shuffled across the three ratios
units = [(data[l][0][data[l][1] == k], data[l][2]) for l in data for k in np.unique(data[l][1])]
labels = np.array([u[1] for u in units])


def slope(vals_lg):
    x = np.concatenate([np.full(len(v), g) for v, g in vals_lg]); y = np.concatenate([v for v, _ in vals_lg])
    return np.polyfit(x, y, 1)[0]


obs = slope(units); c = 0; n = 10000
for _ in range(n):
    lp = rng.permutation(labels); c += abs(slope([(u[0], g) for u, g in zip(units, lp)])) >= abs(obs) - 1e-12
rows.append(dict(test='trend per 10x protein:RNA', value=obs, p=(c + 1) / (n + 1)))
R = pd.DataFrame(rows); R.to_csv('analysis_outputs/vitro_spotsum/area_test.csv', index=False)
pd.set_option('display.width', 200); print(R.round(4).to_string(index=False))
