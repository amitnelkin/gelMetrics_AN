"""GAP_TEST  Is the burst-size gap between two conditions entirely explained by
noise? (RESULTS.md section 28). Null hypothesis H0: the gap between the
conditions' granule bursts equals the gap between their burst-free copies (same
brightness, same noise level, no bursts), i.e. the part beyond noise is 0.
Statistic: beyond = (real_A - copies_A) - (real_B - copies_B), with the median
fraction of the granule shed per burst (main), and the mean of log(fraction)
(sensitivity: uses every burst, less sensitive to the median's jumpiness).
In vivo R2, spot sum, window 13, SP 12, all bursts (no trim), traces >= 180 frames.
p: movie-level permutation (whole movies reassigned, 10000x, seeded), two-sided,
and one-sided in the direction the manuscript claims (FUS > PP7; Qbeta vs PP7:
direction of the observed difference, reported but not pre-specified).
CIs: nested movie -> granule bootstrap, 5000x, each condition resampled independently.
Usage: python checks/gap_test.py
Writes analysis_outputs/gap_test/{per_dataset.csv, tests.csv, tests.md}."""
import os
import numpy as np
import pandas as pd

out = os.path.join('analysis_outputs', 'gap_test'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(29); NPERM, NB = 10000, 5000
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
COMPS = [('4PP7-FUS', '4PP7', True), ('30PP7-FUS', '30PP7', True), ('5Qb', '4PP7', False), ('10Qb', '30PP7', False)]
B = pd.read_csv('analysis_outputs/detection_grid/bursts.csv')
B = B[(B.window == 13) & (B.SP == 12) & (B.levelBefore > 0) & B.dataset.isin(VIVO)].copy(); B['f'] = B.amp.abs() / B.levelBefore
T = pd.read_csv('analysis_outputs/detection_grid/traces.csv')
e = np.zeros(0)
U = {}
for ds in VIVO:
    t = T[T.dataset == ds]; b = B[B.dataset == ds]
    R = {k: g.f.values for k, g in b[b.kind == 'real'].groupby('trace')}; S = {k: g.f.values for k, g in b[b.kind == 'surr'].groupby('trace')}
    U[ds] = {m: [(R.get(k, e), S.get(k, e)) for k in g.trace] for m, g in t.groupby('movie')}
STATS = {'median': lambda x: 100 * np.median(x), 'mean_log': lambda x: np.mean(np.log(100 * x))}


def parts(items, st):
    r = np.concatenate([x[0] for x in items]); c = np.concatenate([x[1] for x in items])
    return STATS[st](r), STATS[st](c)


def flat(G):
    return [x for g in G.values() for x in g]


def boot(G, st):
    keys = list(G); v = []
    for _ in range(NB):
        it = []
        for k in rng.choice(len(keys), len(keys)):
            g = G[keys[k]]; it.extend(g[i] for i in rng.integers(0, len(g), len(g)))
        v.append(parts(it, st))
    return np.array(v)


rows, per = [], []
BS = {(ds, st): boot(U[ds], st) for ds in VIVO for st in STATS}
for ds in VIVO:
    for st in STATS:
        r, c = parts(flat(U[ds]), st); bs = BS[(ds, st)]
        lo, hi = np.percentile(bs[:, 0] - bs[:, 1], [2.5, 97.5])
        per.append(dict(dataset=ds, stat=st, real=r, copies=c, beyond=r - c, beyond_lo=lo, beyond_hi=hi))
for a, b, claimed in COMPS:
    for st in STATS:
        (ra, ca), (rb, cb) = parts(flat(U[a]), st), parts(flat(U[b]), st)
        gap, noise = ra - rb, ca - cb; beyond = gap - noise
        ba, bb = BS[(a, st)], BS[(b, st)]
        d_b = (ba[:, 0] - ba[:, 1]) - (bb[:, 0] - bb[:, 1]); d_g = ba[:, 0] - bb[:, 0]; d_n = ba[:, 1] - bb[:, 1]
        movies = list(U[a].values()) + list(U[b].values()); na = len(U[a]); null = np.empty(NPERM)
        for i in range(NPERM):
            p = rng.permutation(len(movies))
            xa = parts([x for j in p[:na] for x in movies[j]], st); xb = parts([x for j in p[na:] for x in movies[j]], st)
            null[i] = (xa[0] - xa[1]) - (xb[0] - xb[1])
        p2 = (np.sum(np.abs(null) >= abs(beyond) - 1e-12) + 1) / (NPERM + 1)
        p1 = (np.sum(null >= beyond - 1e-12) + 1) / (NPERM + 1) if beyond >= 0 else (np.sum(null <= beyond + 1e-12) + 1) / (NPERM + 1)
        rows.append(dict(compare=f'{a} vs {b}', stat=st, gap=gap, gap_lo=np.percentile(d_g, 2.5), gap_hi=np.percentile(d_g, 97.5),
                         noise_part=noise, noise_lo=np.percentile(d_n, 2.5), noise_hi=np.percentile(d_n, 97.5),
                         beyond=beyond, beyond_lo=np.percentile(d_b, 2.5), beyond_hi=np.percentile(d_b, 97.5),
                         share_beyond=beyond / gap if gap else np.nan, p_two_sided=p2, p_one_sided=p1,
                         one_sided_prespecified=claimed and beyond > 0))
P = pd.DataFrame(per); P.to_csv(os.path.join(out, 'per_dataset.csv'), index=False)
R = pd.DataFrame(rows); R.to_csv(os.path.join(out, 'tests.csv'), index=False)
pd.set_option('display.width', 250)
print(P.round(2).to_string(index=False)); print(R.round(4).to_string(index=False))
u = {'median': '% points', 'mean_log': 'log units'}
with open(os.path.join(out, 'tests.md'), 'w', encoding='utf-8', newline=chr(13) + chr(10)) as fh:
    fh.write('| comparison | statistic | total gap [95% CI] | explained by noise [95% CI] | beyond noise [95% CI] | p, two-sided | p, one-sided |' + chr(10) + '|---|---|---|---|---|---|---|' + chr(10))
    for _, r in R.iterrows():
        f = '.1f' if r.stat == 'median' else '.2f'
        fh.write(f"| {r['compare']} | {'median size (% points)' if r.stat == 'median' else 'mean log size'} | {r.gap:{f}} [{r.gap_lo:{f}}, {r.gap_hi:{f}}] | "
                 f"{r.noise_part:{f}} [{r.noise_lo:{f}}, {r.noise_hi:{f}}] | {r.beyond:{f}} [{r.beyond_lo:{f}}, {r.beyond_hi:{f}}] | {r.p_two_sided:.3f} | "
                 f"{r.p_one_sided:.3f}{' (pre-specified)' if r.one_sided_prespecified else ''} |" + chr(10))
