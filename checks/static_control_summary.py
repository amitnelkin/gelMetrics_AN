"""STATIC_CONTROL_SUMMARY  Apparent dynamics in static synthetic movies
(checks/static_control.m) vs the real movies, per extraction version
(thrWin 0 = per-frame thresholds, published; 31 = steady). RESULTS.md section 23.
Excess turnover (%/h) and false-discovery share at window 13 and 17, SP 12;
noise-free fluctuation G30 = 100*sqrt(mean(D30 - D1)). Nested (movie -> granule)
bootstrap SE, 1000x. Usage: python checks/static_control_summary.py
Writes analysis_outputs/steady/static_control_summary.csv."""
import os
import numpy as np
import pandas as pd

d = os.path.join('analysis_outputs', 'steady')
rng = np.random.default_rng(24)
H = 10 / 3600
S = pd.read_csv(os.path.join(d, 'static_control.csv'))
Dr = pd.read_csv(os.path.join(d, 'detection.csv')); Mr = pd.read_csv(os.path.join(d, 'm2b.csv'))
VNAME = {'v0': 'per-frame (R2)', 'v31': 'steady (R4)'}


def nested(stat, t, nb=1000):
    movs = {m: g.index.values for m, g in t.groupby('movie')}; v = []
    for _ in range(nb):
        pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        v.append(stat(t.loc[pick]))
    return np.std(v)


rows = []
for (ds, ver), t in S.groupby(['dataset', 'version']):
    t = t.reset_index(drop=True)
    r = dict(dataset=ds, version=VNAME[ver], granules=len(t))
    for w in (13, 17):
        exc = lambda x, w=w: 100 * (x[f'fracReal{w}'].sum() - x[f'fracSurr{w}'].sum()) / (x.lengthFrames.sum() * H)
        r[f'static_excess_w{w}'] = exc(t); r[f'static_excess_se_w{w}'] = nested(exc, t)
        r[f'static_FDR_w{w}'] = t[f'nSurr{w}'].sum() / max(t[f'nReal{w}'].sum(), 1)
        real = Dr[(Dr.version == VNAME[ver]) & (Dr.dataset == ds) & (Dr.window == w) & (Dr.SP == 12)].iloc[0]
        r[f'real_excess_w{w}'] = real.excess; r[f'real_FDR_w{w}'] = real.FDR
    tt = t.dropna(subset=['D1', 'D30'])
    g = lambda x: 100 * np.sqrt(max((x.D30 - x.D1).mean(), 0))
    r['static_G30'] = g(tt); r['static_G30_se'] = nested(g, tt.reset_index(drop=True))
    r['real_G30'] = Mr[(Mr.version == VNAME[ver]) & (Mr.dataset == ds)].G30.iloc[0]
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(os.path.join(d, 'static_control_summary.csv'), index=False)
pd.set_option('display.width', 250)
print(R.round(2).T.to_string())
