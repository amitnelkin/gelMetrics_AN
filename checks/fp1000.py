"""1-in-1000 false-burst calibration for 4PP7 vs 4PP7-FUS (RESULTS.md section 37).
Input: analysis_outputs/fp1000/ from checks/fp1000.m (100 burst-free copies per granule).
Rule fixed before looking at the real granules:
  A (run length only): lowest SP at which <= 1/1000 copies show any burst (either sign).
  B (SP 12 + size):    burst size = |amplitude| / granule level (larger of the levels at
                       the burst ends); cutoff = size that <= 1/1000 copies exceed.
Each rule is set per dataset (each condition at its own 1/1000) and in common (the
stricter of the two, the same rule for both). Real granules are then counted at that
setting: granules with >= 1 burst vs the 1/1000 expected from noise, and 4x vs FUS
by movie-level permutation."""
import numpy as np, pandas as pd, matplotlib, os
from scipy.stats import binom
matplotlib.use('Agg'); import matplotlib.pyplot as plt

D = os.path.join(os.path.dirname(__file__), '..', 'analysis_outputs', 'fp1000')
sg = pd.read_csv(os.path.join(D, 'segs.csv')); tr = pd.read_csv(os.path.join(D, 'traces.csv'))
names = pd.read_csv(os.path.join(D, 'datasets.csv'), header=None)[0].tolist()
NC = int(sg['copy'].max()); RATE = 1e-3; rng = np.random.default_rng(1)
den = np.maximum(sg.levelBefore, sg.levelAfter)
sg['size'] = np.where(den > 0, sg.amp.abs() / den, np.inf)
SPs = sorted(sg.SP.unique())

def any_rate(df, n):           # fraction of granules (or copies) with >= 1 burst
    return df.groupby(['copy', 'trace']).ngroups / n

out, curves = [], []
cal = {}
for j, name in enumerate(names, 1):
    s = sg[sg.dataset == j]; n = (tr.dataset == j).sum(); cp = s[s['copy'] > 0]
    # A: false rate per copy at each SP
    rA = {sp: any_rate(cp[cp.SP == sp], n * NC) for sp in SPs}
    spA = next((sp for sp in SPs if rA[sp] <= RATE), None)
    # B: per-copy largest size at SP 12 (0 if no burst) -> 99.9th percentile
    mx = cp[cp.SP == 12].groupby(['copy', 'trace'])['size'].max()
    allmx = np.zeros(n * NC); allmx[:len(mx)] = mx.values
    cutB = np.quantile(allmx, 1 - RATE)
    cal[name] = dict(n=n, spA=spA, cutB=cutB, rA=rA)
    for sp in SPs: curves.append(dict(dataset=name, SP=sp, copy_rate=rA[sp],
                                       real_rate=s[(s['copy'] == 0) & (s.SP == sp)].trace.nunique() / n))
spC = max(c['spA'] for c in cal.values() if c['spA']) if all(c['spA'] for c in cal.values()) else None
cutC = max(c['cutB'] for c in cal.values())

def granules_hit(j, rule, val):
    s = sg[(sg.dataset == j) & (sg['copy'] == 0)]
    if rule == 'A': s = s[s.SP == val]
    else: s = s[(s.SP == 12) & (s['size'] > val)]
    return set(s.trace), s

hits = {}
for j, name in enumerate(names, 1):
    c = cal[name]; n = c['n']
    for lab, rule, val in (('A own', 'A', c['spA']), ('A common', 'A', spC), ('B own', 'B', c['cutB']), ('B common', 'B', cutC)):
        if val is None: continue
        g, s = granules_hit(j, rule, val); hits[(name, lab)] = g
        out.append(dict(dataset=name, rule=lab, setting=(f'SP {val}' if rule == 'A' else f'size > {100*val:.0f}%'),
                        granules=n, granules_with_burst=len(g), expected_from_noise=round(RATE * n, 2),
                        p_excess=binom.sf(len(g) - 1, n, RATE), bursts_neg=int((s.sign < 0).sum()), bursts_pos=int((s.sign > 0).sum()),
                        median_size=float(np.median(s['size'])) if len(s) else np.nan))
res = pd.DataFrame(out)

# 4x vs FUS: share of granules with a burst, movie-level permutation (movies as units)
cmp = []
if len(names) == 2:
    for lab in res.rule.unique():
        frac = []
        for j, name in enumerate(names, 1):
            t = tr[tr.dataset == j].copy(); t['hit'] = t.trace.isin(hits.get((name, lab), set()))
            frac.append(t.groupby('movie').hit.agg(['sum', 'count']))
        a, b = frac; obs = a['sum'].sum() / a['count'].sum() - b['sum'].sum() / b['count'].sum()
        M = pd.concat([a, b]).values; na = len(a); cnt = 0
        for _ in range(20000):
            p = rng.permutation(len(M)); x, y = M[p[:na]], M[p[na:]]
            cnt += abs(x[:, 0].sum() / x[:, 1].sum() - y[:, 0].sum() / y[:, 1].sum()) >= abs(obs) - 1e-12
        cmp.append(dict(rule=lab, share_4x=a['sum'].sum() / a['count'].sum(), share_FUS=b['sum'].sum() / b['count'].sum(),
                        p_movie_perm=(cnt + 1) / 20001))
cmp = pd.DataFrame(cmp)
res.to_csv(os.path.join(D, 'summary.csv'), index=False); cmp.to_csv(os.path.join(D, 'compare.csv'), index=False)
cv = pd.DataFrame(curves); cv.to_csv(os.path.join(D, 'rates_by_SP.csv'), index=False)
pd.set_option('display.width', 200)
print({k: dict(n=v['n'], spA=v['spA'], cutB=round(v['cutB'], 3)) for k, v in cal.items()}, 'common SP', spC, 'common cut', round(cutC, 3))
print(cv.round(4).to_string(index=False)); print(res.round(4).to_string(index=False)); print(cmp.round(4).to_string(index=False))

# figure: share of granules with >= 1 burst vs SP, real vs copies
COL = {'4PP7': '#0072B2', '4PP7-FUS': '#D55E00'}
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for name in names:
    c = cv[cv.dataset == name]
    ax[0].semilogy(c.SP, c.real_rate.clip(lower=1e-4), 'o-', color=COL.get(name), label=f'{name} granules')
    ax[0].semilogy(c.SP, c.copy_rate.clip(lower=1e-4), 's--', color=COL.get(name), mfc='white', label=f'{name} noise copies')
ax[0].axhline(RATE, color='k', lw=.8, ls=':'); ax[0].text(SPs[0], RATE * 1.2, '1 in 1000', fontsize=8)
ax[0].set_xlabel('Stringency SP (run length, frames)'); ax[0].set_ylabel('Granules with ≥1 burst (fraction)'); ax[0].legend(frameon=False, fontsize=7)
ax[0].set_title('A: raising SP (bottom row = none)')
for name in names:
    j = names.index(name) + 1; cp = sg[(sg.dataset == j) & (sg['copy'] > 0) & (sg.SP == 12)]
    rl = sg[(sg.dataset == j) & (sg['copy'] == 0) & (sg.SP == 12)]
    mxc = np.zeros(cal[name]['n'] * NC); v = cp.groupby(['copy', 'trace'])['size'].max().values; mxc[:len(v)] = v
    mxr = np.zeros(cal[name]['n']); v = rl.groupby('trace')['size'].max().values; mxr[:len(v)] = v
    xs = np.linspace(0, 1.5, 301)
    ax[1].semilogy(100 * xs, [(mxr > x).mean() for x in xs], '-', color=COL.get(name), label=f'{name} granules')
    ax[1].semilogy(100 * xs, [(mxc > x).mean() for x in xs], '--', color=COL.get(name), label=f'{name} noise copies')
ax[1].axhline(RATE, color='k', lw=.8, ls=':'); ax[1].set_ylim(1e-4, 1.2); ax[1].set_xlim(0, 60)
ax[1].set_xlabel('Size cutoff (% of granule), SP 12'); ax[1].set_ylabel('Granules with a burst above cutoff'); ax[1].set_title('B: size cutoff')
ax[1].legend(frameon=False, fontsize=7)
for a in ax: a.tick_params(direction='in')
fig.tight_layout(); fig.savefig(os.path.join(D, 'fp1000.png'), dpi=150)
