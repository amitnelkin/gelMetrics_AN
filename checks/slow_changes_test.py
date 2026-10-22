"""SLOW_CHANGES_TEST  Condition comparisons for the slow-change counts of
checks/slow_changes.py (RESULTS.md section 35). Statistic: excess events per
granule-hour (real minus noise-only copies), rises + drops, at each threshold.
p: movie-level permutation (whole movies reassigned between the two conditions,
10000x, seeded), two-sided and one-sided (FUS > PP7 pre-specified).
Usage: python checks/slow_changes_test.py [--absolute]
Writes analysis_outputs/slow_changes/tests.csv."""
import os
import sys
import numpy as np
import pandas as pd

rng = np.random.default_rng(36); NP = 10000
FOLD = 'slow_changes' + ('_absolute' if '--absolute' in sys.argv else '')
P = pd.read_csv(os.path.join('analysis_outputs', FOLD, 'per_trace.csv'))
COMPS = [('4PP7-FUS', '4PP7', True), ('30PP7-FUS', '30PP7', True), ('5Qb', '4PP7', False), ('10Qb', '30PP7', False)]
rows = []
for thr in [0.2, 0.3, 0.4]:
    P['ex'] = P[f'real_up_{thr}'] + P[f'real_dn_{thr}'] - P[f'copy_up_{thr}'] - P[f'copy_dn_{thr}']
    for a, b, claimed in COMPS:
        units = [g for _, g in P[P.dataset == a].groupby('movie')] + [g for _, g in P[P.dataset == b].groupby('movie')]
        na = P[P.dataset == a].movie.nunique()
        ex = [(g.ex.sum(), g.hours.sum()) for g in units]
        rate = lambda idx: sum(ex[i][0] for i in idx) / sum(ex[i][1] for i in idx)
        obs = rate(range(na)) - rate(range(na, len(units))); null = np.empty(NP)
        for k in range(NP):
            p = rng.permutation(len(units)); null[k] = rate(p[:na]) - rate(p[na:])
        rows.append(dict(threshold=int(thr * 100), compare=f'{a} vs {b}', excessA=rate(range(na)), excessB=rate(range(na, len(units))), diff=obs,
                         p_two_sided=(np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (NP + 1),
                         p_one_sided=(np.sum(null >= obs - 1e-12) + 1) / (NP + 1) if obs >= 0 else (np.sum(null <= obs + 1e-12) + 1) / (NP + 1),
                         one_sided_prespecified=claimed and obs > 0))
T = pd.DataFrame(rows); T.to_csv(os.path.join('analysis_outputs', FOLD, 'tests.csv'), index=False)
pd.set_option('display.width', 200); print(T.round(4).to_string(index=False))
