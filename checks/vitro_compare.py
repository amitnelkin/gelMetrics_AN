"""VITRO_COMPARE  Are the titration movies (Raw_Results/Titration experiment) the
source of the stored in vitro data? (RESULTS.md section 38.)
Re-extracted traces: in_vitro_reextracted/after_R1_redchannel/<folder>_data.mat
(fluorescence channel, no offset, as the stored data). Stored: in_vitro_experiments/.
Per stored dataset x new folder:
  trace counts (all kept, >= 180 frames), median length, median spot level, median
  cell level;
  trace matching: a stored trace is matched if some new trace of the same length
  (+-2 frames) has values within 1% (median absolute relative difference), allowing
  a shift of -1/0/+1 frame (the original reader read page f instead of f+1 and
  dropped the first point).
Usage: python checks/vitro_compare.py
Writes analysis_outputs/vitro_compare/{summary.csv, matching.csv}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio

out = os.path.join('analysis_outputs', 'vitro_compare'); os.makedirs(out, exist_ok=True)
STORED = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
NEW = ['vitro_1000ngul_1mgml', 'vitro_100ngul_1mgml', 'vitro_1000ng_10ug', 'vitro_1000ngul_0p1mgml']
R1 = 'in_vitro_reextracted/after_R1_redchannel'


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return [np.ravel(x).astype(float) for x in m[[k for k in m if not k.startswith('__')][0]][1:]]


def desc(S, B, name):
    L = np.array([len(s) for s in S])
    return dict(dataset=name, traces=len(S), traces_ge180=int((L >= 180).sum()), median_length=float(np.median(L)),
                median_spot=float(np.median([np.median(s) for s in S])), median_cell=float(np.median([np.median(b) for b in B])))


def match(s, cand):
    best = np.inf
    for n in cand:
        if abs(len(n) - len(s)) > 2:
            continue
        for sh in (-1, 0, 1):
            a, b = (s[max(0, -sh):], n[max(0, sh):])
            k = min(len(a), len(b))
            if k < 20:
                continue
            best = min(best, np.median(np.abs(a[:k] - b[:k]) / np.maximum(np.abs(a[:k]), 1)))
    return best


rows, mrows = [], []
St = {d: (load(f'in_vitro_experiments/{d}_data.mat'), load(f'in_vitro_experiments/{d}_background.mat')) for d in STORED}
Nw = {d: (load(f'{R1}/{d}_data.mat'), load(f'{R1}/{d}_background.mat')) for d in NEW}
for d, (S, B) in St.items():
    rows.append(dict(source='stored', **desc(S, B, d)))
for d, (S, B) in Nw.items():
    rows.append(dict(source='re-extracted', **desc(S, B, d)))
for d, (S, _) in St.items():
    for e, (N, _) in Nw.items():
        bm = np.array([match(s, N) for s in S])
        mrows.append(dict(stored=d, new=e, matched_1pct=float((bm <= 0.01).mean()), matched_5pct=float((bm <= 0.05).mean()),
                          median_best_diff=float(np.median(bm))))
    print(d, 'done', flush=True)
D = pd.DataFrame(rows); M = pd.DataFrame(mrows)
D.to_csv(os.path.join(out, 'summary.csv'), index=False); M.to_csv(os.path.join(out, 'matching.csv'), index=False)
pd.set_option('display.width', 200)
print(D.round(1).to_string(index=False)); print()
print('share of stored traces with a matching re-extracted trace (within 1%):')
print(M.pivot(index='stored', columns='new', values='matched_1pct').round(2).to_string())
print('\nwithin 5%:'); print(M.pivot(index='stored', columns='new', values='matched_5pct').round(2).to_string())
