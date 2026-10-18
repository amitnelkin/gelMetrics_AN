"""BUILD_LAMBDA_K0_TABLE  Collects analysis_outputs/<stage>/lambda_k0_table.csv
for every stage (in commit order) into one before/after table.

Usage (repo root):  py -3 checks/build_lambda_k0_table.py
Writes analysis_outputs/lambda_k0_before_after.csv (long) and .md (wide:
one row per dataset x sign, one "lambda / k0" column per stage; '*' marks a
k0 at the edge of the grid). The 'paper' column holds the manuscript's
values (negative bursts only; the manuscript reports no positive-burst fits).
"""
import os

import pandas as pd

STAGES = ['baseline', 'after_bugfixes', 'after_B1_zerotrunc', 'after_B2_symbursts',
          'after_B3_bleach', 'after_B4_k0grid', 'after_B5_filloutliers']
PAPER_NEG = {  # MS.pdf pp. 7-10 (lambda, k0 or None, n bursts)
    'protein_RNA_1_1': (1, 39, 91), 'protein_RNA_10_1': (1, 63, 96),
    'protein_RNA_100_1': (1, 139, 103), '4PP7': (1, None, 63), '30PP7': (1, None, 112),
    '5Qb': (4, 90, 134), '10Qb': (4, 110, 181), '4PP7-FUS': (2, None, 210),
    '30PP7-FUS': (1, None, 86), 'SIM(lam=2,k0=60)': (4, None, None)}

frames = []
for st in STAGES:
    p = os.path.join('analysis_outputs', st, 'lambda_k0_table.csv')
    if os.path.exists(p):
        frames.append(pd.read_csv(p))
L = pd.concat(frames, ignore_index=True)
L.to_csv(os.path.join('analysis_outputs', 'lambda_k0_before_after.csv'), index=False)

L['cell'] = (L['lam1'].astype(str) + ' / ' + L['k0_1'].astype(str)
             + L['k0AtGridEdge'].astype(str).str.lower().isin(['1', 'true']).map({True: '*', False: ''}) + ' (n=' + L['nFit'].astype(str) + ')')
order = list(dict.fromkeys(L['dataset']))
W = L.pivot_table(index=['dataset', 'sign'], columns='stage', values='cell', aggfunc='first')
W = W[[s for s in STAGES if s in W.columns]]
# Stages that only reran the simulation (the change did not touch the
# experimental code path) carry the previous value forward, marked '='.
cols = list(W.columns)
for a, b in zip(cols[:-1], cols[1:]):
    miss = W[b].isna() & W[a].notna()
    W.loc[miss, b] = '= ' + W.loc[miss, a].astype(str).str.replace('= ', '', regex=False)
W = W.reindex(sorted(W.index, key=lambda ix: (order.index(ix[0]), ix[1] != 'Neg')))


def paper(ix):
    if ix[1] != 'Neg' or ix[0] not in PAPER_NEG:
        return ''
    lam, k0, n = PAPER_NEG[ix[0]]
    return f"{lam} / {k0 if k0 is not None else '?'}" + (f' (n={n})' if n else '')


W.insert(0, 'paper', [paper(ix) for ix in W.index])
W = W.reset_index()
with open(os.path.join('analysis_outputs', 'lambda_k0_before_after.md'), 'w', encoding='utf-8') as fh:
    fh.write('Rank-1 lambda / k0 (n = amplitudes fitted after the top-10% trim) per stage. '
             '* = k0 at a grid edge. Paper = MS.pdf, negative bursts only. '
             'SIM row: SimulationRunWithAnalysis.m, rng(1), true lambda = 2, k0 = 60.\n\n')
    fh.write(W.to_markdown(index=False))
    fh.write('\n')
print(W.to_string(index=False))
