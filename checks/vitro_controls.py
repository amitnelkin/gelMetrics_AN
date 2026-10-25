"""VITRO_CONTROLS  The two noise controls on the re-extracted in vitro movies
(RESULTS.md section 39), as sections 21 and 23 did in vivo.
1. Camera-noise calibration (checks/noise_from_background.m -> noise_background_vitro/):
   measured frame-to-frame noise of each granule's spot mean and of its dilute-phase
   background vs the camera + photon noise predicted from pixels away from granules.
2. Unchanging granule (checks/static_control.m, measure 'mean', per-frame thresholds
   -> steady/static_control_vitro.csv): a synthetic granule that never changes, made
   from each real track's own image, bleached and with calibrated noise, measured and
   analysed like the real traces (window 13, SP 12, burst-free copies). Its "turnover
   beyond noise" should be ~0 if the copies are a valid noise floor.
Nested movie -> granule bootstrap SE (1000x). Traces >= 180 frames. Section 43: movies dropped
by the authors (checks/vitro_movies.py) are left out; the real-granule columns come from
the primary in vitro data (stored granules, corrected pipeline: analysis_outputs/vitro_same).
Usage: python checks/vitro_controls.py
Writes analysis_outputs/noise_background_vitro/summary.csv,
analysis_outputs/steady/static_control_vitro_summary.csv and vitro_controls.md."""
import os
import numpy as np
import sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vitro_movies import dropped, keep_bursts

rng = np.random.default_rng(40); H = 10 / 3600
LAB = {'vitro_1000ngul_0p1mgml': 'protein_RNA_1_1', 'vitro_1000ng_10ug': 'protein_RNA_10_1',
       'vitro_1000ngul_1mgml': 'protein_RNA_100_1', 'vitro_100ngul_1mgml': 'extra_100ngul_1mgml'}
NAME = {'protein_RNA_1_1': '1:1', 'protein_RNA_10_1': '10:1', 'protein_RNA_100_1': '100:1', 'extra_100ngul_1mgml': 'extra'}
ORDER = list(LAB.values())

d = os.path.join('analysis_outputs', 'noise_background_vitro')
T = pd.read_csv(os.path.join(d, 'traces.csv')); T = T[T.lengthFrames >= 180].copy(); T['dataset'] = T.dataset.map(LAB)
T = T[[not dropped(a, b) for a, b in zip(T.dataset, T.movie)]]
C = pd.read_csv(os.path.join(d, 'calibration.csv')); C['dataset'] = C.dataset.map(LAB)
C = C[[not dropped(a, b) for a, b in zip(C.dataset, C.movie)]]
for m in ['Sum', 'Mean', 'Cell']:
    T[f'ratio{m}'] = T[f'sd{m}_meas'] / T[f'sd{m}_pred']
T['relNoiseSum'] = T.sdSum_meas / T.spotSum
T['relNoiseMean'] = T.sdMean_meas / T.spotMean
q = lambda x: f'{np.median(x):.2f} [{np.percentile(x, 25):.2f}-{np.percentile(x, 75):.2f}]'
NBs = T.groupby('dataset').agg(traces=('trace', 'size'), ratioSum=('ratioSum', q), ratioMean=('ratioMean', q), ratioCell=('ratioCell', q),
                               relNoiseSum=('relNoiseSum', 'median'), relNoiseMean=('relNoiseMean', 'median'),
                               fracAboveCal=('fracSpotPxAboveCal', 'median')).reindex(ORDER)
NBs = NBs.join(C.groupby('dataset').agg(gain_a=('gain_a', 'median'), calR2_min=('R2', 'min')))
NBs.to_csv(os.path.join(d, 'summary.csv'))


def nested(stat, t, nb=1000):
    movs = {m: g.index.values for m, g in t.groupby('movie')}; v = []
    for _ in range(nb):
        pick = np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        v.append(stat(t.loc[pick]))
    return np.std(v)


S = pd.read_csv(os.path.join('analysis_outputs', 'steady', 'static_control_vitro.csv')); S['dataset'] = S.dataset.map(LAB)
S = S[[not dropped(a, b) for a, b in zip(S.dataset, S.movie)]]
vn = pd.read_csv(os.path.join('analysis_outputs', 'vitro_same', 'summary.csv')).set_index('dataset')
G = pd.read_csv(os.path.join('analysis_outputs', 'detection_grid_vitro_same', 'bursts.csv'))
rows = []
for ds in ORDER:
    if ds not in vn.index:
        continue
    t = S[(S.dataset == ds) & (S.version == 'v0')].reset_index(drop=True)
    exc = lambda x: 100 * (x.fracReal13.sum() - x.fracSurr13.sum()) / (x.lengthFrames.sum() * H)
    g = G[G.dataset == ds]
    tt = t.dropna(subset=['D1', 'D30']); g30 = lambda x: 100 * np.sqrt(max((x.D30 - x.D1).mean(), 0))
    rows.append(dict(dataset=ds, version='per-frame (R2)', granules=len(t), static_excess_w13=exc(t), static_excess_se_w13=nested(exc, t),
                     static_FDR_w13=t.nSurr13.sum() / max(t.nReal13.sum(), 1),
                     real_excess_w13=vn.loc[ds, 'turn0'], real_FDR_w13=(g.kind == 'surr').sum() / max((g.kind == 'real').sum(), 1),
                     static_G30=g30(tt), static_G30_se=nested(g30, tt.reset_index(drop=True))))
R = pd.DataFrame(rows); R.to_csv(os.path.join('analysis_outputs', 'steady', 'static_control_vitro_summary.csv'), index=False)

with open(os.path.join('analysis_outputs', 'steady', 'vitro_controls.md'), 'w', encoding='utf-8', newline=chr(13) + chr(10)) as fh:
    fh.write('| condition | granules | measured ÷ camera noise, granule (spot mean) | measured ÷ camera noise, background | granule frame noise (% of level) | turnover beyond noise, unchanging granule (%/h, ±2 SE) | real granules (%/h) | copies ÷ detections, unchanging / real |' + chr(10))
    fh.write('|---|---|---|---|---|---|---|---|' + chr(10))
    for _, r in R.iterrows():
        ds = r.dataset; n = NBs.loc[ds]
        fh.write(f"| {NAME[ds]} | {r.granules} | {n.ratioMean} | {n.ratioCell} | {100 * n.relNoiseMean:.1f} | {r.static_excess_w13:.1f} ± {2 * r.static_excess_se_w13:.1f} | {r.real_excess_w13:.1f} | {r.static_FDR_w13:.2f} / {r.real_FDR_w13:.2f} |" + chr(10))
print(open(os.path.join('analysis_outputs', 'steady', 'vitro_controls.md'), encoding='utf-8').read())
