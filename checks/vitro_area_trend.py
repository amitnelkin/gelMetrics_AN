"""VITRO_AREA_TREND  Why spot mean and spot sum disagree in vitro (RESULTS.md section 53).
spot sum = N x (spot mean - cell mean), so N = spot sum / (spot mean - cell mean) is the
number of granule pixels (the top two multithresh classes of the 15 x 15 crop).
Per granule (>= 180 frames, same_granules_R2): Theil-Sen slope (%/h) of ln(spot mean -
cell), ln(N) and ln(spot sum), all relative to the cell (dilute phase) mean, as the net
change; median N; share of frames with N >= 200 of 225 pixels (granule fills the crop).
Usage: python checks/vitro_area_trend.py -> analysis_outputs/vitro_spotsum/area_trend.csv"""
import numpy as np, scipy.io as sio, pandas as pd
from scipy.stats import theilslopes
D = 'in_vitro_reextracted/same_granules_R2/'
def L(p):
    m = sio.loadmat(p, squeeze_me=True); k = [k for k in m if not k.startswith('__')][0]
    return [np.ravel(x).astype(float) for x in m[k][1:]]
def sl(y):
    ok = np.isfinite(y) & (y > 0)
    if ok.mean() < 0.9: return np.nan
    t = np.arange(len(y)); return theilslopes(np.log(y[ok]), t[ok])[0] * 360 * 100
rows = []
for st, lab in [('same_1_1', '1:1'), ('same_10_1', '10:1'), ('same_100_1', '100:1')]:
    S, B, G = L(D + st + '_data.mat'), L(D + st + '_background.mat'), L(D + st + '_spotsum_data.mat')
    for s, b, g in zip(S, B, G):
        if len(s) < 180: continue
        ex = s - b; N = g / ex
        rows.append(dict(ratio=lab, net_mean=sl(ex / b), net_sum=sl(g / b), area_trend=sl(N), cell_trend=sl(b), areaMedian=np.nanmedian(N), fullCrop=np.nanmean(N >= 200)))
R = pd.DataFrame(rows); R.to_csv('analysis_outputs/vitro_spotsum/area_trend.csv', index=False)
print(R.groupby('ratio', sort=False).agg(granules=('net_mean', 'size'), net_mean=('net_mean', 'mean'), net_sum=('net_sum', 'mean'), area_trend=('area_trend', 'mean'),
      area_median=('areaMedian', 'median'), frames_full_crop=('fullCrop', 'mean')).round(2).to_string())
