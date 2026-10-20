"""NOISE_JITTER_SOURCES  Where the extra frame-to-frame jitter comes from
(RESULTS.md section 21). Per in vivo trace (R2 data, >= 180 frames), robust
relative frame-to-frame noise (1.4826*MAD(diff)/sqrt(2)/median) of the spot
mean, spot sum, cell background mean and multithresh spot area; the rank
correlation between the frame-to-frame changes of the spot mean and the cell
mean; and the spot mean / spot sum divided frame by frame by the cell mean.
Usage: python checks/noise_jitter_sources.py
Writes analysis_outputs/noise_background/jitter_sources.csv."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import spearmanr


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


def rs(x):
    d = np.diff(x)
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) / np.median(x)


R2, SU = 'in_vivo_reextracted/after_R2_offset', 'in_vivo_reextracted/after_R2_sums'
rows = []
for ds in ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']:
    S, B = load(f'{R2}/{ds}_data.mat'), load(f'{R2}/{ds}_background.mat')
    SS, A = load(f'{SU}/{ds}_spotsum_data.mat'), load(f'{SU}/{ds}_spotarea.mat')
    for i, (s, b, ss, a) in enumerate(zip(S, B, SS, A)):
        s, b, ss, a = (np.ravel(v).astype(float) for v in (s, b, ss, a))
        if len(s) < 180 or np.any(b <= 0):
            continue
        rows.append(dict(dataset=ds, trace=i + 1, spotMean=rs(s), spotSum=rs(ss), cellMean=rs(b), spotArea=rs(a),
                         spotMean_overCell=rs(s / b), spotSum_overCell=rs(ss / b),
                         rankCorr_dSpot_dCell=spearmanr(np.diff(s), np.diff(b))[0]))
J = pd.DataFrame(rows)
J.to_csv(os.path.join('analysis_outputs', 'noise_background', 'jitter_sources.csv'), index=False)
pd.set_option('display.width', 250)
print(J.groupby('dataset', sort=False).median(numeric_only=True).drop(columns='trace').round(3).to_string())
