"""FIXEDMASK_MOTION  Are the slow changes measured with the fixed mask (R5) granule
exchange, or the granule moving / changing shape relative to the mask?
(RESULTS.md section 33). Per trace (>= 180 frames): the spot sum and the share
of the granule's light inside the fixed mask (maskFrac), each detrended (cubic)
and smoothed over 13 frames as the detector does; their correlation, and the
variance of the mask share relative to that of the spot sum. Exchange changes the
amount of light but hardly its share inside a fixed mask; movement or reshaping
changes the share.
Usage: python checks/fixedmask_motion.py
Writes analysis_outputs/fixedmask/motion.csv."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.ndimage import uniform_filter1d


def load(p):
    m = sio.loadmat(p, squeeze_me=True)
    return m[[k for k in m if not k.startswith('__')][0]][1:]


R = 'in_vivo_reextracted/after_R5_fixedmask'; rows = []
for ds in ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']:
    for i, (s, f) in enumerate(zip(load(f'{R}/{ds}_spotsum_data.mat'), load(f'{R}/{ds}_maskfrac.mat'))):
        s = np.ravel(s).astype(float); f = np.ravel(f).astype(float)
        if len(s) < 180:
            continue
        t = np.arange(len(s))
        a = uniform_filter1d(s / np.polyval(np.polyfit(t, s, 3), t) - 1, 13)
        c = uniform_filter1d(f / np.polyval(np.polyfit(t, f, 3), t) - 1, 13)
        rows.append(dict(dataset=ds, trace=i + 1, corr=np.corrcoef(a, c)[0, 1], var_ratio=np.var(c) / np.var(a)))
M = pd.DataFrame(rows); M.to_csv(os.path.join('analysis_outputs', 'fixedmask', 'motion.csv'), index=False)
print(M.groupby('dataset', sort=False)[['corr', 'var_ratio']].median().round(2).to_string())
