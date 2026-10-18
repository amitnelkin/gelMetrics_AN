"""CHECK_INTERLEAVE  Data-quality check on the stored traces (read-only).

ExampleTIF.tif is a 2-channel x 361-frame interleaved hyperstack, and
generate_movArray.m reads every TIFF page as one frame. If a movie given to
Collect_intensity_data.m was interleaved like this, its traces alternate
between two levels (channel or illumination state) from frame to frame.

Per trace (>= 60 frames) this reports:
  alt   - fraction of consecutive frame differences that change sign
          (independent noise ~ 0.67, strict alternation = 1.0)
  ratio - mean of odd frames / mean of even frames
A trace is flagged when alt > 0.85 or |log(ratio)| > 0.2.

Usage:  py -3 checks/check_interleave.py   (from the repo root)
Writes analysis_outputs/baseline/data_quality/interleave_per_trace.csv and
interleave_summary.csv.
"""
import glob
import os

import numpy as np
import pandas as pd
import scipy.io as sio

out_dir = os.path.join('analysis_outputs', 'baseline', 'data_quality')
os.makedirs(out_dir, exist_ok=True)
rows = []
for f in sorted(glob.glob(os.path.join('in_*_experiments', '*_data.mat'))):
    m = sio.loadmat(f)
    key = [k for k in m if not k.startswith('__')][0]
    cells = m[key].ravel()[1:]          # first cell is the empty placeholder
    name = os.path.basename(f).replace('_data.mat', '')
    for i, c in enumerate(cells):
        x = np.asarray(c, float).ravel()
        if x.size < 60:
            continue
        s = np.sign(np.diff(x))
        s = s[s != 0]
        alt = float(np.mean(s[1:] != s[:-1]))
        ratio = float(np.mean(x[0::2]) / np.mean(x[1::2]))
        rows.append(dict(dataset=name, trace=i + 1, length=x.size, alt=alt, ratio=ratio,
                         flagged=bool(alt > 0.85 or abs(np.log(ratio)) > 0.2)))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(out_dir, 'interleave_per_trace.csv'), index=False)
L = T[T.length >= 180]
S = (T.groupby('dataset').agg(n60=('trace', 'size'), flagged60=('flagged', 'sum'),
                              median_alt=('alt', 'median'))
     .join(L.groupby('dataset').agg(n180=('trace', 'size'), flagged180=('flagged', 'sum'),
                                    strict_alt180=('alt', lambda a: int((a > 0.95).sum())))))
S.to_csv(os.path.join(out_dir, 'interleave_summary.csv'))
print(S.to_string())
