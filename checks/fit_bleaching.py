"""FIT_BLEACHING  Photobleaching rate of the real background traces, used to
set expConst in SimulationRunWithAnalysis.m (CLAUDE.md [model] change 3:
the original 0.0005 /min gives 3% over 60 min).

Model per trace: B(t) = c + B0 * exp(-k t), t in minutes (10 s frames),
fitted as a straight line to log(B - c). Two values of the camera offset c:
  c = 0   : decay of the stored intensities as analysis.m sees them
  c = 200 : offset of the Andor iXon (ExampleTIF.tif floor ~195-200), i.e.
            the bleaching rate of the fluorescence itself; generateSim.m has
            no offset, so this is the rate that matches the simulation.
Traces >= 180 frames, as in analysis.m.

Usage (repo root):  py -3 checks/fit_bleaching.py <outDir>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio

out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join('analysis_outputs', 'bleaching')
os.makedirs(out_dir, exist_ok=True)
dt_min = 10 / 60
rows = []
for f in sorted(glob.glob(os.path.join('in_*_experiments', '*_background.mat'))):
    m = sio.loadmat(f)
    key = [k for k in m if not k.startswith('__')][0]
    name = os.path.basename(f).replace('_background.mat', '')
    for i, c in enumerate(m[key].ravel()[1:]):
        b = np.asarray(c, float).ravel()
        if b.size < 180:
            continue
        t = np.arange(b.size) * dt_min
        for off in (0, 200):
            y = b - off
            if np.any(y <= 0):
                continue
            k = -np.polyfit(t, np.log(y), 1)[0]
            rows.append(dict(dataset=name, trace=i + 1, length=b.size, offset=off, k_per_min=k))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(out_dir, 'bleach_per_trace.csv'), index=False)
S = (T.groupby(['offset', 'dataset']).k_per_min
     .agg(n='size', median='median', q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75))
     .reset_index())
S['bleach_pct_60min_at_median'] = 100 * (1 - np.exp(-S['median'] * 60))
S.to_csv(os.path.join(out_dir, 'bleach_summary.csv'), index=False)
print(S.to_string(index=False))
for off in (0, 200):
    med = T[T.offset == off].groupby('dataset').k_per_min.median()
    print(f'offset={off}: median over dataset medians = {med.median():.5f} /min '
          f'({100 * (1 - np.exp(-med.median() * 60)):.1f}% over 60 min)')
