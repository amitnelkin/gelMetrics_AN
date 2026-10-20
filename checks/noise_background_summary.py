"""NOISE_BACKGROUND_SUMMARY  Summarises checks/noise_from_background.m: measured
frame-to-frame noise of each trace vs the camera + photon noise predicted from
pixels away from all granules (RESULTS.md section 21). Traces >= 180 frames.
Usage: python checks/noise_background_summary.py
Writes analysis_outputs/noise_background/{summary.csv, noise_measured_vs_predicted.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'noise_background')
T = pd.read_csv(os.path.join(d, 'traces.csv')); T = T[T.lengthFrames >= 180].copy()
C = pd.read_csv(os.path.join(d, 'calibration.csv'))
ORDER = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
for m in ['Sum', 'Mean', 'Cell']:
    T[f'ratio{m}'] = T[f'sd{m}_meas'] / T[f'sd{m}_pred']
T['relNoiseSum'] = T.sdSum_meas / T.spotSum
q = lambda x: f'{np.median(x):.2f} [{np.percentile(x, 25):.2f}-{np.percentile(x, 75):.2f}]'
S = T.groupby('dataset').agg(traces=('trace', 'size'), ratioSum=('ratioSum', q), ratioMean=('ratioMean', q), ratioCell=('ratioCell', q),
                             relNoiseSum=('relNoiseSum', 'median'), fracAboveCal=('fracSpotPxAboveCal', 'median')).reindex(ORDER)
S = S.join(C.groupby('dataset').agg(gain_a=('gain_a', 'median'), calR2_min=('R2', 'min')))
S.to_csv(os.path.join(d, 'summary.csv'))
pd.set_option('display.width', 250)
print('measured / predicted noise, median [IQR] over traces:'); print(S.to_string())
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, m, lab in zip(axs, ['Sum', 'Mean', 'Cell'], ['granule spot sum', 'granule spot mean', 'cell background mean']):
    for ds in ORDER:
        t = T[T.dataset == ds]
        ax.loglog(t[f'sd{m}_pred'], t[f'sd{m}_meas'], '.', ms=3, alpha=0.6, label=ds)
    lo = min(T[f'sd{m}_pred'].min(), T[f'sd{m}_meas'].min()); hi = max(T[f'sd{m}_pred'].max(), T[f'sd{m}_meas'].max())
    ax.plot([lo, hi], [lo, hi], 'k-', lw=1); ax.plot([lo, hi], [2 * lo, 2 * hi], 'k:', lw=1)
    ax.set_xlabel('predicted camera + photon noise\n(from pixels away from granules)'); ax.set_ylabel('measured noise on the trace')
    ax.set_title(lab, fontsize=10, loc='left'); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
axs[0].legend(fontsize=7, markerscale=3); axs[2].text(0.02, 0.97, 'solid: equal; dotted: 2x', transform=axs[2].transAxes, va='top', fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(d, 'noise_measured_vs_predicted.png'), dpi=140); fig.savefig(os.path.join(d, 'noise_measured_vs_predicted.pdf'))
