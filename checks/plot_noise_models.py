"""PLOT_NOISE_MODELS  Summary of checks/noise_model_validation.m (R2 4b, R2 6).
Metrics per model (lambda_true pooled, both signs):
  within25  share of fits with lambda within +-25% of the truth (higher = better)
  rmsLog    RMS of log(lambda_hat / lambda_true) (lower = better)
  medRatio  median lambda_hat / lambda_true (1 = unbiased)
Usage (repo root): python checks/plot_noise_models.py
Writes analysis_outputs/noise_models/validation_summary.csv, noise_models.png/.pdf."""
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'noise_models')
V = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(d, 'validation_k0_*.csv'))])
V['r'] = V.lamHat / V.lamTrue
models = ['fixed', 'free', 'growN', 'prop', 'calib+floor', 'mix', 'mix+free', 'mix+prop', 'floor', 'floor+free', 'mixH+free', 'mixH+prop']
agg = dict(within25=('r', lambda x: np.mean((x > .75) & (x < 1.25))), rmsLog=('r', lambda x: np.sqrt(np.mean(np.log(x) ** 2))),
           medRatio=('r', 'median'), n=('r', 'size'))
S = V.groupby(['model', 'noise']).agg(**agg).reset_index()
P = V.groupby('model').agg(**agg).reset_index(); P['noise'] = 'all'
S = pd.concat([S, P]); S['model'] = pd.Categorical(S.model, models, ordered=True)
S = S.sort_values(['noise', 'model'], key=lambda c: c.astype(str) if c.name == 'noise' else c)
S.to_csv(os.path.join(d, 'validation_summary.csv'), index=False)
pd.set_option('display.width', 200)
print(S.pivot(index='model', columns='noise', values='within25').round(2).to_string()); print()
print(S.pivot(index='model', columns='noise', values='rmsLog').round(2).to_string()); print()
print(S.pivot(index='model', columns='noise', values='medRatio').round(2).to_string())

INK, MUTED = '#1f1e1c', '#898781'
COL = dict(zip(models, ['#d62a2a', '#8c564b', '#b8860b', '#7f7f7f', '#17becf', '#2a78d6', '#1baf7a', '#7a4fd6', '#eb6834', '#e377c2', '#000000', '#9467bd']))
fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
for i, k0 in enumerate(sorted(V.k0True.unique())):
    for j, nz in enumerate(sorted(V.noise.unique())):
        ax = axs[i, j]
        lt = np.array(sorted(V.lamTrue.unique()))
        ax.plot(lt, lt, color=MUTED, ls='--', lw=1)
        for m in ['fixed', 'free', 'mix+free', 'floor+free', 'mixH+free', 'mixH+prop']:
            g = V[(V.model == m) & (V.k0True == k0) & (V.noise == nz)].groupby('lamTrue').lamHat.median().reindex(lt)
            ax.plot(lt, g, 'o-', ms=3.5, lw=1.2, color=COL[m], label=m)
        ax.set_ylim(0, 8)
        ax.set_title(f'k0 {k0}, noise SD {nz}', fontsize=9, loc='left', color=INK)
        ax.set_xlabel('true lambda', fontsize=8)
        if j == 0:
            ax.set_ylabel('estimated lambda (median of 2 signs)', fontsize=8)
        ax.tick_params(labelsize=7, colors=MUTED)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
axs[0, 0].legend(fontsize=7, frameon=False)
fig.suptitle('Lambda recovery by noise / false-burst model (held-out simulations; dashed = truth)', fontsize=11, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'noise_models.{ext}'), dpi=140)
