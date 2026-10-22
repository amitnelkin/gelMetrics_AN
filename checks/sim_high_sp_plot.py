"""SIM_HIGH_SP_PLOT  Plots checks/sim_high_sp.m: what raising SP costs when the
truth is known (RESULTS.md section 34).
Usage: python checks/sim_high_sp_plot.py
Writes analysis_outputs/detection_grid_highSP/sim_high_sp.png/.pdf."""
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'detection_grid_highSP'); S = pd.read_csv(os.path.join(d, 'sim_summary.csv'))
fig, axs = plt.subplots(1, 3, figsize=(18, 5))
for st, c in zip(S.setting.unique(), ['tab:blue', 'tab:red']):
    s = S[S.setting == st]
    axs[0].plot(s.SP, s.detected_per_granule_h, 'o-', color=c, label=f'{st}: with bursts')
    axs[0].plot(s.SP, s.false_per_granule_h, 'x--', color=c, label=f'{st}: burst-free twin')
    axs[0].axhline(s.true_losses_per_granule_h.iloc[0], color=c, ls=':', lw=1)
    axs[1].plot(s.SP, 100 * s.recall, 'o-', color=c, label=st)
    axs[2].plot(s.SP, s.median_detected_size_pct, 'o-', color=c, label=f'{st}: detected')
    axs[2].axhline(s.median_true_size_pct.iloc[0], color=c, ls=':', lw=1.5, label=f'{st}: true')
axs[0].set_yscale('log'); axs[0].set_ylabel('negative bursts per granule-hour (dotted: true losses)')
axs[1].set_ylabel('true losses found by the detector (%)'); axs[2].set_ylabel('median burst size (% of granule)'); axs[2].set_yscale('log')
for ax, t in zip(axs, ['(a) detected bursts: with real bursts vs burst-free', '(b) how many true bursts are found', '(c) size of what is found vs true size']):
    ax.set_xlabel('stringency SP'); ax.set_title(t, fontsize=10, loc='left'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(d, 'sim_high_sp.png'), dpi=130); fig.savefig(os.path.join(d, 'sim_high_sp.pdf'))
