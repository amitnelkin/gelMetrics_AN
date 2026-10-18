"""Step 1 figure: 4PP7 traces from the original pipeline, each point coloured by
the TIFF page it was read from (red-fluorescence page vs DIC page)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

O = r'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\zigzag_scratch\out'
P = pd.read_csv(O + r'\points.csv', header=None,
                names=['trace', 'movie', 'traj', 'frame', 'page', 'x', 'y', 'spot', 'cell', 'nP'])
RED, DIC, INK, MUTED = '#eb6834', '#2a78d6', '#1f1e1c', '#898781'
# (my trace id, stored 4PP7 trace id, title)
panels = [(3, 3, '005 (1 channel): every page is red - clean'),
          (24, 53, '007 (2 ch), Mosaic steps of 1: alternates red/DIC pages - zigzag'),
          (23, 52, '007 (2 ch), even frames only: every point read from a DIC page'),
          (36, 65, '007 (2 ch), odd frames only: detected on DIC plane, read on red page')]
fig, axs = plt.subplots(4, 1, figsize=(9, 10), sharex=True)
for ax, (t, s, title) in zip(axs, panels):
    g = P[P.trace == t]
    dic = (g.movie == 2) & (g.page % 2 == 0)
    i = np.arange(1, len(g) + 1)
    ax.plot(i, g.spot, color=MUTED, lw=0.6, zorder=1)
    ax.scatter(i[~dic.values], g.spot[~dic.values], s=8, color=RED, label='read from red page', zorder=2)
    ax.scatter(i[dic.values], g.spot[dic.values], s=8, color=DIC, label='read from DIC page', zorder=2)
    ax.set_title(f'stored 4PP7 trace {s}  |  {title}', fontsize=9, color=INK, loc='left')
    ax.set_ylabel('spot intensity (a.u.)', fontsize=8, color=INK)
    ax.tick_params(labelsize=7, colors=MUTED)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.grid(axis='y', color='#e6e4df', lw=0.5)
axs[0].legend(fontsize=8, frameon=False, loc='upper right')
axs[-1].set_xlabel('point index in trace (what analysis.m treats as frame number, 10 s)', fontsize=8, color=INK)
fig.suptitle('Original pipeline on 4PP7 movies 005 and 007 (re-extracted traces = stored traces, bit-exact)',
             fontsize=10, color=INK)
fig.tight_layout()
fig.savefig(O + r'\zigzag_mechanism.png', dpi=150)
print('saved')
