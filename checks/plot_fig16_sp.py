"""PLOT_FIG16_SP  Recreates MS Fig. 16 (PP7 red vs PP7-FUS blue, negative burst
amplitude histograms for SP = 10..15, 4x and 30x RNA) without any fitting.
Differences from the MS figure: both conditions share ONE x-axis per panel
(the MS used separate top/bottom axes, red up to ~1600 vs blue up to ~300-500),
bins are common to both conditions, the last bin collects everything above the
pooled 97.5th percentile. Reviewer points: R2 4a (threshold dependence), FUS claim.
Also: median |amplitude| vs SP with granule-bootstrap 95% CIs (traces resampled).
Inputs: analysis_outputs/fig16_sp/amplitudes_sp.csv (checks/export_amplitudes_sp.m).
Usage (repo root): python checks/plot_fig16_sp.py [R2_reextracted|stored_MS] [Neg|Pos]
Writes analysis_outputs/fig16_sp/fig16_<data>_<sign>.png/.pdf and medians_<data>_<sign>.csv."""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = sys.argv[1] if len(sys.argv) > 1 else 'R2_reextracted'
sign = sys.argv[2] if len(sys.argv) > 2 else 'Neg'
d = os.path.join('analysis_outputs', 'fig16_sp')
A = pd.read_csv(os.path.join(d, 'amplitudes_sp.csv'))
A = A[(A.data == data) & (A.sign == sign)]
RED, BLUE, INK, MUTED = '#d62a2a', '#2a5bd6', '#1f1e1c', '#898781'
rng = np.random.default_rng(1)


def boot_median(sub, nb=1000):
    tr = sub.trace.unique()
    g = {t: sub[sub.trace == t].amp.values for t in tr}
    meds = [np.median(np.concatenate([g[t] for t in rng.choice(tr, len(tr))])) for _ in range(nb)]
    return np.percentile(meds, [2.5, 97.5])


fig = plt.figure(figsize=(13, 13))
outer = fig.add_gridspec(2, 2, width_ratios=[2.2, 1], hspace=0.25, wspace=0.25)
rows = []
for b, (base, fus) in enumerate([('4PP7', '4PP7-FUS'), ('30PP7', '30PP7-FUS')]):
    inner = outer[b, 0].subgridspec(3, 2, hspace=0.45, wspace=0.18)
    for k, SP in enumerate(range(10, 16)):
        ax = fig.add_subplot(inner[k // 2, k % 2])
        a1 = A[(A.dataset == base) & (A.SP == SP)]
        a2 = A[(A.dataset == fus) & (A.SP == SP)]
        top = np.percentile(np.concatenate([a1.amp, a2.amp]), 97.5)
        edges = np.linspace(0, top, 13)
        for sub, col, nm in [(a1, RED, base), (a2, BLUE, fus)]:
            x = np.clip(sub.amp.values, 0, top * 0.999)
            w = np.ones_like(x) / len(x)
            ax.hist(x, bins=edges, weights=w, color=col, alpha=0.45)
            ax.axvline(np.median(sub.amp), color=col, lw=1.2, ls='--')
            lo, hi = boot_median(sub)
            rows.append(dict(dataset=nm, SP=SP, n=len(sub), nTraces=sub.trace.nunique(),
                             median=np.median(sub.amp), medLo=lo, medHi=hi, mean=sub.amp.mean()))
        ax.set_title(f'SP={SP}   n {len(a1)} / {len(a2)}   median {np.median(a1.amp):.0f} / {np.median(a2.amp):.0f}',
                     fontsize=8, loc='left', color=INK)
        ax.tick_params(labelsize=7, colors=MUTED)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        if k % 2 == 0:
            ax.set_ylabel('Frequency', fontsize=8)
        if k >= 4:
            ax.set_xlabel('|amplitude| (a.u.; last bin = above)', fontsize=8)
    # median vs SP
    M = pd.DataFrame(rows)
    ax = fig.add_subplot(outer[b, 1])
    for nm, col, off in [(base, RED, -0.08), (fus, BLUE, 0.08)]:
        m = M[M.dataset == nm]
        ax.errorbar(m.SP + off, m['median'], yerr=[m['median'] - m.medLo, m.medHi - m['median']],
                    fmt='o-', color=col, ms=4, lw=1.2, capsize=2, label=nm)
    ax.set_xlabel('SP', fontsize=8); ax.set_ylabel('median |amplitude| (granule bootstrap 95% CI)', fontsize=8)
    ax.legend(fontsize=8, frameon=False)
    ax.tick_params(labelsize=7, colors=MUTED)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title(f'{"4x" if b == 0 else "30x"} RNA', fontsize=10, loc='left', color=INK)
lab = 'R2 re-extraction (red channel, multithresh spot, offset 200)' if data == 'R2_reextracted' else 'stored data used in the MS'
fig.suptitle(f'Fig. 16 recreated, {"negative" if sign == "Neg" else "positive"} bursts, {lab}\n'
             f'red = PP7, blue = PP7-FUS; one common x-axis per panel; dashed = medians; no fitting',
             fontsize=11, color=INK)
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'fig16_{data}_{sign}.{ext}'), dpi=140, bbox_inches='tight')
pd.DataFrame(rows).to_csv(os.path.join(d, f'medians_{data}_{sign}.csv'), index=False)
print(pd.DataFrame(rows).round(0).to_string(index=False))
