"""MS_FIGURES  Figures for the revised manuscript, in the style of the original
(MATLAB look: boxed axes, inward ticks, bold (a) labels, salmon-red histograms,
blue open circles, red dashed lines). Reviewer points: R2 4a/4b/5/6 and the
minor point on Fig. 4 colours (Okabe-Ito for Fig. 4).

Usage (repo root):  py -3 checks/ms_figures.py [--vitro-new]
--vitro-new (RESULTS.md sections 39, 43): the in vitro panels use the stored granules re-read
with the corrected pipeline (vitro_same/, detection_grid_vitro_same*/; controls from
noise_background_vitro/, steady/static_control_vitro_summary.csv) instead of the stored traces,
[--vitro-spotsum: the same granules with the spot sum, the user's choice (section 53):
vitro_spotsum/{summary,bursts,sp_grid}.csv, steady/static_control_vitro_sum_summary.csv;
output figs_vitro_spotsum/]
with movie-level CIs;
only Figs. 6, 18 and 22 (manuscript Figs. 6, 14, 18) are drawn, into figs_vitro_new/.
Writes analysis_outputs/ms_revision/figs/fig*.png (300 dpi) and .pdf.

Inputs (all committed): ms_revision/fig4_*.csv (checks/ms_fig4_sim.m),
ms_grid_w13sp12/ (checks/detection_grid.m, window 13, SP 12), control_table/,
gap_test/, ms_sim_test/, lambda_k0_shape_notrim/, detection_grid*/,
robustness/, noise_background/, steady/, slow_changes*/, halves/, fig16_sp/,
ms_fig15_R2/ (grid fit on the re-extracted data), in vivo/in vitro .mat files.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import chi2_contingency

A = 'analysis_outputs'
OUT = os.path.join(A, 'ms_revision', 'figs')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3,
    'legend.frameon': False, 'legend.fontsize': 7, 'axes.labelsize': 8.5, 'mathtext.default': 'regular',
})
HIST = '#FF6666'          # salmon-red bars, as in the original figures
BLUE = '#0072BD'          # MATLAB blue (open circles)
COPY = '0.25'             # outline for burst-free copies
# condition colours: Okabe-Ito (vermillion / blue / bluish green), as Qb / PP7 / PP7-FUS in old Fig. 10(b)
CQB, CPP, CFUS = '#D55E00', '#0072B2', '#009E73'
OI = ['#0072B2', '#D55E00', '#009E73', '#E69F00', '#56B4E9', '#CC79A7', '#F0E442']

VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
VITRO = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
LAB = {'4PP7': 'PP7\nQCP-5x/PCP-4x', '4PP7-FUS': 'PP7-FUS\nQCP-5x/PCP-4x', '30PP7': 'PP7\nPCP-30x',
       '30PP7-FUS': 'PP7-FUS\nPCP-30x', '5Qb': 'Qβ\nQCP-5x/PCP-4x', '10Qb': 'Qβ\nQCP-10x',
       'protein_RNA_1_1': '1:1', 'protein_RNA_10_1': '10:1', 'protein_RNA_100_1': '100:1'}
SHORT = {'4PP7': 'PP7, 4x', '4PP7-FUS': 'PP7-FUS, 4x', '30PP7': 'PP7, 30x', '30PP7-FUS': 'PP7-FUS, 30x',
         '5Qb': 'Qβ, 5x', '10Qb': 'Qβ, 10x', 'protein_RNA_1_1': 'in vitro 1:1',
         'protein_RNA_10_1': 'in vitro 10:1', 'protein_RNA_100_1': 'in vitro 100:1'}
COL = {'4PP7': CPP, '30PP7': CPP, '4PP7-FUS': CFUS, '30PP7-FUS': CFUS, '5Qb': CQB, '10Qb': CQB}


def label(ax, s, x=-0.16, y=1.04):
    ax.text(x, y, f'({s})', transform=ax.transAxes, fontweight='bold', fontsize=8.5, va='bottom')


def box(ax):
    for sp in ax.spines.values():
        sp.set_visible(True)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300)
    fig.savefig(os.path.join(OUT, name + '.pdf'))
    plt.close(fig)
    print('wrote', name)


CT = pd.read_csv(os.path.join(A, 'control_table', 'control_table.csv')).set_index('dataset')
GB = pd.read_csv(os.path.join(A, 'ms_grid_w13sp12', 'bursts.csv'))
GB = GB[(GB.levelBefore > 0) & (GB.amp < 0)].copy()
GB['f'] = 100 * GB.amp.abs() / GB.levelBefore
import sys
VSUM = '--vitro-spotsum' in sys.argv
VNEW = '--vitro-new' in sys.argv or VSUM
if VSUM:
    OUT = os.path.join(A, 'ms_revision', 'figs_vitro_spotsum'); os.makedirs(OUT, exist_ok=True)
    vs = pd.read_csv(os.path.join(A, 'vitro_spotsum', 'summary.csv')).set_index('dataset').loc[VITRO]
    for c in ['meanBurst', 'turnover', 'net']:
        for suf in ['', '_lo', '_hi']:
            CT.loc[VITRO, c + suf] = vs[c + suf].values
    CT.loc[VITRO, 'sizeReal'] = vs.sizeGranules.values; CT.loc[VITRO, 'sizeCopy'] = vs.sizeCopies.values
    gv = pd.read_csv(os.path.join(A, 'vitro_spotsum', 'bursts.csv')); gv['f'] = 100 * gv.f
    GB = pd.concat([GB[~GB.dataset.isin(VITRO)], gv])
elif VNEW:
    OUT = os.path.join(A, 'ms_revision', 'figs_vitro_new'); os.makedirs(OUT, exist_ok=True)
    vn = pd.read_csv(os.path.join(A, 'vitro_same', 'summary.csv')).set_index('dataset').loc[VITRO]
    for new, old in [('meanBurst', 'mean_burst'), ('turnover', 'turn0'), ('net', 'net0')]:
        for suf in ['', '_lo', '_hi']:
            CT.loc[VITRO, new + suf] = vn[old + suf].values
    CT.loc[VITRO, 'sizeReal'] = vn.size_real.values; CT.loc[VITRO, 'sizeCopy'] = vn.size_copy.values
    gv = pd.read_csv(os.path.join(A, 'detection_grid_vitro_same', 'bursts.csv'))
    gv = gv[(gv.levelBefore > 0) & (gv.amp < 0) & gv.dataset.isin(VITRO)].copy(); gv['f'] = 100 * gv.amp.abs() / gv.levelBefore
    GB = pd.concat([GB[~GB.dataset.isin(VITRO)], gv])


# ---------------------------------------------------------------- Fig. 4
def fig4():
    amp = pd.read_csv(os.path.join(A, 'ms_revision', 'fig4_amplitudes.csv'))
    tr = pd.read_csv(os.path.join(A, 'ms_revision', 'fig4_traces.csv'))
    fig, ax = plt.subplots(2, 3, figsize=(6.5, 4.0), gridspec_kw={'height_ratios': [1, 1.15]})
    cb, cbl, cn = OI[0], OI[1], OI[2]    # base / + bleaching / + noise
    for k, lab_t, lab_h in [(1, 'a', 'b'), (2, 'c', 'd'), (3, 'e', 'f')]:
        t = tr[tr.type == k]
        a = ax[0, k - 1]
        a.plot(t.t, t.noisy, color=cn, lw=0.5, label='+Gaussian noise')
        a.plot(t.t, t.base, color=cb, lw=0.9, label='Base signal')
        a.plot(t.t, t.bleach, color=cbl, lw=0.9, ls='--', label='+Exponential component')
        a.set_xlim(0, 60); a.set_xlabel('Time (min)')
        if k == 1: a.set_ylabel('Intensity (a.u.)')
        label(a, lab_t)
        x = amp[amp.type == k].amp
        h = ax[1, k - 1]
        bins = np.arange(-600, 601, 10)
        h.hist(x[x > 0], bins=bins, color=OI[2], alpha=0.75, label='Positive bursts')
        h.hist(x[x < 0], bins=bins, color=OI[1], alpha=0.75, label='Negative bursts')
        h.set_xlim(-600, 600); h.set_xlabel('Amplitudes (a.u.)')
        if k == 1: h.set_ylabel('Count')
        label(h, lab_h)
    h_, l_ = ax[0, 0].get_legend_handles_labels()
    order = [1, 2, 0]
    fig.legend([h_[i] for i in order], [l_[i] for i in order], loc='upper center', ncol=3, fontsize=7)
    ax[1, 0].legend(loc='upper left', fontsize=6.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94], h_pad=1.5, w_pad=1.2)
    save(fig, 'fig4')
    c = pd.read_csv(os.path.join(A, 'ms_revision', 'fig4_counts.csv')).set_index('type')
    p = chi2_contingency([[c.loc[1, 'nPos'], c.loc[1, 'nNeg']], [c.loc[2, 'nPos'], c.loc[2, 'nNeg']]])[1]
    tot = chi2_contingency([[c.loc[1, 'nPos'] + c.loc[1, 'nNeg'], 1000 * 360], [c.loc[2, 'nPos'] + c.loc[2, 'nNeg'], 1000 * 360]])[1]
    print('Fig 4 counts:', c.to_dict('index'), 'chi2 p (pos/neg split)', p, 'chi2 p (event totals)', tot)


# ---------------------------------------------------------------- Fig. 5
def fig5():
    amp = pd.read_csv(os.path.join(A, 'ms_sim_test', 'amplitudes.csv'))
    fits = pd.read_csv(os.path.join(A, 'ms_sim_test', 'fits.csv'))
    shp = pd.read_csv(os.path.join(A, 'lambda_k0_shape_notrim', 'fits.csv'))
    hs = pd.read_csv(os.path.join(A, 'detection_grid_highSP', 'sim_summary.csv'))
    fig, ax = plt.subplots(2, 3, figsize=(6.5, 3.9))
    for j, lam in enumerate([1, 2, 4]):
        a = ax[0, j]
        x = amp[(amp['lambda'] == lam) & (amp.noise == 40) & (amp.rep == 1)].amp
        a.hist(x, bins=np.arange(0, 1001, 10), color=HIST, density=True)
        a.set_xlim(0, 1000); a.set_ylim(0, 0.012); a.set_yticks([0, 0.005, 0.01]); a.set_xlabel('Amplitude (a.u.)')
        if j == 0: a.set_ylabel('Frequency')
        a.text(0.95, 0.85, f'$\\lambda_{{sim}}$={lam}', transform=a.transAxes, ha='right', fontweight='bold', fontsize=9)
        label(a, 'abc'[j])
    # (d) lambda chosen by the WS grid score and by the likelihood, vs the truth (SNR 1.5)
    a = ax[1, 0]
    d = fits[(fits.source == 'detected') & (fits.noise == 40)]
    for sc, mk, cl, nm in [('WS', 'o', BLUE, 'WS score'), ('NLL', 's', OI[1], 'likelihood')]:
        g = d[d.score == sc].groupby('lambdaTrue').lam
        m = g.median(); lo = g.min(); hi = g.max()
        a.errorbar(m.index, m.values, yerr=[m - lo, hi - m], fmt=mk, mfc='none', color=cl, ms=4, lw=0.8, capsize=2, label=nm)
    a.plot([0, 5.5], [0, 5.5], 'r--', lw=1, label='truth')
    a.set_xlim(0, 5.5); a.set_ylim(0, 5.5); a.set_xlabel('$\\lambda_{sim}$'); a.set_ylabel('$\\lambda_{est}$')
    a.legend(loc='upper left'); label(a, 'd')
    # (e) mean burst size lambda*k0 (shape fit) vs truth 60*lambda, noise 40
    a = ax[1, 1]
    s = shp[shp.kind != 'in vivo'].copy()
    s = s[s.dataset.str.contains('noise 40')]
    s['lamTrue'] = s.dataset.str.extract(r'lambda ([\d.]+)')[0].astype(float)
    tru = 60 * s['lamTrue']
    a.errorbar(tru, s['mean'], yerr=[s['mean'] - s.mean_lo, s.mean_hi - s['mean']], fmt='o', mfc='none', color=BLUE, ms=4, lw=0.8, capsize=2)
    a.plot([0, 330], [0, 330], 'r--', lw=1)
    a.set_xlim(0, 330); a.set_ylim(0, 330)
    a.set_xlabel('true $\\lambda\\cdot k_0$ (a.u.)'); a.set_ylabel('estimated $\\lambda\\cdot k_0$ (a.u.)'); label(a, 'e')
    # (f) detection limit: share of true bursts found vs SP, at manuscript-like and PP7-like noise
    a = ax[1, 2]
    from matplotlib.lines import Line2D
    for nm, cl, mk in [('MS-like', BLUE, 'o'), ('PP7-like', OI[1], 's')]:
        g = hs[hs.setting.str.startswith(nm)].sort_values('SP')
        a.plot(g.SP, 100 * g.recall, '-' + mk, color=cl, mfc='none', ms=4, lw=0.8)
        a.plot(g.SP, 100 * g.false_over_detected.clip(upper=1), ':', color=cl, lw=1)
    a.legend(handles=[Line2D([], [], color=BLUE, lw=1.5, label='original noise'),
                      Line2D([], [], color=OI[1], lw=1.5, label='PP7 noise'),
                      Line2D([], [], color='k', lw=0.8, marker='o', mfc='none', ms=3, label='bursts found'),
                      Line2D([], [], color='k', lw=1, ls=':', label='false share')], loc='upper right', fontsize=5.5, handlelength=1.8)
    a.set_xlabel('Stringency parameter j'); a.set_ylabel('%'); a.set_ylim(0, 105); a.set_yticks([0, 25, 50, 75, 100])
    label(a, 'f')
    fig.tight_layout(h_pad=1.2)
    save(fig, 'fig5')


# ---------------------------------------------------------------- burst-size panels (Figs. 6-9)
def burst_hist(a, ds, xmax=60, ytxt=True):
    g = GB[GB.dataset == ds]
    r, c = g[g.kind == 'real'].f, g[g.kind == 'surr'].f
    bins = np.arange(0, xmax + 0.1, 2.5)
    a.hist(r.clip(upper=xmax), bins=bins, color=HIST, density=True, label='granules')
    a.hist(c.clip(upper=xmax), bins=bins, histtype='step', color=COPY, lw=1, density=True, label='burst-free copies')
    a.set_xlim(0, xmax); a.set_xlabel('Burst size (% of granule)')
    if ytxt: a.set_ylabel('Frequency')
    row = CT.loc[ds]
    a.text(0.97, 0.95, f"median {row.sizeReal:.1f}% vs {row.sizeCopy:.1f}%\n"
           f"turnover beyond noise\n{row.turnover:.1f} [{row.turnover_lo:.1f}, {row.turnover_hi:.1f}] %/h",
           transform=a.transAxes, ha='right', va='top', fontsize=6.8, fontweight='bold')
    return len(r)


def ci_bars(a, keys, col, lo, hi, colors, ylab, zero=True):
    y = CT.loc[keys, col].values
    e = [y - CT.loc[keys, lo].values, CT.loc[keys, hi].values - y]
    x = np.arange(len(keys))
    a.bar(x, y, color=colors, width=0.6, alpha=0.85)
    a.errorbar(x, y, yerr=e, fmt='none', ecolor='k', lw=0.8, capsize=2)
    if zero: a.axhline(0, color='k', lw=0.5)
    a.set_xticks(x); a.set_xticklabels([LAB[k] for k in keys], fontsize=6.5); a.set_ylabel(ylab)
    a.set_xlim(-0.6, len(keys) - 0.4)


def fig6():
    fig = plt.figure(figsize=(6.5, 3.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])
    for j, ds in enumerate(VITRO):
        a = fig.add_subplot(gs[0, j]); n = burst_hist(a, ds, xmax=60, ytxt=(j == 0))
        a.set_title(f'protein:RNA {LAB[ds]}  (n = {n})', fontsize=7.5)
        label(a, 'abc'[j])
        if j == 0: a.legend(loc='center right', fontsize=6)
    cols = [HIST] * 3
    for j, (c, lo, hi, yl) in enumerate([('meanBurst', 'meanBurst_lo', 'meanBurst_hi', 'Mean burst size, $\\lambda\\cdot k_0$ (a.u.)'),
                                         ('turnover', 'turnover_lo', 'turnover_hi', 'Turnover beyond noise (%/h)'),
                                         ('net', 'net_lo', 'net_hi', 'Net change beyond\nbleaching (%/h)')]):
        a = fig.add_subplot(gs[1, j]); ci_bars(a, VITRO, c, lo, hi, cols, yl)
        a.set_xlabel('protein:RNA'); label(a, 'def'[j])
    fig.tight_layout(h_pad=1.0)
    save(fig, 'fig6')


def fig_pair(keys, name, extra=None):
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.3))
    for j, ds in enumerate(keys):
        n = burst_hist(ax[j], ds, xmax=60, ytxt=(j == 0))
        ax[j].set_title(f"{LAB[ds].replace(chr(10), ' + ')}  (n = {n})", fontsize=7.5)
        label(ax[j], 'ab'[j], x=-0.1)
    ax[0].legend(loc='center right', fontsize=6.5)
    fig.tight_layout()
    save(fig, name)


def fig9():
    fig = plt.figure(figsize=(6.5, 4.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1])
    for j, ds in enumerate(['4PP7-FUS', '30PP7-FUS']):
        a = fig.add_subplot(gs[0, j]); n = burst_hist(a, ds, xmax=60, ytxt=(j == 0))
        a.set_title(f"{LAB[ds].replace(chr(10), ' + ')}  (n = {n})", fontsize=7.5); label(a, 'ab'[j], x=-0.1)
        if j == 0: a.legend(loc='center right', fontsize=6.5)
    # (c) granule brightness, 4x and 30x pairs
    a = fig.add_subplot(gs[1, 0])
    keys = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS']
    y = CT.loc[keys, 'brightness'].values / 1000
    e = [y - CT.loc[keys, 'brightness_lo'].values / 1000, CT.loc[keys, 'brightness_hi'].values / 1000 - y]
    x = np.array([0, 1, 2.6, 3.6])
    a.bar(x, y, color=[CPP, CFUS, CPP, CFUS], width=0.7, alpha=0.85)
    a.errorbar(x, y, yerr=e, fmt='none', ecolor='k', lw=0.8, capsize=2)
    a.set_xticks(x); a.set_xticklabels(['PP7', 'PP7-FUS', 'PP7', 'PP7-FUS'], fontsize=7)
    a.text(0.5, -0.22, 'QCP-5x/PCP-4x', transform=a.get_xaxis_transform(), ha='center', fontsize=7)
    a.text(3.1, -0.22, 'PCP-30x', transform=a.get_xaxis_transform(), ha='center', fontsize=7)
    a.set_ylabel('Granule brightness (a.u./1000)')
    a.text(0.5, 0.93, 'p = 0.002', transform=a.get_xaxis_transform(), ha='center', fontsize=7, fontweight='bold')
    a.text(3.1, 0.93, 'p = 0.55', transform=a.get_xaxis_transform(), ha='center', fontsize=7)
    label(a, 'c', x=-0.1)
    # (d) burst-size gap PP7-FUS minus PP7: total, explained by noise, beyond noise
    a = fig.add_subplot(gs[1, 1])
    T = pd.read_csv(os.path.join(A, 'gap_test', 'tests.csv'))
    T = T[T.stat == 'median'].set_index('compare')
    for i, (cmp, lab) in enumerate([('4PP7-FUS vs 4PP7', 'QCP-5x/PCP-4x'), ('30PP7-FUS vs 30PP7', 'PCP-30x')]):
        r = T.loc[cmp]
        for k, (v, lo, hi, cl, nm) in enumerate([(r.gap, r.gap_lo, r.gap_hi, '0.55', 'total'),
                                                   (r.noise_part, r.noise_lo, r.noise_hi, '0.85', 'explained by noise'),
                                                   (r.beyond, r.beyond_lo, r.beyond_hi, CFUS, 'beyond noise')]):
            xx = i * 1.1 + (k - 1) * 0.28
            a.bar(xx, v, width=0.26, color=cl, edgecolor='k', lw=0.4, label=nm if i == 0 else None)
            a.errorbar(xx, v, yerr=[[v - lo], [hi - v]], fmt='none', ecolor='k', lw=0.8, capsize=2)
        ptxt = (f'p = {r.p_one_sided:.3f} (one-sided)\np = {r.p_two_sided:.3f} (two-sided)' if i == 0
                else f'p = {r.p_two_sided:.2f}')
        a.text(i * 1.1, 0.70, ptxt, transform=a.get_xaxis_transform(), ha='center', fontsize=6.5, fontweight='bold' if i == 0 else 'normal')
    a.axhline(0, color='k', lw=0.5)
    a.set_xticks([0, 1.1]); a.set_xticklabels(['QCP-5x/PCP-4x', 'PCP-30x'], fontsize=7)
    a.set_ylabel('PP7-FUS − PP7, median burst\nsize (% points)'); a.set_ylim(-4, 17)
    a.legend(loc='upper center', fontsize=6, ncol=3, columnspacing=0.8); label(a, 'd', x=-0.1)
    fig.tight_layout(h_pad=1.2)
    save(fig, 'fig9')


# ---------------------------------------------------------------- Fig. 10 (b-d): summary of all in vivo conditions
def fig10():
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.2))
    keys = VIVO
    cols = [COL[k] for k in keys]
    for j, (c, lo, hi, yl) in enumerate([('meanBurst', 'meanBurst_lo', 'meanBurst_hi', 'Mean burst size,\n$\\lambda\\cdot k_0$ (a.u.)'),
                                         ('turnover', 'turnover_lo', 'turnover_hi', 'Turnover beyond\nnoise (%/h)'),
                                         ('net', 'net_lo', 'net_hi', 'Net change beyond\nbleaching (%/h)')]):
        a = ax[j]
        y = CT.loc[keys, c].values; e = [y - CT.loc[keys, lo].values, CT.loc[keys, hi].values - y]
        x = np.array([0, 1, 2.4, 3.4, 4.8, 5.8])
        a.bar(x, y, color=cols, width=0.75, alpha=0.85)
        a.errorbar(x, y, yerr=e, fmt='none', ecolor='k', lw=0.8, capsize=2)
        a.axhline(0, color='k', lw=0.5)
        a.set_xticks(x); a.set_xticklabels(['4x', '4x', '30x', '30x', '5x', '10x'], fontsize=7)
        a.set_ylabel(yl); label(a, 'bcd'[j], x=-0.22)
    from matplotlib.patches import Patch
    ax[2].legend(handles=[Patch(color=CPP, label='PP7'), Patch(color=CFUS, label='PP7-FUS'), Patch(color=CQB, label='Qβ')],
                 loc='upper left', fontsize=6.5)
    fig.tight_layout(w_pad=1.5)
    save(fig, 'fig10bcd')


# ---------------------------------------------------------------- Appendix: Fig. 15 (grid method on re-extracted data)
def fig15():
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.4))
    for i, ds in enumerate(VIVO):
        R = sio.loadmat(os.path.join(A, 'ms_fig15_R2', f'{ds}_run.mat'), squeeze_me=True, struct_as_record=False)
        fits = {str(e.title): e for e in np.atleast_1d(R['fits'])}
        e = fits['Negative amplitudes']
        W = np.atleast_2d(e.w_score); K = np.atleast_1d(e.K0Vector); L = np.atleast_1d(getattr(e, 'lambda'))
        ws = W.min(axis=1); k0 = K[W.argmin(axis=1)]
        ls = '-' if ds.startswith('4') or ds == '5Qb' else '--'
        ax[0].plot(L, ws, ls + 'o', color=COL[ds], mfc='none', ms=3.5, lw=0.9, label=SHORT[ds])
        ax[1].plot(L, k0, ls + 'o', color=COL[ds], mfc='none', ms=3.5, lw=0.9)
    ax[0].set_xlabel('$\\lambda_{est}$'); ax[0].set_ylabel('WS'); ax[1].set_xlabel('$\\lambda_{est}$'); ax[1].set_ylabel('$k_0$ (a.u.)')
    ax[0].legend(fontsize=6, ncol=2); label(ax[0], 'a', x=-0.13); label(ax[1], 'b', x=-0.13)
    fig.tight_layout()
    save(fig, 'fig15')


# ---------------------------------------------------------------- Appendix: Fig. 16 (PP7 vs PP7-FUS across SP, re-extracted)
def fig16():
    d = pd.read_csv(os.path.join(A, 'fig16_sp', 'amplitudes_sp.csv'))
    d = d[(d.data == 'R2_reextracted') & (d['sign'] == 'Neg')]
    SPs = [10, 12, 14]
    fig, ax = plt.subplots(2, 3, figsize=(6.5, 3.6), sharex=True)
    for r, (a_, b_, nm) in enumerate([('4PP7', '4PP7-FUS', 'QCP-5x/PCP-4x'), ('30PP7', '30PP7-FUS', 'PCP-30x')]):
        for c, sp in enumerate(SPs):
            a = ax[r, c]
            bins = np.arange(0, 1201, 50)
            for ds, cl, lab in [(a_, '#FF0000', 'PP7'), (b_, '#0000FF', 'PP7-FUS')]:
                x = d[(d.dataset == ds) & (d.SP == sp)].amp
                a.hist(x.clip(upper=1200), bins=bins, color=cl, alpha=0.45, density=True, label=f'{lab} (n={len(x)})')
            a.set_xlim(0, 1200)
            if r == 1: a.set_xlabel('Amplitude (a.u.)')
            if c == 0: a.set_ylabel(f'{nm}\nFrequency')
            if r == 0: a.set_title(f'SP = {sp}', fontsize=8)
            a.legend(fontsize=5.5, loc='upper right')
    fig.tight_layout()
    save(fig, 'fig16')


# ---------------------------------------------------------------- Appendix: Fig. 17 (two-channel movies)
def fig17():
    def cells(p):
        m = sio.loadmat(p); k = [x for x in m if not x.startswith('__')][0]
        return [np.ravel(c).astype(float) for c in m[k].ravel()[1:]]
    old = cells('in_vivo_experiments/4PP7_data.mat')
    new = cells('in_vivo_reextracted/after_R2_offset/4PP7_data.mat')
    t_old = old[48]                                      # stored trace 49 (movie 007, two-channel)
    t_new = max(new, key=len)
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.0), gridspec_kw={'width_ratios': [3, 1.2]})
    ax[0].plot(np.arange(t_old.size) / 6, t_old, color=HIST, lw=0.5, label='original extraction (two-channel movie)')
    ax[0].plot(np.arange(t_new.size) / 6, t_new + 200, color=BLUE, lw=0.6, label='corrected extraction (fluorescence channel only)')
    ax[0].set_xlabel('Time (min)'); ax[0].set_ylabel('Spot intensity (a.u.)'); ax[0].legend(fontsize=6, loc='center right', framealpha=1, frameon=True)
    ax[0].set_xlim(0, 60); label(ax[0], 'a', x=-0.08)
    ax[1].plot(np.arange(16), t_old[:16], 'o-', color=HIST, ms=3, lw=0.8)
    ax[1].set_xlabel('Time point'); ax[1].set_ylabel('Spot intensity (a.u.)'); label(ax[1], 'b', x=-0.3)
    fig.tight_layout()
    save(fig, 'fig17')


# ---------------------------------------------------------------- Appendix: Fig. 18 (noise floor and its source)
def fig18():
    nf = pd.read_csv(os.path.join(A, 'robustness', 'noise_floor_allSP.csv'))
    if VSUM:
        b = pd.read_csv(os.path.join(A, 'vitro_spotsum', 'sp_grid.csv')); b['measure'] = 'mean'; b['FDR'] = b.falseNeg / b.realNeg
        nf = pd.concat([nf[~nf.dataset.isin(VITRO)], b])
    elif VNEW:
        b = pd.read_csv(os.path.join(A, 'detection_grid_vitro_same_sp', 'bursts.csv'))
        b = b[b.dataset.isin(VITRO)].groupby(['dataset', 'SP', 'kind']).size().unstack(fill_value=0).reset_index()
        b = pd.DataFrame(dict(dataset=b.dataset, measure='mean', SP=b.SP, realNeg=b.real, falseNeg=b.surr, FDR=b.surr / b.real))
        nf = pd.concat([nf[~nf.dataset.isin(VITRO)], b])
    nb = pd.read_csv(os.path.join(A, 'noise_background', 'summary.csv')).set_index('dataset')
    st = pd.read_csv(os.path.join(A, 'steady', 'static_control_summary.csv'))
    if VNEW:
        st = pd.concat([st, pd.read_csv(os.path.join(A, 'steady', 'static_control_vitro_sum_summary.csv' if VSUM else 'static_control_vitro_summary.csv'))])
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.3)) if not VNEW else plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={'width_ratios': [1, 1.5, 1.05]})
    a = ax[0]
    for i, ds in enumerate(VIVO + VITRO):
        m = 'spotsum' if ds in VIVO else 'mean'
        g = nf[(nf.dataset == ds) & (nf.measure == m) & (nf.SP <= 20)].sort_values('SP')
        c = COL.get(ds, '0.4'); ls = '-' if ds in VIVO else ':'
        a.plot(g.SP, g.FDR, ls, color=c, lw=0.9, marker='o', mfc='none', ms=2.5, label=SHORT[ds])
    a.axhline(1, color='k', lw=0.5, ls='--')
    a.set_xlabel('Stringency parameter j'); a.set_ylabel('Bursts in copies ÷ bursts in granules'); a.set_ylim(0, 2.6); a.set_xlim(9.5, 20.5)
    a.legend(fontsize=5, ncol=2, loc='upper left'); label(a, 'a')
    a = ax[1]
    NB = VIVO
    if VNEW:
        nb = pd.concat([nb, pd.read_csv(os.path.join(A, 'noise_background_vitro', 'summary.csv')).set_index('dataset')]); NB = VIVO + VITRO
    gcol = 'ratioSum' if VSUM else 'ratioMean'   # --vitro-spotsum: spot-sum ratios for every dataset (the user's choice, section 55)
    r = [float(str(nb.loc[ds, gcol]).split()[0]) for ds in NB]
    rc = [float(str(nb.loc[ds, 'ratioCell']).split()[0]) for ds in NB]
    x = np.arange(len(NB))
    a.bar(x - 0.18, r, 0.36, color=[COL.get(d, HIST) for d in NB], alpha=0.85, label='granule (spot sum)' if VSUM else 'granule (spot mean)')
    a.bar(x + 0.18, rc, 0.36, color='0.75', label='cell background' if not VNEW else 'cell / dilute-phase background')
    a.axhline(1, color='k', lw=0.5, ls='--')
    a.set_xticks(x); a.set_xticklabels([SHORT[d].replace(', ', '\n').replace('in vitro ', 'vitro\n') for d in NB], fontsize=5.5 if not VNEW else 5)
    a.set_ylabel('Measured ÷ camera-predicted noise' if not VNEW else 'Measured ÷ camera-' + chr(10) + 'predicted noise'); a.legend(fontsize=5.5, loc='upper left'); label(a, 'b')
    a = ax[2]
    s = st[st.version.str.startswith('per-frame')].set_index('dataset')
    keys = [k for k in (VIVO + VITRO if VNEW else VIVO) if k in s.index]
    x = np.arange(len(keys))
    a.bar(x - 0.18, s.loc[keys, 'real_excess_w13'], 0.36, color=[COL.get(k, HIST) for k in keys], alpha=0.85, label='granules')
    a.bar(x + 0.18, s.loc[keys, 'static_excess_w13'], 0.36, color='0.75', label='static synthetic granule')
    a.errorbar(x + 0.18, s.loc[keys, 'static_excess_w13'], yerr=1.96 * s.loc[keys, 'static_excess_se_w13'], fmt='none', ecolor='k', lw=0.7, capsize=2)
    a.axhline(0, color='k', lw=0.5)
    a.set_xticks(x); a.set_xticklabels([SHORT[k].replace(', ', '\n').replace('in vitro ', 'vitro\n') for k in keys], fontsize=5.5 if not VNEW else 5)
    a.set_ylabel('Turnover beyond noise (%/h)'); a.legend(fontsize=5.5, loc='upper left'); label(a, 'c')
    fig.tight_layout()
    save(fig, 'fig18')


# ---------------------------------------------------------------- Appendix: Fig. 19 (detection settings)
def fig19():
    g = pd.read_csv(os.path.join(A, 'detection_grid', 'grid.csv'))
    g = g[g.dataset.isin(VIVO)]
    med = g.groupby(['SP', 'window']).FDR.median().unstack()
    hs = pd.read_csv(os.path.join(A, 'detection_grid_highSP', 'summary.csv'))
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.4), gridspec_kw={'width_ratios': [1, 1.15, 0.9]})
    a = ax[0]
    im = a.imshow(med.values, origin='lower', cmap='viridis', vmin=0.6, vmax=1.0, aspect='auto')
    a.set_xticks(range(med.shape[1])); a.set_xticklabels(med.columns); a.set_yticks(range(med.shape[0])); a.set_yticklabels(med.index)
    a.set_xlabel('Moving-average window (frames)'); a.set_ylabel('Stringency parameter j')
    for i in range(med.shape[0]):
        for j in range(med.shape[1]):
            a.text(j, i, f'{med.values[i, j]:.2f}', ha='center', va='center', fontsize=5.5, color='w')
    fig.colorbar(im, ax=a, label='copies ÷ granules (median, in vivo)'); label(a, 'a')
    a = ax[1]
    print('highSP summary columns:', list(hs.columns))
    sub = hs.copy()
    rc = [c for c in sub.columns if c.lower().startswith('real')]
    cc = [c for c in sub.columns if c.lower().startswith('copy') or c.lower().startswith('surr') or c.lower().startswith('false')]
    for ds in VIVO:
        d = sub[sub.dataset == ds].sort_values('SP')
        if d.empty or not rc or not cc: continue
        a.plot(d.SP, d[rc[0]], '-o', color=COL[ds], mfc='none', ms=3, lw=0.9, label=SHORT[ds])
        a.plot(d.SP, d[cc[0]], ':', color=COL[ds], lw=0.9)
    a.set_yscale('log'); a.set_xticks([10, 12, 14, 16, 18, 20, 22, 25]); a.set_xlabel('Stringency parameter j'); a.set_ylabel('Negative bursts detected')
    a.text(0.97, 0.95, 'solid: granules\ndotted: burst-free copies', transform=a.transAxes, ha='right', va='top', fontsize=6)
    a.legend(fontsize=5, loc='lower left', ncol=1, handlelength=1.5, labelspacing=0.2); label(a, 'b')
    # (c) detection limit (was Fig. 5(f) before section 57): share of true bursts found vs j, at MS-like and PP7-like noise
    a = ax[2]
    from matplotlib.lines import Line2D
    sim = pd.read_csv(os.path.join(A, 'detection_grid_highSP', 'sim_summary.csv'))
    for nm, cl, mk in [('MS-like', BLUE, 'o'), ('PP7-like', OI[1], 's')]:
        d = sim[sim.setting.str.startswith(nm)].sort_values('SP')
        a.plot(d.SP, 100 * d.recall, '-' + mk, color=cl, mfc='none', ms=3, lw=0.8)
        a.plot(d.SP, 100 * d.false_over_detected.clip(upper=1), ':', color=cl, lw=1)
    a.legend(handles=[Line2D([], [], color=BLUE, lw=1.5, label='original noise'),
                      Line2D([], [], color=OI[1], lw=1.5, label='PP7 noise'),
                      Line2D([], [], color='k', lw=0.8, marker='o', mfc='none', ms=3, label='bursts found'),
                      Line2D([], [], color='k', lw=1, ls=':', label='false share')], loc='upper right', fontsize=5, handlelength=1.8)
    a.set_xlabel('Stringency parameter j'); a.set_ylabel('Simulated bursts (%)'); a.set_ylim(0, 105); a.set_yticks([0, 25, 50, 75, 100])
    a.set_xticks([10, 15, 20, 25]); label(a, 'c')
    fig.tight_layout(w_pad=0.6)
    save(fig, 'fig19')


# ---------------------------------------------------------------- Appendix: Fig. 20 (which lambda the fit picks)
def fig20():
    c = pd.read_csv(os.path.join(A, 'ms_sim_test', 'fig13_choice_summary.csv'))
    rows = [r for r in c.row.unique() if 'SNR 1.5' in r]
    scores = list(c.score.unique())
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.4), gridspec_kw={'width_ratios': [1.4, 1]})
    a = ax[0]
    ws = c[c.score == 'WS'].set_index('row').loc[rows]
    bottom = np.zeros(len(rows))
    for k in range(1, 6):
        v = ws[f'P({k})'].values
        a.bar(range(len(rows)), v, bottom=bottom, color=OI[k - 1], width=0.7, label=f'$\\lambda_{{est}}$ = {k}')
        bottom += v
    a.set_xticks(range(len(rows))); a.set_xticklabels([r.split(',')[0].replace('λ = ', '') for r in rows])
    a.set_xlabel('$\\lambda_{sim}$'); a.set_ylabel('Fraction of fits (WS score)'); a.set_ylim(0, 1)
    a.legend(fontsize=5.5, loc='upper left', bbox_to_anchor=(1.0, 1.0)); label(a, 'a', x=-0.12)
    a = ax[1]
    for i, sc in enumerate(scores):
        d = c[c.score == sc].set_index('row').loc[rows]
        mean = sum(k * d[f'P({k})'] for k in range(1, 6))
        lt = [float(r.split(',')[0].replace('λ = ', '')) for r in rows]
        a.plot(lt, mean.values, '-o', color=OI[i % 7], mfc='none', ms=3, lw=0.8, label=sc)
    a.plot([0, 5.5], [0, 5.5], 'r--', lw=1)
    a.set_xlim(0, 5.5); a.set_ylim(0, 5.5); a.set_xlabel('$\\lambda_{sim}$'); a.set_ylabel('Mean $\\lambda_{est}$')
    a.legend(fontsize=5.5, ncol=2, loc='upper left'); label(a, 'b')
    fig.tight_layout()
    save(fig, 'fig20')


# ---------------------------------------------------------------- Appendix: Fig. 21 (large slow changes)
def fig21():
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.3), sharey=False)
    for j, (folder, ttl) in enumerate([('slow_changes', 'relative to the granule'), ('slow_changes_absolute', 'same amount of light')]):
        s = pd.read_csv(os.path.join(A, folder, 'summary.csv'))
        s = s[s.threshold == 30]
        tot = s.groupby('dataset').excess.sum()
        lo = s.groupby('dataset').lo.sum(); hi = s.groupby('dataset').hi.sum()
        a = ax[j]; x = np.arange(len(VIVO))
        y = tot.loc[VIVO].values
        a.bar(x, y, color=[COL[d] for d in VIVO], alpha=0.85, width=0.65)
        a.errorbar(x, y, yerr=[y - lo.loc[VIVO].values, hi.loc[VIVO].values - y], fmt='none', ecolor='k', lw=0.7, capsize=2)
        a.axhline(0, color='k', lw=0.5)
        a.set_xticks(x); a.set_xticklabels([SHORT[d].replace(', ', '\n') for d in VIVO], fontsize=6)
        a.set_ylabel('Excess changes ≥ 30%\nper granule-hour'); a.set_title(ttl, fontsize=7.5)
        t = pd.read_csv(os.path.join(A, folder, 'tests.csv'))
        p = t[(t.threshold == 30) & (t['compare'] == '4PP7-FUS vs 4PP7')].p_two_sided.iloc[0]
        a.text(0.5, 0.92, f'PP7-FUS vs PP7 (4x): p = {p:.3f}', transform=a.get_xaxis_transform(), ha='left', fontsize=6.5)
        label(a, 'ab'[j], x=-0.13)
    fig.tight_layout()
    save(fig, 'fig21')


# ---------------------------------------------------------------- Appendix: Fig. 22 (first vs second half)
def fig22():
    h = pd.read_csv(os.path.join(A, 'halves', 'summary.csv'))
    if VSUM:
        vs = pd.read_csv(os.path.join(A, 'vitro_spotsum', 'summary.csv')).set_index('dataset').loc[VITRO]
        add = [dict(dataset=d, measure=m, first=vs.loc[d, k + '1'], second=vs.loc[d, k + '2'], p=vs.loc[d, 'pHalves' + P]) for d in VITRO for m, k, P in (('turn', 'turnover', 'Turnover'), ('net', 'net', 'Net'))]
        h = pd.concat([h[~h.dataset.isin(VITRO)], pd.DataFrame(add)])
    elif VNEW:
        vn = pd.read_csv(os.path.join(A, 'vitro_same', 'summary.csv')).set_index('dataset').loc[VITRO]
        add = [dict(dataset=d, measure=m, first=vn.loc[d, m + '1'], second=vn.loc[d, m + '2'], p=vn.loc[d, f'halves_{m}_p']) for d in VITRO for m in ('turn', 'net')]
        h = pd.concat([h[~h.dataset.isin(VITRO)], pd.DataFrame(add)])
    print('halves measures:', h.measure.unique())
    keys = VIVO + VITRO
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.4))
    for j, (meas, yl) in enumerate([('turn', 'Turnover beyond noise (%/h)'), ('net', 'Net change beyond bleaching (%/h)')]):
        d = h[h.measure == meas].set_index('dataset')
        a = ax[j]; x = np.arange(len(keys))
        a.plot(x - 0.12, d.loc[keys, 'first'], 'o', color=BLUE, mfc='none', ms=4, label='first 30 min')
        a.plot(x + 0.12, d.loc[keys, 'second'], 's', color=OI[1], mfc='none', ms=4, label='second 30 min')
        for i, k in enumerate(keys):
            a.plot([i - 0.12, i + 0.12], [d.loc[k, 'first'], d.loc[k, 'second']], '-', color='0.6', lw=0.6)
            if d.loc[k, 'p'] < 0.05:
                a.text(i, max(d.loc[k, 'first'], d.loc[k, 'second']) + 4, '*', ha='center', fontsize=9)
        a.axhline(0, color='k', lw=0.5)
        a.set_xticks(x); a.set_xticklabels([SHORT[k].replace(', ', '\n').replace('in vitro ', 'in vitro\n') for k in keys], fontsize=5.3)
        a.set_ylabel(yl); label(a, 'ab'[j], x=-0.13)
    ax[0].legend(fontsize=6, loc='upper left')
    fig.tight_layout()
    save(fig, 'fig22')


if __name__ == '__main__' and VNEW:
    for f in [fig6, fig18, fig22]:
        f()
elif __name__ == '__main__':
    for f in [fig4, fig5, fig6, lambda: fig_pair(['4PP7', '30PP7'], 'fig7'), lambda: fig_pair(['5Qb', '10Qb'], 'fig8'),
              fig9, fig10, fig15, fig16, fig17, fig18, fig19, fig20, fig21, fig22]:
        try:
            f()
        except Exception as ex:   # keep going, report at the end
            import traceback; traceback.print_exc()
