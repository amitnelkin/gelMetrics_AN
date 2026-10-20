"""REPORT_PARAMS  Proposed report parameters per dataset (replacing lambda and k0),
for the simulated (checks/report_params_sim.m) and in vivo datasets, as a table
and as mock-ups of manuscript Figs. 7, 10(b) and 15. RESULTS.md section 24.
Parameters (in vivo: R2 extraction = published per-frame multithresh, spot sum;
window 13, SP 12):
  turnover beyond noise (%/h)   sum of fraction shed per negative burst on the
                                traces minus on their burst-free surrogates,
                                per granule-hour (section 18-20 'excess turnover')
  fluctuation beyond noise (%)  G30 = 100*sqrt(mean(D30 - D1)) of the relative
                                granule intensity (white noise cancels). Check:
                                the pixel-level static-granule control (section
                                23) gives mean(D30 - D1) consistent with 0 in all
                                six in vivo datasets, so nothing is subtracted
                                (subtracting a noisy zero only adds noise); its
                                value is shown as 'fluct_static'.
  net change (%/h)              slope of ln(granule / cell background)
                                (section 22; Theil-Sen per granule, mean)
  noise-explained share         turnover on the surrogates / turnover on the traces
                                (burst counts mislead at high SNR: many tiny noise bursts)
95% CIs: nested movie -> granule bootstrap in vivo, granule bootstrap in the
simulations, 1000x (seeded). True values for the simulations from the
noise-free twin.
Usage: python checks/report_params.py
Writes analysis_outputs/report_params/{params.csv, params.md, fig7_mock.png,
fig10b_mock.png, fig15_mock.png} (+ .pdf)."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import theilslopes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'report_params')
rng = np.random.default_rng(24)
H = 10 / 3600; NB = 1000
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
LABEL = {'4PP7': 'PP7 + QCP-5x/PCP-4x', '4PP7-FUS': 'PP7-FUS + QCP-5x/PCP-4x', '30PP7': 'PP7 + PCP-30x',
         '30PP7-FUS': 'PP7-FUS + PCP-30x', '5Qb': 'Qβ + QCP-5x/PCP-4x', '10Qb': 'Qβ + QCP-10x'}
SIMS = ['no bursts', 'lambda 1', 'lambda 2', 'lambda 4', 'lambda 2, SNR 5', 'lambda 2, SNR 0.5']
TAUS = np.array([1, 2, 3, 6, 12, 18, 30, 45, 60])


def boot_groups(t):
    """index arrays of nested (movie -> granule) bootstrap resamples"""
    movs = {m: np.where(t.movie.values == m)[0] for m in np.unique(t.movie.values)}
    for _ in range(NB):
        if len(movs) > 1:
            yield np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))])
        else:
            yield rng.choice(len(t), len(t))


def ci(v):
    return np.percentile(v, [2.5, 97.5])


def turnover(t):
    """t: per-trace frame with fracReal, fracSurr, lengthFrames, movie"""
    f = lambda x: 100 * (x.fracReal.sum() - x.fracSurr.sum()) / (x.lengthFrames.sum() * H)
    return f(t), ci([f(t.iloc[p]) for p in boot_groups(t)])


def g2(D):
    return np.mean(D[:, 1] - D[:, 0])


def gnet(t, ctrl):
    """t, ctrl: per-trace frames with D1, D30, movie; G beyond the control (in %)"""
    a = t[['D1', 'D30']].values; c = ctrl[['D1', 'D30']].values if ctrl is not None else None
    f = lambda a_, c_: 100 * np.sqrt(max(g2(a_) - (g2(c_) if c_ is not None else 0), 0))
    bs = []
    gc = list(boot_groups(ctrl)) if ctrl is not None else None
    for k, p in enumerate(boot_groups(t)):
        bs.append(f(a[p], c[gc[k]] if c is not None else None))
    return f(a, c), ci(bs)


def meanci(t, col):
    t = t.dropna(subset=[col]).reset_index(drop=True)
    return t[col].mean(), ci([t[col].values[p].mean() for p in boot_groups(t)])


rows = []
# ---------- simulations
ST = pd.read_csv(os.path.join(out, 'sim_traces.csv')); ST['movie'] = 0
SG = sio.loadmat(os.path.join(out, 'sim_signals.mat'), squeeze_me=True, struct_as_record=False)
ble = np.exp(-0.0026 * np.linspace(0, 60, 361))
for j, name in enumerate(SIMS):
    t = ST[ST.dataset == name].reset_index(drop=True); s = SG[f's{j + 1}']
    t['net'] = [theilslopes(np.log(a / b), np.arange(len(a)))[0] * 36000 if np.all(a > 0) and np.all(b > 0) else np.nan for a, b in zip(s.S, s.B)]
    t['trueNet'] = [theilslopes(np.log(a / ble), np.arange(len(a)))[0] * 36000 for a in s.S0]
    to, toc = turnover(t); g, gc = gnet(t, None); n, nc = meanci(t, 'net')
    tt = t.assign(fracSurr=0, fracReal=t.trueFrac)
    rows.append(dict(kind='simulated', dataset=name, granules=len(t), movies=np.nan, burstsReal=int(t.nReal.sum()), burstsNoise=int(t.nSurr.sum()),
                     noiseShare=t.fracSurr.sum() / t.fracReal.sum(), turnover=to, turnover_lo=toc[0], turnover_hi=toc[1],
                     fluct=g, fluct_lo=gc[0], fluct_hi=gc[1], net=n, net_lo=nc[0], net_hi=nc[1],
                     true_turnover=100 * tt.fracReal.sum() / (tt.lengthFrames.sum() * H),
                     true_fluct=100 * np.sqrt(max(t.trueD30minusD1.mean(), 0)), true_net=t.trueNet.mean()))
# ---------- in vivo (R2)
GB = pd.read_csv('analysis_outputs/detection_grid/bursts.csv')
GB = GB[(GB.window == 13) & (GB.SP == 12) & (GB.levelBefore > 0) & GB.dataset.isin(VIVO)].copy(); GB['frac'] = GB.amp.abs() / GB.levelBefore
GT = pd.read_csv('analysis_outputs/detection_grid/traces.csv')
SF = pd.read_csv('analysis_outputs/noise_methods/real_sf.csv'); SF = SF[SF.measure == 'spotsum']
SC = pd.concat([pd.read_csv('analysis_outputs/steady/static_control.csv'), pd.read_csv('analysis_outputs/steady/static_control_more.csv')])
SC = SC[SC.version == 'v0']
NL = pd.read_csv('analysis_outputs/net_loss/per_trace.csv')
for ds in VIVO:
    t = GT[GT.dataset == ds].reset_index(drop=True); b = GB[GB.dataset == ds]
    t['fracReal'] = t.trace.map(b[b.kind == 'real'].groupby('trace').frac.sum()).fillna(0)
    t['fracSurr'] = t.trace.map(b[b.kind == 'surr'].groupby('trace').frac.sum()).fillna(0)
    to, toc = turnover(t)
    sf = SF[SF.dataset == ds]
    D = np.array([[float(v) for v in x.strip('[]').split()] for x in sf.Dreal])
    tf = pd.DataFrame(dict(D1=D[:, 0], D30=D[:, 6], movie=sf.movie.values))
    ctrl = SC[SC.dataset == ds].dropna(subset=['D1', 'D30']).reset_index(drop=True)
    g, gc = gnet(tf, None); n, nc = meanci(NL[NL.dataset == ds].reset_index(drop=True), 'sumNet')
    nr, ns = (b.kind == 'real').sum(), (b.kind == 'surr').sum()
    rows.append(dict(kind='in vivo', dataset=ds, granules=len(t), movies=t.movie.nunique(), burstsReal=nr, burstsNoise=ns, noiseShare=t.fracSurr.sum() / t.fracReal.sum(),
                     turnover=to, turnover_lo=toc[0], turnover_hi=toc[1], fluct=g, fluct_lo=gc[0], fluct_hi=gc[1],
                     net=n, net_lo=nc[0], net_hi=nc[1], fluct_raw=100 * np.sqrt(max(g2(tf[['D1', 'D30']].values), 0)),
                     fluct_static=100 * np.sqrt(max(g2(ctrl[['D1', 'D30']].values), 0))))
P = pd.DataFrame(rows); P.to_csv(os.path.join(out, 'params.csv'), index=False)
fmt = lambda r, c: f"{r[c]:.1f} [{r[c + '_lo']:.1f}, {r[c + '_hi']:.1f}]"
with open(os.path.join(out, 'params.md'), 'w', encoding='utf-8', newline='\r\n') as fh:
    fh.write('| dataset | turnover beyond noise (%/h) | fluctuation beyond noise, 5 min (%) | net change (%/h) | share of turnover explained by noise | bursts (real / noise) | granules (movies) | true turnover / fluct. / net |\n|---|---|---|---|---|---|---|---|\n')
    for _, r in P.iterrows():
        name = LABEL.get(r.dataset, 'sim: ' + r.dataset)
        tr = f"{r.true_turnover:.1f} / {r.true_fluct:.1f} / {r.true_net:.1f}" if r.kind == 'simulated' else '–'
        mv = f" ({int(r.movies)})" if r.kind == 'in vivo' else ''
        fh.write(f"| {name} | {fmt(r, 'turnover')} | {fmt(r, 'fluct')} | {fmt(r, 'net')} | {r.noiseShare:.2f} | {r.burstsReal} / {r.burstsNoise} | {r.granules}{mv} | {tr} |\n")
print(open(os.path.join(out, 'params.md'), encoding='utf-8').read())

# ---------- Fig. 7 mock: burst-size distribution real vs noise; G(tau) real vs noise
COL = {'PP7': 'tab:blue', 'PP7-FUS': 'tab:green', 'Qβ': 'tab:red'}
prot = lambda ds: 'Qβ' if 'Qb' in ds else ('PP7-FUS' if 'FUS' in ds else 'PP7')
fig, axs = plt.subplots(2, 6, figsize=(21, 7))
bins = np.linspace(0, 0.8, 17)
for k, ds in enumerate(VIVO):
    b = GB[GB.dataset == ds]; hours = GT[GT.dataset == ds].lengthFrames.sum() * H; r = P[P.dataset == ds].iloc[0]
    ax = axs[0, k]
    hr, _ = np.histogram(b[b.kind == 'real'].frac.clip(upper=0.79), bins); hs, _ = np.histogram(b[b.kind == 'surr'].frac.clip(upper=0.79), bins)
    ax.bar(bins[:-1], hr / hours, width=np.diff(bins), align='edge', color=COL[prot(ds)], alpha=0.55, edgecolor='k', lw=0.5, label='granules')
    ax.step(bins[:-1], hs / hours, where='post', color='k', lw=1.5, label='burst-free copies (noise)')
    ax.set_title(LABEL[ds], fontsize=9, loc='left'); ax.set_xlabel('fraction of granule shed per burst')
    ax.text(0.97, 0.95, f"turnover beyond noise\n{fmt(r, 'turnover')} %/h\nnoise explains {100 * r.noiseShare:.0f}% of turnover",
            transform=ax.transAxes, ha='right', va='top', fontsize=7.5)
    ax = axs[1, k]
    sf = SF[SF.dataset == ds]; D = np.array([[float(v) for v in x.strip('[]').split()] for x in sf.Dreal])
    Dsur = np.array([[float(v) for v in x.strip('[]').split()] for x in sf.Dsurr])
    G = 100 * np.sqrt(np.maximum((D - D[:, :1]).mean(0), 0)); Gs = 100 * np.sqrt(np.maximum((Dsur - Dsur[:, :1]).mean(0), 0))
    ax.plot(TAUS / 6, G, 'o-', color=COL[prot(ds)], label='granules'); ax.plot(TAUS / 6, Gs, 's--', color='0.4', ms=4, label='white-noise copies')
    ax.errorbar([5], [r.fluct_static], fmt='kD', ms=6, label='static-granule control')
    ax.set_xlabel('time lag (min)'); ax.set_ylim(bottom=0)
    ax.text(0.03, 0.95, f"beyond noise at 5 min:\n{fmt(r, 'fluct')} %", transform=ax.transAxes, va='top', fontsize=7.5)
axs[0, 0].set_ylabel('negative bursts per granule-hour'); axs[1, 0].set_ylabel('fluctuation G (% of granule)')
axs[0, 0].legend(fontsize=7, loc='center right'); axs[1, 0].legend(fontsize=7, loc='lower right')
fig.suptitle('Mock-up of Fig. 7 / 9 (in vivo): what replaces the amplitude histogram + Poisson fit', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig7_mock.png'), dpi=130); fig.savefig(os.path.join(out, 'fig7_mock.pdf'))

# same for the simulations
SB = None
fig, axs = plt.subplots(1, 6, figsize=(21, 3.8))
for k, name in enumerate(SIMS):
    r = P[P.dataset == name].iloc[0]; t = ST[ST.dataset == name]
    ax = axs[k]
    ax.bar([0, 1, 2], [r.turnover, r.fluct, r.net], yerr=[[r.turnover - r.turnover_lo, r.fluct - r.fluct_lo, r.net - r.net_lo],
           [r.turnover_hi - r.turnover, r.fluct_hi - r.fluct, r.net_hi - r.net]], color=['0.6', '0.75', '0.9'], edgecolor='k', capsize=3, label='estimated')
    ax.plot([0, 1, 2], [r.true_turnover, r.true_fluct, r.true_net], 'r*', ms=12, label='true')
    ax.axhline(0, color='k', lw=0.7)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['turnover\n%/h', 'fluct. 5 min\n%', 'net\n%/h'], fontsize=8)
    ax.set_title(f'sim: {name}\n(noise explains {100 * r.noiseShare:.0f}% of turnover)', fontsize=9, loc='left')
axs[0].legend(fontsize=8)
fig.suptitle('Simulated datasets (as Fig. 13): estimated parameters (95% CI) vs true values', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig7_mock_sim.png'), dpi=130); fig.savefig(os.path.join(out, 'fig7_mock_sim.pdf'))

# ---------- Fig. 10(b) mock
MK = {'4PP7': ('o', '-'), '4PP7-FUS': ('o', '-'), '5Qb': ('o', '-'), '30PP7': ('s', '--'), '30PP7-FUS': ('s', '--'), '10Qb': ('^', ':')}
fig, axs = plt.subplots(1, 3, figsize=(17, 5))
V = P[P.kind == 'in vivo']
for _, r in V.iterrows():
    m, ls = MK[r.dataset]; c = COL[prot(r.dataset)]
    e = axs[0].errorbar(r.turnover, r.fluct, xerr=[[r.turnover - r.turnover_lo], [r.turnover_hi - r.turnover]],
                        yerr=[[r.fluct - r.fluct_lo], [r.fluct_hi - r.fluct]], fmt=m, color=c, ms=8, capsize=3, label=LABEL[r.dataset])
    e[-1][0].set_linestyle(ls); e[-1][1].set_linestyle(ls)
axs[0].axvline(0, color='k', lw=0.7); axs[0].axhline(0, color='k', lw=0.7)
axs[0].set_xlabel('turnover beyond noise (%/h)'); axs[0].set_ylabel('fluctuation beyond noise over 5 min (%)')
axs[0].set_title('(b) in vivo: dynamics beyond noise', fontsize=10, loc='left'); axs[0].legend(fontsize=7)
for k, (_, r) in enumerate(V.iterrows()):
    axs[1].errorbar(k, r.net, yerr=[[r.net - r.net_lo], [r.net_hi - r.net]], fmt=MK[r.dataset][0], color=COL[prot(r.dataset)], ms=8, capsize=4)
axs[1].axhline(0, color='k', lw=0.7); axs[1].set_xticks(range(len(V))); axs[1].set_xticklabels([LABEL[d] for d in V.dataset], rotation=35, ha='right', fontsize=7)
axs[1].set_ylabel('net change beyond bleaching (%/h)'); axs[1].set_title('(c) in vivo: net gain / loss', fontsize=10, loc='left')
S_ = P[P.kind == 'simulated']
for c_, (col, lab) in enumerate([('turnover', 'turnover %/h'), ('fluct', 'fluctuation 5 min %')]):
    axs[2].errorbar(S_[f'true_{col}'], S_[col], yerr=[S_[col] - S_[f'{col}_lo'], S_[f'{col}_hi'] - S_[col]], fmt='os'[c_], ms=7, capsize=3, label=lab)
    for _, r in S_.iterrows():
        axs[2].annotate(r.dataset, (r[f'true_{col}'], r[col]), fontsize=6, xytext=(3, 3), textcoords='offset points')
lim = [0, max(S_.true_turnover.max(), S_.turnover_hi.max()) * 1.1]
axs[2].plot(lim, lim, 'k-', lw=0.7); axs[2].set_xlabel('true value (simulation)'); axs[2].set_ylabel('estimated')
axs[2].set_title('(d) simulations: estimated vs true', fontsize=10, loc='left'); axs[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig10b_mock.png'), dpi=130); fig.savefig(os.path.join(out, 'fig10b_mock.pdf'))

# ---------- Fig. 15 mock: datasets x parameters
cols = [('turnover', 'turnover beyond\nnoise (%/h)'), ('fluct', 'fluctuation beyond\nnoise, 5 min (%)'), ('net', 'net change\n(%/h)'), ('noiseShare', 'share of turnover\nexplained by noise')]
names = [('sim: ' + d) for d in SIMS] + [LABEL[d] for d in VIVO]
fig, ax = plt.subplots(figsize=(13, 7))
Mx = np.zeros((len(P), len(cols)))
for j, (c, _) in enumerate(cols):
    v = P[c].values.astype(float); Mx[:, j] = (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v) + 1e-12)
ax.imshow(Mx, cmap='viridis', aspect='auto')
for i, (_, r) in enumerate(P.iterrows()):
    for j, (c, _) in enumerate(cols):
        txt = f"{r[c]:.2f}" if c == 'noiseShare' else f"{r[c]:.1f}\n[{r[c + '_lo']:.1f}, {r[c + '_hi']:.1f}]"
        if r.kind == 'simulated' and c != 'noiseShare':
            txt += f"\ntrue {r['true_' + c]:.1f}"
        ax.text(j, i, txt, ha='center', va='center', fontsize=7, color='w' if Mx[i, j] < 0.6 else 'k')
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols], fontsize=9); ax.xaxis.tick_top()
ax.set_yticks(range(len(P))); ax.set_yticklabels(names, fontsize=8); ax.axhline(len(SIMS) - 0.5, color='w', lw=3)
ax.set_title('Mock-up of Figs. 13-15: proposed parameters per dataset (95% CI; colour scaled per column)', fontsize=10, pad=40)
fig.tight_layout(); fig.savefig(os.path.join(out, 'fig15_mock.png'), dpi=130); fig.savefig(os.path.join(out, 'fig15_mock.pdf'))
