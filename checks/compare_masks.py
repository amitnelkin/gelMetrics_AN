"""COMPARE_MASKS  Three extractions side by side, each judged against its own
unchanging-granule control (RESULTS.md section 33):
  R2  per-frame multithresh thresholds (published)
  R4  thresholds = running median over 31 frames (section 23)
  R5  one fixed mask per granule from its average crop, smoothed track (section 33)
Per dataset: frame noise sqrt(D1/2) (real and control), false/real negative
bursts (window 13, SP 12), turnover beyond noise (%/h; real minus burst-free
copies) for the real granules and for the control, and real minus control.
R5 only: share of the granule's light inside the fixed mask, first vs last 30
frames of each trace (does the granule outgrow the mask?).
Comparisons (4PP7-FUS vs 4PP7, 30PP7-FUS vs 30PP7, 5Qb vs 4PP7, 10Qb vs 30PP7):
  turnover beyond its control: z from the bootstrap SEs of real and control in both
  conditions (nested movie -> granule, 1000x);
  burst size beyond noise copies (median fraction shed, real - copies): movie-level
  permutation p (5000x), two-sided and one-sided (FUS > PP7 pre-specified).
Usage: python checks/compare_masks.py
Writes analysis_outputs/fixedmask/{compare.csv, tests.csv, compare_masks.png/.pdf}."""
import os
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = os.path.join('analysis_outputs', 'fixedmask'); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(33); H = 10 / 3600; NB, NPERM = 1000, 5000
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
COMPS = [('4PP7-FUS', '4PP7', True), ('30PP7-FUS', '30PP7', True), ('5Qb', '4PP7', False), ('10Qb', '30PP7', False)]
S = os.path.join('analysis_outputs', 'steady')
VER = {'R2 per-frame (published)': dict(grid='analysis_outputs/detection_grid', sf='analysis_outputs/noise_methods/real_sf.csv', ver='v0',
                                        static=[f'{S}/static_control.csv', f'{S}/static_control_more.csv'], col='0.5'),
       'R4 steady thresholds': dict(grid=f'{S}/grid', sf=f'{S}/real_sf.csv', ver='v31', static=[f'{S}/static_control.csv', f'{S}/static_control_more.csv'], col='tab:green'),
       'R5 fixed mask': dict(grid=f'{out}/grid', sf=f'{out}/real_sf.csv', ver='vfix', static=[f'{S}/static_control_fixed.csv'], col='tab:purple')}


def nested_idx(mov):
    movs = {m: np.where(mov == m)[0] for m in np.unique(mov)}
    return [np.concatenate([rng.choice(movs[m], len(movs[m])) for m in rng.choice(list(movs), len(movs))]) for _ in range(NB)]


def excess(fr, fs, L, idx=None):
    if idx is None:
        idx = slice(None)
    return 100 * (fr[idx].sum() - fs[idx].sum()) / (L[idx].sum() * H)


rows, items = [], {}
for vname, v in VER.items():
    B = pd.read_csv(os.path.join(v['grid'], 'bursts.csv')); B = B[(B.window == 13) & (B.SP == 12) & (B.levelBefore > 0)].copy()
    B['f'] = B.amp.abs() / B.levelBefore
    T = pd.read_csv(os.path.join(v['grid'], 'traces.csv'))
    SF = pd.read_csv(v['sf']); SF = SF[SF.measure == 'spotsum']
    SC = pd.concat([pd.read_csv(p) for p in v['static']]); SC = SC[SC.version == v['ver']]
    for ds in VIVO:
        t = T[T.dataset == ds].reset_index(drop=True); b = B[B.dataset == ds]
        fr = t.trace.map(b[b.kind == 'real'].groupby('trace').f.sum()).fillna(0).values
        fs = t.trace.map(b[b.kind == 'surr'].groupby('trace').f.sum()).fillna(0).values
        L = t.lengthFrames.values.astype(float); mov = t.movie.values
        bi = nested_idx(mov); er = excess(fr, fs, L); er_b = np.array([excess(fr, fs, L, i) for i in bi])
        s = SC[SC.dataset == ds].reset_index(drop=True)
        ss = (s.fracReal13.values, s.fracSurr13.values, s.lengthFrames.values.astype(float))
        si = nested_idx(s.movie.values); es = excess(*ss); es_b = np.array([excess(*ss, i) for i in si])
        D1 = np.array([float(x.strip('[]').split()[0]) for x in SF[SF.dataset == ds].Dreal])
        R = {k: g.f.values for k, g in b[b.kind == 'real'].groupby('trace')}; Su = {k: g.f.values for k, g in b[b.kind == 'surr'].groupby('trace')}
        e0 = np.zeros(0); items[(vname, ds)] = {m: [(R.get(k, e0), Su.get(k, e0)) for k in g.trace] for m, g in t.groupby('movie')}
        rows.append(dict(version=vname, dataset=ds, noise_real=100 * np.sqrt(np.median(D1) / 2), noise_static=100 * np.sqrt(s.D1.dropna().median() / 2),
                         FDR_real=(b.kind == 'surr').sum() / max((b.kind == 'real').sum(), 1), FDR_static=s.nSurr13.sum() / max(s.nReal13.sum(), 1),
                         excess_real=er, excess_real_lo=np.percentile(er_b, 2.5), excess_real_hi=np.percentile(er_b, 97.5), excess_real_se=er_b.std(),
                         excess_static=es, excess_static_lo=np.percentile(es_b, 2.5), excess_static_hi=np.percentile(es_b, 97.5), excess_static_se=es_b.std(),
                         beyond_control=er - es, beyond_control_se=np.hypot(er_b.std(), es_b.std())))
    print(vname, 'done', flush=True)
C = pd.DataFrame(rows)
# mask diagnostic for R5
mf = []
for ds in VIVO:
    m = sio.loadmat(f'in_vivo_reextracted/after_R5_fixedmask/{ds}_maskfrac.mat', squeeze_me=True)['maskFrac'][1:]
    m = [np.ravel(x) for x in m if np.size(x) >= 180]
    mf.append(dict(dataset=ds, maskFrac_first30=np.median([np.nanmean(x[:30]) for x in m]), maskFrac_last30=np.median([np.nanmean(x[-30:]) for x in m])))
C = C.merge(pd.DataFrame(mf).assign(version='R5 fixed mask'), on=['version', 'dataset'], how='left')
C.to_csv(os.path.join(out, 'compare.csv'), index=False)
pd.set_option('display.width', 260); print(C.drop(columns=[c for c in C if c.endswith('_se') or c.endswith('_lo') or c.endswith('_hi')]).round(2).to_string(index=False))


def med_beyond(it):
    r = np.concatenate([x[0] for x in it]); c = np.concatenate([x[1] for x in it])
    return 100 * (np.median(r) - np.median(c))


tests = []
for vname in VER:
    c = C[C.version == vname].set_index('dataset')
    for a, b_, claimed in COMPS:
        z = (c.loc[a, 'beyond_control'] - c.loc[b_, 'beyond_control']) / np.hypot(c.loc[a, 'beyond_control_se'], c.loc[b_, 'beyond_control_se'])
        Ga, Gb = items[(vname, a)], items[(vname, b_)]
        obs = med_beyond([x for g in Ga.values() for x in g]) - med_beyond([x for g in Gb.values() for x in g])
        mv = list(Ga.values()) + list(Gb.values()); na = len(Ga); null = np.empty(NPERM)
        for i in range(NPERM):
            p = rng.permutation(len(mv))
            null[i] = med_beyond([x for j in p[:na] for x in mv[j]]) - med_beyond([x for j in p[na:] for x in mv[j]])
        p2 = (np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (NPERM + 1)
        p1 = (np.sum(null >= obs - 1e-12) + 1) / (NPERM + 1) if obs >= 0 else (np.sum(null <= obs + 1e-12) + 1) / (NPERM + 1)
        tests.append(dict(version=vname, compare=f'{a} vs {b_}', turnover_beyond_control_A=c.loc[a, 'beyond_control'], turnover_beyond_control_B=c.loc[b_, 'beyond_control'],
                          z_turnover_vs_control=z, size_beyond_noise_diff=obs, p_size_two_sided=p2, p_size_one_sided=p1, one_sided_prespecified=claimed and obs > 0))
Tt = pd.DataFrame(tests); Tt.to_csv(os.path.join(out, 'tests.csv'), index=False); print(Tt.round(3).to_string(index=False))

fig, axs = plt.subplots(1, 4, figsize=(24, 5.4))
x = np.arange(len(VIVO)); wd = 0.13
for k, (vname, v) in enumerate(VER.items()):
    t = C[C.version == vname].set_index('dataset').reindex(VIVO)
    for q, (key, hatch, alpha, lab2) in enumerate([('real', None, 1.0, 'real granules'), ('static', '///', 0.35, 'unchanging-granule control')]):
        off = (2 * k + q - 2.5) * wd
        axs[0].bar(x + off, t[f'noise_{key}'], wd, color=v['col'], alpha=alpha, hatch=hatch, edgecolor='k', lw=0.4, label=f'{vname}: {lab2}')
        axs[1].bar(x + off, t[f'FDR_{key}'], wd, color=v['col'], alpha=alpha, hatch=hatch, edgecolor='k', lw=0.4)
        m = t[f'excess_{key}']
        axs[2].errorbar(x + off, m, yerr=[m - t[f'excess_{key}_lo'], t[f'excess_{key}_hi'] - m], fmt='o', color=v['col'], mfc=v['col'] if key == 'real' else 'w', ms=6, capsize=2)
    axs[3].errorbar(x + (k - 1) * 0.22, t.beyond_control, yerr=1.96 * t.beyond_control_se, fmt='s', color=v['col'], ms=7, capsize=3)
axs[0].set_ylabel('frame-to-frame noise (% of granule)'); axs[1].set_ylabel('false / real negative bursts'); axs[1].axhline(1, color='k', ls=':', lw=0.8)
axs[2].set_ylabel('turnover beyond noise (%/h, 95% CI)'); axs[2].axhline(0, color='k', lw=0.8)
axs[3].set_ylabel('real minus control (%/h, 95% CI)'); axs[3].axhline(0, color='k', lw=0.8)
for ax, ttl in zip(axs, ['(1) noise: real vs control', '(2) share of bursts noise also produces', '(3) turnover beyond noise: filled real, open control',
                         '(4) what is left after the control']):
    ax.set_xticks(x); ax.set_xticklabels(VIVO, rotation=30, fontsize=9); ax.set_title(ttl, fontsize=10, loc='left')
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=6, fontsize=9, frameon=False)
fig.tight_layout(rect=[0, 0.07, 1, 1]); fig.savefig(os.path.join(out, 'compare_masks.png'), dpi=130); fig.savefig(os.path.join(out, 'compare_masks.pdf'))
