"""NOISE_METHODS  Three noise-robust ways to score granule dynamics, judged on a
simulation bench with known dynamics (checks/noise_methods_sim.m) and applied
to the real data (checks/score_export2.m 'base'/'surrpoly' bursts;
checks/noise_methods_real.m fluctuations). RESULTS.md section 19.
Reviewer points R2 6 (noise floor), R2 4a, R2 5, "comparison with other methods".
Real data: M1/M3 are also computed with a second surrogate whose noise level
is the plain std of the frame differences (upper estimate; suffix _sd), to
bracket the noise correction (the default uses a robust MAD estimate).
Methods (negative bursts; frac = |amplitude| / granule level before the burst):
  M0 uncorrected  turnover = 100 * sum(frac) / hours          (reference; section 17)
  M1 corrected    excess turnover = 100 * (sum frac_real - sum frac_surr) / hours,
                  each trace paired with its own burst-free surrogate
  M2 fluctuation  F(tau) = 100 * sqrt(mean_traces(D_real - D_surr)) = typical
                  relative change of the granule over lag tau beyond noise (%);
                  tau = 30 frames (5 min) main, 6 and 60 frames reported
  M2b fluctuation G(tau) = 100 * sqrt(mean_traces(D_real(tau) - D_real(1 frame)));
                  white noise adds the same 2 sigma^2 at every lag, so subtracting
                  the 1-frame lag removes it exactly, whatever its level or shape
                  (no surrogate, no noise estimate); tau = 30 frames main
  M3 recalibrated keep bursts with frac >= c, c = the smallest cut at which the
                  surrogate gives <= 10% as many bursts as the real data (per
                  dataset); turnover of the kept bursts
Usage: python checks/noise_methods.py sim | real
Writes analysis_outputs/noise_methods/{sim_scores.csv, sim_eval.csv, real_scores.csv,
real_comparisons.csv, real_stability.csv, sim_methods.png, real_methods.png}."""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'noise_methods')
rob = os.path.join('analysis_outputs', 'robustness')
rng = np.random.default_rng(1)
H = 10 / 3600
CUTS = np.round(np.arange(0.0, 1.01, 0.02), 2)
ORDER = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7'),
         ('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1')]
SC = ['M0_turnover', 'M1_excess_turnover', 'M1sd_excess_turnover', 'M2_F30', 'M2b_G30', 'M3_turnover', 'M3sd_turnover']


def srho(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float); ok = np.isfinite(a) & np.isfinite(b)
    return stats.spearmanr(a[ok], b[ok]).statistic if ok.sum() >= 3 else np.nan


def parse(s):
    return np.array([float(v) for v in s.strip('[]').split()])


def choose_cut(fr, fs):
    for c in CUTS:
        nr, ns = (fr >= c).sum(), (fs >= c).sum()
        if nr >= 5 and ns <= 0.10 * nr:
            return c
    return np.nan


def methods(real, surr, sfr, sfs, hours, cut=None):
    """real/surr: arrays of negative-burst fractions; sfr/sfs: per-trace D arrays (n x taus)."""
    out = dict(M0_turnover=100 * real.sum() / hours, nReal=len(real), nSurr=len(surr))
    out['M1_excess_turnover'] = 100 * (real.sum() - surr.sum()) / hours
    out['M1_excess_rate'] = (len(real) - len(surr)) / hours
    if sfr is not None and len(sfr):
        dex = np.mean(sfr - sfs, axis=0)
        for k, tau in [(4, 6), (6, 30), (8, 60)]:
            out[f'M2_F{tau}'] = 100 * np.sqrt(max(dex[k], 0))
            g = np.mean(sfr[:, k] - sfr[:, 0])
            out[f'M2b_G{tau}'] = 100 * np.sqrt(max(g, 0))
            out[f'M2b_Gsurr{tau}'] = 100 * np.sqrt(max(np.mean(sfs[:, k] - sfs[:, 0]), 0))
    c = choose_cut(real, surr) if cut is None else cut
    out['M3_cut'] = c
    keep = real[real >= c] if np.isfinite(c) else real[:0]
    out['M3_turnover'] = 100 * keep.sum() / hours if np.isfinite(c) else np.nan
    out['M3_n'] = len(keep)
    out['M3_fmean'] = np.mean(keep) if len(keep) else np.nan
    return out


def bootstrap(trs, R, S, SFr, SFs, Lh, cut, nb=300):
    """trs: trace ids; R/S: dict trace->array of fracs; SFr/SFs: dict trace->D."""
    bs = []
    empty = np.zeros(0)
    for _ in range(nb):
        pick = rng.choice(trs, len(trs))
        real = np.concatenate([R.get(t, empty) for t in pick]); surr = np.concatenate([S.get(t, empty) for t in pick])
        sr = np.array([SFr[t] for t in pick if t in SFr]); ss = np.array([SFs[t] for t in pick if t in SFs])
        bs.append(methods(real, surr, sr if len(sr) else None, ss, sum(Lh[t] for t in pick), cut))
    return pd.DataFrame(bs)


def run_sim():
    B = pd.read_csv(os.path.join(d, 'sim_bursts.csv')); B = B[(B.sign < 0) & (B.levelBefore > 0)].copy()
    B['frac'] = B.amp.abs() / B.levelBefore
    T = pd.read_csv(os.path.join(d, 'sim_traces.csv')); SF = pd.read_csv(os.path.join(d, 'sim_sf.csv'))
    rows = []
    for st, g in T.groupby('set'):
        s, r, n = g.s.iloc[0], g.r.iloc[0], g.n.iloc[0]; hours = g.lengthFrames.sum() * H
        sfr = np.array([parse(x) for x in SF[(SF.set == st) & (SF.kind == 'real')].D])
        sfs = np.array([parse(x) for x in SF[(SF.set == st) & (SF.kind == 'surr')].D])
        for sp in [10, 12, 14]:
            b = B[(B.set == st) & (B.SP == sp)]
            m = methods(b[b.kind == 'real'].frac.values, b[b.kind == 'surr'].frac.values, sfr, sfs, hours)
            rows.append(dict(set=st, s=s, r=r, n=n, SP=sp, trueTurnover=100 * (r / 2) * s, truePower=r * s ** 2, **m))
    Q = pd.DataFrame(rows); Q.to_csv(os.path.join(d, 'sim_scores.csv'), index=False)
    q = Q[Q.SP == 12]
    ev = []
    for k, truth in [('M0_turnover', 'trueTurnover'), ('M1_excess_turnover', 'trueTurnover'), ('M3_turnover', 'trueTurnover'),
                     ('M2_F30', 'truePower'), ('M2b_G30', 'truePower')]:
        sq = k.startswith('M2')
        x = q[k] ** 2 if sq else q[k]
        zero = q[q.r == 0][k]
        dyn = q[q.r > 0]
        rho_all = srho(q[truth], x)
        rho_fixn = np.nanmean([srho(g[truth], g[k] ** 2 if sq else g[k])
                               for _, g in q.groupby('n')])
        # noise dependence: score at the highest noise / lowest noise, same dynamics, r > 0
        piv = dyn.groupby(['s', 'r', 'n'])[k].mean().unstack('n')
        ratio = np.nanmedian(piv[0.02] / piv[0.005])
        ev.append(dict(method=k, truth=truth, shareOfSetsScored=float(np.mean(np.isfinite(q[k]))), scoreWhenNoBursts_median=np.nanmedian(zero), scoreWhenNoBursts_max=np.nanmax(zero),
                       rho_vs_truth_all=rho_all, rho_vs_truth_fixedNoise=rho_fixn, highNoise_over_lowNoise=ratio,
                       spStability_CV=Q[Q.r > 0].groupby('set')[k].agg(lambda v: np.nanstd(v) / abs(np.nanmean(v))).median()))
    E = pd.DataFrame(ev); E.to_csv(os.path.join(d, 'sim_eval.csv'), index=False)
    pd.set_option('display.width', 220)
    print(E.round(3).to_string(index=False)); print()
    print(q.groupby(['n', 's', 'r'])[['trueTurnover', 'M0_turnover', 'M1_excess_turnover', 'M3_turnover', 'M3_n', 'M2_F30', 'M2b_G30']].mean().round(2).to_string())
    fig, axs = plt.subplots(1, 5, figsize=(23, 4.5))
    cols = {0.005: '#2a78d6', 0.01: '#eb6834', 0.02: '#d62a2a'}
    for ax, (k, truth, lab) in zip(axs, [('M0_turnover', 'trueTurnover', 'uncorrected turnover (%/h)'),
                                         ('M1_excess_turnover', 'trueTurnover', 'M1 excess turnover (%/h)'),
                                         ('M3_turnover', 'trueTurnover', 'M3 recalibrated turnover (%/h)'),
                                         ('M2_F30', 'truePower', 'M2 fluctuation F(5 min) (%)'),
                                         ('M2b_G30', 'truePower', 'M2b noise-free fluctuation G(5 min) (%)')]):
        for n, g in q.groupby('n'):
            ax.scatter(g[truth], g[k], color=cols[n], label=f'noise {100 * n:.1f}%', s=20)
        lim = q[truth].max()
        if truth == 'trueTurnover':
            ax.plot([0, lim], [0, lim], color='grey', ls='--', lw=1)
        ax.set_xlabel('true turnover (%/h)' if truth == 'trueTurnover' else 'true fluctuation power (rate x size^2)', fontsize=8)
        ax.set_ylabel(lab, fontsize=8); ax.set_title(lab, fontsize=9, loc='left')
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle('Simulation bench (SP 12): a clean score follows the truth and does not move with noise (colours)', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(d, 'sim_methods.png'), dpi=140); fig.savefig(os.path.join(d, 'sim_methods.pdf'))


def run_real():
    def ld(tag):
        B = pd.read_csv(os.path.join(rob, f'{tag}_bursts.csv')); B = B[(B.sign < 0) & (B.levelBefore > 0)].copy()
        B['frac'] = B.amp.abs() / B.levelBefore; return B
    Br, Bs, Bs2 = ld('base'), ld('surrpoly'), ld('surrpolysd')
    T = pd.read_csv(os.path.join(rob, 'base_traces.csv'))
    SF = pd.read_csv(os.path.join(d, 'real_sf.csv'))
    rows, boots = [], {}
    for (ds, ms), Tg in T.groupby(['dataset', 'measure']):
        sfg = SF[(SF.dataset == ds) & (SF.measure == ms)]
        SFr = {t: parse(x) for t, x in zip(sfg.trace, sfg.Dreal)}; SFs = {t: parse(x) for t, x in zip(sfg.trace, sfg.Dsurr)}
        Lh = {t: l * H for t, l in zip(Tg.trace, Tg.lengthFrames)}; hours = sum(Lh.values())
        for sp in [10, 12, 14]:
            r = Br[(Br.dataset == ds) & (Br.measure == ms) & (Br.SP == sp)]; s = Bs[(Bs.dataset == ds) & (Bs.measure == ms) & (Bs.SP == sp)]
            sr = np.array(list(SFr.values())); ss = np.array([SFs[t] for t in SFr])
            m = methods(r.frac.values, s.frac.values, sr, ss, hours)
            s2 = Bs2[(Bs2.dataset == ds) & (Bs2.measure == ms) & (Bs2.SP == sp)]
            m2 = methods(r.frac.values, s2.frac.values, None, None, hours)
            m['M1sd_excess_turnover'] = m2['M1_excess_turnover']; m['M3sd_turnover'] = m2['M3_turnover']
            m['M3sd_cut'] = m2['M3_cut']; m['nSurrSd'] = m2['nSurr']
            row = dict(dataset=ds, measure=ms, SP=sp, **m)
            if sp == 12:
                R = {t: g.frac.values for t, g in r.groupby('trace')}; S = {t: g.frac.values for t, g in s.groupby('trace')}
                bs = bootstrap(list(Tg.trace), R, S, SFr, SFs, Lh, m['M3_cut'])
                S2 = {t: g.frac.values for t, g in s2.groupby('trace')}
                bs2 = bootstrap(list(Tg.trace), R, S2, SFr, SFs, Lh, m['M3sd_cut'])
                bs['M1sd_excess_turnover'] = bs2['M1_excess_turnover']; bs['M3sd_turnover'] = bs2['M3_turnover']
                for k in SC + ['M2_F6', 'M2_F60', 'M2b_G6', 'M2b_G60', 'M2b_Gsurr30', 'M1_excess_rate']:
                    row[k + '_lo'], row[k + '_hi'] = np.nanpercentile(bs[k], [2.5, 97.5]); row[k + '_se'] = np.nanstd(bs[k])
            rows.append(row)
        print(ds, ms, 'done', flush=True)
    R = pd.DataFrame(rows); R.to_csv(os.path.join(d, 'real_scores.csv'), index=False)
    main = lambda ds: 'mean' if ds.startswith('protein') else 'spotsum'
    R12 = R[R.SP == 12]
    M = R12[R12.apply(lambda x: x.measure == main(x.dataset), axis=1)].set_index('dataset').reindex(ORDER)
    pd.set_option('display.width', 250)
    show = pd.DataFrame({k: [f"{M.loc[i, k]:.3g} [{M.loc[i, k + '_lo']:.3g}-{M.loc[i, k + '_hi']:.3g}]" for i in M.index]
                         for k in SC + ['M2b_G6', 'M2b_G60', 'M2b_Gsurr30']}, index=M.index)
    show['M3_cut'] = M.M3_cut; show['M3_n'] = M.M3_n
    show['nReal/nSurr/nSurrSd'] = (M.nReal.astype(int).astype(str) + '/' + M.nSurr.astype(int).astype(str) + '/' + M.nSurrSd.astype(int).astype(str))
    print('\n== SP 12, main measure: score [95% granule CI]'); print(show.to_string())
    C = []
    for a, b in COMPS:
        row = dict(compare=f'{a} vs {b}')
        for k in SC:
            row[k] = (M.loc[a, k] - M.loc[b, k]) / np.hypot(M.loc[a, k + '_se'], M.loc[b, k + '_se'])
        C.append(row)
    C = pd.DataFrame(C); C.to_csv(os.path.join(d, 'real_comparisons.csv'), index=False)
    print('\n== comparisons, z (|z| > 2 ~ p < 0.05)'); print(C.round(1).to_string(index=False))
    stab = []
    for k in SC:
        Mm = R[R.apply(lambda x: x.measure == main(x.dataset), axis=1)]
        cv = Mm.groupby('dataset')[k].agg(lambda v: np.nanstd(v) / abs(np.nanmean(v)) if np.nanmean(v) else np.nan).median()
        a = R12[R12.measure == 'mean'].set_index('dataset')[k]; b = R12[R12.measure == 'spotsum'].set_index('dataset')[k]
        com = b.index.intersection(a.index)
        rho = srho(a[com], b[com])
        relw = ((M[k + '_hi'] - M[k + '_lo']) / M[k].abs()).median()
        stab.append(dict(method=k, medianCV_SP10_14=cv, rankAgreement_mean_vs_spotsum=rho, medianRelCIwidth=relw))
    St = pd.DataFrame(stab); St.to_csv(os.path.join(d, 'real_stability.csv'), index=False)
    print('\n== stability'); print(St.round(3).to_string(index=False))
    COL = ['#d62a2a', '#2a5bd6', '#f08a8a', '#8aa8f0', '#1b8a4a', '#7fcf9f', '#7a4fd6', '#a98be6', '#d3c3f5']
    fig, axs = plt.subplots(1, len(SC), figsize=(4.4 * len(SC), 4.6))
    for ax, k in zip(axs, SC):
        for i, ds in enumerate(ORDER):
            v, lo, hi = M.loc[ds, k], M.loc[ds, k + '_lo'], M.loc[ds, k + '_hi']
            if np.isfinite(v):
                ax.errorbar(i, v, yerr=[[max(v - lo, 0)], [max(hi - v, 0)]], fmt='o', color=COL[i], capsize=3)
        ax.axhline(0, color='grey', lw=0.8)
        ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([o.replace('protein_RNA_', 'iv ') for o in ORDER], rotation=45, fontsize=7)
        ax.set_title(k, fontsize=10, loc='left')
    fig.suptitle('Real data, SP 12 (spot sum in vivo, mean in vitro), 95% granule-bootstrap CIs', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(d, 'real_methods.png'), dpi=140); fig.savefig(os.path.join(d, 'real_methods.pdf'))


if __name__ == '__main__':
    {'sim': run_sim, 'real': run_real}[sys.argv[1]]()
