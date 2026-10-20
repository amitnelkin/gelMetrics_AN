"""MS_FIG13  Manuscript Fig. 13 recreated for every scoring method, and a test
of whether the true lambda was a close runner-up (RESULTS.md section 31).
Data: checks/ms_sim_test.m (amplitudes.csv: lambda 0.5-5 at noise 40 = SNR 1.5
and noise 4; amplitudes_snr.csv: lambda 2 at noise 12 and 120 = SNR 5 and 0.5,
k0 = 60, as the rows of MS Fig. 13). Negative bursts, nothing removed.
Model: Eq. (1.9) as written (checks/ms_eq19_fit.py). Candidate lambda = 1..5 as
in the manuscript; for each, the score is minimised over k0 (2-600, log grid).
Scores: MSE (10-bin density), NLL, KS, CvM, QQ RMSE, WS (the manuscript's rank
score, ranks over the integer-lambda grid).
Heatmaps: per score, (a) the score at each candidate lambda (median over the 5
datasets; NLL shown as NLL - row best) and (b) the best k0 (median), like MS
Fig. 13. Box: true lambda; star: winner.
Close runner-up test: 100 granule (trace) bootstrap resamples per dataset
(5 datasets x 100 = 500 per row), refitted: how often each candidate lambda
wins. For NLL also the likelihood-ratio statistic 2*(NLL(true) - NLL(best)).
Usage: python checks/ms_fig13.py
Writes analysis_outputs/ms_sim_test/{fig13_scores.csv, fig13_choice.csv, fig13_heatmaps.png/.pdf, fig13_choice.png/.pdf}."""
import os
import numpy as np
import pandas as pd
from scipy.special import gammaln
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'ms_sim_test')
rng = np.random.default_rng(31); NBOOT = 100
A = pd.concat([pd.read_csv(os.path.join(d, 'amplitudes.csv')), pd.read_csv(os.path.join(d, 'amplitudes_snr.csv'))])
LAMC = np.array([1, 2, 3, 4, 5.]); K0S = np.exp(np.linspace(np.log(2), np.log(600), 140))
SC = ['MSE', 'NLL', 'KS', 'CvM', 'QQ', 'WS']
ROWS = [(0.5, 40, 'λ = 0.5, SNR 1.5'), (1, 40, 'λ = 1, SNR 1.5'), (2, 40, 'λ = 2, SNR 1.5'), (3, 40, 'λ = 3, SNR 1.5'),
        (4, 40, 'λ = 4, SNR 1.5'), (5, 40, 'λ = 5, SNR 1.5'), (2, 12, 'λ = 2, SNR 5'), (2, 120, 'λ = 2, SNR 0.5')]
XG = {}
for lam in LAMC:
    xg = np.linspace(0, lam + 15 * np.sqrt(lam) + 30, 6000)
    pdf = np.exp(xg * np.log(lam) - lam - gammaln(xg + 1)); Z = np.trapezoid(pdf, xg)
    cdf = np.concatenate([[0], np.cumsum((pdf[1:] + pdf[:-1]) / 2 * np.diff(xg))]) / Z
    XG[lam] = (xg, pdf / Z, cdf)


def scores(a):
    """score[s] : len(LAMC) x len(K0S)"""
    a = np.sort(a); n = len(a); pq = (np.arange(1, n + 1) - 0.5) / n
    cnt, edges = np.histogram(a, bins=10); ctr = (edges[:-1] + edges[1:]) / 2; dens = cnt / (n * np.diff(edges))
    S = {s: np.zeros((len(LAMC), len(K0S))) for s in ['MSE', 'MSEpk', 'NLL', 'KS', 'CvM', 'QQ']}
    hi_, lo_ = np.arange(1, n + 1) / n, np.arange(0, n) / n
    for i, lam in enumerate(LAMC):
        xg, pdf, cdf = XG[lam]; X = a[None, :] / K0S[:, None]
        F = np.interp(X, xg, cdf, right=1.0); f = np.interp(X, xg, pdf, right=0.0) / K0S[:, None]
        S['NLL'][i] = -np.log(f + 1e-300).sum(1)
        S['KS'][i] = np.maximum(np.abs(F - hi_).max(1), np.abs(F - lo_).max(1))
        S['CvM'][i] = 1 / (12 * n) + ((F - (2 * np.arange(1, n + 1) - 1) / (2 * n)) ** 2).sum(1)
        qx = np.interp(pq, cdf, xg); S['QQ'][i] = np.sqrt(((a[None, :] - K0S[:, None] * qx[None, :]) ** 2).mean(1))
        m = np.interp(ctr[None, :] / K0S[:, None], xg, pdf, right=0.0) / K0S[:, None]
        S['MSE'][i] = ((m - dens[None, :]) ** 2).mean(1)
        mp = m / m.max(1, keepdims=True).clip(1e-300) * cnt.max(); S['MSEpk'][i] = ((mp - cnt[None, :]) ** 2).mean(1)
    rk = lambda M: (np.argsort(np.argsort(M.ravel())).reshape(M.shape) + 1) / M.size
    S['WS'] = 0.5 * rk(S['MSEpk']) + 0.25 * rk(S['QQ']) + 0.25 * rk(S['KS'])
    return S


rec, choice = [], []
for lam, nz, lab in ROWS:
    g = A[(A['lambda'] == lam) & (A.noise == nz)]
    for rep, gr in g.groupby('rep'):
        S = scores(gr.amp.values)
        for s in SC:
            for i, lc in enumerate(LAMC):
                j = np.argmin(S[s][i]); rec.append(dict(row=lab, lambdaTrue=lam, noise=nz, rep=rep, score=s, lamCand=lc, value=S[s][i, j], k0=K0S[j]))
        tr = [t.amp.values for _, t in gr.groupby('trace')]
        for b in range(NBOOT):
            a = np.concatenate([tr[k] for k in rng.integers(0, len(tr), len(tr))]); Sb = scores(a)
            for s in SC:
                choice.append(dict(row=lab, lambdaTrue=lam, noise=nz, rep=rep, boot=b, score=s, chosen=LAMC[np.argmin(Sb[s].min(1))]))
    print(lab, 'done', flush=True)
R = pd.DataFrame(rec); R.to_csv(os.path.join(d, 'fig13_scores.csv'), index=False)
C = pd.DataFrame(choice); C.to_csv(os.path.join(d, 'fig13_choice.csv'), index=False)

labs = [r[2] for r in ROWS]
med = R.groupby(['row', 'score', 'lamCand']).agg(value=('value', 'median'), k0=('k0', 'median')).reset_index()
nll = R[R.score == 'NLL'].copy(); nll['d'] = nll.value - nll.groupby(['row', 'rep']).value.transform('min')
dn = nll.groupby(['row', 'lamCand']).d.median()
win = R.loc[R.groupby(['row', 'rep', 'score']).value.idxmin()].groupby(['row', 'score']).lamCand.agg(lambda x: x.mode().iloc[0])
pch = C.groupby(['row', 'score']).apply(lambda x: pd.Series({f'P({int(l)})': (x.chosen == l).mean() for l in LAMC}), include_groups=False)
ptrue = C.assign(ok=C.chosen == C.lambdaTrue).groupby(['row', 'score']).ok.mean().unstack().reindex(labs)[SC]
pd.set_option('display.width', 250)
print('\nP(true lambda chosen), 500 bootstrap fits per row:'); print(ptrue.round(2).to_string())
print('\nmost often chosen lambda:'); print(C.groupby(['row', 'score']).chosen.agg(lambda x: x.mode().iloc[0]).unstack().reindex(labs)[SC].to_string())
print('\nNLL: 2*(NLL(lambda) - NLL(best)), median over datasets:'); print((2 * dn).unstack().reindex(labs).round(1).to_string())
pch.to_csv(os.path.join(d, 'fig13_choice_summary.csv'))

fig, axs = plt.subplots(len(SC), 2, figsize=(15, 4.0 * len(SC)))
for r_, s in enumerate(SC):
    M = med[med.score == s].pivot(index='row', columns='lamCand', values='value').reindex(labs)
    K = med[med.score == s].pivot(index='row', columns='lamCand', values='k0').reindex(labs)
    show = (2 * dn).unstack().reindex(labs) if s == 'NLL' else M
    rel = show.values - show.values.min(1, keepdims=True) if s == 'NLL' else show.values / show.values.min(1, keepdims=True)
    for c_, (Z, lab_, txt) in enumerate([(rel, f'(a) {s}' + (' : 2·ΔNLL vs row best' if s == 'NLL' else ' (colour: relative to row best)'), show.values),
                                          (K.values, f'(b) k0 at best fit, {s}', K.values)]):
        ax = axs[r_, c_]
        im = ax.imshow(np.log10(np.maximum(Z if c_ == 1 else (Z + (1 if s == 'NLL' else 0)), 1e-12)), cmap='viridis_r' if c_ == 0 else 'viridis', aspect='auto')
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                v = txt[i, j]; t = f'{v:.0f}' if (c_ == 1 or s == 'NLL') else f'{v:.3g}'
                ax.text(j, i, t, ha='center', va='center', fontsize=8, color='w' if (c_ == 0 and rel[i, j] > np.median(rel)) or (c_ == 1 and Z[i, j] < np.median(Z)) else 'k')
            lt = ROWS[i][0]
            if lt in LAMC:
                ax.add_patch(plt.Rectangle((list(LAMC).index(lt) - 0.5, i - 0.5), 1, 1, fill=False, ec='red', lw=2.5))
            jb = int(np.argmin(rel[i])); ax.plot(jb + 0.38, i - 0.3, marker='*', color='gold', ms=11, mec='k')
        ax.set_xticks(range(5)); ax.set_xticklabels([f'λest={int(l)}' for l in LAMC]); ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8)
        ax.set_title(lab_, fontsize=10, loc='left')
fig.suptitle('MS Fig. 13 recreated for each scoring method (median of 5 simulated datasets per row). Red box: true λ; star: best', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.985]); fig.savefig(os.path.join(d, 'fig13_heatmaps.png'), dpi=120); fig.savefig(os.path.join(d, 'fig13_heatmaps.pdf'))

fig, axs = plt.subplots(1, len(SC), figsize=(22, 4.6), sharey=True)
for k, s in enumerate(SC):
    P = np.array([[((C.row == lb) & (C.score == s) & (C.chosen == l)).sum() / ((C.row == lb) & (C.score == s)).sum() for l in LAMC] for lb in labs])
    axs[k].imshow(P, cmap='Blues', vmin=0, vmax=1, aspect='auto')
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            axs[k].text(j, i, f'{P[i, j]:.2f}', ha='center', va='center', fontsize=8, color='w' if P[i, j] > 0.5 else 'k')
        lt = ROWS[i][0]
        if lt in LAMC:
            axs[k].add_patch(plt.Rectangle((list(LAMC).index(lt) - 0.5, i - 0.5), 1, 1, fill=False, ec='red', lw=2.5))
    axs[k].set_xticks(range(5)); axs[k].set_xticklabels([str(int(l)) for l in LAMC]); axs[k].set_xlabel('chosen λ'); axs[k].set_title(s, fontsize=11)
axs[0].set_yticks(range(len(labs))); axs[0].set_yticklabels(labs, fontsize=9)
fig.suptitle('How often each λ is chosen when the granules are resampled (500 fits per row). Red box: the true λ', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(d, 'fig13_choice.png'), dpi=130); fig.savefig(os.path.join(d, 'fig13_choice.pdf'))
