"""PLOT_ESTIMATORS  Summaries and figures for checks/estimator_compare_sim.m and
estimator_compare_real.m (R1 2 biased estimator, R2 4a weights/trims, R2 4b).
Simulation metrics per estimator x noise (both signs pooled, lambda_true pooled):
  medRatio  median lambda_hat / lambda_true
  rmseLog   RMS of log(lambda_hat / lambda_true)   (error in factor terms)
  within25  fraction of fits with lambda_hat within +-25% of the truth
  k0Err     median |k0_hat - 60| / 60
  lk0Err    median |lambda_hat*k0_hat - lambda_true*60| / (lambda_true*60)
Usage (repo root): python checks/plot_estimators.py
Writes analysis_outputs/estimators/sim_summary.csv, sim_lambda.png/.pdf,
sim_metrics.png/.pdf, real_summary.md."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = os.path.join('analysis_outputs', 'estimators')
E = pd.read_csv(os.path.join(d, 'sim_estimates.csv'))
K0 = 60
E['ratio'] = E.lamHat / E.lamTrue
E['logerr'] = np.log(E.ratio)
E['k0err'] = (E.k0Hat - K0).abs() / K0
E['lk0err'] = (E.lamHat * E.k0Hat - E.lamTrue * K0).abs() / (E.lamTrue * K0)
order = ['MLE', 'MAP_s1', 'MAP_s05', 'WS_paper', 'WS_equal', 'WS_MSE', 'WS_QQ', 'WS_KS', 'WS_QQheavy', 'WS_KSheavy']
S = (E.groupby(['estimator', 'noise'])
     .agg(medRatio=('ratio', 'median'), rmseLog=('logerr', lambda x: np.sqrt(np.mean(x ** 2))),
          within25=('ratio', lambda x: np.mean((x > 0.75) & (x < 1.25))),
          k0Err=('k0err', 'median'), lk0Err=('lk0err', 'median'), n=('ratio', 'size'))
     .reset_index())
S['estimator'] = pd.Categorical(S.estimator, order, ordered=True)
S = S.sort_values(['noise', 'estimator'])
S.to_csv(os.path.join(d, 'sim_summary.csv'), index=False)
print(S.round(3).to_string(index=False))

INK, MUTED = '#1f1e1c', '#898781'
COL = {'MLE': '#2a78d6', 'MAP_s1': '#1baf7a', 'MAP_s05': '#7a4fd6', 'WS_paper': '#d62a2a', 'WS_equal': '#eb6834',
       'WS_MSE': '#b8860b', 'WS_QQ': '#8c564b', 'WS_KS': '#e377c2', 'WS_QQheavy': '#7f7f7f', 'WS_KSheavy': '#17becf'}
noises = sorted(E.noise.unique())
fig, axs = plt.subplots(2, len(noises), figsize=(5 * len(noises), 9), sharey='row')
for j, nz in enumerate(noises):
    for i, group in enumerate([['MLE', 'MAP_s1', 'MAP_s05', 'WS_paper'],
                               ['WS_paper', 'WS_equal', 'WS_MSE', 'WS_QQ', 'WS_KS', 'WS_QQheavy', 'WS_KSheavy']]):
        ax = axs[i, j]
        lt = np.array(sorted(E.lamTrue.unique()))
        ax.plot(lt, lt, color=MUTED, ls='--', lw=1, label='truth')
        for k, est in enumerate(group):
            g = E[(E.estimator == est) & (E.noise == nz)].groupby('lamTrue').lamHat
            off = (k - len(group) / 2) * 0.04
            ax.errorbar(lt + off, g.median().reindex(lt), yerr=[g.median().reindex(lt) - g.min().reindex(lt),
                        g.max().reindex(lt) - g.median().reindex(lt)], fmt='o-', ms=3.5, lw=1, capsize=1.5,
                        color=COL[est], label=est)
        ax.set_title(f'noise SD {nz}: ' + ('likelihood (MLE / MAP) vs paper WS' if i == 0 else 'WS weight sets'),
                     fontsize=9, loc='left', color=INK)
        ax.set_xlabel('true lambda', fontsize=8)
        ax.set_ylabel('estimated lambda (median, range over seeds x signs)', fontsize=8)
        ax.set_ylim(0, 7.5)
        ax.tick_params(labelsize=7, colors=MUTED)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        if j == 0:
            ax.legend(fontsize=7, frameon=False, loc='upper left')
fig.suptitle('Lambda recovery by estimator (true k0 = 60; generateSim -> analysis.m; 2 seeds x 2 signs per point)',
             fontsize=11, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'sim_lambda.{ext}'), dpi=140)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, (m, lab) in zip(axs, [('rmseLog', 'RMS log error of lambda (lower = better)'),
                              ('within25', 'fraction of lambda within +-25% (higher = better)'),
                              ('lk0Err', 'median rel. error of lambda*k0 (lower = better)')]):
    for k, est in enumerate(order):
        v = S[S.estimator == est].set_index('noise')[m].reindex(noises)
        ax.bar(np.arange(len(noises)) + (k - len(order) / 2) * 0.08, v.values, width=0.08, color=COL[est], label=est)
    ax.set_xticks(range(len(noises))); ax.set_xticklabels([f'noise {n}' for n in noises], fontsize=8)
    ax.set_title(lab, fontsize=9, loc='left', color=INK)
    ax.tick_params(labelsize=7, colors=MUTED)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
axs[0].legend(fontsize=6.5, frameon=False, ncol=2)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(d, f'sim_metrics.{ext}'), dpi=140)

# real data tables
out = []
if os.path.exists(os.path.join(d, 'real_map.csv')):
    M = pd.read_csv(os.path.join(d, 'real_map.csv'))
    L = pd.read_csv(os.path.join('analysis_outputs', 'mle_real', 'mle_real.csv'))
    Wt = pd.read_csv(os.path.join(d, 'real_weights.csv'))
    Tr = pd.read_csv(os.path.join(d, 'real_trims.csv'))
    C = pd.read_csv(os.path.join(d, 'real_counts.csv'))
    neg = lambda X: X[X.sign == 'Neg']
    t = neg(L)[['dataset', 'lambda', 'lamLoBoot', 'lamHiBoot', 'k0', 'k0LoBoot', 'k0HiBoot']].copy()
    t['est'] = 'MLE'
    m2 = neg(M)[['dataset', 'estimator', 'lambda', 'lamLoBoot', 'lamHiBoot', 'k0', 'k0LoBoot', 'k0HiBoot']].rename(columns={'estimator': 'est'})
    T = pd.concat([t, m2])
    T['lambda [boot]'] = T.apply(lambda r: f"{r['lambda']:.2f} [{r.lamLoBoot:.2f}-{r.lamHiBoot:.2f}]", axis=1)
    T['k0 [boot]'] = T.apply(lambda r: f"{r.k0:.0f} [{r.k0LoBoot:.0f}-{r.k0HiBoot:.0f}]", axis=1)
    out.append('### Negative bursts: MLE vs biased (MAP) estimators\n\n' +
               T.pivot(index='dataset', columns='est', values='lambda [boot]').to_markdown() + '\n\n' +
               T.pivot(index='dataset', columns='est', values='k0 [boot]').to_markdown())
    Wt['v'] = Wt['lambda'].astype(str) + '/' + Wt.k0.astype(str)
    out.append('### Negative bursts: WS weight sets (lambda/k0)\n\n' +
               neg(Wt).pivot(index='dataset', columns='weights', values='v').to_markdown())
    Tr['set'] = 'top10%: ' + Tr.topTrimPct.gt(0).map({True: 'on', False: 'off'}) + ', left trunc: ' + Tr.leftTruncPct.gt(0).map({True: 'paper', False: 'off'})
    Tr['v'] = Tr.mleLambda.round(2).astype(str) + '/' + Tr.mleK0.round(0).astype(int).astype(str) + \
        ' | ' + Tr.gridLambda.astype(str) + '/' + Tr.gridK0.astype(str)
    out.append('### Negative bursts: trims on/off (MLE lambda/k0 | WS lambda/k0)\n\n' +
               neg(Tr).pivot(index='dataset', columns='set', values='v').to_markdown())
    C['usedFrac'] = (C.afterLeftTrunc / C.bursts).round(2)
    out.append('### Data used\n\n' + C.to_markdown(index=False))
    with open(os.path.join(d, 'real_summary.md'), 'w', encoding='utf-8') as fh:
        fh.write('\n\n'.join(out) + '\n')
    print('\n\n'.join(out))
