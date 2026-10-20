"""LAMBDA_K0_SHAPE  Shape-only fit of lambda and k0 (RESULTS.md section 26).
Model: the manuscript's model curve exactly as modifiedPoissonFit.m builds it
(after the B1 zero truncation): the Poisson pmf at whole numbers of units n,
P(0) = 0, joined by straight lines, evaluated at A/k0; no noise term. It is
turned into a probability density by normalising over [0, A90] (A90 = 90th
percentile of the dataset's amplitudes). Bounds: lambda 0.05-60, k0 A90/300 to
3*A90 (a fit at a bound is reported as such). (A continuous Gamma-function form
was tried first: it degenerates to an exponential, lambda -> 0 and k0 -> inf,
for every distribution that falls from small to large amplitudes.)
Fitted to every negative-burst amplitude except the largest 10% (the
manuscript's top-10% trim; no lower truncation), by maximum likelihood: the
normalisation removes the overall height, so only the shape of the distribution
is fitted (no binning, no WS ranks).
Uncertainty: granule bootstrap, nested in movies in vivo (200x, refit each);
p values for condition differences: movie-level permutation (1000x, the trim
recomputed in every permuted group, refit each).
Simulations with known lambda (k0 = 60) at three noise levels show what the
shape-fitted lambda measures.
Data: analysis_outputs/lambda_k0_shape/amplitudes.csv (checks/export_amps_shape.m).
Usage: python checks/lambda_k0_shape.py [--plot] [--notrim]   (--notrim: no top-10% trim; output lambda_k0_shape_notrim/)   (--plot: redraw the figure only)
Writes analysis_outputs/lambda_k0_shape/{fits.csv, fits.md, comparisons.csv, shape_fit.png/.pdf}."""
import os
import sys
import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.optimize import minimize
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NOTRIM = '--notrim' in sys.argv   # fit every burst: no top-10% trim either (section 27)
src = os.path.join('analysis_outputs', 'lambda_k0_shape')
out = src + ('_notrim' if NOTRIM else ''); os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(27)
A = pd.read_csv(os.path.join(src, 'amplitudes.csv'))
VIVO = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb']
VITRO = ['protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1']
COMPS = [('4PP7-FUS', '4PP7'), ('30PP7-FUS', '30PP7'), ('5Qb', '4PP7'), ('10Qb', '30PP7'),
         ('protein_RNA_10_1', 'protein_RNA_1_1'), ('protein_RNA_100_1', 'protein_RNA_1_1')]
LAB = {'4PP7': 'PP7, 4x', '4PP7-FUS': 'PP7-FUS, 4x', '30PP7': 'PP7, 30x', '30PP7-FUS': 'PP7-FUS, 30x', '5Qb': 'Qβ, 5x', '10Qb': 'Qβ, 10x',
       'protein_RNA_1_1': 'in vitro 1:1', 'protein_RNA_10_1': 'in vitro 10:1', 'protein_RNA_100_1': 'in vitro 100:1'}
COL = {'4PP7': 'tab:blue', '30PP7': 'tab:blue', '4PP7-FUS': 'tab:green', '30PP7-FUS': 'tab:green', '5Qb': 'tab:red', '10Qb': 'tab:red',
       'protein_RNA_1_1': '0.75', 'protein_RNA_10_1': '0.45', 'protein_RNA_100_1': '0.1'}
NB, NPERM = 200, 1000


NN = np.arange(0, 400)


def shape_x(x, lam):
    """the manuscript's model curve (modifiedPoissonFit.m after B1): Poisson pmf at
    whole numbers of units n, P(0) = 0 (zero truncation), joined by straight lines
    and evaluated at x = A/k0"""
    pm = np.exp(NN * np.log(lam) - lam - gammaln(NN + 1)); pm[0] = 0
    return np.interp(x, NN, pm, right=0)


def logpdf(a, lam, k0, hi):
    g = np.linspace(0, hi / k0, 2000)
    z = np.trapezoid(shape_x(g, lam), g) * k0          # integral over A of shape_x(A/k0)
    return np.log(shape_x(a / k0, lam) + 1e-300) - np.log(z)


def trim(a):
    if NOTRIM:
        return a, a.max() * 1.001
    hi = np.percentile(a, 90); return a[a <= hi], hi


def fit(a_all):
    a, hi = trim(a_all)
    lo_k, hi_k = hi / 300, 3 * hi
    def nll(th):
        lam, k0 = np.exp(th)
        if not (0.05 <= lam <= 60 and lo_k <= k0 <= hi_k):
            return 1e12
        return -logpdf(a, lam, k0, hi).sum()
    best = None
    for l0 in [0.3, 1, 3, 10]:
        for f in [0.7, 1.5]:
            r = minimize(nll, [np.log(l0), np.log(f * np.mean(a) / max(l0, 1))], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-4, maxiter=3000))
            if best is None or r.fun < best.fun:
                best = r
    return np.exp(best.x[0]), np.exp(best.x[1]), -best.fun, hi


def groups(df):
    return {m: [g2.amp.values for _, g2 in g.groupby('trace')] for m, g in df.groupby('movie')}


def boot(G):
    keys = list(G); res = []
    for _ in range(NB):
        a = []
        for k in rng.choice(len(keys), len(keys)):
            tr = G[keys[k]]; a.extend(tr[i] for i in rng.integers(0, len(tr), len(tr)))
        l, k0, _, _ = fit(np.concatenate(a)); res.append((l, k0, l * k0))
    return np.percentile(np.array(res), [2.5, 97.5], axis=0)


PLOT_ONLY = '--plot' in sys.argv   # redraw the figure from fits.csv / comparisons.csv
if not PLOT_ONLY:
    rows = []
    for name in VIVO + VITRO + sorted(A[A.kind == 'simulated'].dataset.unique(), key=lambda s: (float(s.split()[1].rstrip(',')), int(s.split()[-1]))):
        d = A[A.dataset == name]
        lam, k0, ll, hi = fit(d.amp.values); ci = boot(groups(d))
        rows.append(dict(kind=d.kind.iloc[0], dataset=name, lambdaTrue=d.lambdaTrue.iloc[0], noise=d.noise.iloc[0], bursts=len(d), burstsFitted=int((d.amp <= hi).sum()),
                         movies=d.movie.nunique(), lam=lam, lam_lo=ci[0, 0], lam_hi=ci[1, 0], k0=k0, k0_lo=ci[0, 1], k0_hi=ci[1, 1],
                         mean=lam * k0, mean_lo=ci[0, 2], mean_hi=ci[1, 2], cv_data=np.std(trim(d.amp.values)[0]) / np.mean(trim(d.amp.values)[0])))
        print(name, f'lambda {lam:.2f} [{ci[0, 0]:.2f}, {ci[1, 0]:.2f}]  k0 {k0:.1f} [{ci[0, 1]:.1f}, {ci[1, 1]:.1f}]', flush=True)
    F = pd.DataFrame(rows); F.to_csv(os.path.join(out, 'fits.csv'), index=False)

    comp = []
    for a_, b_ in COMPS:
        Ga, Gb = groups(A[A.dataset == a_]), groups(A[A.dataset == b_])
        fa, fb = F.set_index('dataset').loc[a_], F.set_index('dataset').loc[b_]
        obs = np.array([fa.lam - fb.lam, np.log(fa.k0 / fb.k0), np.log(fa.lam * fa.k0 / (fb.lam * fb.k0))])
        movies = list(Ga.values()) + list(Gb.values()); na = len(Ga)
        cnt = np.zeros(3)
        if len(Ga) > 1 and len(Gb) > 1:   # in vivo: whole movies
            units = movies
        else:                             # in vitro: no movie ids, permute granules (traces)
            units = [[t] for g in movies for t in g]; na = sum(len(g) for g in Ga.values())
        for _ in range(NPERM):
            p = rng.permutation(len(units))
            xa = np.concatenate([t for i in p[:na] for t in units[i]]); xb = np.concatenate([t for i in p[na:] for t in units[i]])
            la, ka, _, _ = fit(xa); lb, kb, _, _ = fit(xb)
            cnt += np.abs([la - lb, np.log(ka / kb), np.log(la * ka / (lb * kb))]) >= np.abs(obs) - 1e-12
        pl, pk, pm = (cnt + 1) / (NPERM + 1)
        comp.append(dict(compare=f'{a_} vs {b_}', lamA=fa.lam, lamB=fb.lam, p_lambda=pl, k0A=fa.k0, k0B=fb.k0, p_k0=pk, lamk0A=fa.lam * fa.k0, lamk0B=fb.lam * fb.k0, p_lamk0=pm,
                         unit='movies' if len(Ga) > 1 and len(Gb) > 1 else 'granules'))
        print(comp[-1], flush=True)
    Cp = pd.DataFrame(comp); Cp.to_csv(os.path.join(out, 'comparisons.csv'), index=False)
    with open(os.path.join(out, 'fits.md'), 'w', encoding='utf-8', newline=chr(13) + chr(10)) as fh:
        fh.write('| dataset | bursts (fitted) | λ [95% CI] | k0 [95% CI] | mean burst λ·k0 [95% CI] |' + chr(10) + '|---|---|---|---|---|' + chr(10))
        for _, r in F.iterrows():
            fh.write(f"| {LAB.get(r.dataset, 'sim ' + r.dataset)} | {r.bursts} ({r.burstsFitted}) | {r.lam:.2f} [{r.lam_lo:.2f}, {r.lam_hi:.2f}] | "
                     f"{r.k0:.1f} [{r.k0_lo:.1f}, {r.k0_hi:.1f}] | {r['mean']:.0f} [{r.mean_lo:.0f}, {r.mean_hi:.0f}] |" + chr(10))
else:
    F = pd.read_csv(os.path.join(out, 'fits.csv')); Cp = pd.read_csv(os.path.join(out, 'comparisons.csv'))

# ---------------- figure
fig, axs = plt.subplots(1, 4, figsize=(23, 5.4))
ax = axs[0]
for ds, c in [('4PP7', 'tab:blue'), ('4PP7-FUS', 'tab:green')]:
    a, hi = trim(A[A.dataset == ds].amp.values); r = F.set_index('dataset').loc[ds]
    ax.hist(a, bins=np.linspace(0, hi, 13), density=True, color=c, alpha=0.35, edgecolor=c)
    xx = np.linspace(0, hi, 300); ax.plot(xx, np.exp(logpdf(xx, r.lam, r.k0, hi)), color=c, lw=2.5,
                                          label=f'{LAB[ds]}: λ {r.lam:.2f} [{r.lam_lo:.2f}-{r.lam_hi:.2f}], k0 {r.k0:.0f} [{r.k0_lo:.0f}-{r.k0_hi:.0f}]')
pc = Cp[Cp['compare'] == '4PP7-FUS vs 4PP7'].iloc[0]
ax.set_title(f"(a) shape fits, {'all bursts' if NOTRIM else 'top 10% removed'}" + chr(10) + f"FUS vs PP7: λ p = {pc.p_lambda:.2g}, k0 p = {pc.p_k0:.2g}" + (f", λ·k0 p = {pc.p_lamk0:.2g}" if 'p_lamk0' in pc else '') + ' (movie-level)', fontsize=10, loc='left')
ax.set_xlabel('burst amplitude (a.u., spot mean)'); ax.set_ylabel('probability density'); ax.legend(fontsize=7.5)
ax = axs[1]
a, hi = trim(A[A.dataset == '4PP7-FUS'].amp.values)
r = F.set_index('dataset').loc['4PP7-FUS']
LAMS = np.exp(np.linspace(np.log(r.lam / 8), np.log(r.lam * 8), 70)); K0S = np.exp(np.linspace(np.log(r.k0 / 8), np.log(r.k0 * 8), 90))
L = np.array([[logpdf(a, l, k, hi).sum() for k in K0S] for l in LAMS]); D = L - L.max()
cs = ax.contourf(K0S, LAMS, np.clip(D, -20, 0), levels=np.linspace(-20, 0, 21), cmap='viridis'); fig.colorbar(cs, ax=ax, label='log-likelihood vs best')
ax.contour(K0S, LAMS, D, levels=[-3.0], colors='w', linewidths=2); ax.plot(r.k0, r.lam, 'r*', ms=14)
ax.set_yscale('log'); ax.set_xscale('log'); ax.set_xlabel('k0'); ax.set_ylabel('λ'); ax.set_title('(b) PP7-FUS 4x: likelihood surface (white: 95% region)', fontsize=10, loc='left')
ax = axs[2]
S = F[F.kind == 'simulated']
for nz, c, lab in [(12, 'tab:purple', 'SNR 5'), (40, 'tab:orange', 'SNR 1.5'), (120, 'tab:brown', 'SNR 0.5')]:
    s = S[S.noise == nz]
    if NOTRIM:   # lambda*k0 recovery
        ax.errorbar(60 * s.lambdaTrue * (1 + (nz - 40) / 1000), s['mean'], yerr=[s['mean'] - s.mean_lo, s.mean_hi - s['mean']], fmt='o-', color=c, capsize=3, label=f'{lab} (noise {nz})')
    else:
        ax.errorbar(s.lambdaTrue * (1 + (nz - 40) / 1000), s.lam, yerr=[s.lam - s.lam_lo, s.lam_hi - s.lam], fmt='o-', color=c, capsize=3, label=f'{lab} (noise {nz})')
lim = [24, 330] if NOTRIM else [0.4, 5.5]
ax.plot(lim, lim, 'k:', lw=1, label='perfect recovery')
ax.set_xlabel('true λ·k0 (simulation, k0 = 60)' if NOTRIM else 'true λ (simulation, k0 = 60)'); ax.set_ylabel('fitted λ·k0 (shape fit)' if NOTRIM else 'fitted λ (shape fit)'); ax.legend(fontsize=8)
ax.set_title(('(c) does the shape-fitted λ·k0 track the truth?' if NOTRIM else '(c) does the shape-fitted λ track the true λ?') + chr(10) + 'same detection pipeline as the real data', fontsize=10, loc='left')
ax = axs[3]
for _, r in F[F.kind != 'simulated'].iterrows():
    ax.errorbar(r.k0, r.lam, xerr=[[r.k0 - r.k0_lo], [r.k0_hi - r.k0]], yerr=[[r.lam - r.lam_lo], [r.lam_hi - r.lam]], fmt='s' if '30' in r.dataset else 'o',
                color=COL[r.dataset], ms=8, capsize=3, label=LAB[r.dataset])
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('k0 (a.u.)'); ax.set_ylabel('λ'); ax.legend(fontsize=7, ncol=2)
ax.set_title('(d) all datasets (as Fig. 10b), 95% CI', fontsize=10, loc='left')
for ax in axs:
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.tight_layout(); fig.savefig(os.path.join(out, 'shape_fit.png'), dpi=130); fig.savefig(os.path.join(out, 'shape_fit.pdf'))
