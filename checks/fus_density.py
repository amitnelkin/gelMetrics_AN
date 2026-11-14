# Is PP7-FUS dimmer because of quenching (denser) or less protein? Per granule (track medians):
# spot mean (= brightness per pixel, the 'level' of Fig. 9(d)), cell background around it, enrichment
# (spot mean / cell), spot area (px) and spot sum (total light above the cell). Movie-level permutation p.
import numpy as np, scipy.io as sio
R = 'in_vivo_experiments/'; A = 'in_vivo_reextracted/after_R2_sums/'
def load(c):
    g = lambda f, k: [np.ravel(x).astype(float) for x in sio.loadmat(f)[k][0] if np.size(x)]
    I, B = g(R + c + '/data.mat', 'intensity'), g(R + c + '/background.mat', 'background')
    S, Ar = g(R + c + '/spotsum.mat', 'spotSum'), g(A + c + '_spotarea.mat', 'spotArea')
    mv = sio.loadmat(R + c + '/movies.mat')['movieId'].ravel()
    assert len(I) == len(B) == len(S) == len(Ar) == len(mv), c
    keep = np.array([len(x) >= 180 for x in I])          # traces of at least 30 min, as the analysis
    I, B, S, Ar = [[x for x, k in zip(L, keep) if k] for L in (I, B, S, Ar)]; mv = mv[keep]
    med = lambda L: np.array([np.nanmedian(x) for x in L])
    i, b = med(I), med(B)
    return dict(mean=i, cell=b, enrich=i / b, area=med(Ar), sum=med(S)), mv
def perm(x1, m1, x2, m2, n=20000, rng=np.random.default_rng(1)):
    # movie medians, shuffle movie labels
    a = np.array([np.median(x1[m1 == u]) for u in np.unique(m1)]); b = np.array([np.median(x2[m2 == u]) for u in np.unique(m2)])
    obs = np.median(a) - np.median(b); allv = np.r_[a, b]; k = len(a); c = 0
    for _ in range(n):
        p = rng.permutation(allv); c += abs(np.median(p[:k]) - np.median(p[k:])) >= abs(obs) - 1e-12
    return (c + 1) / (n + 1)
keys = [('mean', 'brightness per pixel (a.u.)'), ('cell', 'cell around it (a.u.)'), ('enrich', 'enrichment (granule / cell)'),
        ('area', 'granule area (px)'), ('sum', 'total light (a.u.)')]
for r in ['4', '30']:
    (P, mp), (F, mf) = load(r + 'PP7'), load(r + 'PP7-FUS')
    print(f'\n{r}x  (granules: PP7 {len(mp)}, FUS {len(mf)}; movies {len(set(mp))}, {len(set(mf))})')
    for k, lab in keys:
        print(f'  {lab:30s} PP7 {np.median(P[k]):9.2f}   FUS {np.median(F[k]):9.2f}   ratio {np.median(F[k]) / np.median(P[k]):5.2f}   p(movies) {perm(F[k], mf, P[k], mp):.3f}')
