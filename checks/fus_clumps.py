# Is PP7-FUS "stickier" (small clumps everywhere in the cell)? Two checks on the raw red movies (lab PC):
# (1) texture of the cytoplasm: fine-grained static structure, measured as the covariance between the
#     high-passed means of the even and of the odd frames (photon/camera noise is independent between
#     frames and cancels; static grain stays), relative to the cytoplasm brightness; glass as control.
# (2) granules per cell area: tracked granules (Mosaic tracks present in the first 40 frames) per um^2 of
#     fluorescent cell area (cells cannot be outlined one by one: the FUS 4x movies have no DIC channel).
# Output: analysis_outputs/fus_clumps/movies.csv, summary.txt, masks.png. Section 87.
import numpy as np, os, glob, scipy.io as sio
from PIL import Image
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

ROOT = '../zigzag_scratch/reextract/'
OUT = 'analysis_outputs/fus_clumps/'
NF, OFFSET, PX = 40, 200, 0.16            # frames used, camera offset, pixel size (um)
SETS = ['4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS']


def frames(f):
    im = Image.open(f); out = []
    for i in range(NF):
        im.seek(i); out.append(np.array(im, dtype=float))
    return np.array(out) - OFFSET


def otsu(v):
    h, e = np.histogram(v, 256); c = (e[:-1] + e[1:]) / 2; w = np.cumsum(h); m = np.cumsum(h * c)
    b = (m[-1] * w - m * w[-1]) ** 2 / np.maximum(w * (w[-1] - w), 1)
    return c[np.argmax(b)]


def movie(ds, name):
    S = frames(ROOT + ds + '/track/' + name + '.tif')
    M = S.mean(0); A = S[0::2].mean(0); B = S[1::2].mean(0)
    sm = ndi.gaussian_filter(M, 2)
    # granules: tracked positions in these frames, and anything very bright; with a margin for the blur
    xy = sio.loadmat(ROOT + ds + '/work/' + name + '_xy.mat')
    G = np.zeros(M.shape, bool); nTr = 0
    for fr, x, y in zip(xy['frame'][0], xy['x'][0], xy['y'][0]):
        fr, x, y = fr.ravel(), x.ravel(), y.ravel()
        k = fr < NF
        if k.any():
            nTr += 1
            G[np.clip(np.round(y[k]).astype(int), 0, M.shape[0] - 1), np.clip(np.round(x[k]).astype(int), 0, M.shape[1] - 1)] = True
    # cells vs glass: Otsu on the log of the smoothed image, granules left out
    lv = np.log(np.maximum(sm, 1))
    t = otsu(lv[~ndi.binary_dilation(G, iterations=6)])
    cell = ndi.binary_fill_holes(ndi.binary_opening(lv > t, iterations=2))
    cyto0 = cell & ~ndi.binary_dilation(G, iterations=6)
    bright = sm > np.median(sm[cyto0]) + 4 * np.std(sm[cyto0])
    gran = ndi.binary_dilation(G | bright, iterations=6)
    cyto = ndi.binary_erosion(cell, iterations=3) & ~gran          # away from the cell edge and the granules
    glass = ndi.binary_erosion(~ndi.binary_dilation(cell, iterations=8), iterations=2)
    HA = A - ndi.gaussian_filter(A, 2); HB = B - ndi.gaussian_filter(B, 2)
    cov = lambda m: np.mean((HA[m] - HA[m].mean()) * (HB[m] - HB[m].mean()))
    noise = lambda m: np.mean((HA[m] - HB[m]) ** 2) / 2                    # frame-to-frame part, for scale
    r = dict(dataset=ds, movie=name, cyto_mean=M[cyto].mean(), glass_mean=M[glass].mean() if glass.sum() > 500 else np.nan,
             cyto_px=int(cyto.sum()), cell_um2=cell.sum() * PX ** 2, tracks=nTr,
             cov_cyto=cov(cyto), cov_glass=cov(glass) if glass.sum() > 500 else np.nan, noise_cyto=noise(cyto))
    r['grain_cv'] = np.sqrt(max(r['cov_cyto'] - (r['cov_glass'] if np.isfinite(r['cov_glass']) else 0), 0)) / r['cyto_mean']
    r['granules_per_100um2'] = 100 * nTr / r['cell_um2']
    return r, (M, cell, cyto, gran)


def perm(a, b, n=20000, rng=np.random.default_rng(1)):
    obs = np.median(a) - np.median(b); v = np.r_[a, b]; k = len(a)
    c = sum(abs(np.median(p[:k]) - np.median(p[k:])) >= abs(obs) - 1e-12 for p in (rng.permutation(v) for _ in range(n)))
    return (c + 1) / (n + 1)


if __name__ == '__main__':
    import pandas as pd
    rows, show = [], {}
    for ds in SETS:
        for f in sorted(glob.glob(ROOT + ds + '/work/*_xy.mat')):
            name = os.path.basename(f)[:-7]
            r, ims = movie(ds, name); rows.append(r)
            if ds not in show: show[ds] = (name, ims)
            print(ds, name, round(r['grain_cv'], 3), round(r['granules_per_100um2'], 1), flush=True)
    T = pd.DataFrame(rows); T.to_csv(OUT + 'movies.csv', index=False)
    L = []
    for r in ['4', '30']:
        P, F = T[T.dataset == r + 'PP7'], T[T.dataset == r + 'PP7-FUS']
        L.append(f'{r}x (movies: PP7 {len(P)}, FUS {len(F)}); medians over movies, movie permutation p')
        for k, lab in [('cyto_mean', 'cytoplasm brightness (a.u.)'), ('glass_mean', 'glass (a.u.)'),
                       ('grain_cv', 'cytoplasm grain (CV, beyond glass)'), ('noise_cyto', 'frame noise var (a.u.^2)'),
                       ('cov_glass', 'glass static var (a.u.^2)'), ('cell_um2', 'cell area per field (um^2)'),
                       ('tracks', 'granules per field'), ('granules_per_100um2', 'granules per 100 um^2 cell')]:
            a, b = P[k].dropna().values, F[k].dropna().values
            L.append(f'  {lab:36s} PP7 {np.median(a):9.3f}  FUS {np.median(b):9.3f}  ratio {np.median(b) / np.median(a):5.2f}  p {perm(b, a):.3f}')
    open(OUT + 'summary.txt', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.4))
    for a, (ds, (name, (M, cell, cyto, gran))) in zip(ax, show.items()):
        a.imshow(M, cmap='gray', vmin=np.percentile(M, 2), vmax=np.percentile(M, 97))
        a.contour(cell, [0.5], colors='c', linewidths=0.6); a.contour(gran, [0.5], colors='r', linewidths=0.6)
        a.set_title(f'{ds} {name}: cells (cyan), granules excluded (red)', fontsize=8); a.axis('off')
    plt.tight_layout(); plt.savefig(OUT + 'masks.png', dpi=110)
