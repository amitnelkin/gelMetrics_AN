"""VITRO_MOVIES  In vitro movies left out of the analysis (RESULTS.md section 43).
The original authors dropped these movies on purpose (the user, 2026-10-25; reason to
be added), and the stored in_vitro_experiments traces contain none of their granules
(section 42). Movie ids = file order (001, 002, ...) as in the extraction_info CSVs."""
DROP = {'protein_RNA_10_1': [5], 'protein_RNA_100_1': [5, 7]}
FOLDER = {'protein_RNA_1_1': 'vitro_1000ngul_0p1mgml', 'protein_RNA_10_1': 'vitro_1000ng_10ug',
          'protein_RNA_100_1': 'vitro_1000ngul_1mgml', 'extra_100ngul_1mgml': 'vitro_100ngul_1mgml'}


def dropped(dataset, movie):
    """True if this movie (int id, or file name such as '005') is left out."""
    return int(movie) in DROP.get(dataset, [])


def keep_bursts(bursts, traces):
    """Remove bursts of granules in dropped movies (detection_grid bursts.csv + traces.csv)."""
    t = traces[[not dropped(d, m) for d, m in zip(traces.dataset, traces.movie)]]
    key = set(zip(t.dataset, t.trace))
    return bursts[[(d, r) in key for d, r in zip(bursts.dataset, bursts.trace)]]
