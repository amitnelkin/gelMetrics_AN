# GelMetrics revision: results log

Branch `revision`. Each section says what was run, the key numbers, what they mean, and which reviewer point they bear on. The per-stage outputs are in `analysis_outputs/<stage>/`. The cumulative λ/k0 table is `analysis_outputs/lambda_k0_before_after.md` (and `.csv`).

## How things were run (applies to every section)
- MATLAB R2025b, run headless (`matlab -batch`).
- **MATLAB figure rendering hangs on this machine** in batch and in normal mode (`exportgraphics`, `print`, `drawnow`, `heatmap`). Figures are therefore drawn in Python from saved data:
  - `checks/gm_noplot.m` runs a copy of `modifiedPoissonFit.m` with only its plotting section removed.
  - `checks/gm_strip_plot_lines.m` removes the runners' single-line plot calls.
  - The computation is unchanged.
- Runners are run without editing them:
  - `checks/run_runners.m` drives `run_experimental.m` / `run_simulation.m`, which run the script text with the dataset load injected.
  - Settings are the per-dataset paper values from `ExperimentalAnalysis.m` line 12 (SP / GOF):

    | Datasets | SP | GOF |
    |---|---|---|
    | in vitro 1:1 | 12 | 15 |
    | in vitro 10:1, 100:1 | 12 | 10 |
    | 4PP7, 4PP7-FUS | 10 | 15 |
    | 30PP7, 30PP7-FUS | 12 | 10 |
    | 5Qb, 10Qb | 14 | 15 |

  - The simulation runner uses `rng(1)`.
- Check scripts t1, t3 and t4 take `GM_SCORER=orig|patched`:
  - `orig` = `checks/gm_mpfit.m`, the original scoring.
  - `patched` = the repo's current `modifiedPoissonFit.m`.

## 1. Baseline reproduces the manuscript (untouched code)
`analysis_outputs/baseline/`. The wrappers reproduce every negative-burst number in MS.pdf. The burst counts and λ match exactly. k0 matches to within 2 (paper 63 vs 62 for 10:1, 110 vs 112 for 10Qb).

| Neg. bursts | paper λ / k0 (n) | baseline λ / k0 (n bursts) |
|---|---|---|
| in vitro 1:1 | 1 / 39 (91) | 1 / 39 (91) |
| in vitro 10:1 | 1 / 63 (96) | 1 / 62 (96) |
| in vitro 100:1 | 1 / 139 (103) | 1 / 139 (103) |
| 4PP7 | 1 / – (63) | 1 / 188 (63) |
| 4PP7-FUS | 2 / – (210) | 2 / 51 (210) |
| 30PP7 | 1 / – (112) | 1 / 192 (112) |
| 30PP7-FUS | 1 / – (86) | 1 / 97 (86) |
| 5Qb | 4 / 90 (134) | 4 / 90 (134) |
| 10Qb | 4 / 110 (181) | 4 / 112 (181) |

The manuscript fits **negative bursts only**. Positive-burst fits are in the table for completeness. They often disagree with the negative ones; for example, 5Qb is λ=1 / k0=270 for positive bursts and λ=4 / k0=90 for negative ones.

The baseline check scripts (t1–t4, original scoring) confirm the indicative Python-port numbers in CLAUDE.md:
- **Noise floor (t2a).** False bursts per flat trace were 4.64, 3.09, 1.34 and 0.93 at SP = 10, 12, 14 and 15. A real λ=2 trace gives 4.55 at SP=15.
- **Bleaching asymmetry (t2b).** With the original simulation at expConst 0.005, symmetric bursts gave mean positive amplitude 645 vs mean negative 107. With the fixed simulation: 172 vs 164.
- **Degeneracy (t1).** At true λ=2 and k0=60, recovered pairs ranged from (1, 116) to (4, 42) across seeds and signs.
- **Granule bootstrap (t4, 40 resamples).**
  - 30PP7-FUS: λ=1 in 35/40, k0 95% interval 61–168.
  - 30PP7: λ=1 in 35/40, k0 95% interval 63–400.
  - The intervals overlap.

## 2. Data-quality findings (raw-data facts)
`analysis_outputs/baseline/tiff_check/`, `analysis_outputs/baseline/data_quality/`.

- **Camera dark offset.** ExampleTIF.tif is an ImageJ hyperstack with 2 channels × 361 frames. The camera metadata was not carried over by bfconvert.
  - Channel 1, over 25 sampled frames:
    - per-frame minimum 193–200;
    - 0.1th percentile 222–230;
    - no pixel below 193.
  - This is consistent with an offset of **about 200** and rules out 100. Read noise puts a few pixels just below the offset, so the true value is about 195–200. A closed-shutter dark frame would pin it exactly.
  - Channel 2 never drops below about 1400, so it has no dark pixels.
- **Frame count.** Movies have 361 frames, not 360. Stored traces are at most 360 long, consistent with Mosaic frames 0–360 with the first point dropped by `Collect_intensity_data.m`.
- **Frame alternation in PP7/PCP data.** `generate_movArray.m` reads every TIFF page as a frame. A 2-channel interleaved TIFF would therefore make traces alternate between channels from frame to frame.
  - Per trace, the fraction of consecutive differences that flip sign is about 0.66 for noise and about 1.0 for strict alternation.

  | Dataset | traces ≥180 frames | flagged | of which alt > 0.95 |
  |---|---|---|---|
  | 4PP7 | 67 | 38 | 27 |
  | 30PP7 | 100 | 29 | 21 |
  | 30PP7-FUS | 88 | 22 | 16 |
  | 4PP7-FUS, 5Qb, 10Qb, all in vitro | – | 0 | 0 |

  - The cause cannot be established without the TIFFs: a different stack layout for those movies, alternating illumination, or something else. These traces enter every PP7 result. **This needs your input** (see the report).
- **Bleaching** (`checks/fit_bleaching.py`, `analysis_outputs/after_B3_bleach/bleach_summary.csv`). Single exponential fit to background traces ≥180 frames, offset 200 subtracted:
  - Qβ: 0.0087/min (about 40% over 60 min).
  - PP7 in vivo: about 0.001/min (about 5%).
  - in vitro: 0.0026–0.011/min.
  - Median of the dataset medians: **0.0026/min (14.5%)**, used for the simulation.
  - The original simulation used 0.0005 (3%).
