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

## 3. Bug fixes and model changes: effect on λ and k0
**What was run.**
- The runners were rerun after every commit (per-dataset paper settings). The simulation runner uses `rng(1)`, true λ=2, k0=60.
- Changes B2–B4 touch only the simulation, so for those stages only the simulation runner was rerun. Experimental values carry forward unchanged.
- Check scripts were rerun where affected, with both scorers. Logs: `analysis_outputs/<stage>/checks_<scorer>/`.
- Full table: `lambda_k0_before_after.md`.
- Stages:
  - A = bug fixes.
  - B1 = zero-truncated Poisson.
  - B2 = symmetric simulated bursts.
  - B3 = expConst 0.0026.
  - B4 = symmetric k0 grid.
  - B5 = movmedian filloutliers.

Rank-1 λ / k0 (`*` = k0 at a grid edge):

| dataset | sign | paper | base | A | B1 | B2–B4 | B5 (final) | n base → final | λ·k0 base → final |
|---|---|---|---|---|---|---|---|---|---|
| in vitro 1:1 | Neg | 1/39 | 1/39 | 1/35 | 1/31 | = | 1/28 | 91 → 94 | 39 → 28 |
| | Pos | – | 1/53 | 1/39 | 1/35 | = | 1/33 | 71 → 77 | 53 → 33 |
| in vitro 10:1 | Neg | 1/63 | 1/62 | 1/62 | **2/39** | = | 2/42 | 96 → 93 | 62 → 84 |
| | Pos | – | 1/139 | 1/88 | 1/74 | = | 1/79 | 78 → 74 | 139 → 79 |
| in vitro 100:1 | Neg | 1/139 | 1/139 | 1/110 | **2/59** | = | 2/59 | 103 → 105 | 139 → 118 |
| | Pos | – | 1/145 | 1/121 | 2/60 | = | 2/60 | 63 → 63 | 145 → 120 |
| 4PP7 | Neg | 1/– | 1/188 | 1/199 | 1/180 | = | 1/196 | 63 → **155** | 188 → 196 |
| | Pos | – | 1/379 | 1/212 | 1/176 | = | 2/126 | 66 → **161** | 379 → 252 |
| 4PP7-FUS | Neg | 2/– | 2/51 | 2/53 | 2/53 | = | **1/83** | 210 → 209 | 102 → 83 |
| | Pos | – | 1/73 | 1/79 | 1/79 | = | 1/77 | 207 → 203 | 73 → 77 |
| 30PP7 | Neg | 1/– | 1/192 | 1/204 | 2/102 | = | 2/91 | 112 → **177** | 192 → 182 |
| | Pos | – | 1/138 | 1/135 | 2/85 | = | 1/238 | 87 → **173** | 138 → 238 |
| 30PP7-FUS | Neg | 1/– | 1/97 | 1/118 | 1/115 | = | 2/96 | 86 → 124 | 97 → 192 |
| | Pos | – | 1/117 | 1/163 | 1/144 | = | 2/108 | 86 → 124 | 117 → 216 |
| 5Qb | Neg | 4/90 | 4/90 | 3/116 | 3/116 | = | **2/155** | 134 → 134 | 360 → 310 |
| | Pos | – | 1/270 | 1/322 | 1/316 | = | 1/325 | 155 → 158 | 270 → 325 |
| 10Qb | Neg | 4/110 | 4/112 | 5/103 | **1/383** | = | 3/153 | 181 → 184 | 448 → 459 |
| | Pos | – | 3/181 | 4/148 | 4/148 | = | 1/438 | 177 → 172 | 543 → 438 |
| SIM λ=2, k0=60 | Neg | 4/– (λsim 2) | 2/63 | 4/41* | 4/41* | 3/51 → 3/52 → 4/39 | 4/39 | | 126 → 156 |
| | Pos | – | 2/72 | 1/107 | 2/64 | 2/72 → 2/74 → 2/74 | 2/72 | | 144 → 144 |

(n = negative/positive bursts found before the top-10% trim.)

**Checks, final code (B5).**
- **t1** (fixed simulation, true k0=60, 3 seeds; λ_true counts as recovered when both signs agree):
  - Original scoring recovers λ_true in 2/18 fits at λ_true ≥ 3.
  - Patched scoring:
    - λ_true = 4 and 5: correct in 11/12 fits, with k0 biased high (69–85).
    - λ_true = 2: wrong in 5/6 fits (estimates 1–4).
    - λ_true = 1: estimate 2 in 3/6 fits.
- **t2 noise floor:** false bursts per flat trace 4.73, 3.12, 1.43 and 0.98 at SP = 10, 12, 14 and 15. A real λ=2 trace gives 4.20 at SP=15. So at SP=10–12, most detected "bursts" in a λ=2-like trace could be noise.
- **t4 granule bootstrap** (40 resamples, 30PP7-FUS vs 30PP7):

  | Scoring | 30PP7-FUS: λ=1 wins | 30PP7-FUS: k0 95% | 30PP7: λ=1 wins | 30PP7: k0 95% |
  |---|---|---|---|---|
  | original | 30/40 | 76–500 | 25/40 | 280–500 |
  | patched | 8/40 | 95–245 | 30/40 | 122–358 |

  Upper ends at 500 are the grid edge. Intervals are wide and overlap.

**Interpretation.**
- The individual λ and k0 values are not stable under the fixes. λ changes for 8 of the 9 datasets on at least one sign (all except in vitro 1:1), and k0 moves by factors of up to 3.
- λ·k0 (about the mean amplitude) is much more stable. For example, 10Qb negative goes 448 → 459 while λ goes 4 → 5 → 1 → 3. This is the λ–k0 degeneracy, and it is visible directly in every score surface (`figs/*.png`, bottom rows).
- The simulation still returns λ=4 for true λ=2 (k0 39 vs 60) after all fixes. The reviewer's 2-vs-4 point is therefore not a bug artefact: the grid-score estimator cannot resolve λ=2 at this SNR.

**Manuscript claims at the final code** (negative bursts, no error bars yet):

| Claim | Status |
|---|---|
| "λ_est = 1 for all three in vitro cases" | **Fails.** 10:1 and 100:1 give λ=2 from the zero-truncation step onwards. |
| k0 increases with protein ratio (39 → 63 → 139) | k0 alone is 28 / 42 / 59, but λ differs, so the k0s are not comparable. λ·k0 still increases: 28 / 84 / 118. |
| Qβ (λ=4) differs from PP7 (λ=1) | **Much weaker.** Qβ 2–3 vs PP7 1–2. Only λ·k0 separates them clearly (Qβ ~310–460 vs PP7 ~180–200). |
| FUS shifts λ up and k0 down for the 4x RNA | **λ part fails** (4PP7 1/196 vs 4PP7-FUS 1/83). The k0 decrease holds as a point estimate. |
| No FUS effect for 30x | Holds as point estimates (2/91 vs 2/96). |

Every one of these needs the likelihood fit with CIs (analysis 2) before it can be stated.

**Reviewer points.**
- **R2 4a (thresholds and error bars) and 4b (integer λ, 2 vs 4):** supported. The current estimator is unstable and degenerate.
- **R2 6 (noise floor):** supported. The false-burst rate is comparable to the real-burst rate at SP ≤ 12.

**B5 caveat.** Local outlier filling roughly doubles the burst counts in 4PP7, 30PP7 and 30PP7-FUS, and only there. These are the same datasets flagged for frame-to-frame alternation (section 2). The global 'mean' criterion was probably flattening those traces. Their results should be treated as provisional until the alternation is explained.

## 4. Cause of the PP7 zigzag traces (lab PC, .nd2 re-tracking)
`analysis_outputs/zigzag_diagnosis/` (figure `zigzag_mechanism.png`, `per_trace.csv`, `page_stats.csv`, `nd2_dims.csv`, plus the driver scripts). The working files are in `..\zigzag_scratch\`, outside the repo.

**What was run.**
- 4PP7 movies 005 and 007 (`Raw_Results\Genescript-PCP4x-QCP5x_K133-tdPP7-mCherry\New`):
  - converted with plain `bfconvert` (as the user did);
  - tracked with the repo `BatchParticleTracker.ijm`, unchanged except for the input folder and a final exit;
  - run through the ORIGINAL `read_tables.m` + `Collect_intensity_data.m` (f7eded0, unchanged; a `uigetdir` shim supplies the folder).
- The channel count of every .nd2 in the four PP7 datasets was read with `showinf`.

**Numbers.**
- **Channels per movie:**
  - 4PP7: 005 = 1, 006-013 = 2.
  - 4PP7-FUS: all 12 = 1.
  - 30PP7: all 8 = 2.
  - 30PP7-FUS: all 8 = 2.
  - Naor 5Qb / 10Qb / 4pp7 001: 2.
- **Mosaic frame numbering.**
  - On 007 (2 channels, 722 planes), Mosaic reports frames 0-360, and track steps are 2 (4792) or 1 (1933).
  - Frame f is raw plane f+1 (interleaved red, DIC, ...). This covers only time points 0-180.
  - 4130 of the 6901 detections are on odd frames, which are DIC planes.
  - On 005 (1 channel), steps are 1 (2581) or 2 (71), with no parity structure.
- **TIFF pages:** in 007, odd pages are red (median ~605, min ~200) and even pages are DIC (median ~6900). `Collect_intensity_data.m` reads page f.
- **Three kinds of trace in a 2-channel movie (007, 40 traces):**
  - even frames only: every point read from a DIC page (median spot ~7500-9400);
  - odd frames only: detected on a DIC plane, read on a red page;
  - mixed parity: alternating red/DIC, the zigzag; 8/40 are flagged, all of this kind.
- **005: 0/19 traces flagged.**
- **Match to stored data:** the 59 re-extracted traces equal stored 4PP7 traces 1-19 (005) and 49-88 (007) **bit-exactly**.

**Interpretation.** The zigzag is a frame-to-page mapping artefact of 2-channel movies. It is not biology, and the cause is not alternating illumination. It is worse than the alternation flag suggested. In every 2-channel movie:
- no trace measures a red condensate correctly;
- even "clean" traces are pure DIC, or red at DIC-detected positions;
- every trace covers only the first ~30 min.

That makes 30PP7 and 30PP7-FUS invalid in full, and 4PP7 apart from movie 005 (19 of 360 traces). 4PP7-FUS is unaffected. Every λ/k0 number for 4PP7, 30PP7 and 30PP7-FUS in sections 1-3, and the manuscript's PP7 comparisons that use them, must be redone from re-extracted data. Open: the Qβ movies checked are also 2-channel, yet their traces show no zigzag. Either they were tracked differently or the even-only/odd-only failure applies to them silently. This must be checked before any Qβ number is trusted.

**Reviewer points.** General data validity, and R2 5 (replicates/variance): the PP7 results cannot be defended until re-extraction.
