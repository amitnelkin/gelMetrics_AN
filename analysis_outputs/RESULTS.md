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

## 5. Likelihood fitter (analysis 2): validation on simulations
`analysis/gelmetricsMLE.m` (+ `gm_mle_nll.m`, `gm_fit_dataset.m`, `gm_amplitudes_by_trace.m`, `gm_sigma_from_ampvar.m`). Checks: `checks/mle_unit_test.m`, `checks/mle_sim_validation.m`, `checks/plot_mle_validation.py`. Outputs: `analysis_outputs/mle_validation/`.

**Model.**
- Amplitude A ~ Σ_{n≥1} ZTP(n | λ) · N(n·k0, σ²), with continuous λ and k0.
- σ = robust SD (1.4826·MAD) of the quiescent-segment amplitudes (`ampVar`).
- It fits the same amplitudes as the grid method (|A| after the top-10% trim, left-truncated at the 15th percentile); both cuts are modelled in the normaliser.
- Optimisation: a 40×40 seed grid, then fmincon in log space.
- CIs come from the profile likelihood (2ΔNLL ≤ 3.84) and from a granule bootstrap that resamples traces, never bursts. σ is re-estimated in every resample.
- The grid-score method (current `modifiedPoissonFit.m` scoring) is run alongside.

**Numerical fix found during validation.** The truncation-window mass Φ(z_max) − Φ(z_min) cancelled to 0 when both z were large, which made the NLL spuriously very negative at λ→0, k0→0. It is now computed from upper tails.

**Unit test (amplitudes drawn from the model itself, 300 per fit, 5 seeds; `unit_test.csv`).**
- At k0/σ = 4: median λ̂ is within 5% of the truth for λ = 0.5 to 5 (including 1.5), and k0 is 59.6–60.4. The k0 profile CI covers in 35/35 fits.
- At k0/σ = 2: λ < 1.5 is not identified (true 0.5 gives λ̂ = 2.9, k0 = 20, along the λ·k0 ridge), but the CIs still cover it.

**Full-pipeline validation** (generateSim → analysis.m, one trace at a time → fit). 150 tracks of 361 frames; noise 40; k0 = 60; λ_true ∈ {0.5, 1, 1.5, 2, 3, 4, 5}; 3 seeds; 100 bootstrap resamples on negative bursts only. Same settings as the experimental runner, with SP = 15. Figure: `sim_validation.png`.
- **λ is resolved by the MLE and not by the grid method.**
  - The grid method returns λ = 1 in 17/21 negative and 14/21 positive fits, whatever the true λ.
  - The MLE λ rises with the truth, but is biased high: median λ̂/λ = 1.29 (negative) and 1.22 (positive).
  - At λ_true = 2: MLE 2.58, 2.37, 2.03 (negative); grid 1, 1, 3.
- **k0.** MLE median 58 (IQR 51–72) for negative bursts and 70 for positive bursts. The grid k0 absorbs λ (k0 ≈ λ·k0 when λ = 1), up to 280.
- **Coverage of the 95% CIs** (all λ):

  | | λ, profile | λ, bootstrap | k0, profile | k0, bootstrap |
  |---|---|---|---|---|
  | negative bursts | 0.38 | 0.62 | 0.24 | 0.81 |
  | positive bursts | 0.57 | – | 0.29 | – |

**Interpretation.** The likelihood fit answers the integer-λ objection: λ is estimated continuously, and it follows the true λ where the grid score cannot. But the detection pipeline does not deliver amplitudes of the form n·k0 + N(0, σ_quiescent²):
- smoothing and per-segment line fits distort amplitudes;
- σ from `ampVar` underestimates the noise on a burst amplitude.

So λ comes out about 25% high, and the profile CIs are too narrow. The granule-bootstrap CIs are closer to nominal (k0 0.81) and should be the ones reported, with the pipeline bias stated. This supports R2 4b (integer λ, 2 vs 4) for the MLE over the grid method, but it does not give exact recovery. Next to test: σ fitted as a free parameter, more seeds, and an SNR grid (analysis 3).
## 6. Re-extraction of the PP7 datasets (Step 2): effect on λ and k0
**What was run.**
- All 37 PP7 movies were processed:
  - 4PP7: 9 movies;
  - 4PP7-FUS: 12;
  - 30PP7: 8;
  - 30PP7-FUS: 8.
- Each movie was converted with plain bfconvert, re-tracked with the R1 `BatchParticleTracker.ijm` (fluorescence channel only), and re-extracted by `checks/reextract_pp7.m` at three stages, one commit each:
  - R1: fluorescence channel only, page = frame+1, xy+1;
  - R2: + offset 200;
  - R3: + 5 px aperture.
- Output: `in_vivo_reextracted/<stage>/` (new files; the repo .mat files are untouched).
- The runners were rerun with `run_runners(stage, true, false, dataFolders)` with unchanged per-dataset paper settings. The simulation and the non-PP7 datasets are unaffected, so their B5 values carry forward ('=' in the table).
- Qβ was not re-extracted: no source .nd2 (see CLAUDE.md).

**Data checks.**
- Alternating (zigzag) traces ≥ 180 frames, old → new:
  - 4PP7: 38 → 0;
  - 30PP7: 29 → 0;
  - 30PP7-FUS: 22 → 0;
  - 4PP7-FUS: 0 → 0.
- Traces ≥ 180 frames, old → new: 4PP7 67 → 76, 4PP7-FUS 90 → 94, 30PP7 100 → 145, 30PP7-FUS 88 → 111.
- Median spot level (R1, before offset):
  - 4PP7: 6903 → 5725;
  - 30PP7: 7043 → 4198;
  - 30PP7-FUS: 6782 → 3636. The old values were inflated by DIC pages.
  - 4PP7-FUS: 3242 → 3236. This is the 1-channel control, where only the ±1 mapping changed.
- **Movie 4PP7 012-recovered-after-corruption has an empty dark 2nd channel, not DIC.** The first channel rule picked the dark channel (7290 junk tracks, 0 traces). It was replaced by a contrast rule (commit b2345cc); that movie now gives 54 traces, and the channel choice is unchanged for the other 36 movies.

**λ / k0 (grid-score method, rank 1; λ·k0 in brackets; n = fitted amplitudes).**

| dataset | sign | B5 (stored data) | R1 red channel | R2 offset | R3 aperture |
|---|---|---|---|---|---|
| 4PP7 | Neg | 1/196 [196] (140) | 1/148 [148] (160) | 1/148 [148] (158) | 2/30 [60] (192) |
| | Pos | 2/126 [252] (145) | 1/101 [101] (149) | 1/103 [103] (150) | 1/48 [48] (200) |
| 4PP7-FUS | Neg | 1/83 [83] (188) | 2/49 [98] (189) | 1/80 [80] (188) | 2/25 [50] (259) |
| | Pos | 1/77 [77] (183) | 3/38 [114] (184) | 2/53 [106] (184) | 1/33 [33] (234) |
| 30PP7 | Neg | 2/91 [182] (159) | 1/140 [140] (203) | 1/141 [141] (203) | 1/44 [44] (275) |
| | Pos | 1/238 [238] (156) | 1/155 [155] (192) | 1/156 [156] (191) | 1/38 [38] (248) |
| 30PP7-FUS | Neg | 2/96 [192] (112) | 1/147 [147] (145) | 1/146 [146] (140) | 1/33 [33] (206) |
| | Pos | 2/108 [216] (112) | 1/134 [134] (141) | 1/130 [130] (146) | 2/19 [38] (214) |

**Interpretation.**
- **R1 is the real correction.** In 4PP7, 30PP7 and 30PP7-FUS, λ·k0 for negative bursts falls by 25% or more (196 → 148, 182 → 140, 192 → 147), and the number of fitted bursts rises by 15-30%. The B5 values for these three datasets came largely from DIC pages and should not be used.
- **R2 (offset) barely moves anything in vivo** (λ·k0 within 2%), as predicted. The 4PP7-FUS λ flip (2 → 1) is the grid score's instability, not an effect of the offset.
- **R3 (aperture) rescales the amplitudes about 3-4×** (λ·k0 140 → 44 for 30PP7 negative bursts), and the burst counts rise by about 35%.
  - The aperture mean includes dimmer edge pixels, so S/B − 1 is smaller.
  - Absolute k0 values from R3 are therefore not comparable with the stored data or with Qβ/in-vitro numbers extracted the old way. Only within-extraction comparisons are valid.
- **Manuscript PP7 claims, negative bursts, R3:**
  - "FUS shifts λ up and k0 down for 4x": λ 2 vs 2 (no shift), k0 30 vs 25, λ·k0 60 vs 50.
  - "No FUS effect for 30x": 1/44 vs 1/33, λ·k0 44 vs 33.

  These are grid-score point estimates without CIs. The grid method returns λ = 1 almost regardless of the true λ on simulations (section 5), so its λ comparisons carry no information. The MLE with granule-bootstrap CIs on the R3 data is the next step.
- **PP7 vs Qβ (λ 4 vs 1 in the manuscript) can no longer be compared** until Qβ is re-extracted the same way.
## 7. Likelihood fit on the real data (analysis 2, real data)
`checks/mle_real_data.m`, `checks/plot_mle_real.py`.
- Outputs (`analysis_outputs/mle_real/`):
  - `mle_real.csv`;
  - `mle_real_estimates.png/pdf` (estimates and CIs);
  - `mle_real_surfaces.png/pdf` (λ-k0 surfaces);
  - `boot.mat`, `surfaces.mat`.
- PP7 uses the R3 re-extraction; in vitro and Qβ use the stored data.
- Settings: the per-dataset paper SP and truncation, 13-frame window, 10% top trim.
- 200 granule-bootstrap resamples per dataset and sign, seeded. A rerun reproduces every number exactly.

Negative bursts. Bootstrap CIs are the ones to report (section 5: profile CIs undercover; λ is ~25% high on simulations).

| dataset | n fitted | σ | MLE λ [boot 95%] | MLE k0 [boot 95%] | λ·k0 [boot 95%] | grid λ/k0 |
|---|---|---|---|---|---|---|
| in vitro 1:1 | 72 | 15 | 1.53 [1.07–2.21] | 29 [20–40] | 45 [23–72] | 1/28 |
| in vitro 10:1 | 76 | 26 | 1.44 [0.78–2.50] | 84 [51–107] | 120 [42–236] | 2/42 |
| in vitro 100:1 | 86 | 18 | 3.47 [1.83–6.28] | 64 [35–106] | 223 [107–460] | 2/59 |
| 4PP7 (R3) | 163 | 30 | 0.84 [0.18–4.97] | 46 [11–71] | 39 [10–64] | 2/30 |
| 4PP7-FUS (R3) | 220 | 23 | 0.96 [0.68–1.45] | 37 [26–49] | 35 [28–52] | 2/25 |
| 30PP7 (R3) | 248 | 34 | 1.21 [1.01–1.47] | 38 [26–50] | 46 [33–59] | 1/44 |
| 30PP7-FUS (R3) | 185 | 25 | 1.38 [0.63–2.03] | 22 [16–36] | 30 [18–39] | 1/33 |
| 5Qb | 103 | 262 | 0.93 [0.02–4.05] | 137 [44–294] | 128 [3–240] | 2/155 |
| 10Qb | 141 | 321 | 0.02 [0.02–7.62] | 369 [50–486] | 7 [6–355] | 3/153 |

**Interpretation**
- **Identifiability.**
  - In vitro: the λ-k0 surfaces are closed and compact (`mle_real_surfaces.png`).
  - PP7: the surfaces are elongated along λ·k0 = const, closed for 4PP7-FUS and 30PP7, nearly open for 4PP7 and 30PP7-FUS.
  - Qβ: the surfaces are open bands, so λ is not identified. The quiescent noise σ (262–321) is about 10× that of every other dataset and comparable to k0.
- **Manuscript claims, tested with bootstrap CIs:**
  - **"λ = 1 for all in vitro cases": fails.** 100:1 gives λ = 3.47 [1.83–6.28], which excludes 1, and 1:1 gives 1.53 [1.07–2.21]. Even allowing for the ~25% upward bias found on simulations, 100:1 is not 1.
  - **"k0 increases with protein:RNA ratio" (39 → 63 → 139): not supported as stated.** k0 is 29 → 84 → 64, and the 10:1 and 100:1 CIs overlap. 1:1 is lower than both.
  - **"Qβ λ = 4 vs PP7 λ = 1": not testable.** Qβ λ is not identified (CIs span the whole search range), and Qβ was extracted the old way.
  - **"FUS shifts λ up and k0 down for 4x": not supported.** λ is 0.84 [0.18–4.97] vs 0.96 [0.68–1.45], and k0 is 46 [11–71] vs 37 [26–49]; both overlap fully.
  - **"No FUS effect for 30x": consistent for λ** (1.21 [1.01–1.47] vs 1.38 [0.63–2.03]). k0 is 38 [26–50] vs 22 [16–36]: the CIs overlap only marginally, so a k0 decrease with FUS cannot be excluded.
- **Grid vs MLE.** The grid λ (1-3) is often inside the MLE CIs, but on simulations the grid does not track λ (section 5). Its agreement here is not evidence.
- **Reviewer points.**
  - R2 4a (error bars): now provided.
  - R2 4b: λ is continuous, and the pipeline bias is quantified.
  - R2 5: the granule bootstrap is in place; FOV/replicate nesting is not yet done.
  - These results undermine several of the manuscript's difference claims.
## 8. Goodness of fit: amplitude histograms with the fitted model
`checks/export_amplitudes.m` (amplitudes via the same code path) and `checks/plot_mle_histograms.py`.
- Outputs (`analysis_outputs/mle_real/`): `mle_hist_Neg.png/pdf`, `mle_hist_Pos.png/pdf`, `mle_hist_gof.csv`, `amplitudes.csv`.
- Each panel shows:
  - the histogram of the fitted amplitudes (blue), with the amplitudes outside the fit window in grey;
  - the fitted truncated mixture (orange) and its n = 1-4 components (dashed).
- The fit window was rebuilt exactly; n matches `mle_real.csv` in all 18 fits.
- KS test against the window-truncated model CDF. The parameters were estimated from the same data, so p is optimistic.

KS p, negative / positive bursts:

| | 1:1 | 10:1 | 100:1 | 4PP7 | 4PP7-FUS | 30PP7 | 30PP7-FUS | 5Qb | 10Qb |
|---|---|---|---|---|---|---|---|---|---|
| Neg | 0.001 | <0.001 | <0.001 | 0.63 | 0.36 | 0.011 | 0.77 | 0.98 | 0.80 |
| Pos | 0.014 | <0.001 | <0.001 | 0.90 | 0.11 | 0.006 | 0.63 | 0.49 | 0.97 |

**Interpretation**
- **PP7 (R3 re-extraction).** The model describes the data: a unimodal, right-skewed amplitude distribution, which the fit explains with a few overlapping n·k0 components. 30PP7 is borderline.
  - No dataset shows resolved peaks at multiples of k0. λ and k0 are therefore set by the skew of one broad hump, which is why the λ-k0 surfaces are ridges (section 7).
- **Qβ.** The fits pass only because σ (262-321) is as wide as the whole distribution, so the model is one broad Gaussian and λ carries no information.
- **In vitro: the model is rejected.**
  - σ from `ampVar` is small (15-26), so the model predicts a comb of narrow peaks at n·k0.
  - The data are smooth and heavy-tailed: after the 10% top trim, 10:1 and 100:1 still keep amplitudes up to 870 and 1480, where the model density is ~0.
  - **This qualifies section 7:** the in-vitro MLE values, including 100:1 λ = 3.47, come from a misspecified model and should not be quoted as estimates.
  - The manuscript's in-vitro λ = 1 is not supported either, because the grid method does not track λ (section 5). The data do not identify λ in vitro under this model.
- **What this points to** (a model change; not applied):
  - amplitude noise that grows with n (σ² = σ0² + n·σ1²), or σ fitted from the data;
  - a heavy-tail or outlier component.
- **Reviewer point.** R2 4a/4b: "the model is a poor fit for the in vitro data" is a fair criticism, now shown directly.

## 9. Revision of sections 6-8 (2026-10-19): aperture reverted, fit visuals, 4x FUS test
**This section supersedes the PP7 parts of sections 6-8.** Those sections used R3 (5 px aperture) data.

### 9a. Aperture spot measure reverted (user decision)
- The R3 commit (e450505) was reverted. The spot is again the mean of the top-2 multithresh classes, the published method.
- Kept, because they are bug fixes, not method changes:
  - R1: red channel only; Mosaic frame/xy +1 mapping;
  - R2: camera offset 200.
- All PP7 analyses now use `in_vivo_reextracted/after_R2_offset/`. `reextract_pp7.m` no longer makes an R3 stage.
- The R3 outputs (`analysis_outputs/after_R3_aperture/`, `in_vivo_reextracted/after_R3_aperture/`) are left on disk for the record and dropped from the λ/k0 table.

Likelihood fit, R2 PP7 data, negative bursts (`analysis_outputs/mle_real/`, rerun). In vitro and Qβ are unchanged from section 7.

| dataset | n | σ | MLE λ [boot 95%] | MLE k0 [boot 95%] | grid λ/k0 | KS p (model) |
|---|---|---|---|---|---|---|
| 4PP7 | 134 | 109 | 0.02 [0.02-20.2] (λ at bound) | 159 [9-198] | 1/148 | 0.75 |
| 4PP7-FUS | 160 | 71 | 0.72 [0.02-2.16] | 64 [33-106] | 1/80 | 0.91 |
| 30PP7 | 183 | 103 | 1.06 [0.55-1.50] | 121 [91-160] | 1/141 | 0.19 |
| 30PP7-FUS | 126 | 110 | 0.90 [0.05-1.64] | 127 [69-189] | 1/146 | 0.33 |

With multithresh, σ (71-110) is comparable to k0, so λ is weakly identified for 4x. The model describes all four PP7 distributions (KS p 0.19-0.91).

### 9b. Why the fits looked worse than the manuscript figures
`checks/plot_ms_style.py` → `analysis_outputs/mle_real/ms_style_Neg.png`. Same data, drawn as `modifiedPoissonFit.m` draws MS Figs. 6-9.

Three reasons, in order of importance:
1. **Rendering.**
   - The MS curve is the Poisson pmf at n·k0, linearly interpolated, evaluated only at the 10 bin centres, and **rescaled so its peak equals the tallest bin** (`yq*max(binCounts)./max(yq)`).
   - It compares shape only, at 10 points, with a free vertical scale, and it has no measurement noise. A broad hump is matched by almost any (λ, k0) on the λ·k0 ridge.
   - The section 8 plots showed the likelihood model itself: a fixed-normalisation density with Gaussian noise σ, on 20 bins. That is a stricter view.
2. **The in-vitro comb.** In vitro, σ from `ampVar` (15-26) is small relative to k0, so the noise model predicts narrow peaks at n·k0 that the data do not show (KS p ≤ 0.014). The MS-style curve has no σ and so cannot show this misfit. This is a real model problem in vitro, not a plotting artefact.
3. **Data.** Section 8 used the R3 aperture data, about 3-4× smaller amplitudes with a different shape. With R2 (multithresh) the PP7 histograms again resemble the MS:
   - 4PP7-FUS: 20-280 with mode ~110, as in MS Fig. 9a;
   - 30PP7: mode ~150 with a tail to ~600;
   - drawn MS-style, the grid curves sit on the points as in the MS.

### 9c. 4x PP7 vs PP7-FUS: the difference is real; marginal CIs were the wrong test
`checks/fus_test_4x.m`, `checks/plot_fus_test_4x.py` → `analysis_outputs/fus_test_4x/` (csv, png/pdf).
- Negative bursts, SP 10, window 13.
- Granule-level permutation: whole traces are shuffled between conditions.
- Joint likelihood-ratio test of a shared (λ, k0) vs separate ones, on a common amplitude window and with each condition's own σ.

| data | bursts 4PP7/FUS | median 4PP7 vs FUS (perm p) | mean (perm p) | KS p | joint LRT perm p | separate fits λ/k0 |
|---|---|---|---|---|---|---|
| R2 re-extraction | 176/209 | 172 vs 119 (<0.0005) | 179 vs 159 (0.29) | 3e-8 | **0.016** | 4PP7 0.02/152, FUS 0.87/58 |
| stored (MS) data | 155/209 | 242 vs 119 (<0.0005) | 637 vs 156 (<0.0005) | 1e-12 | <0.002 | 4PP7 1.45/125, FUS 0.84/98 |

**Interpretation.**
- **The histograms differ, and the statistics agree once the right test is used.** The FUS distribution is shifted to lower amplitudes: median 119 vs 172, 0 of 2000 permutations as extreme. The means do not differ, because FUS also has a heavier upper tail.
- **Why the earlier result said "no detectable difference":**
  1. It compared λ and k0 one at a time. Their 95% CIs overlap because each parameter alone is poorly identified: the λ-k0 ridge lets either trade against the other. The difference lies across the ridge. The joint LRT, which tests both together, detects it (perm p 0.016).
  2. It used the R3 aperture data.
  3. For 4PP7, λ is at its lower bound: the amplitudes form one hump at ~160 with σ ≈ 109, which the model reads as single units of k0 ≈ 150.
- **Direction.**
  - Under the model, FUS has a smaller unit (k0 58 vs 152) and more units per burst (λ 0.87 vs ~0). This is the direction the manuscript claims (λ up, k0 down).
  - But the individual λ values are not well determined, so it should be stated as a shift of the amplitude distribution toward smaller, more numerous units, with the joint test as the evidence. Not as "λ = 1 → 2".
- **The manuscript's difference was exaggerated by the stored 4PP7 data.**
  - The stored 4PP7 traces were mostly DIC-contaminated (section 4), so the mean amplitude was 637 vs 179 after re-extraction, with a tail to 1400 (MS Fig. 7a).
  - A clear but smaller difference survives re-extraction, carried by the bulk of the distribution (median 172 vs 119), not the tail.
- **Reviewer points.**
  - R2 4a: error bars and tests are now provided.
  - R1/R2 on the FUS claim: supported as a distribution shift; not supported as a precise λ change.
