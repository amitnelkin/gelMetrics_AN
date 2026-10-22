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

## 10. MS Fig. 16 recreated: 4x and 30x, PP7 vs PP7-FUS across SP (no fitting)
`checks/export_amplitudes_sp.m`, `checks/plot_fig16_sp.py`.
- Outputs (`analysis_outputs/fig16_sp/`):
  - `fig16_R2_reextracted_Neg.png/pdf` (main);
  - `fig16_stored_MS_Neg.png/pdf` (same plot on the stored MS data, for reference);
  - `medians_*.csv`, `amplitudes_sp.csv`.
- Negative bursts at SP 10-15, R2 re-extraction, same code path as the runners.
- Each panel uses one common x-axis. The MS figure drew PP7 and PP7-FUS on separate top/bottom axes with different ranges (red to ~1600, blue to ~300-500), so its visual peak shift mixed two scales.
- Medians come with granule-bootstrap 95% CIs (traces resampled, 1000×).

| SP | 4PP7 median [95%] (n) | 4PP7-FUS median [95%] (n) | 30PP7 median (n) | 30PP7-FUS median (n) |
|---|---|---|---|---|
| 10 | 172 [161-187] (176) | 119 [109-135] (209) | 178 (351) | 192 (235) |
| 11 | 175 [162-189] (151) | 124 [111-145] (173) | 198 (278) | 201 (186) |
| 12 | 187 [167-198] (123) | 132 [115-155] (144) | 210 (225) | 209 (156) |
| 13 | 181 [162-201] (85) | 134 [116-158] (103) | 235 (155) | 211 (113) |
| 14 | 188 [167-242] (43) | 151 [124-209] (69) | 255 (96) | 226 (76) |
| 15 | 234 [174-278] (23) | 156 [131-268] (45) | 282 (54) | 238 (53) |

**Interpretation.**
- **The Fig. 16 claim holds on the re-extracted data, on a common axis.**
  - 4x: FUS shifts the negative-burst amplitude distribution down at every SP. The median is 25-35% lower, and the 95% CIs do not overlap for SP 10-13. At SP 14-15 there are too few bursts (23-69) for the CIs to separate.
  - 30x: no shift at any SP (medians within ~10%, CIs overlapping).
- The FUS distribution also has a tail of large amplitudes that 4PP7 lacks, so the means are similar while the bulk is shifted.
- Higher SP keeps fewer, longer bursts, which raises all medians. The between-condition comparison is stable.
- This answers the threshold question for this claim (R2 4a): the effect does not depend on the choice of SP.

## 11. Biased estimator, WS weights, and the data left out of the fit
`analysis/gelmetricsMLE.m` (`opts.lamPrior`), `analysis/gm_ws_weights.m`, `checks/estimator_compare_sim.m`, `checks/estimator_compare_real.m`, `checks/plot_estimators.py`.
- Outputs (`analysis_outputs/estimators/`): `sim_estimates.csv`, `sim_summary.csv`, `sim_lambda.png`, `sim_metrics.png`, `real_*.csv`, `real_summary.md`.
- PP7 uses the R2 data.

**Estimators compared** (all on the same amplitudes, same code path for simulated and real data):
- **MLE:** `gelmetricsMLE`, no prior.
- **Biased estimators (MAP, R1 point 2):** the same likelihood plus a lognormal prior on λ with median 1, either weak (log-SD 1) or moderate (log-SD 0.5).
- **WS weight sets (grid method):** weights on (MSE, QQ, KS):

  | name | weights |
  |---|---|
  | paper | 0.5 / 0.25 / 0.25 |
  | equal | 1/3 each |
  | MSE only, QQ only, KS only | 1 on one component |
  | QQ-heavy | 0.25 / 0.5 / 0.25 |
  | KS-heavy | 0.25 / 0.25 / 0.5 |

**Simulation.**
- Setup: λ_true ∈ {0.5, 1, 1.5, 2, 3, 4, 5} × noise SD {20, 40, 80} × 2 seeds × 2 signs; k0 = 60; 150 tracks.
- Metrics: RMS log error of λ, and the fraction of fits with λ within ±25% of the truth.

| estimator | noise 20 | noise 40 | noise 80 |
|---|---|---|---|
| MLE | 0.81 / 21% | 0.34 / 64% | 0.46 / 43% |
| MAP weak (s=1) | 0.76 / 21% | 0.33 / 68% | **0.35** / 43% |
| MAP moderate (s=0.5) | 0.71 / 21% | **0.31 / 75%** | 0.42 / 29% |
| WS paper | 0.80 / 29% | 1.02 / 18% | 0.84 / 29% |
| WS equal | 0.91 / 29% | 1.02 / 18% | 0.84 / 29% |
| WS MSE only | 0.86 / 25% | 1.06 / 18% | 0.93 / 21% |
| WS QQ only | 0.81 / **50%** | 0.65 / 25% | 0.73 / 32% |
| WS KS only | 0.79 / 46% | 0.72 / 25% | 0.85 / 36% |
| WS QQ-heavy | 0.88 / 32% | 0.97 / 21% | 0.81 / 25% |
| WS KS-heavy | 0.86 / 36% | 0.98 / 25% | 0.79 / 29% |

**Real data, negative bursts** (`real_summary.md` has the full tables).
- **MAP.**
  - Where the data identify λ, the prior changes nothing: in vitro and 30PP7, λ and k0 within 2%.
  - Where they don't (Qβ, 4PP7, 4PP7-FUS), λ is pulled to ~0.7-1 and the CIs narrow. That narrowing is the prior, not the data.
  - 4x with the weak prior: 4PP7 λ 0.70 [0.50-1.76], k0 122 [75-188]; FUS λ 0.92 [0.47-1.19], k0 57 [49-85]. λ is similar; k0 is lower with FUS, and the CIs barely touch.
- **WS weights change the answer on real data.**
  - 4PP7: 1/148 (paper) vs 5/37 (QQ or KS only).
  - 10:1: 2/42 (paper) vs 1/116 (QQ only).
  - 1:1: 1/28 (paper) vs 2/21 (KS only).
  - 4PP7-FUS: 1/80 (paper) vs 3/37 (KS only).
- **Data actually used** (per dataset, negative bursts):
  - Traces shorter than 180 frames are dropped by `analysis.m`. Kept: 4PP7 76/254 (30%), 4PP7-FUS 94/379 (25%), 30PP7 145/272, 30PP7-FUS 111/194, 5Qb 189/292, 10Qb 247/431, in vitro 31-45 of 57-90.
  - Of the detected bursts, the top-10% trim and the left truncation (10 or 15% per dataset) together leave 76-82% for the fit.
  - Also outside the fit by construction: segments shorter than the SP run-length threshold, unclassified segments, and the parts of each trace before the first and after the last significant segment. Positive bursts are not used in the manuscript.
- **Effect of the trims.**
  - Removing the top-10% trim raises λ and k0 in vitro: MLE 1:1 1.53 → 2.53; 100:1 3.47 → 5.17; 10:1 1.44 → 2.6. The in-vitro tail carries much of the information.
  - For PP7 (R2) the MLE moves little, except where λ is not identified (4PP7, 4PP7-FUS, Qβ), where it jumps.
  - The grid method moves more: 4PP7 1 → 5 without the trim.

**Interpretation.**
1. **Biased estimator (R1 2).** A weak prior on λ helps where λ is weakly identified. It lowers the λ error at noise 40 and 80 and stabilises the real-data fits for Qβ and 4x PP7, without moving well-identified estimates. It is a reasonable reporting choice if the prior is stated and the unpenalised MLE is shown alongside.
   - It does not fix the low-noise problem. At noise 20 every likelihood estimator overestimates λ about 2× (medians 1.6-7.5 against 0.5-5): the amplitude noise σ from `ampVar` is too small for the pipeline's amplitude errors, so the model splits each hump into sub-units (k0 too small). This is the in-vitro situation (σ 15-26), where the model was rejected (section 8).
2. **Weights (R2 4a).** The paper weights are no better than equal weights and are among the worst in simulation (λ within ±25% in 18-29% of fits). QQ-only or KS-only does somewhat better on λ (25-50%) but changes several real-data answers.
   - The WS λ is not robust to the weighting, which supports the reviewer's "arbitrary" criticism. No weight set makes the grid method competitive with the likelihood at noise 40.
3. **Excluded data (R2 4a).** About 20-25% of detected bursts are excluded by the two trims, and up to 70-75% of traces by the 180-frame rule (PP7 4x). The in-vitro λ depends strongly on the top trim. These cuts must be stated, and their effect shown (this table), in the revision.

## 12. Qβ re-extracted from the Naor .nd2 movies: same behaviour as the stored data
`checks/qb_compare.m`, `checks/plot_qb_compare.py` → `analysis_outputs/qb_compare/` (traces/bursts/fits/amplitudes csv, `qb_compare.png/pdf`).
- **Movies:**
  - 5Qb: `Raw_Results\Naor_Data-All_night_induction\5Qb\250ms`, 10 movies;
  - 10Qb: `...\10Qb`, 9 movies (it has no 250 ms subfolder).
- **Pipeline:** same as PP7. bfconvert → red-channel tracking (R1 macro) → R1/R2 extraction (multithresh spot, offset 200) → `in_vivo_reextracted/after_R2_offset/5Qb_*.mat` and `10Qb_*.mat`.
- **Movie layout:** all 19 movies are 2-channel. The second channel is dim and nearly featureless (median ~1050, max ~1300); the red granules are bright (99th percentile 9,000-35,000).
- **`combined` subfolders.** The movie folders also contain later (Nov 2025) outputs: `combined/spot.mat`, `background.mat`, `*_xy.mat`. These are not the stored data: no trace matches, and the levels differ (median 1065 vs 13431).

| | 5Qb stored | 5Qb re-extracted | 10Qb stored | 10Qb re-extracted |
|---|---|---|---|---|
| traces / ≥180 frames | 292 / 189 | 296 / 191 | 431 / 247 | 517 / 279 |
| median spot / background level | 13431 / 7963 | 13207 / 7777 | 17765 / 10753 | 18264 / 11169 |
| zigzag-flagged traces | 0 | 0 | 0 | 0 |
| neg / pos bursts at SP 14 | 134 / 158 | 145 / 179 | 184 / 172 | 217 / 210 |
| median neg. amplitude, SP 14 | 397 | 426 | 542 | 574 |
| σ (ampVar) | 262 | 306 | 321 | 373 |
| grid λ/k0, negative | 2/155 | 1/298 | 3/153 | 2/250 |
| MLE λ [boot 95%], negative | 0.93 [0.02-40] | 0.44 [0.02-50] | 0.02 [0.02-25] | 50 [0.02-50] |

**Interpretation.**
- **The movies are the source of the stored Qβ data.**
  - Trace counts, lengths, levels and burst counts agree within a few percent. The 5Qb level drop of ~200 is exactly the offset subtraction.
  - The amplitude histograms overlap at every SP.
  - Bursts per trace vs SP are identical. The re-extracted medians are 5-10% higher, as expected from tracking the red channel only.
- **Qβ was not affected by the zigzag artefact.** The stored Qβ traces are clean and at red-granule levels.
- **λ for Qβ is not identified, in either version.**
  - σ (260-370) is as large as the amplitudes, so the likelihood CIs span the whole search range.
  - The grid λ moves between 1 and 3 between versions that differ only slightly.
  - The manuscript's "Qβ λ = 4" is therefore not supported by either method on either version.
  - What the data do show robustly is that the Qβ amplitudes are 2-3× larger than PP7's (median ~400-570 vs ~120-210), with a broad, roughly symmetric hump.

## 13. Linking unit brightness by valency (joint fits)
`analysis/gelmetricsJoint.m`, `checks/valency_sim_validation.m`, `checks/valency_real.m` → `analysis_outputs/valency/` (`sim_validation.csv`, `real_fits.csv`, `real_boot.mat`).

**Idea.** Within one protein, the brightness of one unit (k0) should scale with the number of binding sites on the RNA. Two versions:
- **"full occupancy":** k0 ratio = valency ratio;
- **"at most full":** k0 ratio between 1 and the valency ratio.

Each dataset keeps its own λ.

**Groups:**

| group | datasets | valency (k0 ratio cap) | data | SP |
|---|---|---|---|---|
| Qβ | 5Qb, 10Qb | 2 | R2 | 14 |
| PP7 | 4PP7, 30PP7 | 7.5 | R2 | 12 |
| PP7-FUS | 4PP7-FUS, 30PP7-FUS | 7.5 | R2 | 12 |
| in vitro | 1:1, 10:1, 100:1 (PCP-14x RNA) | at most 1 : 10 : 14 proteins per RNA | stored | 12 |

Truncation 15%, top trim 10%, negative bursts, 50 granule-bootstrap resamples per joint fit.

**Simulation check (known truth): the link does not help yet.**
- Share of fits with λ within ±25% of the truth:

  | noise SD | separate | full occupancy | at most full |
  |---|---|---|---|
  | 40 | 21% | 7% | 11% |
  | 80 | 43% | 50% | 54% |

- The bright member of each pair is already fitted badly on its own: λ is overestimated, e.g. 4.6 for a true 1. This is the known problem that the noise width σ from the quiescent segments is too small (sections 5, 11). Linking passes that error to the partner.
- In the partial-occupancy case (true ratio 1.5, assumed 2), "full occupancy" overestimates λ of the low-valency dataset 2-3×.

**Real data, highlights** (95% bootstrap intervals in brackets; full table in `real_fits.csv`):
- **Qβ:** both versions give 5Qb λ ~0.5 and 10Qb λ ~0 (single units), k0 ~220 vs ~430. "At most full" chooses a ratio of 1.9, close to full occupancy. The intervals still span almost the whole range (5Qb λ 0.02-19).
- **PP7:**
  - Forcing 30PP7 to be up to 7.5× brighter makes 4PP7 consist of ~12 tiny units (k0 17). That is not physically sensible: the data want a ratio of 7.5 or more (it sits at the cap).
  - Read plainly: the 4x and 30x PP7 amplitudes are too similar for the 30x RNA to be much more occupied than the 4x RNA. This agrees with the MS remark that 30x is not fully occupied.
- **PP7-FUS:** "at most full" gives a ratio of 1.9 and λ ~0.8-1.0 for both. "Full occupancy" fits slightly worse and again needs many tiny units for 4x.
- **In vitro:**
  - "Full occupancy" (1:10:14) fits clearly worse than "at most full" (log-likelihood 1507 vs 1436), so full occupancy is rejected.
  - "At most full" gives a brightness ratio of 2.3 for both 10:1 and 100:1 vs 1:1, with λ 1.5 [1.1-2.4], 2.0 [0.5-4.8], 3.5 [1.2-7.0].

**Interpretation.**
- Valency is the right kind of extra information. On the real data it already says something useful: brightness does not scale with valency (30x vs 4x PP7; 10:1/100:1 vs 1:1 in vitro). The high-valency RNAs are far from full.
- It does not yet give reliable λ values. On simulations it does not improve λ recovery, because the per-burst noise model underestimates the noise for bright bursts. That noise model must be fixed first (σ calibrated from simulations of the full pipeline, or σ growing with burst size); then this joint fit should be re-run and re-validated.

## 14. Fixing the noise model (all approaches tried)
`analysis/gm_mle_nll.m`, `analysis/gelmetricsMLE.m` (`opts.noise`, `opts.falseBursts`, `opts.floorK`, `opts.h0K`), `analysis/gelmetricsJoint.m` (`opts.mixFree`, `opts.h0K`), `checks/noise_calibration.m`, `checks/noise_calibration_h0.m`, `checks/noise_model_validation.m`, `checks/noise_model_real.m`, `checks/plot_noise_models.py` → `analysis_outputs/noise_models/`.

**What was wrong (found on simulations with known answers).** The main problem was not the noise width. The detector reports extra small "bursts" (noise wiggles and fragments of real steps, about 20% of detections), and the 15% left cut does not remove them. The fit then explains them by making units small and λ large. Real bursts also come out about 20% larger than the true step; that only rescales k0.

**Approaches tried, all judged on held-out simulations** (7 λ × 3 noise levels × 2 unit sizes × 2 signs; calibration used different seeds):

| model | what it does | λ within ±25% | large-error score (lower = better) |
|---|---|---|---|
| fixed (old) | σ from quiet segments, 15% left cut | 31% | 0.93 |
| free | one burst scatter, fitted | 62% | 0.86 |
| growN / prop | scatter grows with burst size | 49% / 58% | 1.86 / 1.52 |
| calib+floor | scatter calibrated on simulations (1.33σ, +7.5% of size) + noise floor | 58% | 1.10 |
| mix / mix+free / mix+prop | + false-burst component, width fitted | 38% / 65% / 69% | 1.09 / 0.68 / 1.03 |
| floor / floor+free | left cut at a noise floor of 3.1σ (calibrated) | 37% / 64% | 1.07 / 1.05 |
| **mixH+free** | false-burst component with width fixed at 1.9σ (calibrated) + fitted scatter | **63%** | **0.78** |

- Calibration from burst-free simulations: false bursts stay below about 3σ (95%) at every SP, with a typical size of 1.9σ.
- A slow drift over a whole trace is detected as one huge "burst" that no amplitude cut removes.
- **Chosen model: mixH+free.** On simulations it removes the 1.6× λ overestimate (median λ̂/λ ≈ 1.05-1.09) and tracks λ well when units are clearly larger than the noise. Nothing works for small units in very noisy traces (k0 60 at noise 80): a detection limit.

**Real data with the chosen model** (negative bursts, 100 trace-resampling CIs):

| dataset | λ [95%] | k0 [95%] | false-burst share |
|---|---|---|---|
| in vitro 1:1 | 2.3 [0.02-11] | 106 [1-141] | 0.50 (cap) |
| in vitro 10:1 | 5.0 [0.4-15] | 338 [55-441] | 0.50 (cap) |
| in vitro 100:1 | 1.35 [0.8-4.4] | 180 [91-284] | 0.41 |
| 4PP7 | 0.02 [0.02-50] | 171 [3-380] | 0 |
| 4PP7-FUS | 0.6 [0.02-50] | 96 [1-488] | 0.03 |
| 30PP7 | 5.3 [0.02-50] | 207 [20-717] | 0.50 (cap) |
| 30PP7-FUS | 0.9 [0.02-50] | 179 [18-552] | 0.39 |
| 5Qb | 0.8 [0.02-50] | 346 [2-1178] | 0 |
| 10Qb | 0.02 [0.02-50] | 444 [2-1331] | 0 |

- **The λ values change with the noise model.** For example, 30PP7 gives 1.1 (old), 4.3 (floor+free), 1.2 (mix+free) and 5.3 (mixH+free); 10:1 gives 1.4-5.3.
- **The strict noise floor removes almost every PP7 and Qβ burst**: the real quiet-segment σ is large compared with the amplitudes.
- **The valency-linked joint fits are unstable in the same way**; e.g. 30PP7-FUS λ is 1.0 with full occupancy and 16.9 with "at most full". The joint fits' log-likelihoods barely differ between the two versions (1-5 units).

**Plain-English interpretation.**
- The fitting method now works on simulated data whose bursts really are made of equal units.
- On the real data, the amplitude histograms are one smooth hump with no visible steps at 1, 2, 3 units. Many combinations of "how many units" and "how bright a unit" fit equally well. So the number of molecules per burst (λ) cannot be determined from these data, by any of the methods tried, with or without valency.
- The manuscript's λ = 1, 2, 4 values are therefore not supported.
- What the data do support, directly and without a model:
  - the amplitude distributions themselves;
  - the burst size scale (median amplitude);
  - the differences between conditions: Qβ vs PP7 (sections 9-10, 12) and FUS for 4x (sections 9c, 10).
- To get λ, a measured single-unit brightness is needed: single-RNA calibration, or photobleaching steps of one RNA with its proteins. With k0 measured, λ follows from the median burst size.

## 15. Fraction of the granule shed per burst (instead of λ)
`analysis.m` (new optional 7th output `segInfo`: start, end and fitted level of each burst; other outputs unchanged), `collect_intensity.m` (new optional 4th output: integrated intensities), `checks/reextract_sums.m`, `checks/shed_fraction.m`, `checks/plot_shed_fraction.py`.
- Outputs: `analysis_outputs/shed_fraction/` (`fractions.csv`, `medians.csv`, `tests.csv`, `shed_hist.png`, `shed_vs_sp.png`); data in `in_vivo_reextracted/after_R2_sums/`.

**Idea.** For each negative burst, fraction shed = |burst size| ÷ granule intensity just before the burst. It is unitless, so no single-unit brightness (k0) is needed.

**Granule intensity measures:**

| measure | definition | datasets |
|---|---|---|
| mean | current trace: mean of the multithresh spot pixels above the cell level | all, in vitro included |
| spot sum | sum of the multithresh spot pixels above the cell mean (the proposed measure) | PP7, PP7-FUS, Qβ |
| crop sum | sum of the 15×15 crop above its border median (mask-free check) | PP7, PP7-FUS, Qβ |

- The sums were re-extracted from the movies with the same tracking and settings. The mean traces came out identical to R2.
- Burst detection is the same analysis as before. SP 10-15; granule-level permutation tests (2000 shuffles of whole traces); granule-bootstrap CIs.

**Median fraction shed per negative burst, SP 12** [95% CI]:

| dataset | mean | spot sum | crop sum |
|---|---|---|---|
| 4PP7 | 0.066 [0.060-0.072] | 0.110 [0.103-0.125] | 0.043 |
| 4PP7-FUS | 0.096 [0.088-0.110] | 0.172 [0.160-0.187] | 0.062 |
| 30PP7 | 0.092 [0.085-0.101] | 0.153 [0.138-0.168] | 0.058 |
| 30PP7-FUS | 0.114 [0.099-0.126] | 0.172 [0.163-0.181] | 0.058 |
| 5Qb | 0.052 [0.048-0.059] | 0.092 [0.088-0.098] | 0.052 |
| 10Qb | 0.054 [0.051-0.057] | 0.091 [0.087-0.097] | 0.048 |
| in vitro 1:1 / 10:1 / 100:1 | 0.079 / 0.122 / 0.152 | – | – |

**Tests (consistent across SP 10-15 unless stated):**
- **4x FUS vs PP7:** FUS sheds a larger fraction with all three measures (spot sum 17% vs 11%; p ≤ 0.016 at every SP).
- **30x FUS vs PP7:** small increase with mean and spot sum (p 0.002-0.06 at SP 10-14); none with crop sum.
- **Qβ vs PP7:**
  - Qβ sheds a *smaller* fraction per burst with mean and spot sum: 9% vs 11-15%, p < 0.001 at every SP; 5Qb vs 4PP7 (same RNA) p 0.001-0.03.
  - With the crop sum there is no difference (p ≈ 0.6).
  - Qβ granules are 3-10× brighter in total (median spot sum before a burst: 5Qb 433,000, 10Qb 523,000 vs 4PP7 141,000, 30PP7 54,000). Their bursts are larger in absolute terms (38,000-48,000 vs 9,000-17,000) but a smaller share of the granule.
- **In vitro (mean only):** the fraction shed rises with protein:RNA ratio, 8% → 12% → 15%.

**Interpretation (plain).**
- Expressed as a share of the granule, every condition sheds roughly 5-17% of its intensity per burst.
- The FUS effect on the 4x RNA is robust: FUS granules lose a larger share per burst.
- For Qβ vs PP7 the picture reverses compared with absolute burst size. Qβ bursts are bigger, but Qβ granules are much brighter, so each burst is a smaller share (about 9% vs 11-15%).
- The crop-sum measure does not reproduce that difference. A fixed 15×15 box likely cuts off part of the large, bright Qβ granules and includes their edges in the "background" border, so the spot sum (the published segmentation) is the more trustworthy of the two sums here.
- This measure needs no model and no k0, and it can be stated with error bars.

## 16. Why the crop sum disagrees for Qβ (diagnosis)
`collect_intensity.m` (`opts.sumDiag`: extra diagnostic sums; other outputs unchanged), `checks/crop_sum_diagnosis.m`, `checks/plot_crop_diagnosis.py` → `analysis_outputs/crop_diagnosis/` (`radial.csv`, `fractions.csv`, `nearcell.csv`, `summary.csv`, `tests.csv`, `crop_diagnosis.png/pdf`).

**What was checked:**
1. Granule size against the 15×15 box.
2. Radial intensity profiles out to 20 px.
3. The fraction shed recomputed with six intensity measures (SP 12):
   - spot sum;
   - a granule-only "ring sum" (r ≤ 6 px minus the 8-11 px ring);
   - box sums with half-widths 7, 10, 14 and 20 px.
4. Whether the cell light just around a granule (4-8 px ring) rises when the granule sheds.

**Findings:**
- **Qβ granules are not cut off by the box.**
  - Spot masks: Qβ ~56 px, 4PP7 50 px, 30PP7 28 px, all well inside the 225-px box.
  - Half-maximum radius: Qβ 4-5 px, PP7 3 px.
  - Qβ also has a diffuse glow around the granule (~10% of the centre level at 10-20 px), which PP7 lacks.
- **Box sums are mostly not granule light.** Box sum / spot sum ≈ 2.5-3.9 in the 15×15 box and up to ~11 in the 41×41 box. The extra is surrounding cell fluorescence (larger for Qβ). The box-sum fraction is therefore burst ÷ (granule + a variable amount of cell), and that dilution differs between constructs.
- **The shed light does not stay next to the granule.** When the spot sum drops, the 4-8 px ring brightens in only 23-44% of bursts; the median change is slightly negative, because the ring contains the granule's own edge. So the box sum is not insensitive simply because the released material stays inside the box.
- **Granule-level permutation p for each comparison and measure:**

  | measure | Qβ vs PP7 (pooled) | 5Qb vs 4PP7 (same RNA) | 4x FUS vs PP7 |
  |---|---|---|---|
  | spot sum | 0.0005 | 0.0015 | 0.0005 |
  | ring sum | 0.0005 | 0.46 | 0.0005 |
  | box 15×15 | 0.59 | 0.26 | 0.002 |
  | box 21×21 | 0.88 | 0.10 | 0.0005 |
  | box 29×29 | 0.61 | 0.26 | 0.0005 |
  | box 41×41 | 0.0015 | 0.73 | 0.0005 |

**Interpretation (plain).**
- The crop sum disagrees because a fixed box measures the granule *plus* a lot of surrounding cell light. How much extra light it picks up differs between Qβ and PP7 cells, so it blurs the comparison. Missing granule light is not the cause.
- The spot sum (the published multithresh granule mask, minus the local cell level) measures the granule itself. It is the right measure for "fraction of the granule shed".
- **Robust across every measure:** FUS makes 4x granules shed a larger fraction per burst.
- **Qβ shedding a smaller fraction than PP7** is supported by the granule-specific measures (spot sum and ring sum, pooled), but the same-RNA comparison 5Qb vs 4PP7 holds only with the spot sum. It should be stated with that caveat.

## 17. Candidate summary scores for comparing systems
`checks/score_export.m`, `checks/score_compare.py` → `analysis_outputs/scores/` (`scores.csv`, `stability.csv`, `discrimination.csv`, `sim_tracking.csv`, `sim_scores.csv`, `scores.png/pdf`).

**Scores** (negative bursts; fraction = |burst| ÷ granule level just before it):

| score | meaning |
|---|---|
| f_med, f_mean | median / mean fraction of the granule shed per burst |
| A_mean | mean absolute burst size (≈ λ·k0, in intensity units) |
| shape | gamma shape of the fraction distribution: how peaked around a typical size |
| lam_mom | moment estimate of units per burst, mean² ÷ (variance − noise²) |
| rate | bursts per granule per hour |
| turnover | rate × f_mean, % of the granule shed per hour |
| net | % gained per hour minus % shed per hour |

- In vivo: spot sum and mean measures (R2 re-extraction). In vitro: mean measure (stored data).
- SP 10-15; 300 granule-bootstrap resamples at SP 12.
- Simulations: λ 0.5-5 × k0 60/120 × noise 40/80 × 2 seeds.

**How the scores behave:**

| score | variation across SP 10-15 | rank agreement, mean vs spot sum | typical CI width | follows λ in simulations |
|---|---|---|---|---|
| f_med | 9% | 0.94 | ±10% | yes (ρ 0.98) |
| f_mean | 11% | 0.94 | ±8% | yes (ρ 0.98) |
| A_mean | 10% | 0.71 | ±12% | yes (ρ 0.98) |
| shape | 53% | −0.37 | > ±80% | weakly (ρ 0.3-0.8) |
| lam_mom | 13% (often undefined) | – | > ±100% | no (anti-correlated) |
| rate | 51% | 0.83 | ±13% | yes (ρ 0.96-0.98) |
| turnover | 43% | 1.00 | ±17% | yes (ρ 0.98) |
| net | 49% | 0.77 | > ±100% | no |

**Discrimination at SP 12**, z = difference ÷ combined bootstrap SE (|z| > 2 ≈ p < 0.05):

| comparison | f_mean | A_mean | turnover | rate |
|---|---|---|---|---|
| 4PP7-FUS vs 4PP7 (spot sum) | +6.6 | −5.1 | +2.8 | −0.1 |
| 30PP7-FUS vs 30PP7 (spot sum) | +2.1 | −4.6 | +0.2 | −1.1 |
| 5Qb vs 4PP7 (spot sum) | −3.1 | +12.2 | −0.3 | +1.3 |
| 10Qb vs 30PP7 (spot sum) | −9.0 | +23.8 | −3.5 | +2.2 |
| in vitro 10:1 vs 1:1 | +2.1 | +1.9 | +2.1 | +1.5 |
| in vitro 100:1 vs 1:1 | +2.3 | +1.6 | +3.4 | +4.4 |

**Values (SP 12; spot sum in vivo, mean in vitro):**

| dataset | f_mean [95%] | A_mean | turnover %/h [95%] |
|---|---|---|---|
| 4PP7 | 0.106 [0.097-0.116] | 15,000 | 16 [13-20] |
| 4PP7-FUS | 0.158 [0.148-0.172] | 9,000 | 24 [20-28] |
| 30PP7 | 0.145 [0.135-0.158] | 9,600 | 22 [19-26] |
| 30PP7-FUS | 0.162 [0.152-0.174] | 5,900 | 23 [20-26] |
| 5Qb | 0.090 [0.084-0.096] | 36,800 | 16 [14-18] |
| 10Qb | 0.087 [0.083-0.092] | 45,600 | 16 [14-17] |
| in vitro 1:1 | 0.089 [0.071-0.121] | 49 | 20 [15-28] |
| in vitro 10:1 | 0.178 [0.118-0.271] | 130 | 49 [31-82] |
| in vitro 100:1 | 0.167 [0.125-0.247] | 188 | 61 [45-91] |

**Interpretation (plain).**
- **Best single score: the fraction of the granule shed per burst (f_mean or f_med).**
  - It is the "λ × k0" idea made unitless: mean burst size ÷ granule intensity. So it combines both and does not need either on its own.
  - It is stable across detection settings (±10%) and agrees between the two intensity measures. It is precise (±8-10%) and follows the true λ in simulations.
  - It separates FUS vs PP7 for 4x (z 6.6) and Qβ vs PP7 (z −3 to −9), and in vitro 1:1 vs higher protein (z ≈ 2-4).
  - It mixes "more units" and "brighter units" (it rises with both λ and k0), so it is a *relative burst size*, not a molecule count.
- **Useful companions:**
  - **Absolute burst size (A_mean)** separates Qβ from PP7 strongly, in the opposite direction (Qβ bursts larger, z +12 to +24). It is in intensity units, so it is only comparable under identical imaging.
  - **Turnover (% of granule per hour)** is the most biologically readable number and agrees perfectly between measures. Its value moves with SP (the burst count depends on the threshold), so it must be quoted at one fixed SP.
- **Not useful:** the shape-based "number of units" scores (gamma shape, moment λ) and net exchange. They are unstable, imprecise and do not follow λ, which is consistent with sections 13-14: the distribution shape does not carry recoverable unit information in these data.
- **Suggested reporting:** fraction shed per burst as the main parameter, with absolute burst size and turnover at a fixed SP as secondary parameters, all with granule-bootstrap CIs. The noise floor / single-SP choice and the sensitivity tables should then be run on these scores.

## 18. Noise floor for the chosen scores: most detected bursts are indistinguishable from noise
`checks/score_export2.m` (`opts.surrogate`), `checks/robustness.py noise` → `analysis_outputs/robustness/` (`noise_floor_allSP.csv`, `noise_floor.png/pdf`, `*_bursts.csv`, `*_traces.csv`).

**Burst-free control.** Each real trace is replaced by its slow cubic trend (the bleaching curve, as `analysis.m` fits the background) plus fresh white noise at that trace's own noise level. The noise level is 1.4826·MAD of the frame-to-frame differences / √2, which is robust to occasional steps. Everything detected in these surrogates is a false burst.
- The white-noise assumption holds: the lag-1 autocorrelation of the real frame differences is −0.47 to −0.50 (white noise: −0.5).
- A first surrogate (61-frame moving-mean trend + shuffled residuals) was rejected: its trend still contains the real bursts, smeared out.

**False-discovery share** (false ÷ real negative bursts; main measure: spot sum in vivo, mean in vitro):

| dataset | SP 10 | SP 12 | SP 14 | SP 16 | SP 18 | SP 20 | real bursts at SP 12 / 18 |
|---|---|---|---|---|---|---|---|
| 4PP7 | 0.98 | 0.92 | 0.93 | 0.95 | 2.2 | 3.5 | 106 / 5 |
| 4PP7-FUS | 0.97 | 0.92 | 0.75 | 0.79 | 1.1 | 1.7 | 124 / 9 |
| 30PP7 | 1.01 | 0.97 | 0.92 | 0.68 | 0.46 | 1.0 | 207 / 13 |
| 30PP7-FUS | 1.00 | 0.99 | 0.94 | 1.15 | 1.5 | 0.67 | 149 / 10 |
| 5Qb | 0.79 | 0.79 | 0.70 | 0.59 | 0.47 | 0.39 | 313 / 38 |
| 10Qb | 0.80 | 0.78 | 0.67 | 0.69 | 0.63 | 0.64 | 472 / 35 |
| in vitro 1:1 | 0.97 | 0.88 | 0.83 | 1.04 | 1.1 | 1.07 | 94 / 18 |
| in vitro 10:1 | 0.77 | 0.71 | 0.41 | 0.39 | 0.27 | 0.25 | 93 / 33 |
| in vitro 100:1 | 0.91 | 0.72 | 0.79 | 0.76 | 0.74 | 0.81 | 105 / 38 |

- **No SP reaches ≤ 10% false bursts while keeping a usable number of bursts.**
- A size cut does not rescue PP7. Among bursts ≥ 20% of the granule, 33-83% are false in PP7; for Qβ it is 11-15% (9-13 bursts) and for in vitro 10:1 ≥ 30% it is 5% (21 bursts).
- False bursts are similar in size to real ones (median fraction shed, SP 14: e.g. 4PP7 0.13 real vs 0.12 false; 5Qb 0.10 vs 0.08). Their size follows each system's noise-to-brightness ratio, which is lower for the bright Qβ granules.

**Interpretation (plain).**
- At the SP values used in the manuscript (10-14), roughly 70-100% of the detected bursts in every condition could have been produced by measurement noise plus bleaching alone. Only in vitro 10:1 and, partly, Qβ at high SP show a clear excess of real bursts.
- The cause is the known detection-threshold problem: the run-length threshold assumes independent frame-to-frame changes, but it is applied after a 13-frame moving average, which turns plain noise into long monotone runs.
- The burst-based scores (fraction shed, burst size, turnover) and their differences between conditions therefore partly reflect each system's noise level relative to granule brightness, not only granule dynamics.
- **This must be resolved before the sensitivity, nested-bootstrap and window analyses are meaningful.** Those were deliberately not run.
- **Options:**
  1. Correct every score for the expected false bursts from matched surrogates ("excess over noise").
  2. Replace burst detection with a detection-free fluctuation measure: how much more the granule intensity fluctuates than noise alone, at each timescale.
  3. Re-calibrate the detector so its false-discovery share is controlled (e.g. a threshold set from surrogates per dataset), and accept far fewer bursts.

## 19. Three noise-robust approaches compared (simulation bench + real data)
- Code: `checks/noise_methods_sim.m` (bench), `checks/noise_methods_real.m` (fluctuations), `checks/score_export2.m` (surrogates `poly`, `polysd`, `bleach`), `checks/noise_methods.py`.
- Outputs: `analysis_outputs/noise_methods/` (`sim_*.csv`, `sim_eval.csv`, `real_scores.csv`, `real_comparisons.csv`, `real_stability.csv`, `sim_methods.png`, `real_methods.png`) and `analysis_outputs/robustness/noise_floor_sharedbleach.csv`.

**Bench.** Known dynamics with the noise varied independently:
- relative burst size 3, 6 or 12% of the granule, at 0, 3 or 10 bursts per hour (half losses, half gains);
- noise 0.5, 1 or 2% of the granule; bleaching on; 100 traces × 361 frames; 2 seeds.

Exactly the real-data code path: detection with `analysis.m`, the same surrogates, the same fluctuation measure.

**Methods:**
- M0: uncorrected turnover (reference; section 17).
- M1: excess turnover, real minus its own burst-free surrogate (each trace paired with its own surrogate).
- M2: excess fluctuation, F(5 min) = how much the granule's relative intensity changes over 5 min beyond the surrogate's.
- M2b: noise-free fluctuation, G(5 min) = the same, but with the 1-frame lag subtracted instead of a surrogate. White noise adds equally at every lag, so this removes it exactly, whatever its level, with no noise estimate and no burst detection.
- M3: recalibrated detector, keeping only bursts above a per-dataset size cut that leaves ≤ 10% false bursts.

**Bench results (SP 12):**

| method | score with no bursts | follows the truth (Spearman) | high-noise / low-noise score, same dynamics | depends on SP |
|---|---|---|---|---|
| M0 uncorrected | 0.6-2.3 %/h | 0.98 | 1.01 | yes |
| M1 excess turnover | ≈ 0 | 0.95, but ~50% of the true value | 1.29 | yes |
| M2 excess fluctuation | ≈ 0 | 0.98 | 0.96 | no |
| **M2b noise-free fluctuation** | **≈ 0 (max 0.15%)** | **0.98** | **0.97** | **no** |
| M3 recalibrated detector | – | no valid cut in 78% of settings | – | yes |

- M1 under-corrects because each trace's own cubic trend bends with real gains and losses, so its surrogate inherits part of the dynamics.
- A surrogate built on one shared bleaching shape per dataset avoids this. On the real data it changes the false-discovery share only slightly (SP 12: PP7 0.82-1.03, Qβ 0.59-0.68, in vitro 0.66-0.78; previously 0.92-0.99, 0.78-0.79 and 0.71-0.88), so section 18 stands.

**Real data, SP 12 (spot sum in vivo, mean in vitro), 95% granule-bootstrap CI:**

| dataset | M0 turnover %/h | M1 excess turnover %/h | M2b G(5 min) % |
|---|---|---|---|
| 4PP7 | 17.8 [14.6-21.4] | 1.0 [−3.1-5.3] | 2.4 [0-4.0] |
| 4PP7-FUS | 25.9 [21.7-30.0] | 4.2 [−2.5-10.6] | 5.0 [2.3-6.6] |
| 30PP7 | 24.4 [21.6-27.8] | 1.4 [−2.5-5.5] | 3.8 [0-5.3] |
| 30PP7-FUS | 24.8 [21.3-28.7] | 1.6 [−3.5-7.1] | 2.8 [0-5.6] |
| 5Qb | 17.5 [15.4-19.5] | 5.8 [3.2-8.0] | 4.0 [3.2-4.7] |
| 10Qb | 17.6 [15.8-19.4] | 6.4 [4.3-8.6] | 4.1 [3.5-4.6] |
| in vitro 1:1 | 29.6 [21.5-40.7] | 5.5 [−2.0-14.0] | 6.9 [2.9-10.0] |
| in vitro 10:1 | 71.6 [49.9-98.7] | 50.4 [28.9-77.1] | 13.1 [9.5-16.4] |
| in vitro 100:1 | 95.1 [62.2-136] | 31.9 [2.6-57.8] | 11.1 [6.8-15.0] |

- M2b on the surrogates themselves is ≈ 0 (≤ 1%) everywhere.
- M3 finds a valid cut only for 10Qb (6 bursts) and in vitro 10:1 (24 bursts).

**Comparisons, z (|z| > 2 ≈ p < 0.05):**

| comparison | M0 | M1 | M2b |
|---|---|---|---|
| 4PP7-FUS vs 4PP7 | 2.9 | 0.8 | 1.6 |
| 30PP7-FUS vs 30PP7 | 0.1 | 0.0 | −0.5 |
| 5Qb vs 4PP7 | −0.1 | 1.9 | 1.2 |
| 10Qb vs 30PP7 | −3.6 | 2.0 | 0.3 |
| in vitro 10:1 vs 1:1 | 3.1 | 3.4 | 2.4 |
| in vitro 100:1 vs 1:1 | 3.3 | 1.9 | 1.5 |

**Interpretation (plain).**
- **Cleanest method: M2b, the noise-free fluctuation measure.** On the bench it reads zero without dynamics, follows the true dynamics and does not move with the noise level. It needs no burst detection, no SP and no noise estimate. It answers "how much does the granule change beyond noise over 5 minutes?" (it scales with burst rate × size²).
- **What it says about the real data:**
  - In vivo, all six conditions show only small slow changes beyond noise (about 2-5% over 5 min), with overlapping CIs. With noise removed, FUS vs non-FUS and Qβ vs PP7 are not distinguishable (|z| ≤ 1.6).
  - In vitro granules change clearly more (7-13%), and 10:1 is distinguishable from 1:1 (z 2.4).
  - Noise-corrected burst counting (M1) agrees qualitatively: excess turnover is indistinguishable from 0 for every PP7 dataset. It is positive for Qβ (6 %/h) and large in vitro (32-50 %/h).
- **The differences between in-vivo conditions seen with the uncorrected burst scores (sections 15-17) are therefore not supported once noise is accounted for.** They mainly reflect how noisy each system's traces are relative to granule brightness.
- **What remains robust:**
  1. in vitro, higher protein:RNA ratios give clearly more dynamic granules;
  2. Qβ granules show a small but significant excess of real bursts over noise (M1, CI excludes 0), while PP7 granules do not;
  3. in vivo granules are, beyond noise, largely static on the 1-10 min timescale at this imaging noise level.
- **Recommendation:** report M2b with granule CIs as the primary dynamics measure. Report M1 as the burst-level cross-check. Present the burst-based GelMetrics scores only together with their noise-floor control.

## 20. 4PP7 vs 4PP7-FUS: more replicates, or a different detection setting?

Scripts: `checks/power_fus.py` (outputs in `analysis_outputs/power/`) and `checks/detection_grid.m` + `checks/grid_eval.py` (outputs in `analysis_outputs/detection_grid/`; `bursts.csv`, 33 MB, is not committed and can be regenerated). Reviewer points: R2 5, R2 6, R2 2b, R2 4a.

### 20a. Would more replicates make it significant? (same effect, same spread)

The SE shrinks as 1/√(data), so z grows as √k when every condition has k times more data. Nested = movies resampled, then granules within each movie.

| score | 4PP7-FUS vs 4PP7 | z (granule) | z (nested) | × data for p < 0.05 | × data for 80% power | movies needed (80% power) |
|---|---|---|---|---|---|---|
| M2b, noise-free fluctuation (5 min) | 4.98 vs 2.36 % | 1.53 | 1.20 | 2.7 | 5.4 | ~49 / ~66 (now 9 / 12) |
| M1, excess turnover (SP 12, w 13) | 5.8 vs 2.8 %/h | 0.77 | 0.56 | 12 | 25 | not feasible |

- Movies add their own variability (nested SE 1.28× the granule SE). More **movies** are needed, not only more granules per movie.
- Other in vivo comparisons with M2b need 18-87× more data (30PP7-FUS vs 30PP7: no trend at all, z −0.4).
- The FUS z for M2b depends on the lag, from 0.7 to 2.0 (2.0 at 1 min). The 5-min lag was fixed in advance (section 19). Choosing the lag afterwards would be selecting on the outcome.

### 20b. Lower SP or a different smoothing window

Grid: SP {3 4 5 6 8 10 12} × window {5 7 9 13 17 21} frames, negative bursts, main measure, traces ≥ 180 frames. Each setting is run on the real traces and on burst-free surrogates (shared bleaching shape + white noise at each trace's own level). False-discovery share = surrogate bursts ÷ real bursts. Excess turnover = (real − surrogate) fraction shed per hour, with a nested bootstrap SE (500×).

**Rule fixed before looking at any comparison:** pick the setting with the lowest median in vivo false-discovery share among those that keep ≥ 50 real bursts per in vivo dataset. Result: **window 17, SP 12** (median 0.78, worst 1.09).

Median in vivo false-discovery share:

| SP \ window | 5 | 7 | 9 | 13 | 17 | 21 |
|---|---|---|---|---|---|---|
| 3 | 0.97 | 0.97 | 0.98 | 0.95 | 0.95 | 0.95 |
| 5 | 0.96 | 0.95 | 0.96 | 0.94 | 0.94 | 0.96 |
| 8 | 0.96 | 0.91 | 0.92 | 0.98 | 0.88 | 0.94 |
| 12 | 0.89* | 0.90* | 0.80* | 0.91 | **0.78** | 0.88 |

\* fewer than 50 real bursts in at least one dataset.

- **Lowering SP does not help.** It adds real and false bursts in the same proportion: at SP 3-5, 94-98% of the detected bursts are matched by noise.
- **The window does not help either.** No window separates real bursts from noise better than another.
- **4PP7 alone:** the false-discovery share is ≈ 1.0 at every one of the 42 settings (0.70-1.53). At this noise level, no burst in 4PP7 can be told apart from noise.

At the chosen setting (excess turnover, %/h; z nested):

| comparison | A | B | z |
|---|---|---|---|
| 4PP7-FUS vs 4PP7 | 6.3 | −0.5 | 1.5 |
| 30PP7-FUS vs 30PP7 | 3.5 | 2.3 | 0.3 |
| 5Qb vs 4PP7 | 7.2 | −0.5 | 2.2 |
| 10Qb vs 30PP7 | 6.7 | 2.3 | 1.3 |
| in vitro 10:1 vs 1:1 | 51.1 | 8.2 | 3.2 |
| in vitro 100:1 vs 1:1 | 39.5 | 8.2 | 1.9 |

Across all 42 settings (z median [min-max]; settings with z > 1.96):

| comparison | z | settings significant |
|---|---|---|
| 4PP7-FUS vs 4PP7 | 1.4 [−0.6-2.4] | 8 / 42 |
| 30PP7-FUS vs 30PP7 | 0.0 [−1.6-1.2] | 0 / 42 |
| 5Qb vs 4PP7 | 2.3 [0.1-3.3] | 31 / 42 |
| 10Qb vs 30PP7 | 2.7 [0.3-4.3] | 32 / 42 |
| in vitro 10:1 vs 1:1 | 2.6 [1.2-3.5] | 39 / 42 |
| in vitro 100:1 vs 1:1 | 1.6 [1.1-2.2] | 4 / 42 |

**Interpretation (plain).**
- **The 4x FUS trend is consistent but weak.** FUS is ahead at almost every setting (excess turnover median 12 vs 2 %/h), but it reaches p < 0.05 at only 8 of 42 settings. When 42 settings are tried, some will cross p 0.05 by chance, so reporting one of those 8 would be cherry-picking. At the setting chosen by the noise rule, z = 1.5.
- **Changing detection does not rescue it.** The limit is the noise in the 4PP7 traces, not the detector settings: 4PP7 bursts look like noise at every SP and window.
- **More data would.** If the effect is real and the spread stays the same, about 5× more movies per condition (~50-65 movies) would give 80% power with M2b; ~3× would give an even chance of p < 0.05. Brighter or less noisy imaging of 4PP7 would help more than extra movies.
- **More robust than the FUS effect** (significant at most settings): Qβ shows more excess bursts than PP7 (5Qb vs 4PP7, 10Qb vs 30PP7), and in vitro 10:1 is more dynamic than 1:1. For Qβ this agrees with M1 in section 19, while M2b (section 19) did not separate Qβ from PP7. That is worth stating in the response.
- **30PP7-FUS vs 30PP7:** no difference at any setting.

## 21. Where the trace noise comes from: an independent check from the movies

Scripts: `checks/noise_from_background.m` → `checks/noise_background_summary.py`; `checks/noise_jitter_sources.py`; examples from `checks/example_traces.m` → `checks/plot_example_traces.py`. Outputs in `analysis_outputs/noise_background/`. Reviewer point: R2 6. In vivo only (no in vitro movies).

**Question.** The burst-free surrogates use the noise measured on each trace itself. Could that "noise" be fast real granule changes?

**Method.**
- **Prediction from pixels away from all granules.** For each movie, take every pixel > 10 px from any tracked point. Its frame-to-frame variance against its mean intensity is a straight line (R² ≥ 0.99), the EMCCD photon-noise law.
- **Predicted vs measured noise.** Feed each trace's spot and cell pixels (same crop and multithresh masks) through that line. This gives the expected camera + photon noise of the spot sum, the spot mean and the cell mean. Compare it with the noise actually measured on each trace.
- **Caveats on the prediction:**
  - Granule pixels are brighter than any calibration pixel (all spot pixels lie above the calibrated range), so the line is extrapolated 5-10×.
  - The intercept is negative.
  - The prediction is therefore approximate. The conclusions below do not depend on it.

**Measured ÷ predicted noise (median over traces ≥ 180 frames):**

| dataset | spot sum | spot mean | cell mean |
|---|---|---|---|
| 4PP7 | 4.1 | 2.4 | 4.3 |
| 4PP7-FUS | 3.7 | 2.2 | 3.8 |
| 30PP7 | 4.1 | 2.4 | 4.2 |
| 30PP7-FUS | 3.6 | 2.2 | 4.0 |
| 5Qb | 6.1 | 3.5 | 5.9 |
| 10Qb | 6.5 | 3.9 | 6.2 |

**Frame-to-frame jitter, relative to the level (robust):**

| dataset | spot mean | spot sum | cell mean | spot area (pixels in mask) | rank corr. of spot vs cell changes |
|---|---|---|---|---|---|
| 4PP7 | 3.8% | 16% | 6.7% | 16% | 0.79 |
| 4PP7-FUS | 5.3% | 22% | 7.6% | 24% | 0.82 |
| 30PP7 | 5.4% | 21% | 7.2% | 24% | 0.76 |
| 30PP7-FUS | 7.1% | 23% | 9.1% | 28% | 0.72 |
| 5Qb | 2.3% | 11% | 3.8% | 11% | 0.80 |
| 10Qb | 2.3% | 11% | 3.5% | 12% | 0.81 |

**Interpretation (plain).**
- **The jitter is mostly not camera/photon noise.** It is 2-6× larger than photon noise predicts.
- **The extra jitter is a measurement effect, not granule dynamics:**
  1. The cell background, which has no granule bursts, shows the same excess (4-6×).
  2. The spot and cell levels jump up and down together (correlation ~0.8).
  3. The number of pixels the threshold assigns to the granule changes by 11-28% from one frame to the next, i.e. every 10 s. A real granule does not change size that fast and back.
- **Cause.** The multithresh thresholds are recomputed in every frame, and each noisy frame moves the class boundaries. That changes which pixels count as "spot" and "cell" and shifts both class means together.
- **This answers the concern about the surrogates.** The noise level they use is not inflated by real granule dynamics. It is real measurement jitter, so the false-burst rates (sections 18-20) stand.
- **The spot sum (the "fraction shed" measure) is hit hardest.** Its jitter is 11-23% per frame, against 2-7% for the spot mean, because it multiplies by the jittering mask area. Dividing by the cell level frame by frame does not help (the two jitters are not proportional).
- **Implication.** This jitter is avoidable. Freezing the mask, or the thresholds, over time for each granule (e.g. thresholds from a running median) would remove most of it while keeping multithresh. If the jitter fell to the spot-mean level, the noise on the spot sum would drop about 3-4×. That is the main lever on the noise floor and on the 4x FUS power problem (section 20). It is an extraction change, not yet made.
- **Example traces** (`example_traces.png`): a 4PP7-FUS trace and its burst-free copy look alike raw and after smoothing, and the detector calls bursts in both. In the in vitro 10:1 example, the granule fades much faster than the dataset's typical bleaching curve. The detector calls that faster fade negative bursts, while the burst-free copy (typical curve) fades less.
- **Net loss vs bleaching.**
  - The detector divides by the cell-background trend, so a granule that fades faster than its surroundings keeps that net loss in the signal.
  - The surrogates use the dataset's typical granule curve, so a net loss shared by most granules cancels out of the "real minus noise" burst comparison.
  - Net loss is better measured directly: granule fade rate minus cell-background fade rate, per granule. Not yet done.

## 22. Net loss of granule material: granule fade rate vs cell-background fade rate

Script: `checks/net_loss.py` → `analysis_outputs/net_loss/` (`per_trace.csv`, `summary.csv`, `comparisons.csv`, `net_loss.png/pdf`). No burst detection, no noise model.

**Idea.** The cell background is the same fluorescent protein under the same light, so it bleaches at the same rate as the granule. A granule that fades faster than its surroundings is losing material; one that fades slower is gaining.

**Method.**
- Per granule (traces ≥ 180 frames), I fit robust straight-line slopes (Theil-Sen) to ln(intensity) over time, in %/h.
- Net change = slope of ln(granule ÷ cell).
- Granule measured two ways: the spot sum (in vivo only) and the spot mean above the cell.
- Errors: nested movie → granule bootstrap, 2000×. In vitro has no movie ids, so the bootstrap there is over granules.

**Rates, %/h, mean [95% CI]:**

| dataset | cell fade | net, spot sum | net, mean above cell | granules with net < 0 |
|---|---|---|---|---|
| 4PP7 | −7.3 [−12.1, −3.4] | +0.7 [−3.4, 4.2] | −1.3 [−3.8, 1.5] | 43-62% |
| 4PP7-FUS | −7.8 [−12.1, −4.2] | −6.6 [−14.8, 1.5] | −0.9 [−2.9, 1.0] | 54% |
| 30PP7 | +1.1 [−2.4, 5.4] | −5.1 [−11.0, 0.7] | +1.3 [−3.0, 5.0] | 41-59% |
| 30PP7-FUS | −10.5 [−17.6, −5.5] | −1.3 [−4.6, 2.0] | −1.3 [−4.7, 2.0] | 51-52% |
| 5Qb | −54.5 [−60.0, −49.7] | +3.5 [−0.2, 7.6] | +1.2 [−2.8, 5.4] | 41-43% |
| 10Qb | −55.4 [−62.7, −50.1] | +8.9 [5.3, 13.5] | −3.5 [−9.8, 1.3] | 27-56% |
| in vitro 1:1 | −46.9 [−67.7, −27.5] | – | −21.1 [−34.9, −8.5] | 62% |
| in vitro 10:1 | −19.8 [−25.5, −14.6] | – | −39.7 [−62.4, −20.6] | 84% |
| in vitro 100:1 | −102.7 [−136, −72.5] | – | −70.8 [−87.9, −53.9] | 94% |

**Comparisons (z):**
- 4PP7-FUS vs 4PP7: −1.6 (spot sum), +0.2 (mean above cell). The two measures disagree, so there is no effect.
- 30PP7-FUS vs 30PP7: +1.1 / −1.0.
- 10Qb vs 30PP7: +3.8 / −1.4. The measures disagree in sign, so this is not robust.
- In vitro 100:1 vs 1:1: −4.5. In vitro 10:1 vs 1:1: −1.5.

**Checks (in vitro).**
- The camera offset is assumed 200 for the in vitro data (not confirmed). With offset 0 or 100, net loss stays clearly negative in all three ratios (1:1 −29 to −35, 10:1 −41 to −43, 100:1 −85 to −94 %/h), and 100:1 remains the largest.
- The loss is faster in the second half of the movie (e.g. 100:1: −31 vs −87 %/h at offset 200). Part of this is arithmetic: the same absolute loss is a larger percentage once little granule is left above the background.

**Interpretation (plain).**
- **In vivo granules are in net balance.** In every in vivo condition they fade at the same rate as their surroundings, within a few %/h, i.e. by bleaching only. Any shedding is matched by uptake. About half the granules drift up and half down. No FUS effect, and no PP7 vs Qβ difference that both granule measures agree on.
- **In vitro granules lose material steadily**, beyond bleaching. The loss grows with the protein:RNA ratio: 1:1 ≈ −21, 10:1 ≈ −40, 100:1 ≈ −71 %/h. At 100:1, 94% of the granules lose material. This is the clearest, most robust between-condition difference in the whole analysis. It needs no burst detection and is independent of the noise problem.
- It fits the burst results: in vitro shows a real excess of negative bursts over noise (sections 18-20), and in vivo does not.
- **Caveat.** A granule drifting out of focus would also fade faster than the diffuse background. For in vitro (granules settling or moving on the glass) this cannot be excluded from intensity alone. Checking the spot area or size over time, or the focus, would help.
- **Bleaching itself:** PP7 movies bleach slowly (about 7-10 %/h). Qβ movies bleach about 55 %/h, probably reflecting different imaging settings (not checked).

## 23. Steady thresholds: less noise, but a slow artefact; checked with static-granule controls

Scripts:
- `collect_intensity.m` (`opts.thrWin`, default 0 = published behaviour unchanged);
- `checks/reextract_steady.m` → `in_vivo_reextracted/after_R4_steadythr/`;
- `checks/detection_grid.m` and `checks/noise_methods_real.m` (now take a data folder);
- `checks/steady_compare.py`, `checks/static_control.m`, `checks/static_control_summary.py`.

Outputs are in `analysis_outputs/steady/`. In vivo only: the in vitro movies are not on the lab PC (`InVitro_PCP14-MCP15` holds a single 2 MB 005.nd2), so in vitro could not be re-extracted.

**Change tested (R4).** Same tracking, crop, offset and multithresh classes as R2. But each track's four thresholds are replaced by their running median over 31 frames (~5 min, fixed in advance), instead of being recomputed on every frame. The masks then follow bleaching but no longer jump with each frame's noise. All six datasets keep exactly the same traces and lengths as R2.

**Noise drop (relative frame-to-frame jitter, median over traces):**

| dataset | spot sum R2 → R4 | spot mean | cell mean | spot area |
|---|---|---|---|---|
| 4PP7 | 16.4 → 7.1% | 3.8 → 1.5% | 6.7 → 2.0% | 16 → 6% |
| 4PP7-FUS | 21.5 → 10.4% | 5.3 → 2.0% | 7.6 → 2.1% | 24 → 9% |
| 30PP7 | 20.6 → 9.6% | 5.4 → 2.1% | 7.2 → 2.0% | 24 → 9% |
| 30PP7-FUS | 23.2 → 12.0% | 7.1 → 2.9% | 9.1 → 2.6% | 28 → 12% |
| 5Qb | 10.9 → 4.3% | 2.3 → 0.8% | 3.8 → 1.2% | 11 → 4% |
| 10Qb | 11.0 → 4.1% | 2.3 → 0.8% | 3.5 → 1.1% | 12 → 4% |

**Taken at face value (R4 vs R2, window 13, SP 12):**
- Burst counts go up 1.4-2×.
- The false-discovery share falls from 0.6-1.05 to 0.28-0.64.
- Excess turnover rises in every dataset (e.g. 4PP7 1.6 → 12.3 %/h).
- 4PP7-FUS vs 4PP7 becomes significant: z 1.7-3.1 on excess turnover over 6 settings, and z 5.3 on the noise-free fluctuation G30 (was 1.2).

**Negative control: static granules.**
- For each real 4PP7 and 4PP7-FUS track (≥ 24 points), I built a synthetic granule that never changes:
  - the median of its first 31 real crops;
  - kept centred along the track;
  - bleached with the movie's own background trend;
  - with fresh camera + photon noise per pixel and frame (variance a·I + c, calibrated per movie as in section 21).
- The same measurement (multithresh, per-frame or steady thresholds) and the same downstream treatment were then applied: surrogates, detection, structure function.

| | 4PP7 R2 | 4PP7 R4 | 4PP7-FUS R2 | 4PP7-FUS R4 |
|---|---|---|---|---|
| frame noise, real / static | 16.9 / 16.7% | 7.2 / 6.6% | 22.2 / 21.1% | 10.6 / 9.2% |
| excess turnover w13, real / static (%/h) | 1.6 / 1.6 | 12.3 / 13.0 | 3.0 / 2.3 | 22.7 / 14.2 |
| excess turnover w17, real / static | −0.5 / 0.7 | 12.0 / 11.6 | 6.3 / −2.1 | 23.1 / 13.3 |
| false-discovery share w13, real / static | 0.87 / 1.00 | 0.54 / 0.45 | 0.98 / 0.94 | 0.57 / 0.53 |
| G30, real / static (%) | 2.4 / 2.8 | 4.9 / 4.4 | 5.0 / 3.9 | 8.6 / 6.2 |

**What the control shows.**
1. **Section 21 is confirmed.** A static granule with nothing but photon noise reproduces the real frame noise almost exactly (16.7 vs 16.9% in 4PP7). The real jitter is photon noise amplified by per-frame threshold jitter.
2. **With per-frame thresholds (published), the white-noise surrogates are valid.** A static granule gives excess ≈ 0 and a false-discovery share ≈ 1. The noise-floor results of sections 18-20 stand. Real 4PP7 is indistinguishable from a static granule.
3. **Steady thresholds create a slow artefact.** A static granule reads 12-14 %/h of "excess turnover" and G30 4.4-6.2%. The running-median thresholds drift slowly, so the measured noise is no longer white and the white-noise surrogates no longer apply. Most of the face-value gains above are this artefact: real 4PP7 under R4 equals its static control.
4. **Measured against matched static controls** (real minus static, propagated bootstrap SEs):

| 4PP7-FUS vs 4PP7 | per-frame (R2) | steady (R4) |
|---|---|---|
| 4PP7, beyond static: turnover w13 / w17, G30² | 0.0 / −1.2 %/h, −2.3 | −0.7 / 0.4 %/h, 4.1 |
| 4PP7-FUS, beyond static | 0.8 / 8.3 %/h, 9.9 | 8.5 / 9.8 %/h, 36.0 (G ≈ 6%) |
| z, turnover w13 / w17 | 0.1 / 1.7 | 1.3 / 1.8 |
| z, G30² | 0.6 | 2.5 |
| × data for 80% power (best measure) | ≈ 3 | ≈ 1.3-2.3 |

**Interpretation (plain).**
- Holding the thresholds steady lowers the noise 2-3×, but it is not a free win: it adds a slow fake signal of its own. Results from R4 are only meaningful against a matched static-granule control, which is now built (`static_control.m`).
- **Against that control:**
  - 4PP7 granules look static: no dynamics beyond what a non-changing granule shows, with either extraction.
  - 4PP7-FUS granules show a modest excess with steady thresholds: about 9 %/h more turnover and ~6% fluctuation over 5 min beyond static. FUS vs 4PP7 reaches z 2.5 on the fluctuation measure and 1.3-1.8 on turnover. That is borderline: it is significant on one measure only, and it is better than R2 (z ≤ 1.7).
- **The trend is the same in both extractions (FUS > 4PP7).** A firm claim would still need roughly 1.5-2× more movies, or an extraction whose steady masks do not drift. Possible fix: hold one mask per granule over the whole trace and correct for bleaching separately. Not tried.
- **The static control is the strongest noise control in the revision.** It is a pixel-level negative control with the real track geometry and calibrated camera noise. For the response letter it should replace the trace-level white-noise surrogates.
- **Limits of the control:**
  - it ignores tracking jitter (the granule is kept centred);
  - it ignores real slow focus drift;
  - photon noise is extrapolated above the calibrated range (section 21).
  - It is therefore a lower bound on non-biological fluctuation.

## 24. Proposed report parameters per dataset (replacing λ and k0), with mock-ups of Figs. 7, 10(b) and 15

Scripts: `checks/report_params_sim.m` (simulations with known truth), `checks/static_control.m` (now all six in vivo datasets; `static_control_more.csv`), and `checks/report_params.py`. Outputs in `analysis_outputs/report_params/`: `params.csv`, `params.md`, `fig7_mock.png` (in vivo), `fig7_mock_sim.png`, `fig10b_mock.png`, `fig15_mock.png`. In vitro is not included: no raw movies, so no static control and no re-extraction.

**Parameters.** In vivo uses R2 (the published per-frame multithresh, spot sum), window 13, SP 12. 95% CIs: nested movie → granule bootstrap in vivo, granule bootstrap in the simulations.
1. **Turnover beyond noise (%/h).** Fraction of the granule shed per detected negative burst, summed per granule-hour, minus the same on burst-free copies of each trace.
2. **Fluctuation beyond noise over 5 min (%).** G30 = 100·√mean(D30 − D1) of the relative granule intensity. White noise cancels exactly. The pixel-level static-granule control gives mean(D30 − D1) consistent with 0 in all six in vivo datasets (within ~2 SE), so nothing is subtracted.
3. **Net change beyond bleaching (%/h).** Granule vs cell-background fade rate (section 22).
4. **Quality check: share of turnover explained by noise.** Turnover on the burst-free copies ÷ turnover on the traces. Size-weighted, because at high SNR noise gives many tiny bursts.

**Simulations (as manuscript Fig. 13).** generateSim type 3, k0 = 60, 200 traces each, seeded. True values come from the noise-free twin (same seed, noise 0).

| dataset | turnover beyond noise (%/h) | fluctuation beyond noise, 5 min (%) | net change (%/h) | share of turnover explained by noise | bursts (real / noise) | granules (movies) | true turnover / fluct. / net |
|---|---|---|---|---|---|---|---|
| sim: no bursts | 0.2 [-0.1, 0.5] | 0.0 [0.0, 0.1] | -0.0 [-0.1, 0.1] | 0.87 | 317 / 308 | 200 | 0.0 / 0.0 / 0.0 |
| sim: lambda 1 | 5.6 [4.6, 6.7] | 1.3 [1.1, 1.4] | -0.3 [-1.1, 0.4] | 0.21 | 541 / 304 | 200 | 4.0 / 1.3 / -0.3 |
| sim: lambda 2 | 10.5 [8.9, 12.0] | 1.9 [1.7, 2.1] | -0.9 [-1.9, 0.2] | 0.13 | 642 / 313 | 200 | 7.1 / 1.9 / -0.6 |
| sim: lambda 4 | 25.4 [21.0, 30.4] | 4.0 [3.6, 4.4] | 0.4 [-1.9, 3.0] | 0.05 | 729 / 260 | 200 | 15.4 / 4.0 / 0.2 |
| sim: lambda 2, SNR 5 | 10.8 [9.1, 12.9] | 2.0 [1.8, 2.2] | -0.8 [-1.9, 0.3] | 0.07 | 686 / 477 | 200 | 7.0 / 2.0 / -0.5 |
| sim: lambda 2, SNR 0.5 | 10.0 [7.6, 12.5] | 2.3 [2.0, 2.5] | -1.2 [-2.6, 0.2] | 0.32 | 460 / 275 | 200 | 8.3 / 2.3 / -1.2 |
| PP7 + QCP-5x/PCP-4x | 1.6 [-5.0, 7.6] | 2.4 [0.0, 4.7] | 0.7 [-3.0, 4.4] | 0.91 | 106 / 92 | 76 (9) | – |
| PP7-FUS + QCP-5x/PCP-4x | 3.0 [-5.5, 9.9] | 5.0 [0.0, 7.4] | -6.6 [-14.7, 1.6] | 0.88 | 124 / 122 | 94 (12) | – |
| PP7 + PCP-30x | 1.5 [-3.4, 7.5] | 3.7 [0.0, 5.6] | -5.1 [-11.2, 0.7] | 0.94 | 207 / 217 | 145 (8) | – |
| PP7-FUS + PCP-30x | 2.2 [-5.5, 10.4] | 2.8 [0.0, 6.8] | -1.3 [-4.4, 2.0] | 0.91 | 149 / 141 | 111 (8) | – |
| Qβ + QCP-5x/PCP-4x | 9.3 [5.0, 13.9] | 4.0 [3.0, 5.1] | 3.5 [-0.3, 7.4] | 0.47 | 313 / 204 | 191 (10) | – |
| Qβ + QCP-10x | 9.7 [6.5, 13.0] | 4.1 [3.1, 4.9] | 8.9 [5.0, 13.9] | 0.45 | 472 / 286 | 279 (9) | – |

**What the parameters do on the simulations.**
- **Fluctuation (G30) recovers the truth almost exactly** in every simulation, from SNR 0.5 to 5, and reads 0.0 without bursts.
- **Turnover beyond noise:**
  - reads ≈ 0 without bursts;
  - rises with λ at fixed k0 (5.6 → 10.5 → 25.4 %/h for λ 1, 2, 4);
  - does not change with noise level at fixed λ (10.0-10.8);
  - but over-reads the true turnover by 1.2-1.65×, because each detected run also picks up noise and smoothing.
  - It is therefore good for comparing conditions, not as an absolute rate.
- **Net change** reads the true value (≈ 0) within its CI.
- **None of the three separates λ from k0.** All three scale with burst rate × burst size, the identifiable combination (sections 5-9). They replace λ and k0, they do not estimate them.

**What the parameters say about the in vivo datasets.**
- **PP7 and PP7-FUS (both RNAs): no turnover beyond noise.** Noise explains 88-94% of the detected turnover, and the CIs include 0. Fluctuation is small (3-5% over 5 min) with CIs reaching 0. There is no net change beyond bleaching.
- **Qβ (both RNAs): clear turnover beyond noise**, about 9-10 %/h (CI excludes 0), with noise explaining only ~45%. Fluctuation is 4% [3-5] over 5 min. 10Qb shows a net gain beyond bleaching: +8.9 %/h [5, 14] by the spot sum, but section 22 notes the mean-above-cell measure disagrees.
- **Recommendation for the figures:**
  - Fig. 7/9 → `fig7_mock.png`: per dataset, the burst-size distribution of the granules over that of their burst-free copies (so the excess over noise is visible), plus the fluctuation-vs-lag curve.
  - Fig. 10(b) → `fig10b_mock.png`: turnover beyond noise × fluctuation, with CIs, colour = protein, marker = RNA. Also net change, and estimated vs true for the simulations.
  - Figs. 13-15 → `fig15_mock.png`: datasets × parameters, with CIs and (for the simulations) true values.
- **Claims that change:**
  - "Qβ granules are more dynamic than PP7 granules" is supported (turnover beyond noise).
  - "PP7-FUS increases λ for low-valency RNA" is not supported: no PP7 condition shows dynamics distinguishable from noise with this extraction (see section 23 for the steady-threshold test).

## 25. Why 4x PP7 vs PP7-FUS differ consistently: brightness vs dynamics; and why λ and k0 are not reported

Scripts: `checks/fus_brightness.py` → `analysis_outputs/fus_brightness/` (`pvalues.csv/.md`, `fus_brightness.png/pdf`); `checks/lambda_k0_degeneracy.py` → `analysis_outputs/lambda_k0_degeneracy/` (`degeneracy.png/pdf`, `summary.csv`).

### 25a. Brightness vs dynamics, with p values

- In vivo R2 (published extraction), spot sum, window 13, SP 12, traces ≥ 180 frames.
- **p values:** movie-level permutation (whole movies reassigned between the two conditions, 10000×), because movies are the replicates. For the raw burst sizes, the naive burst-level Mann-Whitney p (every burst counted as independent) is also shown.
- **CIs:** nested movie → granule bootstrap.

| measure | 4PP7-FUS vs 4PP7 | 30PP7-FUS vs 30PP7 | 5Qb vs 4PP7 | 10Qb vs 30PP7 |
|---|---|---|---|---|
| granule brightness (spot mean) | 3.04e+03 vs 5.64e+03, p = 0.0016 | 3.44e+03 vs 4e+03, p = 0.55 | 1.32e+04 vs 5.64e+03, p = 0.0002 | 1.75e+04 vs 4e+03, p = 0.0005 |
| frame noise (% of granule) | 22.2 vs 16.4, p = 0.002 | 23.3 vs 20.7, p = 0.23 | 11 vs 16.4, p = 0.0002 | 11 vs 20.7, p = 0.0002 |
| median burst size, granules (%) | 17.2 vs 11, p = 0.0006 (naive burst-level p = 8e-11) | 17.2 vs 15.3, p = 0.1 (naive burst-level p = 0.02) | 9.25 vs 11, p = 0.0006 (naive burst-level p = 0.001) | 9.14 vs 15.3, p = 0.0001 (naive burst-level p = 5e-31) |
| median burst size, burst-free copies (%) | 14.9 vs 12.3, p = 0.01 | 16.1 vs 14, p = 0.19 | 7.18 vs 12.3, p = 0.0001 | 6.55 vs 14, p = 0.0002 |
| burst size beyond noise (granules - copies, % points) | 2.26 vs -1.23, p = 0.062 | 1.05 vs 1.31, p = 0.78 | 2.07 vs -1.23, p = 0.0003 | 2.58 vs 1.31, p = 0.22 |
| turnover beyond noise (%/h) | 3.04 vs 1.61, p = 0.73 | 2.24 vs 1.52, p = 0.84 | 9.33 vs 1.61, p = 0.017 | 9.73 vs 1.52, p = 0.0019 |
| fluctuation over 5 min (%) | 4.98 vs 2.36, p = 0.17 | 2.79 vs 3.75, p = 0.67 | 3.96 vs 2.36, p = 0.052 | 4.06 vs 3.75, p = 0.66 |
| net change beyond bleaching (%/h) | -6.58 vs 0.72, p = 0.14 | -1.31 vs -5.14, p = 0.26 | 3.55 vs 0.72, p = 0.34 | 8.86 vs -5.14, p = 0.0001 |

**Reading (plain).**
- **4x pair:**
  - FUS granules are about half as bright (p = 0.002), so they are noisier relative to their size (22% vs 16% per frame, p = 0.002).
  - Their bursts are larger: 17.2% vs 11.0% of the granule. That holds even with movies as replicates (p = 0.0006); the naive burst-level test (p = 8e-11) greatly overstates the evidence.
  - Of the 6.1-point gap, 2.6 points are reproduced by burst-free copies with the same brightness and noise (p = 0.01). The remaining 3.5 points are beyond noise (p = 0.06, a trend).
  - Turnover beyond noise, fluctuation and net change do not differ (p 0.14-0.73).
- **30x pair:** similar brightness (p = 0.55) and similar noise (p = 0.23); no difference in any measure (p ≥ 0.1).
- **The consistency you noticed is therefore real, but mostly a brightness difference:**
  - FUS at low valency makes dimmer, smaller granules. Their apparent bursts are larger partly for that reason, plus a borderline extra.
  - At high valency the two proteins make granules of similar brightness, and nothing differs.
  - The brightness difference is itself a robust, reportable biological result. It fits the manuscript's lower k0 for PP7-FUS at 4x.
- **Qβ vs PP7:** brighter granules and less noise (p ≤ 0.0005), and clearly more turnover beyond noise (p = 0.017 and 0.002). 5Qb bursts exceed their noise copies (p = 0.0003 vs 4PP7).

### 25b. Why λ and k0 are not fitted

- **Model:** that of `gelmetricsMLE.m` in basic form: amplitude = n·k0 + Gaussian noise (σ = 71, from the quiescent segments), n ≥ 1 zero-truncated Poisson(λ).
- **Data:** 4PP7-FUS negative-burst amplitudes (R2, SP 10), between the 15th percentile (truncation) and the 90th percentile (top trim); n = 156.
- **Results:**
  - (a) λ = 0.3/k0 = 85, λ = 1/k0 = 56, λ = 2/k0 = 38 and λ = 4/k0 = 23 fit the same histogram equally well. The log-likelihood is within 0.6 of the best for all four, and KS p is 0.51-0.83. Only the mean burst size (≈ 87-99) is fixed by the data.
  - (b) The likelihood surface is a long ridge along the line of constant mean burst size. The 95% region spans λ from < 0.05 to > 8.
  - (c) Data simulated from three of these truths (λ 0.3, 1, 4) at the same sample size look alike: pairwise KS p 0.12, 0.31, 0.03.
  - (d) 200 datasets simulated from a known truth (λ = 1, k0 = 56), each fitted: the fitted λ is within 2× of the truth in only 44%. The 90% range of fitted λ is 0.05-5.7, and the estimates scatter along the ridge (k0 17-106).
- **This is why the revision reports the dynamics parameters of section 24 instead of λ and k0.** The manuscript's own Fig. 15(b) shows the same ridge: the best k0 falls steadily as λ is raised from 1 to 5, at a nearly unchanged WS.

## 26. Top-10% trim for the burst-size pipeline; and a shape-only fit of λ and k0

### 26a. Section 25a repeated with the largest 10% of bursts removed
`python checks/fus_brightness.py --trim 10` → `analysis_outputs/fus_brightness_trim10/`. Per dataset, the largest 10% of burst sizes are removed, separately for the granules and their burst-free copies (302 of 2941 bursts). Everything else is as in section 25a.

| measure | 4PP7-FUS vs 4PP7 | 30PP7-FUS vs 30PP7 | 5Qb vs 4PP7 | 10Qb vs 30PP7 |
|---|---|---|---|---|
| granule brightness (spot mean) | 3.04e+03 vs 5.64e+03, p = 0.0016 | 3.44e+03 vs 4e+03, p = 0.55 | 1.32e+04 vs 5.64e+03, p = 0.0002 | 1.75e+04 vs 4e+03, p = 0.0005 |
| frame noise (% of granule) | 22.2 vs 16.4, p = 0.002 | 23.3 vs 20.7, p = 0.23 | 11 vs 16.4, p = 0.0002 | 11 vs 20.7, p = 0.0002 |
| median burst size, granules (%) | 16.5 vs 10.6, p = 0.0009 (naive burst-level p = 1e-11) | 16.7 vs 14.5, p = 0.056 (naive burst-level p = 0.007) | 8.82 vs 10.6, p = 0.0002 (naive burst-level p = 0.0003) | 8.68 vs 14.5, p = 0.0001 (naive burst-level p = 6e-34) |
| median burst size, burst-free copies (%) | 14.5 vs 11.8, p = 0.0034 | 15.6 vs 13.3, p = 0.071 | 6.55 vs 11.8, p = 0.0001 | 6.19 vs 13.3, p = 0.0002 |
| burst size beyond noise (granules - copies, % points) | 2 vs -1.24, p = 0.073 | 1.09 vs 1.2, p = 0.88 | 2.27 vs -1.24, p = 0.0003 | 2.49 vs 1.2, p = 0.24 |
| turnover beyond noise (%/h) | 2.83 vs 1.28, p = 0.66 | 1.64 vs 1.29, p = 0.91 | 7.53 vs 1.28, p = 0.018 | 7.77 vs 1.29, p = 0.0049 |
| fluctuation over 5 min (%) | 4.98 vs 2.36, p = 0.17 | 2.79 vs 3.75, p = 0.67 | 3.96 vs 2.36, p = 0.052 | 4.06 vs 3.75, p = 0.66 |
| net change beyond bleaching (%/h) | -6.58 vs 0.72, p = 0.14 | -1.31 vs -5.14, p = 0.26 | 3.55 vs 0.72, p = 0.34 | 8.86 vs -5.14, p = 0.0001 |

The trim changes nothing qualitatively:
- 4x FUS vs PP7: burst size 16.5 vs 10.6% (p = 0.0009); the noise-only copies differ too (p = 0.003); beyond noise, a trend (p = 0.07); no difference in dynamics beyond noise.
- 30x: nothing (p ≥ 0.056).
- Qβ vs PP7: more turnover beyond noise (p = 0.018 and 0.005).

### 26b. Shape-only fit of λ and k0
- **How section 25b fitted (not WS):** the likelihood model with a Gaussian noise term, on the 15-90% window.
- **Here:** the manuscript's own model curve, exactly as `modifiedPoissonFit.m` builds it after B1:
  - the zero-truncated Poisson pmf at whole numbers of units, joined by straight lines, evaluated at A/k0;
  - no noise term;
  - normalised into a density over [0, A90];
  - fitted to every negative burst except the largest 10% (no lower 15% cut), by maximum likelihood.
- The normalisation removes the overall height, so only the shape is fitted: no binning, no WS ranks. Bounds: λ 0.05-60.
- A continuous Gamma-function version was tried first. It degenerates to an exponential (λ → 0, k0 → ∞) for any distribution that falls from small to large amplitudes, so it was dropped.
- Data: spot mean (manuscript units), R2 in vivo, stored in vitro, window 13, SP 12 (`checks/export_amps_shape.m`).
- CIs: 200× granule bootstrap, nested in movies in vivo. p values: 1000× movie-level permutation, with the trim recomputed and the model refitted in every permuted group (granules for in vitro: no movie ids).
- Scripts: `checks/lambda_k0_shape.py` → `analysis_outputs/lambda_k0_shape/` (`fits.csv/.md`, `comparisons.csv`, `shape_fit.png/pdf`).

| dataset | bursts (fitted) | λ [95% CI] | k0 [95% CI] | mean burst λ·k0 [95% CI] |
|---|---|---|---|---|
| PP7, 4x | 123 (110) | 0.55 [0.05, 8.60] | 174.5 [24.5, 196.8] | 95 [9, 550] |
| PP7-FUS, 4x | 144 (129) | 0.44 [0.19, 8.01] | 114.5 [14.4, 133.0] | 50 [21, 214] |
| PP7, 30x | 225 (202) | 0.84 [0.34, 2.81] | 155.9 [57.3, 220.4] | 130 [51, 275] |
| PP7-FUS, 30x | 156 (140) | 0.72 [0.33, 6.16] | 176.3 [35.8, 207.7] | 126 [52, 265] |
| Qβ, 5x | 336 (302) | 0.57 [0.23, 5.67] | 331.9 [70.7, 367.7] | 191 [82, 425] |
| Qβ, 10x | 497 (447) | 5.46 [0.06, 6.50] | 98.2 [76.9, 515.2] | 536 [34, 602] |
| in vitro 1:1 | 94 (84) | 1.53 [1.26, 1.99] | 26.4 [16.9, 37.8] | 41 [26, 66] |
| in vitro 10:1 | 93 (83) | 1.09 [0.74, 1.64] | 87.7 [52.0, 131.3] | 96 [49, 212] |
| in vitro 100:1 | 105 (94) | 1.45 [1.08, 1.97] | 115.1 [47.4, 272.5] | 166 [87, 433] |
| sim lambda 0.5, noise 12 | 483 (434) | 2.10 [1.89, 2.26] | 16.9 [16.1, 17.5] | 35 [31, 39] |
| sim lambda 0.5, noise 40 | 429 (386) | 0.99 [0.84, 1.23] | 35.7 [33.3, 37.8] | 35 [30, 42] |
| sim lambda 0.5, noise 120 | 357 (321) | 5.13 [0.71, 6.83] | 22.8 [16.3, 96.3] | 117 [68, 129] |
| sim lambda 1, noise 12 | 615 (553) | 2.40 [2.06, 2.44] | 23.9 [23.0, 28.7] | 57 [52, 61] |
| sim lambda 1, noise 40 | 511 (460) | 1.80 [1.61, 2.02] | 39.4 [37.9, 41.2] | 71 [63, 79] |
| sim lambda 1, noise 120 | 370 (333) | 0.86 [0.44, 5.54] | 94.0 [22.6, 111.1] | 81 [49, 128] |
| sim lambda 2, noise 12 | 680 (612) | 2.18 [1.50, 2.42] | 42.7 [37.4, 58.3] | 93 [85, 99] |
| sim lambda 2, noise 40 | 654 (588) | 2.16 [1.91, 2.33] | 51.3 [47.9, 56.6] | 111 [101, 121] |
| sim lambda 2, noise 120 | 467 (420) | 1.50 [0.87, 1.82] | 97.7 [92.2, 125.7] | 146 [107, 171] |
| sim lambda 3, noise 12 | 737 (663) | 2.26 [2.09, 2.39] | 64.3 [63.4, 68.5] | 145 [134, 157] |
| sim lambda 3, noise 40 | 680 (612) | 2.27 [2.14, 2.41] | 65.8 [62.7, 69.0] | 149 [140, 161] |
| sim lambda 3, noise 120 | 544 (489) | 1.78 [1.55, 2.20] | 105.5 [90.8, 118.7] | 188 [173, 209] |
| sim lambda 4, noise 12 | 740 (666) | 2.51 [2.29, 2.63] | 71.6 [69.7, 84.7] | 180 [169, 195] |
| sim lambda 4, noise 40 | 763 (686) | 2.48 [2.31, 2.57] | 78.7 [75.4, 85.1] | 195 [180, 209] |
| sim lambda 4, noise 120 | 635 (571) | 2.05 [1.84, 2.22] | 110.3 [104.9, 120.1] | 226 [207, 245] |
| sim lambda 5, noise 12 | 786 (707) | 2.38 [2.22, 2.50] | 98.8 [96.2, 103.8] | 235 [217, 257] |
| sim lambda 5, noise 40 | 735 (661) | 2.60 [2.45, 2.78] | 107.2 [102.0, 112.2] | 279 [260, 300] |
| sim lambda 5, noise 120 | 705 (634) | 2.28 [2.11, 2.49] | 121.5 [117.2, 127.3] | 277 [257, 304] |

**Comparisons (p, permutation):**

| comparison | λ | p (λ) | k0 | p (k0) |
|---|---|---|---|---|
| 4PP7-FUS vs 4PP7 | 0.44 vs 0.55 | 0.98 | 114 vs 175 | 0.61 |
| 30PP7-FUS vs 30PP7 | 0.72 vs 0.84 | 0.78 | 176 vs 156 | 0.50 |
| 5Qb vs 4PP7 | 0.58 vs 0.55 | 0.98 | 332 vs 175 | 0.058 |
| 10Qb vs 30PP7 | 5.46 vs 0.84 | 0.001 | 98 vs 156 | 0.048 |
| in vitro 10:1 vs 1:1 | 1.09 vs 1.53 | 0.33 | 88 vs 26 | 0.008 |
| in vitro 100:1 vs 1:1 | 1.45 vs 1.53 | 0.90 | 115 vs 26 | 0.007 |

**What the shape-fitted λ measures (simulations, k0 = 60, same detection pipeline).**
- At SNR 5 and 1.5, the fitted λ is 2.1-2.6 for every true λ from 1 to 5, and 1-2 at true λ 0.5.
- At SNR 0.5 it is erratic (0.9-5.1).
- The fitted k0 is not 60 either: it grows with the true λ, from 17 to 121. Only the mean burst size λ·k0 grows steadily with the true λ.
- The shape of the detected-burst distribution is set mainly by the detector, which removes small bursts and merges close ones, and by the noise. It is not set by the number of units per burst.

**Interpretation (plain).**
- Fitting the manuscript's curve to the shape alone gives clean-looking fits and decent CIs in vitro (λ 1.1-1.5), but the simulations show this λ does not track the true number of units.
- In vivo, λ's CIs span 0.05-8 in most datasets: the shape does not pin it down there either.
- Differences in k0 between in vitro ratios (p ≈ 0.007) are differences in typical burst size, i.e. the mean λ·k0. That agrees with the burst-size and turnover results.
- The one λ difference (10Qb, 5.5 vs 0.8, p = 0.001) is not trustworthy: the 10Qb λ CI is 0.06-6.5 (two near-equal solutions).
- **This supports reporting burst size and turnover with their noise controls (sections 24-25), not λ and k0.**

## 27. Mean burst size (the identifiable λ·k0) and a minimal control table

Script: `checks/control_table.py` → `analysis_outputs/control_table/` (`control_table.csv/.md/.png`).

- **Mean burst size.** The mean negative-burst amplitude over **all** bursts (no top-10% trim, no 15th-percentile cut), in spot-mean units, window 13, SP 12.
  - It is the combination of λ and k0 that the data pin down. For the zero-truncated Poisson model it equals λ·k0/(1 − e^−λ), i.e. ≈ λ·k0 for λ ≳ 2 and ≈ k0 for small λ. λ and k0 separately are not identified (sections 25-26).
  - 95% CI: nested movie → granule bootstrap (in vitro: granules), 2000×.
- **Minimal control table**, one row per dataset. Each parameter is given with its noise control: frame noise, and the share of turnover that burst-free copies reproduce.

| dataset | granules (movies) | granule brightness (a.u.) | frame noise (% of granule) | mean burst size, λ·k0 (a.u.) | turnover beyond noise (%/h) | share of turnover explained by noise | net change beyond bleaching (%/h) |
|---|---|---|---|---|---|---|---|
| PP7 + QCP-5x/PCP-4x | 76 (9) | 5644 [4645, 6267] | 16.4 | 190 [170, 215] (n = 123) | 1.6 [-5.1, 7.9] | 0.91 | 0.7 [-3.5, 4.3] |
| PP7-FUS + QCP-5x/PCP-4x | 94 (12) | 3036 [1843, 3735] | 22.2 | 182 [113, 269] (n = 144) | 3.0 [-5.6, 10.1] | 0.88 | -6.6 [-15.0, 1.5] |
| PP7 + PCP-30x | 145 (8) | 3998 [2305, 9095] | 20.7 | 295 [166, 433] (n = 225) | 1.5 [-3.4, 7.0] | 0.94 | -5.1 [-11.4, 0.7] |
| PP7-FUS + PCP-30x | 111 (8) | 3436 [2777, 3932] | 23.3 | 253 [195, 304] (n = 156) | 2.2 [-5.4, 10.0] | 0.91 | -1.3 [-4.6, 1.9] |
| Qβ + QCP-5x/PCP-4x | 191 (10) | 13197 [11838, 14888] | 11.0 | 461 [395, 548] (n = 336) | 9.3 [5.2, 13.7] | 0.47 | 3.5 [-0.2, 7.7] |
| Qβ + QCP-10x | 279 (9) | 17534 [16083, 19160] | 11.0 | 574 [501, 660] (n = 497) | 9.7 [6.2, 13.3] | 0.45 | 8.9 [5.1, 13.6] |
| in vitro, protein:RNA 1:1 | 45 | 914 [660, 1238] | 3.0 | 129 [51, 261] (n = 94) | 8.8 [-2.5, 21.2] | 0.70 | -21.1 [-35.3, -8.3] |
| in vitro, protein:RNA 10:1 | 37 | 1715 [1527, 1796] | 3.0 | 234 [149, 318] (n = 93) | 50.1 [25.1, 76.2] | 0.30 | -39.7 [-63.0, -20.4] |
| in vitro, protein:RNA 100:1 | 31 | 1726 [809, 2714] | 2.6 | 432 [205, 728] (n = 105) | 46.2 [17.5, 79.3] | 0.51 | -70.8 [-88.9, -55.1] |

**Notes.**
- Brightness and mean burst size are in spot-mean units (as the manuscript's k0). Turnover is from the spot sum, as a fraction of the granule.
- A noise share near 1 means the detected turnover is about what noise alone produces, so it should not be interpreted as dynamics (all PP7 rows).
- In vitro has no movie ids, so its CIs are granule-level only.

## 28. Can noise explain the whole FUS-PP7 gap? (hypothesis test); λ·k0 by shape fit on all bursts; updated control table

### 28a. Hypothesis test: is the burst-size gap entirely explained by noise?
- Script: `checks/gap_test.py` → `analysis_outputs/gap_test/` (`tests.csv/.md`, `per_dataset.csv`).
- **H0:** the gap between the conditions' granule bursts equals the gap between their burst-free copies, which have the same brightness and noise but no bursts. In other words, the part beyond noise is 0.
- **Statistic:** the difference in differences of the median fraction shed per burst. Sensitivity check: the mean of log(fraction), which uses every burst.
- In vivo R2, all bursts, window 13, SP 12.
- **p:** movie-level permutation, 10000×. One-sided only where the manuscript states the direction (FUS > PP7).
- **CIs:** nested bootstrap, 5000×.

| comparison | statistic | total gap [95% CI] | explained by noise [95% CI] | beyond noise [95% CI] | p, two-sided | p, one-sided |
|---|---|---|---|---|---|---|
| 4PP7-FUS vs 4PP7 | median size (% points) | 6.1 [2.9, 8.8] | 2.6 [0.1, 5.2] | 3.5 [-0.0, 6.5] | 0.061 | 0.031 (pre-specified) |
| 4PP7-FUS vs 4PP7 | mean log size | 0.40 [0.20, 0.60] | 0.24 [-0.01, 0.51] | 0.16 [-0.11, 0.43] | 0.174 | 0.088 (pre-specified) |
| 30PP7-FUS vs 30PP7 | median size (% points) | 1.9 [-1.1, 4.4] | 2.2 [-0.5, 4.5] | -0.3 [-3.2, 2.2] | 0.789 | 0.388 |
| 30PP7-FUS vs 30PP7 | mean log size | 0.14 [-0.05, 0.33] | 0.25 [0.08, 0.42] | -0.10 [-0.33, 0.10] | 0.305 | 0.151 |
| 5Qb vs 4PP7 | median size (% points) | -1.8 [-4.0, -0.3] | -5.1 [-7.1, -3.5] | 3.3 [0.8, 5.4] | 0.000 | 0.000 |
| 5Qb vs 4PP7 | mean log size | -0.16 [-0.31, -0.00] | -0.58 [-0.82, -0.32] | 0.42 [0.14, 0.69] | 0.001 | 0.001 |
| 10Qb vs 30PP7 | median size (% points) | -6.1 [-9.0, -3.9] | -7.4 [-9.9, -5.5] | 1.3 [-1.2, 3.6] | 0.205 | 0.106 |
| 10Qb vs 30PP7 | mean log size | -0.49 [-0.68, -0.32] | -0.68 [-0.86, -0.50] | 0.19 [-0.02, 0.38] | 0.034 | 0.016 |

**Reading.**
- **4x:** 3.5 of the 6.1 points (57%) are beyond noise.
  - With the manuscript's direction pre-specified, H0 is rejected (one-sided p = 0.031; two-sided p = 0.061).
  - With the every-burst statistic (mean log size), H0 is not rejected (one-sided p = 0.088).
  - So the beyond-noise FUS effect at 4x is supported at the 5% level by one reasonable statistic, but not robustly across statistics.
- **30x:** noise explains the whole gap.
- **5Qb vs 4PP7:** the part beyond noise is clear (p ≤ 0.001 on both statistics).

### 28b. λ·k0 by shape fitting, all bursts (no trim, no truncation)
- Script: `python checks/lambda_k0_shape.py --notrim` → `analysis_outputs/lambda_k0_shape_notrim/`. The model is as in section 26b, normalised over [0, max amplitude].

| dataset | bursts (fitted) | λ [95% CI] | k0 [95% CI] | mean burst λ·k0 [95% CI] |
|---|---|---|---|---|
| PP7, 4x | 123 (123) | 0.14 [0.05, 7.95] | 174.5 [24.6, 187.8] | 24 [9, 219] |
| PP7-FUS, 4x | 144 (144) | 0.91 [0.13, 1.86] | 123.8 [85.7, 143.2] | 112 [15, 259] |
| PP7, 30x | 225 (225) | 1.21 [0.34, 1.59] | 177.5 [117.8, 289.3] | 215 [52, 335] |
| PP7-FUS, 30x | 156 (156) | 0.75 [0.34, 1.70] | 179.1 [97.5, 212.7] | 134 [70, 206] |
| Qβ, 5x | 336 (336) | 0.62 [0.40, 1.11] | 350.9 [256.8, 386.6] | 217 [139, 380] |
| Qβ, 10x | 497 (497) | 0.36 [0.10, 5.22] | 478.9 [103.6, 573.0] | 174 [51, 548] |
| in vitro 1:1 | 94 (94) | 0.88 [0.84, 1.61] | 136.1 [33.0, 282.1] | 119 [42, 297] |
| in vitro 10:1 | 93 (93) | 1.45 [1.20, 1.73] | 152.1 [101.2, 186.0] | 220 [130, 316] |
| in vitro 100:1 | 105 (105) | 1.32 [1.10, 1.61] | 327.8 [164.1, 468.1] | 434 [188, 754] |

- **Simulations:** the shape-fitted λ·k0 tracks the true λ·k0 = 60λ (panel c of `shape_fit.png`). At SNR 1.5-5 it is within 0.8-1.35× for true λ·k0 60-300, and it over-reads at true 30 (1.5-1.6×). λ alone stays at about 2 regardless of the truth (as in section 26).
- **In vivo, λ is mostly below 1.** There, the product λ·k0 is unstable, e.g. 4PP7 24 [9, 219]: almost every burst is one "unit", so the product understates the burst size, and its CI is wide. The model's mean burst size, λ·k0/(1 − e^−λ), and the model-free mean burst size (190 for 4PP7) are the stable quantities.
- **Movie-level p for λ·k0:**
  - 4PP7-FUS vs 4PP7: 112 vs 24, p = 0.095;
  - 30x pair: p = 0.54;
  - 5Qb vs 4PP7: p = 0.001;
  - 10Qb vs 30PP7: p = 0.33;
  - in vitro 10:1 and 100:1 vs 1:1: p = 0.29 and 0.071.

### 28c. Control table, updated
`analysis_outputs/control_table/control_table.png/.md/.csv`. Added:
- the median burst size of the granules vs their noise-only copies;
- the burst size beyond noise, with 95% CI;
- λ·k0 from the shape fit.

| dataset | granules (movies) | granule brightness (a.u.) | frame noise (% of granule) | burst size, granules / noise copies (% of granule) | burst size beyond noise (% points) | mean burst size (a.u.) | λ·k0, shape fit (a.u.) | turnover beyond noise (%/h) | net change beyond bleaching (%/h) |
|---|---|---|---|---|---|---|---|---|---|
| PP7 + QCP-5x/PCP-4x | 76 (9) | 5644 [4645, 6267] | 16.4 | 11.0 / 12.3 | -1.2 [-3.0, 1.0] | 190 [170, 215] | 24 [9, 219] | 1.6 [-5.1, 7.9] | 0.7 [-3.5, 4.4] |
| PP7-FUS + QCP-5x/PCP-4x | 94 (12) | 3036 [1861, 3735] | 22.2 | 17.2 / 14.9 | 2.3 [-0.1, 4.9] | 182 [113, 269] | 112 [15, 259] | 3.0 [-5.9, 10.4] | -6.6 [-14.9, 1.5] |
| PP7 + PCP-30x | 145 (8) | 3998 [2305, 8978] | 20.7 | 15.3 / 14.0 | 1.3 [-0.7, 3.5] | 295 [165, 436] | 215 [52, 335] | 1.5 [-3.9, 7.5] | -5.1 [-11.5, 0.7] |
| PP7-FUS + PCP-30x | 111 (8) | 3436 [2645, 3932] | 23.3 | 17.2 / 16.1 | 1.0 [-1.0, 2.8] | 253 [196, 305] | 134 [70, 206] | 2.2 [-5.9, 10.8] | -1.3 [-4.5, 1.9] |
| Qβ + QCP-5x/PCP-4x | 191 (10) | 13197 [11812, 14982] | 11.0 | 9.2 / 7.2 | 2.1 [1.0, 3.2] | 461 [394, 545] | 217 [139, 380] | 9.3 [5.3, 13.5] | 3.5 [-0.3, 7.5] |
| Qβ + QCP-10x | 279 (9) | 17534 [16069, 18904] | 11.0 | 9.1 / 6.6 | 2.6 [1.6, 3.8] | 574 [502, 658] | 174 [51, 548] | 9.7 [6.4, 13.3] | 8.9 [5.4, 13.5] |
| in vitro, protein:RNA 1:1 | 45 | 914 [660, 1238] | 3.0 | 7.9 / 10.4 | -2.5 [-6.4, 1.3] | 129 [52, 267] | 119 [42, 297] | 8.8 [-2.8, 21.5] | -21.1 [-35.1, -8.8] |
| in vitro, protein:RNA 10:1 | 37 | 1715 [1527, 1796] | 3.0 | 12.2 / 8.8 | 3.4 [0.7, 7.3] | 234 [151, 321] | 220 [130, 316] | 50.1 [27.9, 74.2] | -39.7 [-62.5, -20.2] |
| in vitro, protein:RNA 100:1 | 31 | 1726 [809, 2714] | 2.6 | 15.2 / 10.0 | 5.1 [0.3, 7.9] | 432 [200, 725] | 434 [188, 754] | 46.2 [17.0, 78.7] | -70.8 [-88.1, -54.3] |

## 29. Median burst size added; median vs mean; summary deck

- `analysis_outputs/control_table/` now also lists the median burst size [95% CI] (spot-mean units, all bursts).
- **Median vs mean.** The burst-size distributions are right-skewed: a few very large bursts, possibly merged events. Skewness is 0.5-8.9; e.g. in vitro 1:1 has median 36, mean 129, max 5272.
  - **Mean:** what λ·k0 represents, and what adds up to total turnover (rate × mean). It is dominated by the rare largest bursts, so its CIs are wide and it is sensitive to merged events.
  - **Median:** the typical burst. Robust to the tail, tighter CIs; use it for comparing conditions, as in all tests in sections 25-28.
  - **Report the mean only when relating to λ·k0 or total turnover.**
  - Noteworthy: 4PP7-FUS has a smaller median (132 vs 187) but a similar mean (182 vs 190), with skewness 2.2 vs 0.5. FUS has more small bursts and a longer tail of large ones.
- **Summary deck for the PI meeting:** `analysis_outputs/deck/GelMetrics_revision_summary.pptx`, 12 slides with speaker notes, built by `checks/make_deck.py` from the figures above.

## 30. The manuscript's approach revisited: Eq. (1.9) on simulated data, six ways of scoring the fit

Scripts: `checks/ms_sim_test.m` (simulation + detection) and `checks/ms_eq19_fit.py` (fits). Outputs in `analysis_outputs/ms_sim_test/`: `traces.png`, `fits_examples.png`, `recovery.png`, `fits.csv`, `fits_summary.csv`.

**Setup (as MS section III).**
- Simulation: generateSim type 3, 1000 signals × 361 frames, k0 = 60, white noise 40 a.u., bleaching 0.0026/min.
- λ ∈ {0.5, 1, 2, 3, 4, 5}, 5 independent datasets each.
- Detection: analysis.m, window 13, SP 15; negative bursts; **nothing removed** (no 15th-percentile cut, no top trim).
- Controls:
  - the same simulation at noise 4 (almost noise-free);
  - the true jump sizes of the noise-free twin, i.e. what Eq. (1.9) is meant to describe, without noise or detection.
- **Model:** Eq. (1.9) exactly as written, p(I) = (1/k0)·λ^(I/k0)·e^(−λ)/(I/k0)!, with no noise term.
- **Scores** (same grid, λ 0.05-20 plus 1..5, k0 2-600): MSE (10-bin histogram), negative log-likelihood, Kolmogorov-Smirnov, Cramér-von Mises, QQ RMSE, the manuscript's WS, and WS restricted to integer λ 1-5 (as in the MS).

**Fitted λ (median of 5 datasets):**

| true λ | 0.5 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| detected, noise 40 (MS), all scores | 2.4-3.0 | 2.2-3.0 | 2.0-2.4 | 2.0-2.6 | 2.4-3.0 | 2.4-3.2 |
| detected, noise 4, all scores | 0.1-1.0 | 0.7-3.5 | 0.1-2.4 | 0.1-1.8 | 0.6-2.2 | 0.9-3.0 |
| true jumps, NLL / CvM / KS | 6.9 / 4.0 / 0.2 | 4.4 / 3.0 / 3.5 | 3.5 / 3.0 / 3.2 | 3.8 / 3.2 / 3.2 | 4.4 / 4.0 / 4.0 | 5.1 / 5.0 / 5.0 |

**Fitted k0 (true 60), detected, noise 40:** 26-34 at λ 0.5, rising to 99-137 at λ 5, for every score.

**Fitted λ·k0 (true 60λ), detected, noise 40:**
- 76-81 for true 30, 102-108 for 60, 154-158 for 120, 205-210 for 180, 256-276 for 240, 313-335 for 300.
- It tracks the truth with an upward offset that shrinks from +50 to +20.
- At noise 4 it is close to the truth for the likelihood, CvM and QQ scores.

**Goodness of fit:** the best (likelihood) curve is rejected by a KS test for most settings (p < 0.01 at noise 40; p < 0.001 at noise 4).

**Mean burst size vs truth:**

| true λ | 0.5 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| true mean jump | 76 | 94 | 141 | 190 | 245 | 302 |
| detected, noise 4 | 62 | 88 | 138 | 188 | 247 | 305 |
| detected, noise 40 | 83 | 111 | 164 | 213 | 268 | 325 |

**Why the manuscript's approach does not work (plain).**
1. **The detector reshapes the burst distribution.**
   - Smoothing over 13 frames turns each step into a ramp. Steps closer together than the window merge into one larger "burst", and small steps are missed (example traces: a trace with two clear drops gives 0 detected bursts).
   - The detected distribution therefore always has the same smooth hump shape, whatever the true number of units. Fitted λ ends up at about 2-3 for every true λ from 0.5 to 5, with every scoring method. This is the simulated-vs-estimated λ discrepancy the referee pointed out (2 vs 4).
2. **The model does not match how bursts are generated, even on perfect data.**
   - Bursts carry 1, 2, 3, ... whole units, and a burst with 0 units is invisible. Eq. (1.9) is a smooth curve that does not account for either.
   - On the true jumps it recovers λ only for λ ≥ 3, and overestimates it below (λ 0.5 → 4-7; λ 1 → 3-4.4).
3. **The answer depends on how the fit is scored.**
   - With little noise, the six scores give λ from 0.1 to 3.5 for the same data.
   - The integer λ grid of the manuscript hides this by forcing one of five values.
4. **What survives:** the product λ·k0, i.e. the mean burst size. It follows the truth for every score and every noise level, with a modest upward bias at the MS noise level. Separately λ and k0 do not.

**Story for the response (suggested wording).**
"We first estimated λ and k0 by fitting Eq. (1.9) to the burst-amplitude histogram. On simulated data with known parameters this does not recover λ or k0, with any of six scoring methods (MSE, likelihood, KS, Cramér-von Mises, QQ, WS). The burst detector smooths and merges events, so the detected-burst distribution has nearly the same shape for any true λ (fitted λ ≈ 2-3 for true λ 0.5-5). Only the product λ·k0, i.e. the mean burst size, is recovered reliably. We therefore report burst size (median and mean) with a noise control, rather than λ and k0 separately."

## 31. MS Fig. 13 recreated for every scoring method; was the true λ a close runner-up?

Script: `checks/ms_fig13.py` (simulated rows from `checks/ms_sim_test.m`, plus `ms_sim_test(2, [12 120], '_snr')` for the Fig. 13 rows λ = 2 at SNR 5 and 0.5). Outputs in `analysis_outputs/ms_sim_test/`: `fig13_heatmaps.png` (per score: score and k0 at each λest = 1-5), `fig13_choice.png`, `fig13_scores.csv`, `fig13_choice.csv`, `fig13_choice_summary.csv`.

**Setup.**
- Eq. (1.9) as written, all negative bursts, k0 = 60.
- Candidate λ = 1-5 as in the manuscript, k0 optimised for each.
- Rows: true λ 0.5-5 at SNR 1.5 (noise 40, the MS default), and λ 2 at SNR 5 and 0.5. 5 datasets per row.
- **"Close runner-up" test:**
  - 100 granule (trace) bootstrap resamples per dataset (500 per row), refitted: how often does each λ win?
  - For the likelihood, 2·ΔNLL between the true λ and the winner. Above about 4, the true λ is rejected.

**How often the true λ is chosen (500 resampled fits per row):**

| row | MSE | NLL | KS | CvM | QQ | WS |
|---|---|---|---|---|---|---|
| λ = 1, SNR 1.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| λ = 2, SNR 1.5 | 0.87 | 0.59 | 0.75 | 0.77 | 0.72 | 0.63 |
| λ = 3, SNR 1.5 | 0.15 | 0.50 | 0.35 | 0.59 | 0.69 | 0.31 |
| λ = 4, SNR 1.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| λ = 5, SNR 1.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| λ = 2, SNR 5 | 0.68 | 0.16 | 0.09 | 0.28 | 0.91 | 0.97 |
| λ = 2, SNR 0.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

(True λ 0.5 is not among the candidates; it is fitted as 2-3.)

**Likelihood, 2·(NLL(λest) − NLL(best)), median of 5 datasets:**

| row | λest 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| λ = 1 | 124 | 10 | **0** | 68 | 195 |
| λ = 2 | 143 | **0** | 4 | 117 | 311 |
| λ = 3 | 171 | **0** | 1 | 130 | 355 |
| λ = 4 | 191 | 3 | **0** | 136 | 373 |
| λ = 5 | 239 | 25 | **0** | 115 | 334 |
| λ = 2, SNR 5 | **0** | 22 | 228 | 571 | 1018 |
| λ = 2, SNR 0.5 | 304 | 128 | 33 | **0** | 11 |

**k0 if λ were fixed at the truth** (true 60): 63-77 for λ ≥ 2 at SNR ≥ 1.5, but 84-99 at λ 1 and 99-112 at SNR 0.5.

**Reading (plain).**
- **The true λ is not a close runner-up**, except when it happens to be 2 or 3.
  - At true λ 4 the true value is never chosen (0 of 500 resamples, any score), and the likelihood rejects it decisively (2·ΔNLL = 136).
  - Same at λ 5 (334) and λ 1 (124).
- **The method has a preferred answer of λ ≈ 2-3**, set by the detector, and it is "right" only when the truth falls there. A stopped clock is right twice a day.
- **The chosen λ follows the noise level, not the number of units.** The same true λ = 2 is fitted as 1-2 at SNR 5, 2 at SNR 1.5, and 4 at SNR 0.5, with every score. This reproduces the referee's simulated-vs-estimated discrepancy (true 2, estimated 4) at low SNR.
- **The score choice matters at the margins** (e.g. λ 2 at SNR 5: NLL picks 1, QQ and WS pick 2), but no score recovers λ across the range.
- **Each k0 heatmap shows the ridge:** as λest increases, the best k0 falls in proportion. Only λ·k0 is fixed by the data.

**Addendum: "success" = true λ among the top 3 of the 5 candidates.** Rank of the true λ per dataset (5 datasets × 7 rows with true λ 1-5), from `fig13_scores.csv`:

| criterion | all scores | best score | chance (random pick of 5) |
|---|---|---|---|
| true λ ranked 1st | 24% | 37% (QQ) | 20% |
| top 2 | 43% | 46% (WS) | 40% |
| top 3 | 62% | 71% (QQ) | 60% |

- Top-3 success is 100% for true λ 2-4 and λ 2 at SNR 5. It is 0-20% for λ 5, 0-80% for λ 1 (score-dependent), and 0% for λ 2 at SNR 0.5.
- The fixed guess "λ ∈ {2, 3, 4}", which ignores the data, succeeds in 71% of these datasets.
- So the method is no better than chance at ranking λ. Its top-3 hits come from the middle candidates it always favours.

## 32. Two clarifications: the static control on the steady-threshold plot; what "set by noise" means

**32a. `steady_with_static.png`** (`checks/plot_steady_static.py`, values in `analysis_outputs/steady/steady_with_static.csv`). This is the section 23 three-panel comparison, now with the unchanging-granule control drawn next to every real value, for all six datasets:
- frame noise (real vs control);
- share of bursts that noise also produces;
- turnover beyond noise, filled = real, open = control.

**Turnover beyond noise, %/h, real vs control:**

| dataset | per-frame (R2) real | R2 control | steady (R4) real | R4 control |
|---|---|---|---|---|
| 4PP7 | 1.6 | 1.6 | 12.3 | 13.0 |
| 4PP7-FUS | 3.0 | 2.3 | 22.7 | 14.2 |
| 30PP7 | 1.5 | 3.1 | 15.8 | 10.8 |
| 30PP7-FUS | 2.2 | 0.3 | 15.0 | 14.1 |
| 5Qb | 9.3 | 3.0 | 22.4 | 14.0 |
| 10Qb | 9.7 | 1.6 | 21.3 | 13.0 |

**Reading.**
- **Noise:** the control reproduces the real noise in every dataset and both extractions (within 1-2 points), so the trace noise is measurement noise.
- **Per-frame thresholds:** the control shows ≈ 0 turnover, so the white-noise copies are a valid noise floor. PP7 real ≈ control, while Qβ is clearly above its control.
- **Steady thresholds:** the control itself shows 11-14 %/h, which is the artefact. Only 4PP7-FUS and Qβ rise clearly above their controls.

**32b. Why the "set by brightness / noise" part of the 4x gap is significant (p = 0.01).**
- That bar is the difference between the noise-only copies of FUS and of PP7: how much bigger the bursts are that pure noise produces in FUS traces.
- It is significant because the FUS traces are noisier **relative to the granule**, not because FUS adds noise:

| | integrated granule intensity | frame noise, absolute | frame noise, relative |
|---|---|---|---|
| 4PP7 | 123,800 | 19,600 | 16% |
| 4PP7-FUS | 50,000 | 10,700 | 22% |

- In absolute terms FUS traces are *less* noisy (10.7k vs 19.6k), because photon noise grows only with the square root of the signal. But the granule is 2.5× dimmer, so the same kind of noise is a larger fraction of it. Bursts are measured as a fraction of the granule, so noise alone produces fractionally larger "bursts" in FUS.
- At 30x the two proteins give similar brightness and relative noise, so this part is small and not significant (p = 0.19).

## 33. A fixed granule mask (R5): less noise, but it measures granule movement

Scripts:
- `collect_intensity.m` (`opts.fixedMask`, `opts.posWin`; default off = published behaviour unchanged);
- `checks/reextract_fixed.m` → `in_vivo_reextracted/after_R5_fixedmask/`;
- `checks/static_control.m` (new version −1 = fixed mask) → `analysis_outputs/steady/static_control_fixed.csv`;
- `checks/compare_masks.py` and `checks/fixedmask_motion.py` → `analysis_outputs/fixedmask/` (`compare.csv`, `tests.csv`, `motion.csv`, `compare_masks.png/pdf`).

**Design (R5).**
- The track position is smoothed (running median over 15 frames, fixed in advance).
- The 15×15 crops of all of a granule's frames are averaged, and multithresh(4) is run once on that average image.
- The resulting spot pixels (top two classes) and cell pixels (middle two) are used unchanged in every frame, so there is no mask drift.
- Same tracks, crop, offset and classes as R2. All six datasets keep exactly the same traces.
- In vivo only: there are no in vitro movies.

**Results against each extraction's own unchanging-granule control:**

| dataset | noise R2 / R4 / R5 (% per frame) | R5 control noise | R5 turnover beyond noise: real / control (%/h) | R5 real − control | share of light in mask, first / last 30 frames |
|---|---|---|---|---|---|
| 4PP7 | 16.9 / 7.2 / 4.1 | 6.0 | 14.2 / 4.4 | 9.7 | 0.82 / 0.81 |
| 4PP7-FUS | 22.2 / 10.6 / 6.6 | 8.4 | 15.6 / 0.7 | 14.9 | 0.76 / 0.75 |
| 30PP7 | 21.5 / 9.8 / 5.7 | 8.0 | 11.3 / −0.3 | 11.6 | 0.70 / 0.74 |
| 30PP7-FUS | 24.1 / 12.2 / 6.9 | 11.1 | 16.4 / 2.1 | 14.3 | 0.63 / 0.67 |
| 5Qb | 11.8 / 4.5 / 2.2 | 4.8 | 24.7 / 6.5 | 18.1 | 0.81 / 0.81 |
| 10Qb | 11.7 / 4.3 / 2.1 | 4.8 | 16.5 / 0.3 | 16.2 | 0.81 / 0.80 |

**Condition comparisons:**

| comparison | R2: size beyond noise p (one-sided) | R4: turnover vs control z | R5: turnover vs control z | R5: size beyond noise p (one-sided) |
|---|---|---|---|---|
| 4PP7-FUS vs 4PP7 | 0.032 | 1.4 | 0.7 | 0.32 |
| 30PP7-FUS vs 30PP7 | 0.39 | −0.8 | 0.5 | 0.13 |
| 5Qb vs 4PP7 | < 0.001 | 2.5 | 1.1 | 0.009 |
| 10Qb vs 30PP7 | 0.11 | 0.8 | 0.7 | 0.043 |

**Is R5's new signal granule exchange or granule movement?** (`fixedmask_motion.py`)
- Slow changes (13-frame smoothing, as the detector sees them) in the spot sum move together with changes in the share of the granule's light that falls inside the fixed mask.
- Median correlation: 0.76-0.90 for PP7 datasets, 0.57 for Qβ.
- Exchange changes how much light there is, but hardly its share inside a fixed mask. A granule shifting, deforming or drifting in focus relative to the mask changes the share.

**Reading (plain).**
1. The fixed mask works as intended against noise. Frame noise drops about 4× (e.g. 4PP7 16.9% → 4.1%). Its control shows no slow artefact: turnover about 0-7 %/h, and burst-free copies reproduce about 100% of the control's bursts.
2. With noise this low, every dataset shows turnover beyond noise (10-18 %/h). But the check above shows that most of that slow signal in PP7 follows the granule moving or reshaping within the mask. The unchanging-granule control is centred on its track and cannot move, so it does not contain this. **R5's "dynamics" in PP7 are therefore largely granule motion, not exchange.** In Qβ the link is weaker (0.57), consistent with Qβ having genuine exchange (sections 18-24).
3. **The 4x FUS difference does not get stronger; it disappears** (z 0.7, p 0.32). The 30x pair still does not differ. Qβ vs PP7 remains (5Qb vs 4PP7 p = 0.009).
4. The R5 control is noisier than the real data (4PP7: 6.0 vs 4.1%). The pixel-noise model, extrapolated above the calibrated intensity range (section 21), overestimates the photon noise of bright granule pixels by about 1.5×. The controls are therefore somewhat conservative. With per-frame thresholds (R2) the mask jitter dominated, which is why the noise matched there.
- **Net:** no mask design tried so far makes the 4x FUS effect stronger.
  - Per-frame thresholds: borderline (one-sided p 0.03).
  - Steady thresholds: borderline against a large artefact.
  - Fixed mask: absent, and it measures movement.
- The robust results (Qβ more dynamic than PP7; no 30x FUS effect) hold in all three extractions.
- A movement-insensitive measure would need either granule light summed without a mask (but section 16 showed crop sums are dominated by cell light), or a real static sample to calibrate the movement contribution.
