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
