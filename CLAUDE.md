# CLAUDE.md — GelMetrics revision

## Context
This repo (`gelMetrics`, MATLAB) is the code for the manuscript *GelMetrics: an algorithm for analyzing the dynamics of gel-like phase-separated condensates* (Nelkin et al., ref. VM10081E). It is under revision.

The method:
1. Tracks fluorescent condensates over time.
2. Detects "bursts" (runs of monotone change in the smoothed, background-normalised signal).
3. Fits the burst-amplitude histogram with a modified Poisson model to get two parameters:
   - **k0**: intensity of one exchanged unit (one RNA plus its bound fluorescent proteins)
   - **λ**: mean number of units exchanged per burst

Referee 1 is positive with minor points. Referee 2 recommended rejection with resubmission. Referee 2's main concerns:
- no error bars
- arbitrary thresholds that look cherry-picked
- λ restricted to integers, and a simulated-vs-estimated λ discrepancy (2 vs 4)
- no noise floor
- no replicate or variance analysis
- no comparison with other methods
- interpretability of k0 and λ

The goal of the code work is to fix real bugs, then produce the analyses that answer these points. The response letter itself is written elsewhere. Your job is code, analyses, figures and numbers.

## Environment and conventions
- MATLAB R2025b. Toolboxes: Curve Fitting (`fit` in `analysis.m`), Statistics and Machine Learning, Image Processing. `goodnessOfFit` in `modifiedPoissonFit.m` comes from the System Identification Toolbox. If it is missing, replace it with `mean((xref - x).^2)`.
- Run headless with `matlab -batch "scriptName"`.
- Files use CRLF line endings. Keep them.
- Never modify the `.mat` files in `in_vitro_experiments/` or `in_vivo_experiments/`.
- Work on a branch, not `main`. Commit bug fixes and model changes separately (see below), so the effect of each on λ and k0 can be reported to the reviewer.
- Seed every simulation with `rng(...)`.
- No figures pop up in batch runs: save figures to files and write tables to CSV.
- **MATLAB figure rendering hangs on this machine** (batch and normal mode: `exportgraphics`, `print`, `drawnow`, `heatmap` never return). Never render in MATLAB. Save data and plot in Python (matplotlib, scipy, pandas available via `py -3`). Use `checks/gm_noplot.m` (copy of `modifiedPoissonFit.m` without its plotting section) and `checks/gm_strip_plot_lines.m` when running repo code that plots.
- Run the runners without editing them through `checks/run_runners.m` (`run_runners(stage, doExp, doSim)`). It uses the per-dataset paper settings and writes `analysis_outputs/<stage>/lambda_k0_table.csv` + `*_run.mat`. Figures: `py -3 checks/plot_runner_outputs.py <stage>`. Cumulative table: `py -3 checks/build_lambda_k0_table.py`. Checks: `bash checks/run_checks.sh <stage> "<orig|patched>" "<scripts>"`.
- Commits use `git -c user.name=amitnelkin -c user.email=amit.nelkin@gmail.com` (no global git identity is configured).

## Pipeline map
- `BatchParticleTracker.ijm`: FIJI/Mosaic tracking. It produces CSV files with columns traj, frame, x, y.
- `read_tables.m`: converts the CSVs to `*_xy.mat`.
- `Collect_intensity_data.m`: TIFF + xy → `data.mat` and `background.mat`, using a 15×15 crop and multithresh(4), i.e. 5 classes. "Spot" is the top 2 classes, "cell" the middle 2. The "dark" class is computed but never used.
- `analysis.m`: cubic-polynomial background fit, filloutliers, 13-frame moving mean, `data = (S/Bfit - 1)*Bfit(1)`, then `getPosition.m` (sign runs, threshold `-SP/log2(p)`) and a linear fit per segment, which gives the amplitudes.
- `modifiedPoissonFit.m`: grid search over λ = 1:5 and k0 using a weighted rank score (0.5 MSE, 0.25 QQ, 0.25 KS).
- `ExperimentalAnalysis.m` and `SimulationRunWithAnalysis.m`: the runners. `generateSim.m` generates the simulations.

The repo's own comment (`ExperimentalAnalysis.m` line 12) lists per-dataset settings: SP = 10, 12 or 14 and GOF = 10 or 15. Treat this as a known reviewer issue to resolve with one shared setting plus a sensitivity analysis, not as something to preserve.

## Known issues
Line numbers refer to the original files.

### [bug] Corrections (commit 1)
- **`modifiedPoissonFit.m`**:
  - L58–61: QQ and KS compare histogram bar heights, not amplitudes. They are invariant to bin order, and they carry 50% of the score. Fix: compare the sorted amplitude samples with the model CDF in intensity units.
  - L41 and L105: the model is evaluated on `linspace(0, max/k0, nBins)`, but the data sit at bin centres starting at `min + w/2`. Fix: `interp1(xx, yy./max(yy), hist_x_axis./normalizer)`.
  - L15–16: the "15th percentile" is taken over bin centres, not data, and it is hardcoded. Fix: make `percent` an optional 6th argument and use `prctile(data, percent)`.
  - L111: `lambda(j)` should be `lambda(i)`.
- **`ExperimentalAnalysis.m`**:
  - L133 and L150: `min` returns an index, not a k0 value. Map it with `K0Vector(idx)`, and warn when the best k0 is at a grid edge. On real data the fit hits both 1 and 500.
  - L8: `linspace(0, 60, 360)` gives dt = 60/359 min. Use the real frame interval (confirm with the user; probably 10 s).
- **`SimulationRunWithAnalysis.m`**: no `rng` seed. The simulation also skips the experimental pipeline's top-10% trim and truncation. Make both runners use identical fitting.
- **`generateSim.m`**:
  - L57: the baseline is not bleached, but `analysis.m` divides by a bleached background. Use `(base + bursts).*exp(-k t)`.
  - L12 and L16: hardcoded 360 instead of `trackLength`.
  - L44: the comment says 1% but the code uses 2.5%.
- **`analysis.m`**:
  - L98: the `-1` is outside the index. The branch is unreachable, but fix it anyway: the range should end at `startPosValid(tt+1)-2`.
  - L38: the comment says 15 min, but 180 frames is about 30 min.
- **`plotTracksCustom.m`** L19, L21 and L32: segments are drawn one frame late. Use `(x-1)*dt`.
- **`read_tables.m`** L26: the output name uses `name(1:3)`, so files can collide. Use `name(1:end-4)`.
- **`Collect_intensity_data.m`**:
  - L18: `*.mat` also picks up earlier outputs. Use `*_xy.mat`.
  - L39–41: drops the first point of every track and indexes frames 0-based. Mosaic frames start at 0, so add 1. Verify this.
  - L44–45 and L69–70: frames skipped by the tracker (`linkrange = 2`) compress time. Place points on the full frame axis and fill gaps with `fillmissing`.
  - L51 and L53: `length(frame)` should be `size(frame,1)` and `size(frame,2)`.
  - L14: replace the `'\'` path separator with `filesep`.

### [model] Method changes (commit 2+, one at a time, report the effect of each)
- Zero-truncated Poisson: `yy(1) = 0`. Observed bursts have N ≥ 1, and the max-normalisation absorbs the 1/(1-P0) factor. It also removes the λ=1 flat top, where P(0) = P(1).
- Symmetric simulated bursts: 0.5 instead of the 0.6 "positive bias" hack (`generateSim.m` L50).
- Realistic bleaching in the simulation: `expConst` 0.0005 (3% over 60 min) is far too low. Fit it from the real backgrounds.
- Symmetric k0 grid in the simulation runner (±13 SE instead of −8/+18).
- `filloutliers` with `'movmedian'` instead of the global `'mean'`.
- `Collect_intensity_data.m`: subtract the camera dark offset, and use a fixed circular aperture mean instead of the mean of thresholded pixels. The current measure is not additive in molecule number, which undermines k0·N. This needs re-extraction from the TIFFs; ask whether they are available.

### Not fixable by a line edit
- **Detection threshold.** `-SP/log2(p)` assumes independent signs, but it is applied after a 13-point moving average. Calibrate SP from the false-positive rate on burst-free signals instead.
- **λ and k0 are nearly degenerate.** Only λ·k0 (about the mean amplitude) is well identified. Replace the grid search with a likelihood (below).

## State of the revision (branch `revision`)
- Done:
  - `checks/` (Appendix C) created.
  - Baseline run on the untouched code.
  - Appendix A committed ("Bug fixes from review"). It includes a crash fix: `modelCDF.m` computes the Poisson pmf in log space because `poisspdf` underflowed.
  - Appendix B committed one change per commit: zero-truncation, symmetric bursts, bleaching, k0 grid, filloutliers.
  - The appendices were removed from this file after that.
- **Not applied:** the camera-offset + fixed-aperture change in `Collect_intensity_data.m`. The aperture part needs re-extraction and the raw TIFFs are not available. The offset part could be applied to the stored traces (traces are pixel means, so mean(px - c) = mean(px) - c). Ask the user first, since it changes the data pipeline.
- `expConst` in the simulation is 0.0026 /min, the median of per-dataset bleaching fits (`checks/fit_bleaching.py`). It is not the 0.005 placeholder.
- The baseline wrappers reproduce every negative-burst number in MS.pdf (lambda, burst counts, k0 within 2). The manuscript fits negative bursts only.
- Per-stage outputs: `analysis_outputs/<stage>/`. Before/after table: `analysis_outputs/lambda_k0_before_after.md`. Narrative: `analysis_outputs/RESULTS.md`.
- Confirmed in MATLAB, with the original scoring (were Python-port estimates):
  - **Noise floor.** False bursts per flat trace were 4.64, 3.09, 1.34 and 0.93 at SP = 10, 12, 14 and 15.
  - **Bleaching asymmetry.** At expConst 0.005, mean positive vs negative amplitude was 645 vs 107 in the original simulation and 172 vs 164 in the fixed one.
  - **Granule bootstrap (t4, 40 resamples).** k0 95% interval 61-168 for 30PP7-FUS and 63-400 for 30PP7, with lambda = 1 in 35/40 resamples for both.

## Analyses to produce (priority order) and the reviewer point each answers
1. **Apply the [bug] fixes and rerun everything** (both runners, all datasets). Make a before/after table of λ and k0 per dataset. (General)
2. **Likelihood fitter**, as a new function, e.g. `gelmetricsMLE.m`, alongside the old fit. Model:
   - amplitudes are a zero-truncated Poisson mixture of Gaussians: A ~ Σ_{n≥1} P(n | λ, n≥1) · N(n·k0, σ²), with σ from the quiescent segments (`ampVar`);
   - correct for detection truncation at A_min by dividing by P(A > A_min);
   - fit continuous λ and k0 with `fminsearch` or `fmincon` on the negative log-likelihood;
   - get confidence intervals from the profile likelihood and from a granule bootstrap;
   - plot the λ–k0 likelihood surface to show their correlation.

   (R2 4a, 4b, R1 point 2: a weak prior on λ is the "biased estimator")
3. **Simulation validation figure.** Use the fixed simulation and exactly the experimental code path. Make a recovery grid over λ_true ∈ {0.5, 1, 1.5, 2, 3, 4, 5} × SNR (k0/σ), with intervals. Also test closely spaced bursts (refractory period below the smoothing window). (R2 4b)
4. **Noise floor.** Run burst-free simulation types 1 and 2: false bursts per trace and their amplitudes versus SP. Choose one SP from a target false-discovery rate. If the user has a fixed-cell or crosslinked dataset, run it too. (R2 6)
5. **Sensitivity table.** All datasets with one shared parameter set, plus SP × truncation × top-trim × number-of-bins sweeps. Report whether the *ranking* of conditions survives. (R2 4a)
6. **Granule-level bootstrap.** Resample traces, not bursts, ideally nested in FOV or replicate; 200+ resamples. Also compare λ from the first and second 30 minutes, for time-uniformity and aging. Report traces, bursts, FOVs and replicates per condition. (R2 5)
7. **Moving-average window sensitivity** (window 7–21 frames), and state which event durations it suppresses. (R2 2b)

Output each analysis as CSV tables and saved figures in an `analysis_outputs/` folder, with a short `RESULTS.md` summarising the numbers per reviewer point.

## Facts confirmed with the user / from the data
- Frame interval 10 s. Movies are 361 frames: ExampleTIF.tif is 2 channels x 361 frames, ImageJ hyperstack, interleaved. Stored traces are at most 360 long.
- Mosaic frames are 0-based. CSV columns 4/5 are x = image column, y = image row (user). There is no CSV in the repo to verify against a TIFF.
- PSF / spot radius 5 px (BatchParticleTracker.ijm).
- Camera: Andor iXon Ultra EMCCD. Dark offset is about 200: the ExampleTIF channel-1 floor is 193-200 per frame, and no pixel is below 193. A closed-shutter frame would pin it exactly.
- Raw TIFFs are not available for re-extraction (only ExampleTIF.tif).
- Manuscript claims to test: those raised by the reviewers. The list of difference and similarity claims is in RESULTS.md / the MS summary.
- **Open data-quality flag.** In 4PP7, 30PP7 and 30PP7-FUS, 21-40% of the traces with at least 180 frames alternate frame to frame (sign-flip fraction about 1 instead of about 0.67; `checks/check_interleave.py`). This is consistent with an interleaved 2-channel stack being read page-by-page by `generate_movArray.m`. The other datasets are clean. Unresolved; ask the user.
- Still open: is there a fixed-sample or known-liquid control dataset?

## Working style
- Before changing a result-affecting line, say what will change and why.
- After each commit, rerun the affected checks and report how λ and k0 moved.
- If a fix makes a manuscript claim fail, say so plainly. The aim is a defensible paper, not preserving current numbers.
