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
- **MATLAB figure rendering hung on the user's laptop** (batch and normal mode: `exportgraphics`, `print`, `drawnow`, `heatmap` never return). Never render in MATLAB. Save data and plot in Python (matplotlib, scipy, pandas available via `py -3`). Use `checks/gm_noplot.m` (copy of `modifiedPoissonFit.m` without its plotting section) and `checks/gm_strip_plot_lines.m` when running repo code that plots. **On the lab PC the test passes** (2026-10-18, ~40 s), so MATLAB figures may be used there. There, Python is `%LOCALAPPDATA%\Programs\Python\Python312\python.exe` (no `py` launcher), git is `C:\Program Files\Git\cmd\git.exe`, Parallel Computing Toolbox is not installed, Fiji + MOSAIC is in `..\Fiji\Fiji.app`, and bftools is in `..\bftools` (it needs `JAVA_HOME` = Fiji's bundled zulu8 JDK). Fiji must run in GUI mode (`ImageJ-win64.exe -macro`), because `--headless` cannot open .nd2 through the Bio-Formats importer. On another machine, test once with `matlab -batch "f=figure('Visible','off'); plot(1:3); exportgraphics(f, fullfile(tempdir,'t.png'))"` (kill it if it does not return within 1 min); if it works, MATLAB figures may be used there.
- Run the runners without editing them through `checks/run_runners.m` (`run_runners(stage, doExp, doSim)`). It uses the per-dataset paper settings and writes `analysis_outputs/<stage>/lambda_k0_table.csv` + `*_run.mat`. Figures: `py -3 checks/plot_runner_outputs.py <stage>`. Cumulative table: `py -3 checks/build_lambda_k0_table.py`. Checks: `bash checks/run_checks.sh <stage> "<orig|patched>" "<scripts>"`.
- Commits use `git -c user.name=amitnelkin -c user.email=amit.nelkin@gmail.com` (no global git identity is configured).

## Pipeline map
- `BatchParticleTracker.ijm`: FIJI/Mosaic tracking. It produces CSV files with columns traj, frame, x, y.
- `read_tables.m`: converts the CSVs to `*_xy.mat`.
- `Collect_intensity_data.m`: TIFF + xy → `data.mat` and `background.mat`, using a 15×15 crop and multithresh(4), i.e. 5 classes. "Spot" is the top 2 classes, "cell" the middle 2. The "dark" class is computed but never used.
- `analysis.m`: cubic-polynomial background fit, filloutliers, 13-frame moving mean, `data = (S/Bfit - 1)*Bfit(1)`, then `getPosition.m` (sign runs, threshold `-SP/log2(p)`) and a linear fit per segment, which gives the amplitudes.
- `modifiedPoissonFit.m`: grid search over λ = 1:5 and k0 using a weighted rank score (0.5 MSE, 0.25 QQ, 0.25 KS).
- `ExperimentalAnalysis.m` and `SimulationRunWithAnalysis.m`: the runners. `generateSim.m` generates the simulations.

### Extraction pipeline as actually run (from the user, 2026-10-18)
1. Each .nd2 (NIS-Elements; red fluorescence + DIC) was converted to TIFF with bftools, e.g. in a Windows cmd folder of .nd2 files:
   `for %i in (*.nd2) do ..ftoolsfconvert "%i" "%~dpni.tif"`
   The result is a 2-channel interleaved TIFF (red, DIC, red, DIC, ...), like ExampleTIF.tif.
2. Particle tracking with the repo macro `BatchParticleTracker.ijm`. Note that the macro opens the **.nd2** (Bio-Formats, hyperstack), not the TIFF, and tracks it as opened. It has no channel split or selection. The user says tracking was on the red channel. Verify how Mosaic numbered the frames: 0-360 per channel, or 0-721 over all planes.
3. The Mosaic CSVs `Traj_<name>.nd2.csv` were renamed to `<name>.csv` to match the TIFF names. The user's cmd snippet:
   ```
   for %i in (Traj_*.nd2.csv) do (
       set "fname=%~ni"
       setlocal enabledelayedexpansion
       set "newname=!fname:Traj_=!"
       set "newname=!newname:.nd2=!"
       echo "%i" -> "!newname!.csv"
       endlocal
   )
   ```
   (As written it only echoes the new names. The actual rename was presumably a `ren` in place of the `echo`.)
4. `read_tables.m` produced the `*_xy.mat` files, run unchanged. Naming convention: each dataset is its own folder of movies named `001.nd2` … `0NN.nd2` (8-15 per dataset). `name(1:3)` therefore gives unique `001_xy.mat` … names, and directory order pairs `0NN.tif` with `0NN_xy.mat` correctly. No collision or mismatch (also confirmed: no dataset contains a repeated track-length list).
5. `Collect_intensity_data.m` produced `data.mat` / `background.mat`. The user edited the output names locally per dataset run, giving the repo's `<dataset>_data.mat` / `_background.mat`. It reads the interleaved 2-channel TIFF page by page (`generate_movArray.m`), so `movArray(frame)` is a red or a DIC page depending on how Mosaic numbered frames (step 2). This is the prime suspect for the PP7 zigzag traces.

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
- **Noise floor (2026-10-19, RESULTS.md section 18): with a valid burst-free surrogate (cubic trend + white noise at the trace's own level; the real noise is white), 70-100% as many negative bursts are detected as in the real data at SP 10-15 for every condition, with similar sizes. No SP gives <= 10% false bursts with usable counts. The burst-based scores must be corrected for this or replaced before further robustness analyses.**
- **Detection threshold.** `-SP/log2(p)` assumes independent signs, but it is applied after a 13-point moving average. Calibrate SP from the false-positive rate on burst-free signals instead.
- **λ and k0 are nearly degenerate.** Only λ·k0 (about the mean amplitude) is well identified. Replace the grid search with a likelihood (below).

## State of the revision (branch `revision`)
- Done:
  - `checks/` (Appendix C) created.
  - Baseline run on the untouched code.
  - Appendix A committed ("Bug fixes from review"). It includes a crash fix: `modelCDF.m` computes the Poisson pmf in log space because `poisspdf` underflowed.
  - Appendix B committed one change per commit: zero-truncation, symmetric bursts, bleaching, k0 grid, filloutliers.
  - The appendices were removed from this file after that.
- **Extraction fixed and re-run (lab PC, 2026-10-18, user-approved), one commit each:**
  - R1 (c694551): `BatchParticleTracker.ijm` tracks only the fluorescence channel, which it saves as `<name>.tif` in `outputFolder`. `generate_movArray.m` reads only that channel. `Collect_intensity_data.m` is now a wrapper around `collect_intensity.m`, which maps page = frame+1 and column/row = x+1/y+1 (Mosaic is 0-based; spot centroid − Mosaic = +0.87 col / +0.80 row) and pairs TIFF with `_xy.mat` by name.
  - R2 (c259610): camera offset 200 subtracted.
  - R3 (e450505): 5 px circular aperture for the spot. **Reverted 2026-10-19 at the user's request** (the multithresh spot is the published method). Use R2 (fter_R2_offset) for all PP7 analyses.
  - `collect_intensity(path, struct('offset',0|200,'aperture',0|5))` reproduces each stage.
  - Re-extracted data (PP7 only) are in `in_vivo_reextracted/<after_R1_redchannel|after_R2_offset|after_R3_aperture>/`, made by `checks/reextract_pp7.m`. The repo `.mat` files are untouched.
  - Working files are in `..\zigzag_scratch\reextract\<dataset>\` (`tif` = bfconvert output, `track` = Fiji/Mosaic, `work` = hard links + renamed CSVs).
- **Qβ re-extracted (2026-10-19).** The Naor movies (`Raw_Results\Naor_Data-All_night_induction\5Qb\250ms`, 10 movies; `10Qb`, 9 movies) reproduce the stored Qβ data (trace counts, levels, burst counts within a few percent; RESULTS.md section 12), so they are the source. Re-extracted Qβ is in `in_vivo_reextracted/after_R1_redchannel|after_R2_offset/5Qb_*`, `10Qb_*`. Qβ λ is not identified in either version (σ ≈ amplitude).
- **Likelihood fitter (analysis 2) done**: `analysis/gelmetricsMLE.m` and helpers, validated (RESULTS.md section 5). On the full pipeline, λ is ~25% high and the profile CIs undercover; use the granule-bootstrap CIs. `opts.fitSigma` exists but is too slow as implemented (~8 min per fit) and gave open CIs: not used.
- `expConst` in the simulation is 0.0026 /min, the median of per-dataset bleaching fits (`checks/fit_bleaching.py`). It is not the 0.005 placeholder.
- The baseline wrappers reproduce every negative-burst number in MS.pdf (lambda, burst counts, k0 within 2). The manuscript fits negative bursts only.
- Per-stage outputs: `analysis_outputs/<stage>/`. Before/after table: `analysis_outputs/lambda_k0_before_after.md`. Narrative: `analysis_outputs/RESULTS.md`.
- Confirmed in MATLAB, with the original scoring (were Python-port estimates):
  - **Noise floor.** False bursts per flat trace were 4.64, 3.09, 1.34 and 0.93 at SP = 10, 12, 14 and 15.
  - **Bleaching asymmetry.** At expConst 0.005, mean positive vs negative amplitude was 645 vs 107 in the original simulation and 172 vs 164 in the fixed one.
  - **Granule bootstrap (t4, 40 resamples).** k0 95% interval 61-168 for 30PP7-FUS and 63-400 for 30PP7, with lambda = 1 in 35/40 resamples for both.

## How λ and k0 are estimated (two methods, both in the code)
1. **Grid / WS method (the manuscript's method).** `modifiedPoissonFit.m`, called by the runners (`ExperimentalAnalysis.m`, `checks/run_runners.m`). It produces every λ/k0 in `lambda_k0_before_after.md`.
   - **Model:** Eq. 1.9 as implemented in code. The Poisson pmf at integer n is linearly interpolated at I/k0 and divided by its maximum. It is not the gamma-function form of (I/k0)!, and it has no noise term.
   - **Search:** λ ∈ {1..5}, k0 ∈ 1:500.
   - **Score:** WS = 0.5·MSE (10-bin histogram) + 0.25·QQ + 0.25·KS, as ranks.
   - **Revision changes:**
     - bug fixes: the model is evaluated at the data's bin centres; QQ/KS use amplitudes, not bar heights; the truncation percentile is taken over the data;
     - B1: zero-truncation, P(0) = 0.
2. **Likelihood method (new, analysis 2).** `analysis/gelmetricsMLE.m`.
   - **Model:** Eq. 1.9's lattice n·k0 with zero-truncated Poisson weights, each point broadened by Gaussian noise σ (σ from the quiescent-segment amplitudes, `ampVar`). As σ → 0 it reduces to the discrete Poisson on n·k0.
   - **Fit:** continuous λ and k0 by maximum likelihood; no WS.
   - **Truncation:** the same 10% top trim and truncation percentile, modelled in the normaliser.
   - **CIs:** profile likelihood and granule bootstrap.
   - RESULTS.md sections 5, 7-9; the "grid" columns there are method 1 on the same amplitudes.
- **Open decision for the user:** which method the revised manuscript reports. Method 1 does not track λ in simulation (section 5); method 2 is misspecified in vitro (section 8).

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
- Raw TIFFs are not on the laptop, but the original **.nd2 files are on the user's lab PC**, so re-extraction is possible there. The user offered them; the work moves to the lab PC.
- ExampleTIF.tif is one of the TIFFs actually used in the original analysis (red fluorescence + DIC, 2 channels x 361 frames, interleaved red, DIC, red, ...). Channel 1 = red fluorescence (floor ~200), channel 2 = DIC (~6800 everywhere, no dark pixels). The user says all movies were exported the same way, and tracking (FIJI/Mosaic) was run on the red channel only. **Correction (2026-10-18):** the movies differ in channel count, and on 2-channel movies Mosaic tracked raw planes, DIC included (see the zigzag bullet below).
- Each `<name>_data.mat` aggregates ~8-15 movies (`*_xy.mat` files), appended movie by movie, so many traces share a movie.
- Manuscript claims to test: those raised by the reviewers. The list of difference and similarity claims is in RESULTS.md / the MS summary.
- **Zigzag cause found (lab PC, 2026-10-18; scratch work in `..\zigzag_scratch\`, outside the repo).** Supersedes the "open flag" bullet below.
  - Not all movies are 2-channel. `showinf`: 4PP7 (`Raw_Results\Genescript-PCP4x-QCP5x_K133-tdPP7-mCherry\New`, movies 005-013) has 005 = 1 channel and 006-013 = 2 channels. 4PP7-FUS (001-012) is all 1 channel. 30PP7 and 30PP7-FUS (001-008 each) are all 2 channels. Naor 5Qb, 10Qb and 4pp7 001 are 2 channels too, but Qb shows no zigzag: unexplained, check.
  - On a 2-channel movie, Mosaic (via the repo macro) walks the first 361 of the 722 raw planes (red, DIC, red, ...) and calls them frames 0-360. Mosaic frame f = plane f+1, i.e. only the first ~180 time points, with DIC on odd frames. It detects particles on DIC planes too (007: 4130 detections on odd frames vs 2771 on even).
  - `Collect_intensity_data.m` reads TIFF page f (not f+1). So even-only tracks (linked red-to-red, step 2) are read entirely from DIC pages and look clean. Odd-only tracks (DIC detections) are read from red pages. Mixed-parity tracks alternate red/DIC: the zigzag. `check_interleave.py` flags only the last group: in a 2-channel movie, no trace is a valid red-condensate trace.
  - Confirmed: the ORIGINAL f7eded0 `read_tables.m` + `Collect_intensity_data.m` on fresh tracking of 005 and 007 reproduce stored 4PP7 traces 1-19 and 49-88 bit-exactly (59/59).
- **Open data-quality flag (cause unknown; superseded above).** In 4PP7, 30PP7 and 30PP7-FUS, 21-40% of the traces with >= 180 frames jump between two levels on every point for the whole trace. Example: 4PP7 trace 21 goes 4819, 8928, 4340, 8776, ...; the jump is 50-70% of the level, vs ~5% for clean traces (`checks/check_interleave.py`, `analysis_outputs/baseline/data_quality/zigzag_example.png`). Flagged traces come in contiguous index blocks (i.e. particular movies). In vitro, Qb and 4PP7-FUS have none.
  - `generate_movArray.m` reads every TIFF page as a time point, so on a red/DIC interleaved TIFF every other page is DIC. That fits the zigzag (red ~4500 / DIC ~8500), but then why most traces are clean is unexplained.
  - One hypothesis to test, not established: Mosaic run on the 2-channel stack may number all 722 planes as frames. With `linkrange = 2`, a track can link red-to-red every 2 planes (clean trace, max 360 points) or alternate red/DIC on consecutive planes (zigzag).
  - **First thing to do on the lab PC:** re-run tracking on one PP7 movie and on ExampleTIF's movie. Look at the Mosaic CSV frame numbers (0-360 or 0-721? steps of 1 or 2?). Then run `Collect_intensity_data.m` and see which traces zigzag. Do not change the extraction pipeline without the user's OK.
- With re-extraction possible, the deferred [model] change (camera offset ~200 subtraction + fixed circular aperture, radius matched to the PSF of 5 px) can be applied, together with a reader that loads only the red channel. Each of these is a pipeline change: ask first, and report its effect on lambda/k0 like the other model changes. Offset alone, applied to the stored traces, changes amplitudes by < 1% in vivo and 3-10% in vitro (late in the movie).
- Still open: is there a fixed-sample or known-liquid control dataset?

## Next session
1. Done on the lab PC (2026-10-18): zigzag resolved (RESULTS.md section 4), PP7 re-extracted (R1-R3), likelihood fitter built and validated (section 5).
2. Done: likelihood fit on the real data (RESULTS.md sections 7-9; section 9 supersedes the PP7 parts, which now use R2 multithresh data). 4x FUS vs PP7: the amplitude distributions differ (median perm p < 0.0005; joint (λ,k0) LRT perm p 0.016) even though the separate λ and k0 CIs overlap (λ-k0 ridge). Earlier: The model fits PP7 (KS p 0.01-0.9) but is rejected in vitro (KS p ≤ 0.014: a comb of narrow peaks vs smooth heavy-tailed data), so the in vitro λ values (e.g. 100:1 = 3.47) are not trustworthy. Next model step: σ growing with n, or fitted. Qβ λ is not identified (σ ≈ 10× the other datasets). No FUS effect is detectable for 4x.
3. Next: "Analyses to produce" items 3-7 (SNR recovery grid, noise floor, sensitivity, bootstrap/time halves, window sensitivity).
4. Open questions for the user: are the Naor Qβ .nd2 files the source of the stored Qβ data (if so, re-extract them)? Is there a fixed/crosslinked control dataset?

## Working style
- Before changing a result-affecting line, say what will change and why.
- After each commit, rerun the affected checks and report how λ and k0 moved.
- If a fix makes a manuscript claim fail, say so plainly. The aim is a defensible paper, not preserving current numbers.
