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

## Files to recreate (all content is in the appendices of this file)
Nothing from the earlier chat has been added to the repo. Recreate it from the appendices:

- **Appendix A, `p1_bugs.patch`**: every [bug] fix, including the new standalone `modelCDF.m`. Apply as commit 1.
- **Appendix B, `p2_model.patch`**: every [model] change, applied on top of A. Split it into one commit per change (zero-truncation, symmetric simulated bursts, bleaching constant, symmetric k0 grid, filloutliers, camera offset plus aperture).
- **Appendix C, check scripts**: create them in a `checks/` folder. They call the repo's own `analysis.m` and `generateSim.m`, and `addpath` the repo via `repoDir`. Set `repoDir = '..'` when they live in `checks/`.

Notes on applying the patches:
- Both patches were made against the original files with line endings normalised to LF. The repo files are CRLF, so `git apply` may reject them. If it does, apply each hunk with your edit tool by matching content, and keep CRLF endings in the result.
- `modelCDF.m` must be a standalone file on the path, not a local function. An earlier version made it local, and MATLAB reported "Unrecognized function modelCDF".
- `gm_mpfit.m` deliberately reproduces the ORIGINAL scoring (bar-height QQ/KS, unaligned model grid, bin-centre truncation), so the checks give a baseline. After commit 1, add a switch so the checks can also score with the patched `modifiedPoissonFit.m`, and report both.
- The scripts have not yet been run in MATLAB. Fix any errors you hit.

## Suggested order for the first session
1. Create a branch. Create `checks/` from Appendix C.
2. Run `t1`–`t4` on the untouched code and save the outputs as the baseline in `analysis_outputs/baseline/`.
3. Apply Appendix A (commit 1). Rerun the checks and both runners, and tabulate how λ and k0 moved.
4. Apply Appendix B one change per commit, rerunning after each.
5. Then start on the analyses below, beginning with the likelihood fitter.
6. Once Appendices A–C exist as real files and commits, delete the three appendices from this file (keep the rest). They are only needed once, and this file is loaded at the start of every session.

## Preliminary findings (from a Python port, so indicative only; confirm in MATLAB)
- **Noise floor.** On flat, burst-free signals, false bursts per trace were 4.8, 3.2, 1.4 and 0.9 at SP = 10, 12, 14 and 15. Real λ=2 simulations gave 4.7 per trace at SP = 15.
- **Degeneracy.** At true λ=2 and k0=60, recovered (λ, k0) pairs ranged from (1, 111) to (3, 49) across seeds.
- **Bleaching asymmetry.** With the unbleached-baseline simulation at expConst = 0.005, symmetric bursts gave mean positive amplitude 727 vs mean negative 111. With the fixed simulation: 165 vs 155.
- **Sample sizes.** In 4PP7, 67 of 360 traces pass the 180-frame filter, giving 24–66 positive bursts depending on SP.
- **Real-data sweep** (SP × truncation × top trim). k0 for 4PP7 ranged from 74 to 500. λ stayed mostly at 1 with the 10% top trim and flipped to 3–5 without it.
- **Granule bootstrap at paper settings.** For both 30PP7 and 30PP7-FUS, λ = 1 in 36/40 resamples. k0 95% intervals were 42–190 and 102–402, overlapping. Check whether the manuscript claims these systems differ.

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

## Open questions to confirm with the user before relying on results
- Real frame interval (assumed 10 s).
- Are Mosaic frames 0-based in their version? Do CSV columns 4 and 5 map to image columns and rows as `Collect_intensity_data.m` assumes?
- Camera dark offset value, and an aperture radius matching the PSF.
- Are the raw TIFFs available for re-extraction?
- Which differences between conditions does the manuscript claim (so the bootstrap can test exactly those)?
- Is there a fixed-sample or known-liquid control dataset?

## Working style
- Before changing a result-affecting line, say what will change and why.
- After each commit, rerun the affected checks and report how λ and k0 moved.
- If a fix makes a manuscript claim fail, say so plainly. The aim is a defensible paper, not preserving current numbers.

---

# Appendix A: `p1_bugs.patch` (bug fixes, commit 1)

```diff
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/Collect_intensity_data.m bug/Collect_intensity_data.m
--- orig/Collect_intensity_data.m	2026-10-18 08:12:00.110205338 +0000
+++ bug/Collect_intensity_data.m	2026-10-18 08:58:09.925189850 +0000
@@ -11,13 +11,13 @@
 clc;
 
 % Path to folder containing tif and mat files
-file_path = [uigetdir("*", ''), '\'];
+file_path = [uigetdir("*", ''), filesep];
 
 % load all .tiff and all .mat files
 tiff_files = dir([file_path,'*.tif']);
-mat_files = dir([file_path,'*.mat']);
+mat_files = dir([file_path,'*_xy.mat']); % only read_tables outputs (was *.mat: picked up data.mat/background.mat)
 
-% Length of the subframe to
+% Half-width of the subframe: the crop is (2*7+1) = 15 x 15 pixels
 sub_frame_length = 7;
 
 % Initialize variables
@@ -32,25 +32,28 @@
     xdata = xydata.x;
     ydata = xydata.y;
     framedata = xydata.frame;
+    allFrames = cell2mat(framedata(:));
+    frameOffset = double(min(allFrames) == 0); % Mosaic frame numbers start at 0; MATLAB indices at 1
     
     % Start iterating through spots and gather intensity
     for t = 1:length(xdata)
         spot_frames = cell2mat(framedata(t)); 
-        spot_frames = spot_frames(2:end); % remove frame 0 to avoid error
-        spot_x = cell2mat(xdata(t)); spot_x = spot_x(2:end);
-        spot_y = cell2mat(ydata(t)); spot_y = spot_y(2:end);
+        spot_frames = spot_frames + frameOffset; % convert to 1-based (was: dropped first point, kept 0-based index)
+        spot_x = cell2mat(xdata(t));
+        spot_y = cell2mat(ydata(t));
         if length(spot_x)>=24 % continue if spot exists for >24 frames
             % iterate over frames and get intensity data
-            spot_intensity = zeros(length(spot_x),1);
-            background_intensity = zeros(length(spot_x),1);
+            frameIdx = spot_frames - spot_frames(1) + 1;  % position of each detection on the full frame axis
+            spot_intensity = nan(frameIdx(end),1);       % NaN where the tracker skipped frames
+            background_intensity = nan(frameIdx(end),1);
             
             for n = 1:length(spot_x)
                frame = movArray(spot_frames(n)).data; % complete current frame
                % Extract sub frame around the spot coordiantes
                sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
-                                   min(spot_y(n)+sub_frame_length,length(frame)),...
+                                   min(spot_y(n)+sub_frame_length,size(frame,1)),...
                                    max(spot_x(n)-sub_frame_length,1),...
-                                   min(spot_x(n)+sub_frame_length,length(frame))];
+                                   min(spot_x(n)+sub_frame_length,size(frame,2))];
                sub_frame_bounds = round(sub_frame_bounds);
                sub_frame = frame(sub_frame_bounds(1):sub_frame_bounds(2),sub_frame_bounds(3):sub_frame_bounds(4));
                % sub frame is ready from this point
@@ -66,11 +69,13 @@
                
                % frame partition to spot,background is done
                
-               spot_intensity(n,1) = mean(sub_frame(spot_frame==1));
-               background_intensity(n,1) = mean(sub_frame(cell_frame==1));
+               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame==1));
+               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1));
                
             end
                       
+           spot_intensity = fillmissing(spot_intensity,'linear');           % fill skipped frames
+           background_intensity = fillmissing(background_intensity,'linear');
            intensity = [intensity,{spot_intensity}];
            background = [background,{background_intensity}];
         end
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/ExperimentalAnalysis.m bug/ExperimentalAnalysis.m
--- orig/ExperimentalAnalysis.m	2026-10-18 08:12:00.112175745 +0000
+++ bug/ExperimentalAnalysis.m	2026-10-18 08:58:09.918311721 +0000
@@ -5,7 +5,8 @@
 % Algorithm Parameters
 movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
 stringencyParameter = 14;                 % Parameter for burst events detection
-timeVecTotal = linspace(0, 60, 360);      % Time steps vector
+frameInterval_min = 10/60;                % frame interval in minutes - set to your acquisition interval
+timeVecTotal = (0:359) * frameInterval_min; % Time steps vector (was linspace(0,60,360): dt = 60/359 min)
 
 % Load signal and background data
 % StringencyParameter (SP) and truncation percentile (GOF) for replicating the exact plots: 
@@ -124,13 +125,18 @@
 
 lambda = 1:5;
 K0Vector = 1:500; % Experimental data is noisy by nature, using standard errors will lead to unreasonably large grids
+truncPercent = 15; % left-truncation percentile, passed explicitly to modifiedPoissonFit
 
 % Poisson analysis on positive amplitudes
-[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes');
+[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes',truncPercent);
 
 % Final ranking of lambdas
 for ll = 1:length(lambda)
     [min_w_vec{ll}, min_w_k0{ll}] = min(w_score(ll,:));
+    min_w_k0{ll} = K0Vector(min_w_k0{ll}); % convert index to k0 value (only equal when K0Vector = 1:N)
+end
+if any(cell2mat(min_w_k0) == K0Vector(1) | cell2mat(min_w_k0) == K0Vector(end))
+    warning('Best k0 is at the edge of K0Vector - extend the grid');
 end
 
 [min_w_scores, min_w_idx] = sort(cell2mat(min_w_vec));
@@ -143,11 +149,15 @@
 end
 
 % Poisson analysis on negative amplitudes
-[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes');
+[~, ~, ~, ~, ~, w_score] = modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes',truncPercent);
 
 % Final ranking of lambdas
 for ll = 1:length(lambda)
     [min_w_vec{ll}, min_w_k0{ll}] = min(w_score(ll,:));
+    min_w_k0{ll} = K0Vector(min_w_k0{ll}); % convert index to k0 value (only equal when K0Vector = 1:N)
+end
+if any(cell2mat(min_w_k0) == K0Vector(1) | cell2mat(min_w_k0) == K0Vector(end))
+    warning('Best k0 is at the edge of K0Vector - extend the grid');
 end
 
 [min_w_scores, min_w_idx] = sort(cell2mat(min_w_vec));
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/SimulationRunWithAnalysis.m bug/SimulationRunWithAnalysis.m
--- orig/SimulationRunWithAnalysis.m	2026-10-18 08:12:00.113830129 +0000
+++ bug/SimulationRunWithAnalysis.m	2026-10-18 08:58:09.919202952 +0000
@@ -1,5 +1,6 @@
 close all;   % Close all open figure windows
 clear;       % Clear workspace
+rng(1);      % fixed seed so simulation figures are reproducible
 
 % Algorithm Parameters
 movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
@@ -100,7 +101,7 @@
     khat_SEhat_coeff = sqrt(2)/sqrt(length(ampPos)) ;
     khat_SEhat = khat * khat_SEhat_coeff;
     K0Vector = floor(khat - 8 * khat_SEhat):ceil(khat + 18 * khat_SEhat); % The k0 lookup grid is integers spanning >20 standard errors from the k0 estimator, the grid is not symmetrical since we noticed bias towards higher lambdas during simulations
-    modifiedPoissonFit(ampPos,lambda,K0Vector,10,'Positive amplitudes');
+    modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes',15); % same trim and truncation as ExperimentalAnalysis.m
     
     
     % Poisson analysis on negative amplitudes
@@ -108,7 +109,7 @@
     khat_SEhat_coeff = sqrt(2)/sqrt(length(ampNeg)) ;
     khat_SEhat = khat * khat_SEhat_coeff;
     K0Vector = floor(khat - 8 * khat_SEhat):ceil(khat + 18 * khat_SEhat); % The k0 lookup grid is integers spanning >20 standard errors from the k0 estimator, the grid is not symmetrical since we noticed bias towards higher lambdas during simulations
-    modifiedPoissonFit(abs(ampNeg),lambda,K0Vector,10,'Negative amplitudes');
+    modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes',15); % same trim and truncation as ExperimentalAnalysis.m
 
 end
 
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/analysis.m bug/analysis.m
--- orig/analysis.m	2026-10-18 08:12:00.116152865 +0000
+++ bug/analysis.m	2026-10-18 08:58:09.920757248 +0000
@@ -35,7 +35,7 @@
         currentBackground = currentBackground';
     end
     
-    % Discard experimental signals with less than 15 minutes of data
+    % Discard experimental signals with less than 180 frames (30 min at 10 s/frame)
     if length(currentSignal) < 180
         continue;
     end
@@ -95,7 +95,7 @@
         end
         if xVar(end)>length(data)
             xVar = xVar(1:end-1);
-            pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1))-1,1);
+            pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-2),1);
         else
             pVar{tt} = polyfit(xVar,data(startPosValid(tt)+lengthsValid(tt)+1:startPosValid(tt+1)-1),1);
         end
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/generateSim.m bug/generateSim.m
--- orig/generateSim.m	2026-10-18 08:12:00.117762784 +0000
+++ bug/generateSim.m	2026-10-18 08:58:09.920281811 +0000
@@ -9,11 +9,11 @@
 
 function [signals, backgroundSignals] = generateSim(numTracks, trackLength, type, expConst, noiseStd, lambda)
     % Initialize the time vector
-    t = linspace(0, 60, 360)';
+    t = linspace(0, 60, trackLength)';
     
     backgroundSignals = randi([1000,1500],numTracks,1) .* ones(numTracks,trackLength); % Base background signal
     backgroundSignals = backgroundSignals .* exp(-expConst*t'); % Add exponential component
-    backgroundSignals = backgroundSignals + noiseStd*randn(numTracks,360); % Add Gaussian noise component
+    backgroundSignals = backgroundSignals + noiseStd*randn(numTracks,trackLength); % Add Gaussian noise component
 
     % Initialize the signals matrix
     signals = zeros(numTracks, trackLength);
@@ -41,7 +41,7 @@
                 currentLevel = 0;
                 for j = 1:trackLength
                     distFromBurst = abs(j-currentBurstPos);
-                    if rand < 0.025 && distFromBurst > 18 % Assuming a 1% chance of a burst at each point
+                    if rand < 0.025 && distFromBurst > 18 % 2.5% chance of a burst at each point
                         currentBurstPos = j; 
                         change = k * poissrnd(lambda);
                         % We noticed a positive bias in our simulated signals
@@ -54,7 +54,7 @@
                     end
                     signals(i, j) = currentLevel;
                 end
-                signals(i, :) = randi([2000,10000]) + signals(i, :) .* exp(-expConst * t') + noiseStd * randn(1, trackLength);
+                signals(i, :) = (randi([2000,10000]) + signals(i, :)) .* exp(-expConst * t') + noiseStd * randn(1, trackLength); % whole signal bleaches, matching analysis.m
                 
             end
             
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/modelCDF.m bug/modelCDF.m
--- orig/modelCDF.m	1970-01-01 00:00:00.000000000 +0000
+++ bug/modelCDF.m	2026-10-18 08:58:09.916438966 +0000
@@ -0,0 +1,13 @@
+function [a_grid, F_model] = modelCDF(lambda_val, k0, threshold, a_max)
+% MODELCDF  CDF of the fitted model in intensity units, on (threshold, a_max].
+% Same piecewise-linear interpolated Poisson shape as the curve used in
+% modifiedPoissonFit.m.
+a_grid = linspace(threshold, a_max, 2000);
+xx = 0:ceil(a_max/k0)+1;
+yy = poisspdf(xx, lambda_val);
+dens = interp1(xx, yy, a_grid./k0);
+F_model = cumtrapz(a_grid, dens);
+F_model = F_model ./ F_model(end);
+[F_model, iu] = unique(F_model);          % interp1 needs strictly increasing values
+a_grid = a_grid(iu);
+end
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/modifiedPoissonFit.m bug/modifiedPoissonFit.m
--- orig/modifiedPoissonFit.m	2026-10-18 08:12:00.124841925 +0000
+++ bug/modifiedPoissonFit.m	2026-10-18 08:58:09.916304904 +0000
@@ -1,7 +1,7 @@
 % Function modifiedPoissonFit fits the amplitude data to modified Poisson
 % functions, calcualtes the mean squared error between the observed data
 % and the theoretical functions, and plots the results
-function [pvec, MSE, QQ_error, KS_h, KS_dist, w_score] = modifiedPoissonFit(data,lambda,K0Vector,numBins,title_str)
+function [pvec, MSE, QQ_error, KS_h, KS_dist, w_score] = modifiedPoissonFit(data,lambda,K0Vector,numBins,title_str,percent)
 
 % Present the data as 10 bins
 [binCounts,binEdges] = histcounts(data,numBins,'Normalization','pdf');
@@ -12,8 +12,13 @@
 end
 nn_full = binCounts;
 hist_x_axis = abs(mean([binEdges(1:end-1);binEdges(2:end)]));
-percent = 15;
-GOF_threshold = prctile(hist_x_axis, percent); % The pecrentile as a threshold value for goodness-of-fit analysis (we are interested more about the right side of distribution)
+if nargin < 6 || isempty(percent), percent = 15; end % truncation percentile is now an input
+if percent > 0
+    GOF_threshold = prctile(data, percent); % percentile of the AMPLITUDES (was: of the bin centres)
+else
+    GOF_threshold = 0;
+end
+data_thres = data(data > GOF_threshold); % truncated amplitude samples, used for QQ and KS
 % If we want no data to be truncated, GOF_threshold & percent should be set to 0
 % GOF_threshold = 0;percent = 0;
 
@@ -38,7 +43,7 @@
         % Generate a theoretical Poisson distribution using the normalized
         % x axis and the current lambda from the lambda vector
         yy = poisspdf(xx,lambda(j));
-        yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
+        yq = interp1(xx,yy./max(yy),hist_x_axis./normalizer); % model evaluated at the data's own bin centres
         
         % Calcualte goodness of fit between the truncated data and the theoretical function
         yq_thres = yq(hist_x_axis > GOF_threshold);
@@ -53,12 +58,15 @@
         end
         
         MSE(j,i) = goodnessOfFit(yq_thres(~isnan(yq_thres)),binCounts_thres(~isnan(yq_thres))./max(binCounts_thres),'MSE');
-        n = numel(nn_full_thres);
+        % QQ and KS on the amplitude samples vs the model CDF (was: on bar heights)
+        [a_grid, F_model] = modelCDF(lambda(j), normalizer, GOF_threshold, max(data));
+        n = numel(data_thres);
         p = ((1:n) - 0.5) / n;
-        q1 = quantile(nn_full_thres, p);
-        q2 = quantile(yq_thres*max(binCounts_thres)./max(yq_thres), p);
+        q1 = sort(data_thres(:))';                       % sample quantiles
+        q2 = interp1(F_model, a_grid, p, 'linear', 'extrap'); % model quantiles
         QQ_error(j,i) = sqrt(mean((q2 - q1).^2));
-        [~,~,KS_dist(j,i)] = kstest2(yq_thres*max(binCounts_thres)./max(yq_thres), nn_full_thres, 0.01); % Two-sample Kolmogorov-Smirnov test
+        F_emp = (1:n) / n;                                % empirical CDF at sorted samples
+        KS_dist(j,i) = max(abs(F_emp - interp1(a_grid, F_model, q1, 'linear', 1)));
     end
 end
 
@@ -102,13 +110,13 @@
     normalizer = x_axis(1,K_idx(i));
     xx = 0:max(hist_x_axis)/normalizer+1;
     yy = poisspdf(xx,lambda_single);
-    yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
+    yq = interp1(xx,yy./max(yy),hist_x_axis./normalizer); % model at the data's bin centres
     yq_thres = yq(hist_x_axis > GOF_threshold);
     binCounts_thres = binCounts((hist_x_axis > GOF_threshold));
     nn_full_thres = nn_full((hist_x_axis > GOF_threshold));
 
     if all(~yq_thres) || all(~nn_full_thres) || all(~binCounts_thres)
-        warning(['Truncation unsuccesful for lambda=',num2str(lambda(j)),' and k0=',num2str(normalizer),' - working on full data'])
+        warning(['Truncation unsuccesful for lambda=',num2str(lambda(i)),' and k0=',num2str(normalizer),' - working on full data'])
         yq_thres = yq;
         binCounts_thres = binCounts;
         nn_full_thres = nn_full;
@@ -124,14 +132,16 @@
     legend('show');
 
     subplot(1,5,4);
-    qqplot(yq_thres*max(binCounts_thres)./max(yq_thres),nn_full_thres); box on;
+    [a_grid, F_model] = modelCDF(lambda_single, normalizer, GOF_threshold, max(data));
+    p_plot = ((1:numel(data_thres)) - 0.5) / numel(data_thres);
+    plot(interp1(F_model, a_grid, p_plot, 'linear', 'extrap'), sort(data_thres(:))', 'o'); refline(1,0); box on;
     title('Truncated QQ plot');
     ylabel('Quantiles of sample data');
     xlabel('Quantiles of theoretical distribution');
 
     subplot(1,5,5);
-    cdfplot(nn_full_thres); hold on;
-    cdfplot(yq_thres*max(binCounts_thres)./max(yq_thres));
+    cdfplot(data_thres); hold on;
+    plot(a_grid, F_model, 'LineWidth', 1.5);
     title('Truncated CDF plot');
     legend('Experiment','Estimate');
 end
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/plotTracksCustom.m bug/plotTracksCustom.m
--- orig/plotTracksCustom.m	2026-10-18 08:12:00.126527758 +0000
+++ bug/plotTracksCustom.m	2026-10-18 08:58:09.921835171 +0000
@@ -16,9 +16,9 @@
         lengthsValid(tt)+1);
     y = polyval(p{tt},x);
     if y(end)>=y(1)
-        plot(x*(timeVec(2)-timeVec(1)),y,'g','LineWidth',1.2);
+        plot((x-1)*(timeVec(2)-timeVec(1)),y,'g','LineWidth',1.2);
     else
-        plot(x*(timeVec(2)-timeVec(1)),y,'r','LineWidth',1.2);
+        plot((x-1)*(timeVec(2)-timeVec(1)),y,'r','LineWidth',1.2);
     end
 end
 % Plot blue linear line overlayed on quiescent segments
@@ -29,7 +29,7 @@
         continue;
     end
     yVar = polyval(pVar{tt},xVar);
-    plot(xVar*(timeVec(2)-timeVec(1)),yVar,'b','LineWidth',1.2);   
+    plot((xVar-1)*(timeVec(2)-timeVec(1)),yVar,'b','LineWidth',1.2);   
 end
 
 return
diff -ruN '--exclude=.git' '--exclude=*.mat' orig/read_tables.m bug/read_tables.m
--- orig/read_tables.m	2026-10-18 08:12:00.128024520 +0000
+++ bug/read_tables.m	2026-10-18 08:58:09.921965137 +0000
@@ -23,5 +23,5 @@
     x = x(2:end);
     y = y(2:end);
     frame = frame(2:end);
-    save([files(f).name(1:3),'_xy.mat'],'x','y','frame');
+    save([files(f).name(1:end-4),'_xy.mat'],'x','y','frame'); % full name (was first 3 characters: files could overwrite each other)
 end
\ No newline at end of file

```

# Appendix B: `p2_model.patch` (model changes, on top of A)

```diff
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/Collect_intensity_data.m full/Collect_intensity_data.m
--- bug/Collect_intensity_data.m	2026-10-18 08:58:09.925189850 +0000
+++ full/Collect_intensity_data.m	2026-10-18 08:58:09.959422105 +0000
@@ -19,6 +19,8 @@
 
 % Half-width of the subframe: the crop is (2*7+1) = 15 x 15 pixels
 sub_frame_length = 7;
+aperture_radius = 3;  % fixed aperture for spot intensity, in pixels - set to cover the PSF
+camera_offset = 0;    % camera dark offset in A.U. - measure from a dark frame and set here
 
 % Initialize variables
 intensity = cell(1,1);
@@ -48,7 +50,7 @@
             background_intensity = nan(frameIdx(end),1);
             
             for n = 1:length(spot_x)
-               frame = movArray(spot_frames(n)).data; % complete current frame
+               frame = double(movArray(spot_frames(n)).data) - camera_offset; % complete current frame, offset removed
                % Extract sub frame around the spot coordiantes
                sub_frame_bounds = [max(spot_y(n)-sub_frame_length,1),...
                                    min(spot_y(n)+sub_frame_length,size(frame,1)),...
@@ -69,8 +71,13 @@
                
                % frame partition to spot,background is done
                
-               spot_intensity(frameIdx(n),1) = mean(sub_frame(spot_frame==1));
-               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1));
+               % Fixed circular aperture (was: mean of thresholded pixels, not additive in molecule number)
+               [cc, rr] = meshgrid(1:size(sub_frame,2), 1:size(sub_frame,1));
+               cy = spot_y(n) - sub_frame_bounds(1) + 1;
+               cx = spot_x(n) - sub_frame_bounds(3) + 1;
+               aperture = (rr-cy).^2 + (cc-cx).^2 <= aperture_radius^2;
+               spot_intensity(frameIdx(n),1) = mean(sub_frame(aperture));
+               background_intensity(frameIdx(n),1) = mean(sub_frame(cell_frame==1 & ~aperture));
                
             end
                       
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/SimulationRunWithAnalysis.m full/SimulationRunWithAnalysis.m
--- bug/SimulationRunWithAnalysis.m	2026-10-18 08:58:09.919202952 +0000
+++ full/SimulationRunWithAnalysis.m	2026-10-18 08:58:09.958320035 +0000
@@ -12,7 +12,7 @@
 numTracks = 1000;                                % Number of tracks to simulate
 t = linspace(0, 60, 360);                        % Time vector for simulation
 noiseStd = 40;                                   % Standard deviation of noise
-expConst = 0.0005;                               % Exponential constant for photobleaching
+expConst = 0.005;                                % Exponential constant for photobleaching (was 0.0005 = 3% over 60 min; set from your backgrounds)
 lambda_sim = 2;                                  % Poisson parameter
 
 % generateSim function generates simulated trajectories
@@ -100,7 +100,7 @@
     khat = var(abs(ampPos))/mean(abs(ampPos));
     khat_SEhat_coeff = sqrt(2)/sqrt(length(ampPos)) ;
     khat_SEhat = khat * khat_SEhat_coeff;
-    K0Vector = floor(khat - 8 * khat_SEhat):ceil(khat + 18 * khat_SEhat); % The k0 lookup grid is integers spanning >20 standard errors from the k0 estimator, the grid is not symmetrical since we noticed bias towards higher lambdas during simulations
+    K0Vector = max(1, floor(khat - 13 * khat_SEhat)):ceil(khat + 13 * khat_SEhat); % Symmetric k0 grid, +-13 standard errors around the k0 estimator (was -8/+18, tuned to a bias seen in simulation)
     modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes',15); % same trim and truncation as ExperimentalAnalysis.m
     
     
@@ -108,7 +108,7 @@
     khat = var(abs(ampNeg))/mean(abs(ampNeg));
     khat_SEhat_coeff = sqrt(2)/sqrt(length(ampNeg)) ;
     khat_SEhat = khat * khat_SEhat_coeff;
-    K0Vector = floor(khat - 8 * khat_SEhat):ceil(khat + 18 * khat_SEhat); % The k0 lookup grid is integers spanning >20 standard errors from the k0 estimator, the grid is not symmetrical since we noticed bias towards higher lambdas during simulations
+    K0Vector = max(1, floor(khat - 13 * khat_SEhat)):ceil(khat + 13 * khat_SEhat); % Symmetric k0 grid, +-13 standard errors around the k0 estimator (was -8/+18, tuned to a bias seen in simulation)
     modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes',15); % same trim and truncation as ExperimentalAnalysis.m
 
 end
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/analysis.m full/analysis.m
--- bug/analysis.m	2026-10-18 08:58:09.920757248 +0000
+++ full/analysis.m	2026-10-18 08:58:09.958585512 +0000
@@ -50,7 +50,7 @@
     [f,gof,output] = fit(timeVec',currentBackground','poly3');
     currentBackground = feval(f,timeVec);
     % Deal with outliers (mostly results of experimental noise)
-    currentSignal = filloutliers(currentSignal,'nearest','mean');
+    currentSignal = filloutliers(currentSignal,'nearest','movmedian',3*movingAverageWindowSpan); % local outliers only (global mean could clip real plateaus)
     % Apply moving average filter on the signal
     currentSignal = smoothdata(currentSignal, 'movmean', movingAverageWindowSpan);
     % Remove background and denoise from photobleaching component
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/generateSim.m full/generateSim.m
--- bug/generateSim.m	2026-10-18 08:58:09.920281811 +0000
+++ full/generateSim.m	2026-10-18 08:58:09.958468280 +0000
@@ -47,8 +47,8 @@
                         % We noticed a positive bias in our simulated signals
                         % we counter this by giving a higher chance for a
                         % negative burst which reduces signal intensity
-                        if rand < 0.6 
-                            change = -change; % 60% chance of being negative
+                        if rand < 0.5
+                            change = -change; % symmetric: 50% chance of being negative
                         end
                         currentLevel = currentLevel + change;
                     end
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/modelCDF.m full/modelCDF.m
--- bug/modelCDF.m	2026-10-18 08:58:09.916438966 +0000
+++ full/modelCDF.m	2026-10-18 08:58:09.957815996 +0000
@@ -5,6 +5,7 @@
 a_grid = linspace(threshold, a_max, 2000);
 xx = 0:ceil(a_max/k0)+1;
 yy = poisspdf(xx, lambda_val);
+yy(1) = 0;                                % zero-truncated: no 0-unit bursts
 dens = interp1(xx, yy, a_grid./k0);
 F_model = cumtrapz(a_grid, dens);
 F_model = F_model ./ F_model(end);
diff -ruN '--exclude=.git' '--exclude=*.mat' bug/modifiedPoissonFit.m full/modifiedPoissonFit.m
--- bug/modifiedPoissonFit.m	2026-10-18 08:58:09.916304904 +0000
+++ full/modifiedPoissonFit.m	2026-10-18 08:58:09.957602511 +0000
@@ -43,6 +43,7 @@
         % Generate a theoretical Poisson distribution using the normalized
         % x axis and the current lambda from the lambda vector
         yy = poisspdf(xx,lambda(j));
+        yy(1) = 0; % zero-truncated Poisson: a burst of 0 units cannot be detected
         yq = interp1(xx,yy./max(yy),hist_x_axis./normalizer); % model evaluated at the data's own bin centres
         
         % Calcualte goodness of fit between the truncated data and the theoretical function
@@ -110,6 +111,7 @@
     normalizer = x_axis(1,K_idx(i));
     xx = 0:max(hist_x_axis)/normalizer+1;
     yy = poisspdf(xx,lambda_single);
+    yy(1) = 0; % zero-truncated Poisson, as in the fitting loop
     yq = interp1(xx,yy./max(yy),hist_x_axis./normalizer); % model at the data's bin centres
     yq_thres = yq(hist_x_axis > GOF_threshold);
     binCounts_thres = binCounts((hist_x_axis > GOF_threshold));

```

# Appendix C: check scripts (create in `checks/`)

## `checks/gm_mpfit.m`

```matlab
function [lamOrder, kBest, scoreSorted, w_score] = gm_mpfit(data, lambda, K0Vector, numBins, percent)
% GM_MPFIT  Same scoring math as modifiedPoissonFit.m, with three changes
% so it can be used inside parameter sweeps:
%   1. the truncation percentile is an argument (it is hardcoded to 15 in
%      the original function),
%   2. nothing is plotted,
%   3. the lambdas are returned ranked by their best weighted score.
%
% Outputs
%   lamOrder    - lambda values sorted from best to worst score
%   kBest       - best k0 for each lambda, in the same order
%   scoreSorted - best weighted score for each lambda, in the same order
%   w_score     - full weighted-score matrix (lambda x k0)

if nargin < 4 || isempty(numBins), numBins = 10; end
if nargin < 5 || isempty(percent), percent = 15; end
lambda = lambda(:)';
K0Vector = K0Vector(:)';
data = data(:)';

[binCounts, binEdges] = histcounts(data, numBins, 'Normalization', 'pdf');
if binEdges(1) < 0
    binEdges = fliplr(binEdges);
    binCounts = fliplr(binCounts);
end
nn_full = binCounts;
hist_x_axis = abs(mean([binEdges(1:end-1); binEdges(2:end)]));

% As in the original: percentile of the BIN CENTRES, not of the data.
% percent = 0 means no truncation (the original comment says to set both
% GOF_threshold and percent to 0 for that).
if percent > 0
    GOF_threshold = prctile(hist_x_axis, percent);
else
    GOF_threshold = 0;
end

MSE      = zeros(numel(lambda), numel(K0Vector));
QQ_error = zeros(numel(lambda), numel(K0Vector));
KS_dist  = zeros(numel(lambda), numel(K0Vector));

w_MSE = 0.5; w_QQ = 0.25; w_KS = 0.25;

for j = 1:numel(lambda)
    for i = 1:numel(K0Vector)
        normalizer = K0Vector(i);
        xx = 0:max(hist_x_axis)/normalizer+1;
        yy = poisspdf(xx, lambda(j));
        yq = interp1(xx, yy./max(yy), linspace(0, max(hist_x_axis)/normalizer, length(binCounts)));

        sel = hist_x_axis > GOF_threshold;
        yq_thres = yq(sel);
        binCounts_thres = binCounts(sel);
        nn_full_thres = nn_full(sel);

        if all(~yq_thres) || all(~nn_full_thres) || all(~binCounts_thres)
            yq_thres = yq;
            binCounts_thres = binCounts;
            nn_full_thres = nn_full;
        end

        % goodnessOfFit(x, xref, 'MSE') = mean((xref - x).^2); written
        % out here so the System Identification Toolbox is not needed
        ok = ~isnan(yq_thres);
        xref = binCounts_thres(ok) ./ max(binCounts_thres);
        MSE(j,i) = mean((xref - yq_thres(ok)).^2);

        n = numel(nn_full_thres);
        p = ((1:n) - 0.5) / n;
        q1 = quantile(nn_full_thres, p);
        q2 = quantile(yq_thres*max(binCounts_thres)./max(yq_thres), p);
        QQ_error(j,i) = sqrt(mean((q2 - q1).^2));
        [~, ~, KS_dist(j,i)] = kstest2(yq_thres*max(binCounts_thres)./max(yq_thres), nn_full_thres, 0.01);
    end
end

M = numel(MSE);
rankKS  = reshape((tiedrank(KS_dist(:))  - 1) / (M - 1), size(KS_dist));
rankMSE = reshape((tiedrank(MSE(:))      - 1) / (M - 1), size(MSE));
rankQQ  = reshape((tiedrank(QQ_error(:)) - 1) / (M - 1), size(QQ_error));
w_score = w_QQ*rankQQ + w_KS*rankKS + w_MSE*rankMSE;

[bestPerLam, kIdx] = min(w_score, [], 2);
[scoreSorted, ord] = sort(bestPerLam);
lamOrder = lambda(ord);
kBest = K0Vector(kIdx(ord));
end
```

## `checks/gm_generateSim.m`

```matlab
function [signals, backgroundSignals, t] = gm_generateSim(numTracks, type, expConst, noiseStd, lambda, fixBleach, negFrac, k0, refractory, pBurst)
% GM_GENERATESIM  generateSim.m with switches for the issues being tested.
%
%   fixBleach  - false (default): baseline NOT bleached, as in generateSim.m
%                true: whole signal bleached, (base + bursts).*exp(-k t),
%                which matches the model assumed in analysis.m
%   negFrac    - fraction of negative bursts (generateSim.m uses 0.6)
%   k0         - intensity per unit (generateSim.m uses 60)
%   refractory - minimum frames between bursts (generateSim.m uses 18)
%   pBurst     - burst probability per frame (generateSim.m uses 0.025)
%
% type: 1 = constant, 3 = bursts (type 2 is in the original generateSim.m)

if nargin < 6  || isempty(fixBleach),  fixBleach  = false; end
if nargin < 7  || isempty(negFrac),    negFrac    = 0.6;   end
if nargin < 8  || isempty(k0),         k0         = 60;    end
if nargin < 9  || isempty(refractory), refractory = 18;    end
if nargin < 10 || isempty(pBurst),     pBurst     = 0.025; end

T = 360;
t = linspace(0, 60, T);

backgroundSignals = randi([1000,1500], numTracks, 1) .* ones(numTracks, T);
backgroundSignals = backgroundSignals .* exp(-expConst*t);
backgroundSignals = backgroundSignals + noiseStd*randn(numTracks, T);

signals = zeros(numTracks, T);
switch type
    case 1
        for i = 1:numTracks
            signals(i,:) = randi([2000,10000]) .* exp(-expConst*t) + noiseStd*randn(1, T);
        end
    case 3
        for i = 1:numTracks
            pos = -10; level = 0; trace = zeros(1, T);
            for j = 1:T
                if rand < pBurst && abs(j - pos) > refractory
                    pos = j;
                    change = k0 * poissrnd(lambda);
                    if rand < negFrac
                        change = -change;
                    end
                    level = level + change;
                end
                trace(j) = level;
            end
            base = randi([2000,10000]);
            if fixBleach
                signals(i,:) = (base + trace) .* exp(-expConst*t) + noiseStd*randn(1, T);
            else
                signals(i,:) = base + trace .* exp(-expConst*t) + noiseStd*randn(1, T);
            end
        end
    otherwise
        error('type must be 1 or 3');
end
end
```

## `checks/gm_loadDataset.m`

```matlab
function [signals, backgrounds] = gm_loadDataset(repoDir, folder, name, minLen)
% GM_LOADDATASET  Loads a <name>_data.mat / <name>_background.mat pair the
% same way ExperimentalAnalysis.m does (drops the empty first cell).
% If minLen > 0, keeps only traces with at least minLen frames.

if nargin < 4, minLen = 0; end
S = struct2cell(load(fullfile(repoDir, folder, [name '_data.mat'])));
B = struct2cell(load(fullfile(repoDir, folder, [name '_background.mat'])));
S = S{1}; B = B{1};
S = S(2:end); B = B(2:end);
S = S(:); B = B(:);
if minLen > 0
    keep = cellfun(@numel, S) >= minLen;
    S = S(keep); B = B(keep);
end
signals = S; backgrounds = B;
end
```

## `checks/t1_lambda_recovery.m`

```matlab
% T1  Lambda/k0 recovery on simulated data, original vs corrected simulation.
% Uses the repo's own analysis.m; scoring via gm_mpfit (same math as
% modifiedPoissonFit.m, no plots). k0 grid built as in
% SimulationRunWithAnalysis.m (khat - 8 SE : khat + 18 SE).
% True k0 = 60. Watch the lam*k0 column: if it is stable while lam and k0
% swing, the two parameters are not separately identified.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');

nTracks = 300; seeds = 0:2; lamSims = 1:5; lamGrid = 1:5;
win = 13; SP = 15; truncPct = 15; tVec = linspace(0, 60, 360);

rows = {};
for fixBleach = [false true]
    negFrac = 0.6; if fixBleach, negFrac = 0.5; end
    for lamTrue = lamSims
        for s = seeds
            rng(s);
            [S, B] = gm_generateSim(nTracks, 3, 0.0005, 40, lamTrue, fixBleach, negFrac);
            [~, ~, ap, an] = analysis(S, B, tVec, win, SP, 0, 0);
            ap = ap(2:end); an = abs(an(2:end));
            [lp, kp] = gm_mpfit(ap, lamGrid, k0grid(ap), 10, truncPct);
            [ln, kn] = gm_mpfit(an, lamGrid, k0grid(an), 10, truncPct);
            rows(end+1, :) = {fixBleach, lamTrue, s, numel(ap), lp(1), kp(1), lp(1)*kp(1), ...
                              numel(an), ln(1), kn(1), ln(1)*kn(1)}; %#ok<SAGROW>
            fprintf('fix=%d lamTrue=%d seed=%d | pos: lam=%d k0=%d | neg: lam=%d k0=%d\n', ...
                fixBleach, lamTrue, s, lp(1), kp(1), ln(1), kn(1));
        end
    end
end
warning(ws);

T = cell2table(rows, 'VariableNames', {'fixedSim','lamTrue','seed', ...
    'nPos','lamPos','k0Pos','lamK0Pos','nNeg','lamNeg','k0Neg','lamK0Neg'});
disp(T);
writetable(T, 't1_lambda_recovery.csv');

function K = k0grid(a)
    a = abs(a(:));
    kh = var(a) / mean(a);
    se = kh * sqrt(2 / numel(a));
    K = floor(kh - 8*se):ceil(kh + 18*se);
    K = K(K > 0);
end
```

## `checks/t2_false_positives_and_bias.m`

```matlab
% T2  (a) Noise floor: bursts detected in flat, burst-free signals
%         (original generateSim.m, type 1) for several stringency values.
%     (b) Bleaching asymmetry: symmetric bursts (50/50), original vs
%         corrected simulation, at low and realistic bleaching.
%
% Note: analysis.m's eventsPerSignal skips traces with zero events, so
% events per trace here is sum(eventsPerSignal) / nTracks.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');

nTracks = 300; win = 13; tVec = linspace(0, 60, 360);

% (a) Noise floor
rng(1);
[S, B] = generateSim(nTracks, 360, 1, 0.0005, 40, 2);   % original, type 1
fprintf('\n(a) Flat signals, no bursts (noise sd = 40):\n');
for SP = [10 12 14 15]
    [~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, SP, 0, 0);
    ap = ap(2:end); an = an(2:end);
    fprintf('  SP=%2d  events/trace=%.2f  pos=%d  neg=%d  median|amp|=%.1f\n', ...
        SP, sum(ev)/nTracks, numel(ap), numel(an), median(abs([ap; an])));
end

% Reference: real bursts, lambda = 2
rng(1);
[S, B] = generateSim(nTracks, 360, 3, 0.0005, 40, 2);
[~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, 15, 0, 0);
fprintf('  Reference, bursty lambda=2, SP=15: events/trace=%.2f  median pos amp=%.1f\n', ...
    sum(ev)/nTracks, median(ap(2:end)));

% (b) Bleaching asymmetry with symmetric bursts
fprintf('\n(b) Symmetric bursts (50/50), SP = 15:\n');
for fixBleach = [false true]
    for expConst = [0.0005 0.005]
        rng(2);
        [S, B] = gm_generateSim(nTracks, 3, expConst, 40, 2, fixBleach, 0.5);
        [~, ~, ap, an] = analysis(S, B, tVec, win, 15, 0, 0);
        ap = ap(2:end); an = an(2:end);
        fprintf('  fixedSim=%d expConst=%-6g (%.0f%% bleach/60 min): nPos=%d nNeg=%d meanPos=%.0f meanNeg=%.0f\n', ...
            fixBleach, expConst, 100*(1-exp(-expConst*60)), numel(ap), numel(an), mean(ap), mean(abs(an)));
    end
end
warning(ws);
```

## `checks/t3_parameter_sweep.m`

```matlab
% T3  Sensitivity of the top-ranked lambda and its k0 to the stringency
% parameter (SP), the left-truncation percentile (GOF) and the top trim,
% on the experimental data in the repo. k0 grid 1:500 as in
% ExperimentalAnalysis.m. k0 = 1 or 500 means the fit hit a grid edge.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');

datasets = { ...   % name, folder, SP and GOF listed in ExperimentalAnalysis.m
    '4PP7',            'in_vivo_experiments',  10, 15;
    '30PP7-FUS',       'in_vivo_experiments',  12, 10;
    'protein_RNA_1_1', 'in_vitro_experiments', 12, 15};
% Add more rows as needed, e.g. '30PP7', '4PP7-FUS', '5Qb', '10Qb',
% 'protein_RNA_10_1', 'protein_RNA_100_1'.

SPs = [10 12 14]; GOFs = [0 10 15 25]; topTrims = [0 10 20];   % % removed from top
lamGrid = 1:5; K0Vector = 1:500; win = 13; tVec = linspace(0, 60, 360);

rows = {};
for d = 1:size(datasets, 1)
    [S, B] = gm_loadDataset(repoDir, datasets{d,2}, datasets{d,1});
    nLong = sum(cellfun(@numel, S) >= 180);
    fprintf('\n== %s: %d traces, %d with >= 180 frames (repo settings SP=%d, GOF=%d)\n', ...
        datasets{d,1}, numel(S), nLong, datasets{d,3}, datasets{d,4});
    for SP = SPs
        [~, ~, ap] = analysis(S, B, tVec, win, SP, 0, 0);
        ap = ap(2:end);
        fprintf('  SP=%d  nPos=%d |', SP, numel(ap));
        for g = GOFs
            for tt = topTrims
                a = ap;
                if tt > 0, a = rmoutliers(ap, "percentiles", [0, 100 - tt]); end
                [lo, kb] = gm_mpfit(a, lamGrid, K0Vector, 10, g);
                fprintf(' g%d/t%d:lam%d,k%d', g, tt, lo(1), kb(1));
                rows(end+1, :) = {datasets{d,1}, SP, numel(ap), g, tt, lo(1), kb(1)}; %#ok<SAGROW>
            end
        end
        fprintf('\n');
    end
end
warning(ws);

T = cell2table(rows, 'VariableNames', {'dataset','SP','nPos','GOFpct','topTrimPct','lamBest','k0Best'});
writetable(T, 't3_parameter_sweep.csv');
```

## `checks/t4_granule_bootstrap.m`

```matlab
% T4  Granule-level bootstrap at the repo's settings: resample whole
% traces (not bursts) with replacement, rerun detection and fitting, and
% report how often each lambda wins and the spread of k0.
% Bursts from one trace are not independent, so resampling bursts would
% understate the uncertainty.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir);
ws = warning('off', 'all');

datasets = { ...   % name, folder, SP, GOF
    '30PP7-FUS', 'in_vivo_experiments', 12, 10;
    '30PP7',     'in_vivo_experiments', 12, 10};
nBoot = 40;   % increase to 200+ for the paper
lamGrid = 1:5; K0Vector = 1:500; win = 13; tVec = linspace(0, 60, 360);

rng(0);
for d = 1:size(datasets, 1)
    [S, B] = gm_loadDataset(repoDir, datasets{d,2}, datasets{d,1}, 180);
    SP = datasets{d,3}; G = datasets{d,4};
    lamB = zeros(nBoot, 1); kB = zeros(nBoot, 1);
    for b = 1:nBoot
        idx = randi(numel(S), numel(S), 1);
        [~, ~, ap] = analysis(S(idx), B(idx), tVec, win, SP, 0, 0);
        ap = ap(2:end);
        a = rmoutliers(ap, "percentiles", [0, 90]);
        [lo, kb] = gm_mpfit(a, lamGrid, K0Vector, 10, G);
        lamB(b) = lo(1); kB(b) = kb(1);
    end
    fprintf('\n%s (%d traces, SP=%d, GOF=%d, top 10%% trimmed, %d resamples)\n', ...
        datasets{d,1}, numel(S), SP, G, nBoot);
    for l = lamGrid
        fprintf('  lambda=%d wins %d/%d\n', l, sum(lamB == l), nBoot);
    end
    fprintf('  k0 2.5 / 50 / 97.5 percentiles: %s\n', mat2str(round(prctile(kB, [2.5 50 97.5]))));
end
warning(ws);
```
