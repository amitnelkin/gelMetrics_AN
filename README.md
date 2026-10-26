# Instructions for Use
## *GelMetrics: An algorithm for analyzing the dynamics of gel-like phase separated condensates*

### This repository contains the code and data required to recreate the analyses of the paper.
>**Tested on Matlab R2025b**

### Prerequisites
- Curve Fitting Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- ImageJ FIJI with the Mosaic plugin, to detect and track bright spots in the microscopy movies.

### What the analysis reports
For each condition, `GelMetricsReport.m` gives:

| measure | meaning |
|---|---|
| granules (movies) | traces of at least 30 min; movies (fields of view) are the independent replicates |
| brightness | median granule intensity (median per movie, then over movies) |
| frame noise | frame-to-frame noise, % of the granule |
| mean / median burst size (a.u.) | intensity lost per burst; the mean corresponds to λ·k0 of the Poisson model |
| burst size, granules / copies | median burst size as % of the granule, in the granules and in their **burst-free copies** |
| burst size beyond noise | the difference between the two (% points) |
| turnover beyond noise | % of the granule exchanged in bursts per hour, minus the same in the copies |
| net change beyond bleaching | %/h; the granule fading faster (< 0) or slower (> 0) than its own background |
| first vs second 30 min | turnover and net change per half hour, with a test for a change |

Every measure has a 95% confidence interval (movies resampled, then granules within movies).
`GelMetricsCompare.m` tests two conditions by shuffling whole movies between them.

**Noise floor.** Every granule trace is compared with burst-free copies of itself: the same photobleaching trend and the same noise level, but no bursts, analyzed in the identical way. Bursts found in a copy are false. Only the part of a measure beyond its copies is attributed to the granule. The copies are random: with one copy per trace (as used for the numbers in the paper), the copy-based measures change by up to ~2 percentage points of burst size between random draws. The runners therefore use 20 copies per trace (`opts.nCopies`); set `opts.nCopies = 1` to reproduce the paper's procedure.

**λ and k0.** The original method fitted λ (units per burst) and k0 (intensity of one unit) with a weighted score (`modifiedPoissonFit.m`). Only their product, the mean burst size, can be estimated reliably, so the fit is no longer part of the analysis. It is kept in `legacy/` for reference.

### Running the analysis
- **Figures of the paper:** `PaperFigures_revised.m`, one section per figure (Figs. 5-9). Each section starts from the repository data, runs the same functions as below, and saves the figure in `figures/`. It reproduces the manuscript numbers (the noise copies use the paper's seeds); Figure 5 takes about 10 min.
- **Experimental data:** `ExperimentalAnalysis.m`. Analyzes any number of conditions (datasets of the paper or your own folders), tests each against the first one, and draws the burst size distributions (as Figs. 6-9) and a summary of mean burst size, turnover and net change (as Figs. 6 and 10).
- **Simulated data:** `SimulationRunWithAnalysis.m`. Bursty and constant (burst-free) simulated granules with known parameters, through the same code.
- **Datasets of the paper:** `gm_paper_dataset('<name>')` with names `4PP7`, `4PP7-FUS`, `30PP7`, `30PP7-FUS`, `5Qb`, `10Qb` (in vivo) and `vitro_1_1`, `vitro_10_1`, `vitro_100_1` (in vitro, protein:RNA).
  - The data are the corrected extractions in `in_vivo_reextracted/` and `in_vitro_reextracted/same_granules_R2/` (fluorescence channel only, camera offset removed).
  - `in_vivo_experiments/` and `in_vitro_experiments/` hold the original traces, for reference only.

### New data pipeline
1. Acquire time-lapse movies; 8-12 movies (fields of view) per condition, each about 360 time points, is usually sufficient. Movies are the replicates.
2. In ImageJ FIJI, run the Particle Tracker 2D/3D (Mosaic plugin) to get CSV files of granule coordinates. The macro `BatchParticleTracker.ijm` does this for a folder (`Plugins` → `Macros` → `Run...`).
   > For a detailed walkthrough of this step, see https://github.com/naorgk/slncRNA_Analysis/blob/main/Full%20pipeline%20walkthrough.pdf
3. Run `read_tables.m` in the folder of the CSV files; it writes one `<movie>_xy.mat` per movie.
4. Run `Collect_intensity_data.m` and select the folder containing the `.tif` movies and the `_xy.mat` files (one folder per condition).
   - Movies and track files are paired by name. Only the fluorescence channel is read, and the camera offset (200 a.u. for our EMCCD; set `opts.offset` in `collect_intensity.m` for your camera) is subtracted.
   - It writes `data.mat` (spot intensity), `background.mat` (surrounding intensity), `spotsum.mat` (integrated granule intensity) and `movies.mat` (the movie of every trace).
5. In `ExperimentalAnalysis.m`, load each condition with `gm_load_extraction('<folder>')` and run it.

### Main files
| file | role |
|---|---|
| `GelMetricsReport.m` | all measures and confidence intervals for one condition |
| `GelMetricsCompare.m` | movie-level test between two conditions |
| `PaperFigures_revised.m` | the manuscript figures 5-9 |
| `gm_plot_bursthist.m`, `gm_plot_bars.m`, `gm_style.m` | the figure panels (manuscript style) |
| `analysis/gm_eq19_fit.m`, `analysis/gm_shape_fit.m` | the λ and λ·k0 fits of Fig. 5 (Eq. 1.9 and the original model curve) |
| `analysis.m`, `getPosition.m`, `runLengthEncode.m` | burst detection (smoothing, runs of monotone change) |
| `collect_intensity.m`, `Collect_intensity_data.m`, `generate_movArray.m`, `read_tables.m` | extraction of the granule and background intensities from the movies |
| `generateSim.m` | simulated granules |
| `gm_load_extraction.m`, `gm_paper_dataset.m` | load your extraction / the paper's datasets |
| `legacy/` | the original λ/k0 fit and runners, for reference |
| `checks/`, `analysis/`, `analysis_outputs/` | the revision analyses behind every number of the paper (see `analysis_outputs/RESULTS.md`) |
