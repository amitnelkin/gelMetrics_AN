# Instructions for Use
## *GelMetrics: An algorithm for analyzing the dynamics of gel-like phase separated condensates*

### This repository contains all the code and data required to recreate the analyses of the paper.
>**Tested on Matlab R2025b**

### Prerequisites
- Curve Fitting Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- In addition, ImageJ FIJI (specifically, the Mosaic plugin) was used to read microscopy image files, detect bright spots, and convert to the formats required for the scripts.

### What the analysis reports
Every granule trace is compared with burst-free copies of itself: the same photobleaching trend and the same noise level, but no bursts, analyzed in the identical way. Bursts found in a copy are false, so only the part of a measure beyond its copies is attributed to the granule (the noise floor).
For each condition, the analysis reports:
- the number of granules (traces of at least 30 min) and of movies;
- the brightness and the frame-to-frame noise of the granules;
- the mean and median burst size (a.u.); the mean burst size corresponds to the product λ·k0 of the modified Poisson model;
- the median burst size as a fraction of the granule, in the granules and in their burst-free copies, and the difference between the two (burst size beyond noise);
- the turnover beyond noise: the fraction of the granule lost in bursts per hour, minus the same in the copies;
- the net change beyond photobleaching: how much faster (< 0) or slower (> 0) the granule fades than its own background.

Every measure is given with a 95% confidence interval from resampling the replicates (granules, or movies when `moviesAsReplicates = true`). All formulas are in `summarizeGranules.m`.

### Modes
- To run on simulated data, run the script: `SimulationRunWithAnalysis.m`
  You can switch between simulation types using the `simType` variable defined within the script (1: constant, 2: increasing, 3: bursts). The script simulates granules with a known number of units per burst (λ) and intensity of one unit (k0) at different noise levels, analyzes them with the same code as the experimental data, and compares every measure with its true value.

- To run on experimental data, run the script: `ExperimentalAnalysis.m`
  Choose the data folder (`dataFolder`) and the datasets (`conditions`) at the top of the script. Each condition is reported on its own: burst plots, every measure with its confidence interval, the burst size distributions of the granules and of their copies, and a summary of mean burst size, turnover and net change.

- To compare conditions, run the script: `compareConditions.m`
  List pairs of conditions in `pairs`. For each pair, the script prints every measure of both conditions and the p value from randomly reassigning granules (or movies) between them.

### Data
Each dataset folder contains four `.mat` files: `data.mat` (granule intensity), `background.mat` (surrounding cell or dilute phase), `spotsum.mat` (integrated granule intensity) and `movies.mat` (the movie of every trace).

| folder | protein | RNA |
|---|---|---|
| `in_vivo_experiments/PCP_4x` | PCP | QCP-5x/PCP-4x |
| `in_vivo_experiments/PCP-FUS_4x` | PCP-FUS | QCP-5x/PCP-4x |
| `in_vivo_experiments/PCP_30x` | PCP | PCP-30x |
| `in_vivo_experiments/PCP-FUS_30x` | PCP-FUS | PCP-30x |
| `in_vivo_experiments/QCP_5x` | QCP | QCP-5x/PCP-4x |
| `in_vivo_experiments/QCP_10x` | QCP | QCP-10x |
| `in_vitro_experiments/protein_RNA_1_1` | PCP, 1:1 protein:RNA | PCP-14x/MCP-15x |
| `in_vitro_experiments/protein_RNA_10_1` | PCP, 10:1 protein:RNA | PCP-14x/MCP-15x |
| `in_vitro_experiments/protein_RNA_100_1` | PCP, 100:1 protein:RNA | PCP-14x/MCP-15x |

PCP: PP7 coat protein; QCP: Qβ coat protein; PCP-FUS: PCP fused to the disordered region of FUS.

### New data pipeline
  1. Obtain microscopy timelapse images. From our experience, 8-12 images (each containing 360 timepoints) is mostly sufficient.
  2. In ImageJ FIJI, use the Particle Tracker 2D/3D (included in Mosaic plugin) to generate CSV files of spatiotemporal granule coordinates. We provide a macro file (`BatchParticleTracker.ijm`) written to be run through ImageJ FIJI: `Plugins` &#8594; `Macros` &#8594; `Run...`
     Multi-channel movies (e.g., fluorescence + DIC) are reduced to the fluorescence channel before tracking.
      > For a detailed explanation about this part of the analysis, see our previous work in: https://github.com/naorgk/slncRNA_Analysis/blob/main/Full%20pipeline%20walkthrough.pdf

  3. Run script `read_tables.m` on CSV files generated in step 2. The script should be in the same folder as CSV files. The script expects the specific structure of the Mosaic Particle Tracker, and saves a `<movie>_xy.mat` file for every CSV file.
  4. Run script `Collect_intensity_data.m` on the microscopy images and the `_xy.mat` files generated in step 3 (each `.tif` is paired with the `_xy.mat` file of the same name). The camera offset (`offset`, 200 for our camera) is subtracted from all intensities. The script will generate the four `.mat` files described under Data.
  5. Run script `ExperimentalAnalysis.m` (or `compareConditions.m`) on the folder generated in step 4.

### Parameters
- `movingAverageWindowSpan` (13 time points) and `stringencyParameter` (12) are used for all experimental datasets of the paper. Use the same values for all conditions that are compared, and compare the number of detected bursts with that found in the burst-free copies.
- `nCopies` (20) is the number of burst-free copies per granule; `seed` sets the random copies.
- `expConst` in `SimulationRunWithAnalysis.m` is the photobleaching rate of our setup and should be tuned to yours.
