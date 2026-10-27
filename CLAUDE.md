# CLAUDE.md — GelMetrics revision (handoff, 2026-10-25)

This file lets a new chat continue the revision work. The full narrative with every number is `analysis_outputs/RESULTS.md` (sections 1-60). Read the relevant section there before redoing anything.

**Done on the laptop (2026-10-26, section 55): `CLAUDE_laptop_spotsum.md` (finished brief). The in vitro spot-sum results are in the MS, Table I, letter and control file (tracked, backups `* - backup - 261026 before spot sum.docx`). New: `checks/vitro_area_test.py` (area growth 21 / 49 / 85 %/h, 100:1 vs 1:1 p = 0.013, log trend p = 0.031; Flory-Huggins now refers to it); Fig. 14(b) uses the spot sum (`ms_figures.py --vitro-spotsum`).**

**Done on the laptop (2026-10-26, section 59): new Fig. 5 in the MS, letter and control file (tracked; backups `* - backup - 261026 before fig5.docx`); old Fig. 5(f) is now Fig. 15(c).**

**Done on the lab PC (2026-10-26): new Fig. 5 drawn (section 58); drawings of the plain-style scripts checked (section 56); MATLAB figure check (section 52; `figures/Fig5-9`), in vitro spot sums and comparison (section 53).**

**In vitro into the manuscript: done (change list v2, section 44; then the spot sum, section 55). Closed 2026-10-26 by the user: the three dropped in vitro movies had a technical problem with the images and are not mentioned in the MS (the stored granules are the paper's own movies); item A7 (FUS brightness) is already in the MS with both numbers (median over granules 3036 vs 5644, over movies 1865 vs 5694); Fig. 9(c) keeps the per-granule medians. `CLAUDE_laptop_vitro.md` is the finished brief.**

## Context
- Repo: https://github.com/amitnelkin/gelMetrics_AN, branch **`revision`** (never commit to `main`). Code for the manuscript *GelMetrics: an algorithm for analyzing the dynamics of gel-like phase-separated condensates* (Nelkin et al., ref. VM10081E), under revision. `MS.pdf` is in the repo root (not committed; ask the user if missing).
- **The method:**
  1. tracks fluorescent condensates (granules) over 60 min (361 frames, 10 s);
  2. detects "bursts" (runs of monotone change in the smoothed, bleaching-corrected signal);
  3. fits the burst-amplitude histogram with a modified Poisson, Eq. (1.9): p(I) = (1/k0)·λ^(I/k0)·e^−λ/(I/k0)!, giving k0 (one exchanged unit) and λ (units per burst).
- **Referees:** Referee 1 is positive with minor points. Referee 2 recommended rejection with resubmission. Referee 2's concerns:
  - no error bars;
  - arbitrary thresholds that look cherry-picked;
  - integer λ, and a simulated-vs-estimated λ discrepancy (2 vs 4);
  - no noise floor;
  - no replicate analysis;
  - no comparison with other methods;
  - interpretability of k0 and λ.
- **My job:** code, analyses, figures and numbers. The response letter is written elsewhere.

## The user
- Biologist (the biology is "their bread and butter"), basic statistics.
- Wants **plain, concise, short English with no math jargon**: tables, then a few bullets.
- Wants honesty: if a result contradicts a manuscript claim, say so plainly with numbers, and don't tune parameters to rescue a result.

## Rules (from the user; still in force)
- Never commit to `main`. Never modify the `.mat` files in `in_vitro_experiments/` or `in_vivo_experiments/`.
- Ask before deleting files, installing anything, changing the data-extraction pipeline, or running anything expected to take more than ~1 hour.
- Before changing any result-affecting line, say what changes and why.
- Seed every simulation (`rng(...)`), and use the same code path for simulated and experimental data.
- Bootstrap over traces (granules), never over bursts. Prefer nested movie → granule resampling: movies are the replicates.
- If CLAUDE.md is wrong, say so and update it.
- **Commits and line endings:**
  - commit with `git -c user.name=amitnelkin -c user.email=amit.nelkin@gmail.com`;
  - push with `git -c http.sslBackend=schannel push`;
  - keep CRLF line endings. Git Bash `sed -i` and the Write tool produce LF; fix with `sed -i 's/\r*$/\r/' file`, or write via Python with `newline=''`.
- After each piece of work: commit, push, add a RESULTS.md section, and give a short plain report.

## Environment
- **Lab PC** (where all this ran):
  - MATLAB R2025b `"/c/Program Files/MATLAB/R2025b/bin/matlab.exe" -batch "..."` (figures work there, but everything is plotted in Python anyway);
  - Python `%LOCALAPPDATA%\Programs\Python\Python312\python.exe` (no `py` launcher), with numpy/scipy/pandas/matplotlib/PyMuPDF/python-pptx;
  - git `C:\Program Files\Git\cmd\git.exe`;
  - PowerPoint installed (render slides via COM for checks);
  - Fiji+MOSAIC `..\Fiji\Fiji.app` (GUI mode only: `ImageJ-win64.exe -macro`); bftools `..\bftools` (needs `JAVA_HOME` = Fiji's zulu8 JDK).
  - The raw `.nd2` movies are only here (`..\Raw_Results\...`; Qβ in `Raw_Results\Naor_Data-All_night_induction\5Qb\250ms`, `10Qb`). Working copies are in `..\zigzag_scratch\reextract\<dataset>\{track,work}`.
- **Laptop:**
  - MATLAB figure rendering hangs there: never render in MATLAB; use `checks/gm_noplot.m` / `gm_strip_plot_lines.m` when running repo code that plots.
  - Python via `py -3`.
- **Gotchas:**
  - PowerShell `2>&1` on native exes breaks `$?`.
  - A Python string containing a literal `\n` written via a heredoc can turn into a real newline: use `chr(10)` in generated code.
  - pandas: a column named `compare` clashes with `DataFrame.compare`; use `df['compare']`.
  - `PYTHONIOENCODING=utf-8` when printing λ/β.

## What works WITHOUT the .nd2 files (any machine with the repo)
Everything downstream of extraction:
- the re-extracted traces are committed in `in_vivo_reextracted/`;
- the in vitro re-extracted traces are in `in_vitro_reextracted/` (the stored originals are in `in_vitro_experiments/`);
- all CSV outputs and the controls are in `analysis_outputs/`.
- Two large `bursts.csv` files are not committed (`analysis_outputs/detection_grid/`, `steady/grid/`, `fixedmask/grid/`). Regenerate them with `checks/detection_grid.m` (minutes to ~1.5 h), which needs only the repo `.mat` files.

**Needs the .nd2 / TIFF movies (lab PC only):**
- re-extraction (`collect_intensity.m` on the work folders);
- the unchanging-granule control (`checks/static_control.m`);
- the camera-noise calibration (`checks/noise_from_background.m`).
- Their outputs are committed, so they don't need rerunning unless the extraction changes.
- **In vitro:** the raw titration movies were found and re-extracted (sections 38-43). **Primary in vitro data = `in_vitro_reextracted/same_granules_R2/`**: the stored granules (tracks rebuilt by exact point matching), read with R2 (offset already removed: do not subtract 200 again); results in `analysis_outputs/vitro_same/` (`py checks/vitro_new.py --same`). Variants: `after_R2_offset/` (all re-tracked granules; `vitro_new/` without the dropped movies, `vitro_new_all/` with all). Dropped movies: `checks/vitro_movies.py`.

## In vitro raw movies (section 38)
- `..\Raw_Results\Titration experiment\`, 42 single-channel movies. The tracking script is `..\zigzag_scratch\track_all_vitro.ps1`; it uses junctions in `..\zigzag_scratch\vitro_raw` because the folder name contains a space, which Bio-Formats can't handle. Outputs: `..\zigzag_scratch\reextract\vitro_*\{track,work}`.
- **Mapping, from trace matching (it differs from the user's first guess; the user should confirm the ratio names):**

| stored dataset | folder | movies | file stem |
|---|---|---|---|
| protein_RNA_1_1 | 1000ng_ul_RNA_0.1mg_ml_protein | 12 | vitro_1000ngul_0p1mgml |
| protein_RNA_10_1 | 1000ng_RNA_10ug_protein | 5 | vitro_1000ng_10ug |
| protein_RNA_100_1 | 1000ng_ul_RNA_1mg_ml_protein | 10 | vitro_1000ngul_1mgml |
| extra (not in stored data) | 100ng_ul_RNA_1mg_ml_protein | 15 | vitro_100ngul_1mgml |

- **Scripts:** `checks/vitro_compare.py` (verification), `checks/vitro_new.py` (all parameters, movie-level), and `detection_grid.m` with the `vitroFolder`, `vitroJobs` and `doVivo` arguments.
- **Controls on the new data (section 39):** camera noise ratio is about 2.4, and an unchanging granule reads about 0 turnover beyond noise, so the copies are a valid floor in vitro. The in vitro manuscript figures with the new data are made by `py checks/ms_figures.py --vitro-new` → `analysis_outputs/ms_revision/figs_vitro_new/` (Figs. 6, 14, 18). The user applies the text changes by hand.

## Disk: `..\zigzag_scratch\` (about 11 GB after the 2026-10-25 cleanup)
- On 2026-10-25, `tif\` (the 2-channel TIFFs) and the diagnostic leftovers were deleted, about 33 GB in total.
- `work\*.tif` are now hard links to `track\<name>.tif`, the red-only TIFFs. All 56 in vivo movies were checked pixel-identical (`reextract\red_copy_check.csv`).
- Keep `track\` and `work\`.

## Data and extraction stages (`in_vivo_reextracted/`)
- **Datasets in vivo:**
  - 4PP7 (QCP-5x/PCP-4x RNA + PP7), 4PP7-FUS, 30PP7 (PCP-30x), 30PP7-FUS: 8-12 movies each;
  - 5Qb (QCP-5x + Qβ), 10Qb (QCP-10x + Qβ): 10 and 9 movies.
- **In vitro:** protein_RNA_1_1, _10_1, _100_1 (PP7 + PCP-14x RNA).
- **Original bug (section 4):** 2-channel movies (red + DIC interleaved) were tracked and read page by page, so many PP7 traces alternated red/DIC (the "zigzag"). Affected: 4PP7 (8 of 9 movies), 30PP7, 30PP7-FUS. Reproduced bit-exactly (59/59 traces) with the original code.
- **Stages:**
  - **R1:** fluorescence channel only, 0-based frame/xy fix.
  - **R2:** + camera offset 200. **R2 = the published multithresh method on clean data: the primary dataset for everything.**
  - `after_R2_sums`: spot sum (integrated granule intensity above the cell) and spot area.
  - R3 (5 px aperture): reverted at the user's request; kept on disk only.
  - **R4** (`opts.thrWin = 31`): steady thresholds, test only (section 23).
  - **R5** (`opts.fixedMask`): one fixed mask per granule, test only (section 33).
- `collect_intensity(path, struct('offset',200, 'thrWin',0, 'fixedMask',false, ...))`. Defaults reproduce R2.
- `extraction_info*.csv` (in `after_R2_offset`): movie names and kept-track counts, used to assign movie ids to traces.

## Main findings so far (details: RESULTS.md section given)
1. **Baseline** reproduces every manuscript number (section 1). Bug fixes and model changes, with their λ/k0 effects, are in section 3 and `lambda_k0_before_after.md`.
2. **Noise floor (sections 18-21):** burst-free copies of each trace give 70-100% of the detected PP7 bursts. Noise reproduces 88-94% of PP7 turnover and about 45% for Qβ. No SP (3-15) or smoothing window (5-21) fixes this (section 20).
3. **Source of the noise (sections 21, 23, 32):**
   - The traces are 2-6× noisier than camera/photon noise, because the per-frame multithresh mask jumps 11-28% per frame.
   - An unchanging synthetic granule (real image, real track, calibrated noise) reproduces the real noise (16.7 vs 16.9% in 4PP7).
   - With per-frame thresholds that control reads about 0 turnover, so the white-noise burst-free copies are a valid noise floor.
4. **λ and k0 are not identifiable (sections 25, 26, 30, 31):**
   - Only λ·k0 (≈ the mean burst size) is.
   - The manuscript's Eq. (1.9), fitted to simulations exactly as MS section III, returns λ ≈ 2-3 for every true λ from 0.5 to 5, with six scores (MSE, likelihood, KS, Cramér-von Mises, QQ, WS) and with no data removed.
   - The true λ wins only when it is 2-3. Top-3-of-5 success is 62% vs 60% for chance (71% for always guessing 2-4).
   - λ follows the noise level (true 2 → 1 / 2 / 4 at SNR 5 / 1.5 / 0.5).
   - Causes: smoothing merges and misses bursts; the continuous Eq. (1.9) ignores whole units and invisible zero-unit bursts.
   - λ·k0 tracks the truth.
5. **Proposed report parameters (sections 24, 27-29;** table in `analysis_outputs/control_table/`):
   - mean and median burst size (the reliable λ·k0);
   - burst size of granules vs their noise copies, and the part beyond noise;
   - turnover beyond noise (%/h);
   - net change beyond bleaching (%/h; granule fade vs cell-background fade, section 22);
   - granules (movies).
   - Use the median for comparisons (skewed sizes); use the mean when relating to λ·k0 or total turnover.
6. **4x PP7 vs PP7-FUS (sections 20, 25, 26, 28, 32):**
   - FUS granules are about 2× dimmer (p = 0.002), so relatively noisier (22% vs 16%; *absolute* noise is lower), so noise alone inflates their apparent burst size.
   - The burst-size gap is 6.1 points (17.2 vs 11.0%, movie-level p = 0.0006). 2.6 points are reproduced by the noise copies; 3.5 points (57%) are beyond noise.
   - H0 "noise explains the whole gap": rejected one-sided (p = 0.031); p = 0.061 two-sided, p = 0.088 with a mean-log statistic. **A trend, not robust.**
   - Turnover beyond noise does not differ.
   - For 80% power, about 2-5× more movies are needed (about 50 per condition; section 20).
7. **30x:** no FUS effect in any measure. **Qβ vs PP7:** robustly more dynamic (turnover beyond noise 9-10 %/h vs ≈ 0; p 0.002-0.02).
8. **In vitro (primary data, sections 43, 53; stored granules, corrected pipeline, spot sum as in vivo, movie-level):** no net loss (+1 / +14 / +16 %/h) and no difference between protein:RNA ratios in any noise-controlled measure; turnover beyond noise 18 / 15 / 29 %/h; granule area grows +21 / +49 / +85 %/h (spreading or focus drift), which made the earlier spot-mean "net loss" (−21 / −39 / −70 %/h, p = 0.015). Only 100:1 turnover changes between halves (late). Controls hold.
9. **Alternative extractions (tests, not adopted):**
   - R4 steady thresholds: 2-3× less noise but a slow artefact (the control itself shows 11-14 %/h).
   - R5 fixed mask: 4× less noise and a clean control, but its slow signal in PP7 tracks the granule moving relative to the mask (correlation 0.76-0.90). Under R5 the 4x FUS effect vanishes.
   - **No mask design makes FUS stronger. R2 stays primary.**
10. The noise-free fluctuation G30 (section 19) is consistent with all of the above. In vivo conditions are not distinguishable by it.
11. **Raising SP (section 34):** real PP7 traces and their burst-free copies lose bursts together; by SP 18-22 both are nearly empty, while Qβ keeps an excess. In a PP7-like simulation the detector finds 30% of true bursts at SP 12 (81% of detections are noise) and 6% at SP 18, at about 3× their true size. The manuscript's simulation was about 25× less noisy than real PP7.
12. **Large slow changes (section 35):** counted without the run rule, 4x FUS granules show more ≥30% changes than noise, and more than 4x PP7 (p ≈ 0.003), as a fraction of themselves; counted as a fixed amount of light, the FUS excess disappears. Wider smoothing windows (21-31) find more true bursts, but about half the detections are still noise; 4x FUS vs PP7 turnover is not significant at any of 20 settings.
13. **First vs second 30 min (section 36):** in vivo measures are stable between halves. In vitro is not: 1:1 and 100:1 become more dynamic late (100:1 net loss −35 → −90 %/h), 10:1 slows (turnover 76 → 18 %/h). Report in vitro per half, or state the time dependence.

## Status of the original "analyses to produce"
| # | analysis | status |
|---|---|---|
| 1 | bug fixes + before/after λ, k0 | done (section 3) |
| 2 | likelihood fitter | done, `analysis/gelmetricsMLE.m` (sections 5, 7-9); λ not identified |
| 3 | simulation recovery grid λ × SNR | done (sections 5, 26, 30, 31) |
| 4 | noise floor / SP choice | done: no acceptable SP; noise-corrected measures instead (sections 18-21) |
| 5 | sensitivity table | SP × window grid done (section 20), top-trim (26a). The λ-based sweep was dropped along with λ. |
| 6 | granule/movie bootstrap; time halves | movie-level bootstrap/permutation everywhere; first vs second 30-min halves done (section 36) |
| 7 | smoothing-window sensitivity | done in the detection grid (section 20); `checks/window_duration_test.m` (which event durations the window suppresses) **written, not run** |

## Open items / decisions for the user
- Which numbers the paper reports. The recommendation is the control-table parameters, not λ and k0.
- 4x FUS: present as a trend, or collect more movies.
- Is there a fixed / crosslinked granule movie? It would be the direct, model-free noise control. Not available so far.
- Possible next analyses: run `window_duration_test.m`; a movement-insensitive granule measure.
- Summary deck for the PI: `analysis_outputs/deck/GelMetrics_revision_summary.pptx`, 14 slides, built by `checks/make_deck.py`. Speaker notes are in `checks/deck_notes.py`; one entry per slide, so edit both together.

## Manuscript revision documents (done 2026-10-23; local only, NOT in git)
- `revision_docs/` (excluded via .git/info/exclude; contains the unpublished manuscript, do not commit or push):
  - `MS_revised_tracked.docx` (433 tracked changes, 23 pages clean, author 'Amit Nelkin'), `MS_revised_clean.docx` (accepted), `Response_to_reviewers.docx`, `Revision_control.docx` (section-by-section table + 'please check' list).
  - `work/`: generators. `ms_edits.py` (all text edits, inserts, figures, Table I, `DELETE` paragraphs, `FIG_SHIFT` renumbering; keyed by Word paragraph index of `MS_base.docx` = copy of MS.docx; all strings use the OLD figure numbers, build_ops renumbers), `build_ops.py` -> `ops.json` (word-level, anchored), `apply_ops.ps1` (Word COM, Track Changes on, full-width sections for appendix figures/Table I, formatting pass with tracking off), `finalize.ps1` (clean copy + PDFs), `pdf2png.ps1` (Windows PDF renderer -> PNG for checking), `response.py` / `control.py` -> HTML -> `html2docx.ps1`.
  - Rebuild: `cd revision_docs/work; py -3 build_ops.py`, then apply_ops.ps1 (from MS_base.docx), finalize.ps1. **Run Word on copies in %TEMP% (outside OneDrive), then copy the outputs back:** OneDrive AutoSave silently wrote edits back into the input file once (2026-10-23). The scripts now also set `AutoSaveOn = $false`; MS.docx itself was never modified. Paragraph indices refer to MS_base.docx; if MS.docx changes, re-dump with dump.ps1 and re-check the indices.
- Figures: `checks/ms_figures.py` -> `analysis_outputs/ms_revision/figs/` (MS style); data for Fig. 4 from `checks/ms_fig4_sim.m`; Fig. 15 from `run_runners('ms_fig15_R2', ...)`; file names keep the old numbers (fig15/16 = new Figs. 11/12, fig17-22 = new 13-18); burst files for Figs. 6-9 in `analysis_outputs/ms_grid_w13sp12/` (detection_grid, window 13, SP 12).
- Decisions used: report burst size / turnover / net change with noise controls (not lambda, k0); 4x FUS = fractional, brightness-dependent effect (one-sided p = 0.031, two-sided 0.061); fixed-granule control not feasible because condensates cannot be 'paused': like an oil droplet they keep exchanging material to stay in equilibrium (user's reason, 2026-10-23); reference [30] stays (revised code goes to the same repository); FRAP dropped; concede-and-show-data tone.
- **The user now edits the revision_docs/*.docx files by hand. Never rebuild and copy over them** (that wiped the user's edits once, 2026-10-24). Give any further change in chat (where, old text, new text, table values) and let the user apply it. Exception (2026-10-25, on the user's explicit request): edits may go into the user's file as tracked changes, after a dated backup, working on a %TEMP% copy, with a hash check that the original is unchanged before writing (section 49).
- Figures are now 1-18: old appendix Figs. 11-14 (tuning of the original WS fit) were removed at the user's request, Appendices A-B shortened, old 15-22 renumbered to 11-18 in the manuscript, letter and control file.
- Open items for the user are listed at the end of Revision_control.docx (Fig. 10(a) arrow label 'lambda', Eq. 1.1 as Word equation, etc.).

## Where things are (scripts → outputs)
- **User-facing analysis (sections 45, 54; plain style, the user's request):** `ExperimentalAnalysis.m`, `PaperFigures_revised.m` (Figs. 6-9), `SimulationRunWithAnalysis.m` (Fig. 5: WS fit repeated, lambda*k0, all measures on simulations; section 57) (scripts; loading, plotting and the Fig. 5 fits are local functions inside them) and the only function files `burstMeasures.m`, `summarizeGranules.m` (all formulas), `compareConditions.m`. Keep this plain style: no new helper files, no cellfun or anonymous functions. The earlier `GelMetricsReport`/`GelMetricsCompare`/`gm_*` are in `legacy/`. One noise copy per granule (as in the paper) varies by up to ~2 points between seeds; ExperimentalAnalysis uses 20, the figure script 1 (paper seeds 18-23, in vitro 18-20).
- **Manuscript figures:** see the line above; `figures/` holds the PNGs.
- **Pipeline:** `BatchParticleTracker.ijm`, `read_tables.m`, `collect_intensity.m` (+ wrapper `Collect_intensity_data.m`), `generate_movArray.m`, `analysis.m` (optional 7th output `segInfo` per burst), `getPosition.m`, `modifiedPoissonFit.m`, `generateSim.m` (optional k0), runners `ExperimentalAnalysis.m` / `SimulationRunWithAnalysis.m` (run via `checks/run_runners.m`).
- **Helpers:** `checks/gm_loadDataset.m`, `analysis/gm_amplitudes_by_trace.m`, `analysis/gelmetricsMLE.m`, `analysis/gelmetricsJoint.m`.
- **Noise floor and detection:**
  - `checks/score_export2.m`, `checks/robustness.py` (section 18);
  - `checks/detection_grid.m` (args: data folder, out name, SPs, windows, doVitro), `checks/grid_eval.py` (section 20);
  - `checks/power_fus.py` (section 20).
- **Noise source and controls:**
  - `checks/noise_from_background.m`, `noise_background_summary.py`, `noise_jitter_sources.py` (section 21);
  - `checks/static_control.m` (args: work root, datasets, out name, versions 0 / 31 / −1), `static_control_summary.py`, `plot_steady_static.py` (sections 23, 32);
  - `checks/reextract_steady.m`, `steady_compare.py` (section 23); `checks/reextract_fixed.m`, `compare_masks.py`, `fixedmask_motion.py` (section 33).
- **Fluctuation and net change:** `checks/noise_methods_real.m` (args: folders, out file), `noise_methods.py` (section 19); `checks/net_loss.py` (section 22).
- **Report parameters and FUS:** `checks/report_params_sim.m`, `report_params.py` (section 24); `checks/fus_brightness.py [--trim 10]` (sections 25-26); `checks/gap_test.py` (section 28); `checks/control_table.py` (sections 27-29).
- **λ, k0:**
  - `checks/lambda_k0_degeneracy.py` (section 25);
  - `checks/export_amps_shape.m`, `lambda_k0_shape.py [--notrim] [--plot]` (sections 26-28);
  - `checks/ms_sim_test.m` (args: λ list, noises, tag), `ms_eq19_fit.py`, `ms_fig13.py` (sections 30-31).
- **High SP, slow changes, time halves:** `checks/example_traces_sp.m`, `high_sp.py`, `sim_high_sp.m`, `sim_high_sp_plot.py` (section 34); `checks/slow_changes.py [--absolute]`, `slow_changes_test.py`, `win_sp.py` (section 35); `checks/halves.py` (section 36).
- **Deck and summary page:** `checks/make_deck.py`, `checks/deck_notes.py` → `analysis_outputs/deck/`; `analysis_outputs/summary/index.html`.

## Session notes
- Remote Control drops when the app restarts the session (PC sleep, app update). Re-enable it with the ccd_session_mgmt `set_remote_control` tool when the user asks. The app setting "Connect new sessions to Remote Control" is on.
- Another chat may push to `revision` in parallel: fetch before committing, and merge (never force-push or reset).

## Working style that worked
- Fix the decision rule before looking at outcomes (e.g. the setting chosen by noise, not by p value). Report the whole grid.
- Every comparison: movie-level permutation p and nested bootstrap CI, plus the naive burst-level p only for contrast.
- Long MATLAB jobs: run in the background and log to `..\zigzag_scratch\*.log`. Keep scripts fast (vectorise; avoid pandas inside permutation loops).
- Check figures visually before sending, and send key figures with SendUserFile.
