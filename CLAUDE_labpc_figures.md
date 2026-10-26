# CLAUDE_labpc_figures.md: check the MATLAB figure script on the lab PC

Read `CLAUDE.md` first; all its rules hold. Written by the laptop chat on 2026-10-26. The laptop's MATLAB cannot draw figures (it hangs on `plot`), so the drawing code has never been run. The lab PC's MATLAB can draw.

## First
- Run `git -c http.sslBackend=schannel pull` on branch `revision`. You need commit 338071d or later.
- Read `analysis_outputs/RESULTS.md` section 51.

## Task
1. **Run `PaperFigures_revised.m` from the repository root** in MATLAB on the lab PC, the whole script.
   - It takes about 8 min: Figure 5 simulates.
   - It writes `figures/Fig5.png` … `Fig9.png` (300 dpi) and the PDFs.
   - Run it headless if possible: `matlab -batch "cd('<repo>'); PaperFigures_revised"`.
   - Watch out: `run('...')` switches into the script's folder, and the script finds the repository from its own location.
2. **Fix whatever errors.** Only the drawing code is untested; every computed number already matches the manuscript (section 51), so do not change the computations.
   - The likely trouble spots: `legend(nexttile(t, 1), ...)`, the `errorbar`/`bar` calls, the `tiledlayout` sizes, the label positions, and `exportgraphics` with `'ContentType','vector'`.
3. **Compare every figure with the manuscript version** (`analysis_outputs/ms_revision/figs/fig5.png` … `fig9.png`; the in vitro Fig. 6 is in `figs_vitro_new/fig6.png`). Look at each PNG.
   - The panels, data and numbers should be the same. The style should be close: salmon bars, black copy outline, boxed axes, bold (a) labels, Arial 8.
   - Fix overlaps, clipped labels and unreadable text.
   - Small differences in the last digit of CIs and permutation p values are expected (random draws).
4. **Also run `ExperimentalAnalysis.m`** (default: 4PP7 vs 4PP7-FUS, 20 copies) and check its two figures: burst size distributions, and summary bars with CIs.
5. **When done:**
   - commit and push the fixes and the `figures/*.png`;
   - add a RESULTS.md section (52);
   - remove this task from CLAUDE.md;
   - give the user a short plain report with the figures (SendUserFile).

## Do not
- Touch the user's own `PaperFigures.m` (untracked, theirs).
- Touch `revision_docs/` (not on the lab PC anyway).
