# CLAUDE_labpc_figures2.md: check the drawings of the simplified scripts (lab PC)

Read `CLAUDE.md` first; all its rules hold. Written by the laptop chat on 2026-10-26.

The user found the previous code too complex. It was rewritten in a plain style (RESULTS.md section 54); it is the user's code from now on, so keep the plain style when you fix things:
- `ExperimentalAnalysis.m`, `PaperFigures_revised.m`, `SimulationRunWithAnalysis.m`: plain scripts; data loading, plotting and the Fig. 5 fits are local functions at the end of each file;
- `burstMeasures.m`, `summarizeGranules.m`, `compareConditions.m`: the only function files;
- the earlier helpers (`GelMetricsReport.m`, `GelMetricsCompare.m`, `gm_*.m`) are now in `legacy/`.

Every number was checked on the laptop: identical to the previous version, CIs and p values included (section 54). The drawing code was carried over from your section 52 fixes, but it has not been run: the laptop's MATLAB cannot draw.

## Task
1. Pull branch `revision` (`git -c http.sslBackend=schannel pull`).
2. **Run `PaperFigures_revised.m` from the repository root**, the whole script. It takes about 10 min; it writes `figures/Fig5-9.png/.pdf`.
   - The script uses paths relative to the repository folder, so `cd` there first.
   - For example: `matlab -batch "cd('<repo>'); PaperFigures_revised"`.
3. **Compare each PNG** with your checked versions from section 52 (`git show 6971c48:figures/Fig5.png > old5.png`, etc.):
   - they should look the same;
   - Fig. 6 now uses the in vitro spot sum, as it did after section 53.
   - Fix any error or layout problem in the same plain style:
     - no new helper files;
     - no `cellfun`, anonymous functions or option structures;
     - comments in plain words.
4. **Run `ExperimentalAnalysis.m`** (default: 4PP7 vs 4PP7-FUS, 20 copies), with `saveFigures = true` for the check, and look at its two figures.
5. **Run `SimulationRunWithAnalysis.m`** (no figures; just check that it runs).
6. **When done:**
   - commit and push the fixes and `figures/*.png` (set `saveFigures` back to false);
   - add a RESULTS.md section;
   - remove this task from CLAUDE.md;
   - give the user a short plain report, with the figures (SendUserFile).

## Do not
- Touch the user's own `PaperFigures.m` (untracked).
- Change any computation: the numbers are verified.
