# CLAUDE_labpc_vitro_spotsum.md: extract the in vitro spot sums and compare (lab PC)

Read `CLAUDE.md` first; all its rules hold. Written by the laptop chat on 2026-10-26. Do this after the figure check (`CLAUDE_labpc_figures.md`), or whenever that one is done.

## Why
- **In vivo**, burst fractions, turnover, frame noise and net change use the spot sum: the integrated granule light above the background (`after_R2_sums/`).
- **In vitro** only the spot mean was ever saved:
  - the original stored data had nothing else;
  - `checks/vitro_same_granules.m` (line 66) calls `collect_intensity` with 2 outputs only, so `sums.spotSum` was never kept.
- So in vitro fractions are relative to (spot mean − background). That measure does not include how many pixels the granule covers. The user wants the in vitro spot sums, to see whether the in vitro results change.

## Task
1. Pull branch `revision` (`git -c http.sslBackend=schannel pull`).
2. **Extract the spot sums for the primary in vitro data:**
   - The data are the same granules, the same rebuilt tracks and the same R2 settings as `in_vitro_reextracted/same_granules_R2/`.
   - Change `checks/vitro_same_granules.m` so that it also takes the 4th output (`[intensity, background, info, sums] = collect_intensity(...)`).
   - Save `same_<ratio>_spotsum_data.mat` (variable with the same layout as the in vivo `*_spotsum_data.mat`, first cell empty) and `same_<ratio>_spotsum_background.mat` in the same folder.
   - **Do not change** the existing `same_*_data.mat` / `_background.mat`. Check that the re-run reproduces them bit-for-bit.
   - Never modify `in_vitro_experiments/`.
3. **`gm_paper_dataset.m`:** for `vitro_*`, add `D.spotSum` when the spotsum file exists, exactly as for in vivo. Keep a way to get the old version (e.g. `gm_paper_dataset(name, 'spotmean')`), so the paper numbers stay reproducible.
4. **Compare spot mean vs spot sum** with `GelMetricsReport` (paper seeds, `nCopies = 1`; and `nCopies = 20`) and `GelMetricsCompare`:
   - **per ratio:** frame noise, burst size granules / copies, beyond noise, turnover beyond noise, net change, first vs second half;
   - **tests:** 100:1 vs 1:1 and 10:1 vs 1:1 (beyond noise, turnover, net change; two-sided and one-sided p);
   - **the in vitro controls** (`checks/vitro_controls.py`, static granule and noise ratio), if they are cheap to redo with the spot sum.
5. **Report honestly, do not choose:**
   - Give a table of spot mean vs spot sum, and say plainly which conclusions change. The manuscript currently states: net loss −21 / −39 / −70 %/h, 100:1 vs 1:1 p = 0.015; burst size beyond noise is a trend (p ≈ 0.1).
   - The user decides which version the paper uses. Do not change the manuscript numbers or the figure scripts' defaults.
6. **When done:**
   - commit and push (scripts, the new `.mat` files, outputs);
   - add a RESULTS.md section;
   - remove this task from CLAUDE.md;
   - give the user a short plain report.
