# CLAUDE_labpc_fig5.md: draw and check the new Fig. 5 (lab PC)

Read `CLAUDE.md` first; all its rules hold. Written by the laptop chat on 2026-10-26.

## Background
The user asked for a new MS Fig. 5 that shows the same output for simulated data as for experimental data (RESULTS.md section 57):
- **(a-c)** the original WS fit, repeated on 20 simulated experiments per condition:
  - (a) λ;
  - (b) k0 against λ;
  - (c) λ·k0, both from the WS fit and as the mean burst size.
- **(d-f)** burst sizes of the granules vs their burst-free copies.
- **(g-i)** mean burst size, turnover beyond noise and net change, with the true values marked.
- Conditions: λ 1, 2, 4 at SNR 5 and 0.5.

`SimulationRunWithAnalysis.m` makes the whole figure. `PaperFigures_revised.m` no longer makes Fig. 5. The numbers were checked on the laptop (headless). The drawing code has not been run, because the laptop's MATLAB cannot draw.

## Task
1. Pull branch `revision` (`git -c http.sslBackend=schannel pull`).
2. In `SimulationRunWithAnalysis.m`, set `saveFigures = true`. Run it from the repository root (about 15-20 min): `matlab -batch "cd('<repo>'); SimulationRunWithAnalysis"`.
3. **The printed numbers must equal RESULTS.md section 57** (same seeds).
4. **Look at `figures/Fig5.png`.** The layout should match the laptop's Python preview, `analysis_outputs/ms_revision/fig5_sim_preview.png`. Fix any layout problem in the plain style:
   - no new helper files;
   - no `cellfun`, anonymous functions or option structures;
   - comments in plain words;
   - no change to any computation.
5. **When done:**
   - set `saveFigures` back to false;
   - commit and push `figures/Fig5.png/.pdf` and any fixes;
   - add a RESULTS.md section;
   - remove this task from CLAUDE.md;
   - send the user `figures/Fig5.png` (SendUserFile) with a short plain report.

## Do not
- Touch the user's own `PaperFigures.m` (untracked) or the Word files in `revision_docs/`.
