# CLAUDE_laptop_vitro.md: put the re-extracted in vitro results into the manuscript (laptop chat)

Read `CLAUDE.md` first: all its rules still hold. This file adds one task, written by the lab-PC chat on 2026-10-25. The lab PC has the raw movies but not the Word files. The laptop has the Word files but not the movies. All the data work is done and pushed. What is left is the manuscript text, the tables and the figures.

## First
- Run `git -c http.sslBackend=schannel pull` on branch `revision`. You need commits 935ac76 and 66403b9 or later.
- Read `analysis_outputs/RESULTS.md` sections 38 and 39. They hold every number below.
- The new figures are in `analysis_outputs/ms_revision/figs_vitro_new/` (png at 300 dpi, plus pdf). Regenerate them with `py -3 checks/ms_figures.py --vitro-new`; the default run, without the flag, is unchanged.

## What changed (short)
- The raw in vitro titration movies were found (`Raw_Results\Titration experiment`, lab PC only), re-tracked and re-extracted with the same R2 pipeline as in vivo.
  - The data are in `in_vitro_reextracted/after_R2_offset/`, with the offset already removed.
  - All 42 movies are single-channel, so the zigzag bug never affected in vitro.
- **The stored traces came from these movies:** 30-42% of them have an exact re-extracted twin.
- **Mapping found from the data:**

| stored dataset | folder | movies |
|---|---|---|
| 1:1 | 1000ng_ul_RNA_0.1mg_ml_protein | 12 |
| 10:1 | 1000ng_RNA_10ug_protein | 5 |
| 100:1 | 1000ng_ul_RNA_1mg_ml_protein | 10 |

  - A 4th folder, 100ng_ul_RNA_1mg_ml_protein (15 movies), is not in the paper. Leave it out unless the user asks.
- **Movies are now the replicates in vitro:** nested movie → granule bootstrap CIs, and movie-level permutation p values.
- **Controls now include in vitro:**
  - the trace noise is about 2.4× the camera noise (as PP7 in vivo);
  - an unchanging synthetic granule reads −4 to 0 %/h turnover beyond noise, against 13-45 %/h for the real granules. So the burst-free copies are a valid noise floor in vitro too.

## The task: put this into the revised manuscript and the response letter
**Never overwrite the user's files.** The user edits `revision_docs/*.docx` by hand; a rebuild wiped their edits once. Do NOT run `build_ops.py` or `apply_ops.ps1` on them. The paragraph indices in `ms_edits.py` are stale now anyway.

Two allowed ways (ask the user which):
- **A. Change list (default).** Read the current `revision_docs/MS_revised_tracked.docx`, `Response_to_reviewers.docx` and `Revision_control.docx`; read-only is fine. Give the user, for each change: section/paragraph, the exact old text, the exact new text. The user applies them.
- **B. Tracked changes on a new copy.**
  1. Copy the user's current file to %TEMP% (outside OneDrive; set `AutoSaveOn = $false`).
  2. Apply the edits with Word COM, Track Changes on, author 'Amit Nelkin', anchored by finding the text, not by paragraph index.
  3. Save under a NEW name, e.g. `MS_revised_tracked_vitro.docx`, next to the original.
  4. Never touch the original. The user reviews it and accepts or rejects.

### Changes to make
**1. Figures** (replace the images; keep the numbering):
- MS Fig. 6 ← `figs_vitro_new/fig6.png`. Caption: the error bars are now 95% CIs over movies, then granules; the n values in the panel titles are burst counts (114 / 203 / 155).
- MS Fig. 14 (appendix noise floor and its source) ← `figs_vitro_new/fig18.png`. Panels (b) and (c) now include in vitro. Caption: in vitro background = the dilute phase around the granule; the static synthetic granule is made from each in vitro track as in vivo.
- MS Fig. 18 (first vs second 30 min) ← `figs_vitro_new/fig22.png`.

**2. Table I, in vitro rows (old → new; new CIs are over movies):**

| | 1:1 | 10:1 | 100:1 |
|---|---|---|---|
| granules (movies) | 45 → 53 (12) | 37 → 86 (5) | 31 → 59 (10) |
| brightness (a.u.) | 914 → 864 | 1715 → 1456 | 1726 → 741 |
| frame noise (%) | 3.0 → 3.1 | 3.0 → 3.3 | 2.6 → 3.9 |
| burst size, granules / copies (%) | 7.9 / 10.4 → 8.0 / 9.2 | 12.2 / 8.8 → 11.1 / 9.8 | 15.2 / 10.0 → 14.5 / 9.6 |
| beyond noise (points) | −2.5 → −1.2 [−5.2, 5.2] | 3.4 → 1.3 [−1.0, 7.2] | 5.1 → 4.9 [0.6, 9.0] |
| median burst (a.u.) | 36 → 49 [28, 90] | 86 → 57 [45, 86] | 118 → 93 [33, 178] |
| mean burst (a.u.) | 129 → 167 [49, 440] | 234 → 182 [51, 332] | 432 → 417 [114, 1270] |
| turnover beyond noise (%/h) | 8.8 → 13.2 [−4.2, 38.5] | 50.1 → 42.1 [3.5, 105.3] | 46.2 → 45.3 [12.9, 73.6] |
| net change beyond bleaching (%/h) | −21.1 → −24.2 [−54.2, −3.4] | −39.7 → −41.3 [−108.2, −13.5] | −70.8 → −53.4 [−82.4, −34.3] |

- Source files: `analysis_outputs/vitro_new/summary.csv`, `summary.md`, `tests.csv`.
- **λ·k0 shape-fit column:** not redone for in vitro. Ask the user: leave it blank for in vitro, or run `checks/export_amps_shape.m` + `lambda_k0_shape.py` adapted to `in_vitro_reextracted/after_R2_offset` (MATLAB -batch with `checks/gm_noplot.m`; never render figures in MATLAB on the laptop).
- If Table I prints other columns, take them from `summary.csv` (columns `*_lo` / `*_hi` are the CI ends).

**3. Text: claims that change**
- **Methods:** say that the in vitro traces were extracted from the raw movies with the same pipeline as in vivo (fluorescence channel only, camera offset 200 subtracted), and that movies are the replicates: 12 / 5 / 10 movies.
- **Significance of the titration trend:** wherever the text gives it as significant (e.g. "z = −4.5", granule-level), change it.
  - The direction holds: more protein → larger bursts beyond noise, more turnover, more net loss.
  - With movies as replicates it is a trend only. 100:1 vs 1:1: net change p = 0.08, burst size beyond noise p = 0.10, turnover p = 0.18, median burst p = 0.22. 10:1 vs 1:1: all p > 0.4 (only 5 movies at 10:1).
  - Say this plainly; do not soften it. Use the user's "concede and show the data" tone.
- **First vs second half:**

| | turnover beyond noise, first → second half (%/h) | p | net change, first → second half (%/h) | p |
|---|---|---|---|---|
| 1:1 | 6.4 → 21.3 | 0.23 | −19.5 → −31.6 | 0.41 |
| 10:1 | 65.0 → 14.7 | 0.32 | −37.0 → −25.1 | 0.63 |
| 100:1 | −8.4 → 117.1 | <0.001 | −20.7 → −115.9 | 0.001 |

  - Only 100:1 changes significantly. Remove any statement that 1:1 or 10:1 change between halves (the old values were 1:1 turnover p = 0.006; 10:1 76 → 18, p = 0.007).
  - Replace the old 100:1 numbers (−35 → −90 %/h) with −21 → −116 %/h.
  - Keep the caveat: late dissolution and focus drift can't be told apart.
- **Brightness (needs care):** the re-extracted 100:1 granules are dimmer than 10:1 (741 vs 1456 a.u.; the stored traces gave 1726 vs 1715). Find any sentence saying the granules get brighter or denser with more protein (the Flory-Huggins sentence near Fig. 6(d)-(f)), and flag it to the user. Do not rewrite the biology yourself; give the numbers and ask.
- **Appendix noise section:** add one or two sentences: "The same controls on the in vitro movies give a measured ÷ camera noise ratio of 2.4 (2.1-2.6), and a static synthetic granule reads −4 to 0 %/h turnover beyond noise, versus 13-45 %/h in real granules (Fig. 14(b, c))." The numbers are in `analysis_outputs/steady/vitro_controls.md`.

**4. Response letter and Revision_control.docx**
- Referee 2, "no replicate analysis": in vitro now has movie-level statistics too (12 / 5 / 10 movies); give the 100:1 vs 1:1 p values honestly.
- Referee 2, "no noise floor": the in vitro static-granule control.
- Any place where the letter quotes the old in vitro numbers: update them.
- Add the new items to the "please check" list in Revision_control.docx (as a change list for the user, or in the new copy if way B was chosen).

## Questions for the user (ask, don't guess)
1. **Ratio labels.**
   - The three matched folders all have the same RNA (1000 ng), which fits the manuscript's "constant RNA (120 nM), varying protein".
   - But protein 0.1 → 1 mg/ml is only a 10× range for 1:1 → 100:1, and "10 µg" (10:1) has no volume.
   - Ask the user to confirm the 1:1 / 10:1 / 100:1 names and the Methods concentrations.
2. **The 100:1 brightness**, as above.
3. **The λ·k0 column for in vitro:** leave it blank, or recompute?
4. **Way A or way B** for the edits.

## When done
- Write a RESULTS.md section 40: what changed in the MS, and in which version and file.
- Update CLAUDE.md: remove this task, and point to the new file names.
- Commit and push the repo files only. `revision_docs/` stays local and is never committed.
- Give the user a short plain report.
