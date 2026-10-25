# CLAUDE_laptop_vitro.md: update the manuscript to the FINAL in vitro data (laptop chat, v2)

Read `CLAUDE.md` first: all its rules still hold. Written by the lab-PC chat on 2026-10-25. It **replaces** the earlier version of this file.

The earlier in vitro change list (`revision_docs/vitro_change_list.md`, RESULTS.md section 40) used the section 38 numbers: every re-tracked granule and every movie. **Those numbers are superseded.** The user has decided:
- 10:1 movie 5 and 100:1 movies 5 and 7 were dropped on purpose by the original authors. The reason will come from the user; ask for it and add it to the Methods.
- **The in vitro analysis uses the same granules and traces as the original data, read with the corrected pipeline.** That means:
  - fluorescence channel;
  - correct frame and position;
  - camera offset 200 subtracted;
  - per-frame multithresh as published;
  - noise-controlled measures with movie-level statistics.

## First
- Run `git -c http.sslBackend=schannel pull` on branch `revision`.
- Read `analysis_outputs/RESULTS.md` sections 41-43; section 43 is the one that matters.
- **Primary in vitro numbers:**
  - `analysis_outputs/vitro_same/summary.csv`, `summary.md` and `tests.csv`;
  - controls: `analysis_outputs/steady/vitro_controls.md`.
- **Figures** (already drawn from the primary data): `analysis_outputs/ms_revision/figs_vitro_new/`:
  - `fig6` = MS Fig. 6;
  - `fig18` = MS Fig. 14;
  - `fig22` = MS Fig. 18.
- Regenerate the figures with `py -3 checks/ms_figures.py --vitro-new` (`--vitro-new` reads `vitro_same/`).

## How to edit (unchanged rules)
- **Never overwrite or rebuild the user's Word files** (`revision_docs/*.docx`); the user edits them by hand.
- **Default: a change list,** giving section, the exact old text and the exact new text.
  - The old text is whatever is in the user's CURRENT files. The user may already have applied some of `vitro_change_list.md`, so search for both the original numbers and the section-38 numbers.
  - Check that each old text occurs exactly once.
  - Write it to `revision_docs/vitro_change_list_v2.md`. Do not edit the first list; keep it as a record.
- **Only if the user asks:** tracked changes in a NEW copy (for example `MS_revised_tracked_vitro2.docx`).
  - Make it from a copy in %TEMP% with AutoSave off, author 'Amit Nelkin'.
  - Anchor each edit by finding its text, never by paragraph index.

## The numbers (replace every in vitro number in the MS, Table I, the response letter and Revision_control)
**Table I, in vitro rows** (95% CI over movies, then granules; brightness as the user chose = median over movies of the per-movie medians, raw minus offset):

| | 1:1 | 10:1 | 100:1 |
|---|---|---|---|
| granules ≥ 30 min (movies) | 46 (12) | 37 (4) | 31 (8) |
| brightness, median over movies (a.u.) | 929 | 1726 | 1609 |
| brightness, pooled median [CI] | 915 [549, 1503] | 1709 [1476, 1789] | 1694 [290, 3021] |
| frame noise (% of granule) | 3.0 | 2.9 | 2.5 |
| burst size, granules / copies (% of granule) | 7.1 / 6.7 | 11.5 / 9.6 | 15.9 / 8.7 |
| burst size beyond noise (% points) | 0.4 [−4.2, 5.5] | 1.8 [−1.7, 8.5] | 7.2 [0.1, 19.3] |
| median burst size (a.u.) | 49 [30, 71] | 77 [55, 144] | 124 [26, 319] |
| mean burst size (a.u., pooled) | 129 [46, 320] | 233 [67, 396] | 1075 [125, 2993] |
| turnover beyond noise (%/h) | 11.6 [−1.8, 31.6] | 47.8 [0.2, 117.3] | 43.1 [4.6, 103.8] |
| net change beyond bleaching (%/h) | −21.2 [−47.2, −2.9] | −39.2 [−125.1, −7.6] | −70.3 [−101.2, −44.8] |
| share of turnover reproduced by noise | 0.60 | 0.32 | 0.55 |

- **Mean burst at 100:1:** the value is driven by a few very large bursts; the CI runs 125-2993. Say so in the caption, or prefer the median. Ask the user.
- **λ·k0 shape-fit column:** not computed for in vitro. Keep whatever the user decided in the earlier round.

**Comparisons (movies shuffled, 10 000×; `vitro_same/tests.csv`):**

| | burst size beyond noise | turnover | net change | median burst |
|---|---|---|---|---|
| 100:1 vs 1:1 | p = 0.096 | p = 0.20 | **p = 0.015** | p = 0.096 |
| 10:1 vs 1:1 | 0.65 | 0.32 | 0.60 | 0.17 |
| 100:1 vs 10:1 | 0.30 | 0.93 | 0.45 | 0.39 |

**First vs second 30 min (`vitro_same/summary.md`):**

| | turnover (%/h) | p | net change (%/h) | p |
|---|---|---|---|---|
| 1:1 | 5.8 → 18.9 | 0.33 | −18.5 → −30.5 | 0.49 |
| 10:1 | 77.4 → 12.3 | 0.21 | −36.9 → −24.7 | 0.59 |
| 100:1 | −27.7 → 126.5 | < 0.001 | −35.2 → −85.1 | 0.006 |

**Controls (Fig. 14(b, c); `steady/vitro_controls.md`):**
- measured ÷ camera noise for the granule: 2.6 / 2.3 / 2.7;
- the same for the dilute-phase background: 4.9 / 3.6 / 5.1;
- unchanging synthetic granule, turnover beyond noise: −0.5 ± 4.7 / −3.0 ± 5.2 / −2.3 ± 7.1 %/h, against 11.6 / 47.8 / 43.1 for the real granules.

## What the text should now say
1. **Methods, data.**
   - The in vitro granules are the original ones (tracks rebuilt by matching every stored point to the re-tracked movies: 99.2% of points matched exactly).
   - They were read with the corrected extraction: fluorescence channel, correct frame and pixel, camera offset 200 subtracted.
   - Movies: 12 / 4 / 8, plus the user's reason for the dropped movies.
   - The verification can go in the response letter or Appendix: the original code on the raw movies reproduces 27-46% of the stored traces bit-for-bit, with identical bursts (section 41).
2. **Titration result.**
   - Net loss beyond bleaching grows with protein:RNA (−21 / −39 / −70 %/h), and 100:1 vs 1:1 is significant with movies as replicates (p = 0.015).
   - Burst size beyond noise and median burst follow the same order but are only trends (p ≈ 0.1).
   - Turnover beyond noise is higher at 10:1 and 100:1 than at 1:1, but not significantly (p 0.2-0.3).
   - 10:1 has only 4 movies.
   - Remove any granule-level significance claim (e.g. z = −4.5).
3. **Time dependence.** Only 100:1 changes between halves. Its loss and turnover are concentrated in the second half hour (late dissolution or focus drift; can't be told apart). Remove any statement that 1:1 or 10:1 change.
4. **Brightness / Flory-Huggins.**
   - With the original granules, 100:1 is about as bright as 10:1 (1609 vs 1726, median over movies), not dimmer.
   - The earlier "100:1 is dimmer" (741) came from the extra re-tracked granules and the dropped movies; undo any change made for it.
   - Re-check the softened Flory-Huggins sentence with the user.
5. **Appendix noise section.** The in vitro control sentence with the numbers above.
6. **Response letter.**
   - Replicates: in vitro movie-level statistics, with the p values above.
   - Noise floor: the in vitro static-granule control.
   - Data validity: the raw movies were re-analysed; the original granules were kept and re-read with the corrected pipeline.

## Questions for the user
1. The reason the three movies were dropped (for the Methods).
2. The mean burst at 100:1: keep it with a note, or report the median only?
3. Change list (default), or tracked changes in a new copy?

## When done
- Write RESULTS.md section 44, and update CLAUDE.md (task status, new file names).
- Commit and push the repo files only; never `revision_docs/`.
- Give the user a short plain report.
