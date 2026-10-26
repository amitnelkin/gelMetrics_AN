# CLAUDE_laptop_spotsum.md: in vitro switches to the spot sum (laptop chat)

Read `CLAUDE.md` first; all its rules hold. Written by the lab-PC chat on 2026-10-26.

## The user's decision (2026-10-26)
**Use the spot sum for both in vivo and in vitro** for everything that measures material moving:
- burst size (% of the granule);
- burst size beyond noise;
- turnover beyond noise;
- frame noise;
- net change beyond bleaching;
- first vs second half.

**Keep the spot mean** for granule brightness and for the burst amplitudes in a.u. (median / mean burst, λ·k0).

In vivo already uses the spot sum; nothing changes there. **In vitro changes**, and so do its conclusions. Details are in `analysis_outputs/RESULTS.md` section 53; read it first.

## Why (for the text and the letter)
- **What each measure is:**
  - the spot mean is light per pixel;
  - the spot sum is all the granule's light above the dilute phase, so it follows the amount of material.
- **The granules spread:** in vitro granules cover more pixels over the hour (+21 / +49 / +85 %/h at 1:1 / 10:1 / 100:1). The spot mean therefore falls while the total light stays or rises.
- **So the earlier result was the growing area:** the earlier "net loss beyond bleaching, larger at higher protein:RNA" was not lost material.
- **Spreading vs focus drift:** settling and spreading on the glass, or focus drift. The images can't tell these apart; say so.
- **One measure for the whole paper:** the same measure is now used in vivo and in vitro.

## Pull first
`git -c http.sslBackend=schannel pull` on `revision`. The files:
- **Numbers:**
  - `analysis_outputs/vitro_spotsum/summary.csv` (the paper values: paper seeds, nCopies 1);
  - `reports.csv` (nCopies 1 and 20, spot mean and spot sum side by side);
  - `tests.csv` (movie-level p);
  - `area_trend.csv`.
- **Controls:** `analysis_outputs/steady/vitro_controls_sum.md`.
- **Figures** (MS style, Python): `analysis_outputs/ms_revision/figs_vitro_spotsum/`:
  - `fig6` = MS Fig. 6;
  - `fig18` = MS Fig. 14;
  - `fig22` = MS Fig. 18.

  Regenerate with `py -3 checks/ms_figures.py --vitro-spotsum`.
- **MATLAB:** `figures/Fig6.png` (`PaperFigures_revised.m`) now also uses the spot sum, the new default of `gm_paper_dataset`. `gm_paper_dataset(name, 'spotmean')` gives the old version.

## The new in vitro numbers (1:1 / 10:1 / 100:1; 95% CI over movies, then granules)

| | 1:1 | 10:1 | 100:1 |
|---|---|---|---|
| granules ≥ 30 min (movies) | 46 (12) | 37 (4) | 31 (8) |
| brightness, median over movies (a.u., spot mean) | 929 | 1726 | 1609 |
| frame noise (% of granule) | 13.7 | 16.0 | 10.5 |
| burst size, granules / copies (%) | 11.9 / 10.5 | 13.3 / 10.7 | 15.8 / 11.6 |
| burst size beyond noise (% points) | 1.4 [−3.1, 6.1] | 2.5 [−0.2, 21.3] | 4.2 [−1.2, 14.8] |
| median burst (a.u., spot mean) | 49 [27, 71] | 77 [56, 144] | 124 [26, 319] |
| mean burst (a.u., spot mean) | 129 [45, 327] | 233 [67, 402] | 1075 [127, 2930] |
| turnover beyond noise (%/h) | 18.0 [6.2, 31.9] | 14.7 [−4.8, 40.1] | 29.3 [2.6, 63.2] |
| noise share of turnover | 0.47 | 0.53 | 0.54 |
| net change beyond bleaching (%/h) | +0.6 [−18.0, 21.3] | +13.9 [2.9, 41.7] | +16.2 [−11.3, 38.8] |
| turnover, 1st → 2nd half (%/h), p | 22.2 → 12.7, 0.50 | 31.2 → −5.0, 0.15 | −15.9 → 82.5, 0.011 |
| net change, 1st → 2nd half (%/h), p | 1.1 → 4.1, 0.84 | 20.0 → 3.7, 0.56 | 35.1 → 3.1, 0.36 |
| granule area change (%/h) | +21 | +49 | +85 |

**Movie-level tests (10 000 movie shuffles, `vitro_spotsum/tests.csv`, nCopies 1):**

| | beyond noise | turnover | net change | median burst (a.u.) |
|---|---|---|---|---|
| 100:1 vs 1:1 | p = 0.41 | 0.49 | 0.32 | 0.097 |
| 10:1 vs 1:1 | 0.74 | 0.85 | 0.38 | 0.17 |

With nCopies 20, 100:1 vs 1:1 beyond noise gives p = 0.26. Nothing in vitro reaches 0.05 between ratios.

**Controls with the spot sum (Fig. 14(c) in vitro bars):**
- unchanging granule, turnover beyond noise: 1.7 ± 8.8 / −1.3 ± 7.2 / 2.7 ± 5.8 %/h, against 18.0 / 14.7 / 29.3 for the real granules;
- measured ÷ camera noise (spot sum): 4.6 / 3.7 / 4.3, close to PP7 in vivo (3.6-4.1).
- Fig. 14(b) keeps the spot-mean ratios for all datasets (2.6 / 2.3 / 2.7 in vitro), as before. If the user prefers, cite the spot-sum ratios in the text.

## What to change in the manuscript (as a change list for the user; never overwrite their Word files)
Same rules as before:
- read the user's CURRENT `revision_docs/*.docx` (read-only);
- give section, the exact old text and the exact new text, checking that each old text occurs once;
- write it to `revision_docs/vitro_spotsum_change_list.md`;
- tracked changes go in a NEW copy only if the user asks.

Then:
1. **Methods:**
   - in vitro, as in vivo, burst sizes, turnover and net change are computed from the spot sum (integrated granule intensity above the dilute phase);
   - brightness and burst amplitudes in a.u. come from the spot mean;
   - remove any sentence saying in vitro uses spot mean − background.
2. **In vitro results** (Sec. V.A and Fig. 6):
   - replace the net-loss titration result (−21 / −39 / −70 %/h, 100:1 vs 1:1 p = 0.015) with the numbers above;
   - no significant difference between ratios in any noise-controlled measure;
   - turnover beyond noise is positive in all ratios (18 / 15 / 29 %/h);
   - net change is about 0 or slightly positive;
   - the granule area grows with time, most at 100:1 (spreading or focus drift);
   - the median burst amplitude in a.u. rises with protein:RNA (49 / 77 / 124; trend, p = 0.097).
3. **Fig. 6:**
   - replace with `figs_vitro_spotsum/fig6.png`;
   - caption: burst sizes as % of the spot sum; (d) mean burst in a.u. (spot mean); (f) net change from the spot sum.
4. **Table I, in vitro rows:** the table above. Frame noise is now 10-16% (was 2.5-3%).
5. **Time halves (Fig. 18 and text):**
   - replace with `figs_vitro_spotsum/fig22.png`;
   - only 100:1 turnover changes (−16 → 83 %/h, p = 0.011);
   - the in vitro net change no longer differs between halves.
6. **Appendix noise section and Fig. 14:**
   - replace with `figs_vitro_spotsum/fig18.png`;
   - in vitro static granule −1 to 3 %/h with the spot sum.
7. **Discussion / abstract:**
   - remove or soften any claim that in vitro dynamics or net loss depend on protein:RNA;
   - keep the method's validation in vitro: the copies and controls work, and the measures are reported with CIs;
   - Flory-Huggins sentences: re-check with the user, since brightness (spot mean) is unchanged but the "loss" is gone.
8. **Response letter:**
   - Referee 2 (measurement validity, interpretability): explain the switch to the spot sum and the spreading finding, showing the spot-mean numbers once for transparency (`vitro_spotsum/reports.csv`, measure = spotmean);
   - update every in vitro number quoted.

## Questions for the user
1. Show the spot-mean version once (appendix or letter) to explain the change, or only in the letter?
2. Fig. 14(b): keep the spot-mean noise ratios, or switch to spot-sum ratios for all datasets? The in vivo spot-sum ratios exist in `analysis_outputs/noise_background/summary.csv` (ratioSum).
3. Mention spreading vs focus drift in the Results or only in the Discussion?

## When done
- Write a RESULTS.md section (54), and update CLAUDE.md (task status).
- Commit and push the repo files only, never `revision_docs/`.
- Give the user a short plain report.
