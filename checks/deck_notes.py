"""Speaker notes for checks/make_deck.py, one entry per slide, in slide order.
Each content slide: WHAT IT SHOWS / HOW IT WAS MADE / WHAT IT TELLS US.
Details and all numbers: analysis_outputs/RESULTS.md (section given per slide)."""

NOTES = [
# 1 title
"""Purpose of this meeting: agree on what the revised manuscript can claim, and which numbers it should report instead of λ and k0.

Everything in this deck was re-run from the raw .nd2 movies (in vivo) or the stored traces (in vitro, raw movies not available), with the code on the 'revision' branch. Full details and every number: analysis_outputs/RESULTS.md.""",

# 2 bottom line
"""Five take-home messages; each one has its own slide(s) after this.

1. A pipeline bug (slide 3) corrupted the PP7 traces. It is fixed and all PP7 and Qβ data were re-extracted.
2. When we give the burst detector copies of each trace that contain no bursts at all, it still finds about as many "bursts" in PP7 (slides 4-5).
3. λ and k0 cannot be told apart from the burst sizes; the manuscript's own fit returns λ ≈ 2-3 for any true λ (slides 6-9). Only their product, the mean burst size, is reliable.
4. The 4x FUS difference is mostly a brightness difference, with a borderline part beyond noise (slide 11).
5. Two differences are robust: Qβ is more dynamic than PP7, and in vitro granules lose material faster at higher protein:RNA (slide 12).""",

# 3 zigzag
"""WHAT IT SHOWS
Two stored 4PP7 traces, as used in the manuscript, produced by the original code. Top: a trace from a 1-channel movie (005) - every point is fluorescence (orange), a normal trace. Bottom: a trace from a 2-channel movie (007) - points alternate between the red fluorescence page (orange, ~5000) and the DIC brightfield page (blue, ~11000), giving the "zigzag".

HOW IT WAS MADE
I re-ran the original tracking macro, read_tables.m and Collect_intensity_data.m on the raw .nd2 files of movies 005 and 007. The re-created traces match the stored ones exactly (59 of 59 traces, bit for bit), so this is what the original pipeline did.
Mechanism: on 2-channel movies (red + DIC interleaved: red, DIC, red, ...) the tracker stepped through every page as if it were a time point, and the intensity reader then used page = frame, so many traces alternate between red and DIC pages. Even "clean-looking" traces from these movies were read entirely from DIC pages.

WHAT IT TELLS US
- Affected: 4PP7 (8 of 9 movies), all of 30PP7 and 30PP7-FUS. 4PP7-FUS (all 1-channel), Qβ and in vitro were not.
- Fix (user-approved): track and read the fluorescence channel only, correct the 0-based frame numbering, and subtract the camera offset (200). Everything after this slide uses the re-extracted data.
- Part of the published 4x PP7 vs PP7-FUS difference came from DIC contamination of 4PP7 only. (RESULTS.md sections 4, 6, 9)""",

# 4 noise floor
"""WHAT IT SHOWS
Left: one real 4PP7-FUS granule trace. Right: a simulated copy of the same trace that contains no bursts at all. Top row: the raw intensity. Bottom row: what the burst detector sees (smoothed over 13 frames and corrected for bleaching), with the bursts it called (red = losses, blue = gains).

HOW IT WAS MADE
For every real trace we built a burst-free copy: the dataset's typical bleaching curve, scaled to that granule's brightness, plus fresh random noise of the same size as that trace's own frame-to-frame noise (measured with a median, so real bursts barely affect it). Nothing else from the real trace goes in. The copy then goes through exactly the same detector, and anything detected in it is by construction a false burst.
"Share of turnover explained by noise" = the total fraction of the granule shed in bursts detected in the copies, divided by the same total in the real traces.

WHAT IT TELLS US
- The real and the burst-free traces look alike, and the detector calls bursts in both.
- PP7 and PP7-FUS: the copies reproduce 88-94% of the detected turnover - most PP7 "bursts" are indistinguishable from noise.
- Qβ: only ~45% - Qβ has real bursts beyond noise.
- Lowering SP (3-12) or changing the smoothing window (5-21 frames) does not separate real from false bursts in PP7: it adds both in the same proportion. (RESULTS.md sections 18-20)""",

# 5 noise source
"""WHAT IT SHOWS
Each dot is one granule (colour = dataset). x axis: the frame-to-frame noise we would expect from the camera and photon statistics alone. y axis: the noise actually measured on that granule's trace. Solid line = equal, dotted = twice as large. Three panels: the granule spot sum, the granule spot mean, and the cell background around the granule.

HOW IT WAS MADE
Prediction: from pixels far (> 10 px) from any granule we measured how each pixel's frame-to-frame variance grows with its brightness (a straight line, as expected for an EMCCD camera), then applied that line to the actual granule and cell pixels of every trace.
Control: we also built synthetic, unchanging granules - each real granule's average image, kept on its real track, bleached, with fresh camera-like noise every frame - and ran them through the same extraction.

WHAT IT TELLS US
- The traces are 2-6x noisier than camera noise predicts - and the cell background, which has no bursts, shows the same excess.
- Cause: the multithresh thresholds are recomputed every frame, so the set of pixels counted as "granule" changes by 11-28% every 10 s, and granule and cell levels jump together (correlation ~0.8).
- An unchanging synthetic granule reproduces the real noise almost exactly (16.7 vs 16.9% in 4PP7), so the extra noise is measurement, not granule dynamics.
- Holding the thresholds steady cuts the noise 2-3x but creates a slow artefact of its own, so it was not adopted. (RESULTS.md sections 21, 23)""",

# 6 degeneracy
"""WHAT IT SHOWS
Four panels on why λ and k0 cannot be separated:
(a) the burst-size histogram of 4PP7-FUS with four model curves - λ = 0.3, 1, 2 and 4, each with its best k0;
(b) the likelihood map over λ (y) and k0 (x): yellow = good fit; the white contour is the 95% region;
(c) data simulated from three different true (λ, k0) pairs;
(d) 200 simulated datasets from one known truth (λ = 1, k0 = 56), each fitted - every grey dot is one fitted (λ, k0).

HOW IT WAS MADE
Model: bursts carry n = 1, 2, 3 ... units (Poisson, zero units excluded), each unit worth k0, plus measurement noise of the measured size; fitted by maximum likelihood on the manuscript's amplitude window (15th-90th percentile).

WHAT IT TELLS US
- (a) The four very different explanations fit the same histogram equally well (likelihoods within 0.6 of each other, all pass a KS test).
- (b) The good region is a long diagonal ridge along "λ × k0 = constant", i.e. constant mean burst size.
- (c) Data from different truths look alike at this sample size.
- (d) Even with the truth known, the fitted λ lands within 2x of it only 44% of the time; the estimates slide along the ridge.
- So the data fix only the mean burst size, not how it splits into number of units and unit size. The manuscript's integer λ hides this by picking one point on a flat ridge. (RESULTS.md section 25)""",

# 7 manuscript fit on simulations
"""WHAT IT SHOWS
Left: the λ fitted by the manuscript's Eq. (1.9) against the true λ used in the simulation. Right: the fitted λ·k0 against the true λ·k0. Each colour is one scoring method; points are the median and bars the range over 5 independent simulated datasets; the dotted line is the truth.

HOW IT WAS MADE
Simulations exactly as manuscript section III: 1000 signals x 360 frames, bursts of k0 = 60 per unit with a Poisson number of units, white noise 40 a.u., bleaching, SP 15, 13-frame window; true λ = 0.5, 1, 2, 3, 4, 5. Every detected negative burst was kept (no 15% cut, no top trim). Eq. (1.9) was fitted as written, with six ways of scoring the fit: histogram MSE, likelihood, Kolmogorov-Smirnov, Cramér-von Mises (whole-curve distance), QQ error, and the manuscript's WS (also with λ restricted to 1-5 as in the paper). Controls: the same at noise 4 (almost noise-free), and the true jump sizes before noise and detection.

WHAT IT TELLS US
- Every scoring method returns λ ≈ 2-3 for every true λ from 0.5 to 5, and k0 drifts from ~30 to ~130 (true 60).
- λ·k0 (right) does follow the truth for every score, ~20-50 too high at this noise level.
- Why: the 13-frame smoothing turns steps into ramps, merges steps closer than the window and misses small ones, so the detected distribution has the same hump shape for any λ. On top of that, Eq. (1.9) is a smooth curve that ignores that bursts carry whole units and that 0-unit bursts are invisible: even on the true jump sizes it overestimates λ below 3.
- Removing the noise does not rescue it, so the failure is detection + model, not measurement noise. This is the referee's simulated-vs-estimated λ discrepancy. (RESULTS.md section 30)""",

# 8 Fig. 13 revisited
"""WHAT IT SHOWS
The manuscript's Fig. 13, re-done as "how often is each λ chosen". One panel per scoring method; rows = simulated datasets (true λ 0.5-5 at the manuscript's noise, and λ 2 at SNR 5 and 0.5, as in Fig. 13); columns = the λ the fit picked (1-5). Numbers = the fraction of fits that picked that λ; the red box marks the true λ.

HOW IT WAS MADE
Same simulations as the previous slide (5 datasets per row). For each candidate λ = 1-5 the best k0 was found, and the λ with the best score was chosen. To see whether the true λ was a close second, each dataset's granules (traces) were resampled 100 times and refitted: 500 fits per row. Separate heatmaps with the score and k0 at every candidate λ, for every method, are in analysis_outputs/ms_sim_test/fig13_heatmaps.png.

WHAT IT TELLS US
- The true λ is chosen only when it happens to be 2 or 3. At true λ 4 and 5 it is never chosen (0 of 500), and the likelihood rejects it decisively (2·ΔNLL = 136 and 334; > 4 means rejected).
- Relaxed criterion - true λ in the top 3 of the 5 candidates: 62% overall, vs 60% for picking 3 of 5 at random, and 71% for always guessing "2, 3 or 4" without looking at the data.
- The chosen λ follows the noise level, not the number of units: the same true λ = 2 is fitted as 1-2 at SNR 5, 2 at SNR 1.5 and 4 at SNR 0.5.
- If λ were known, k0 would come out about right (63-77 for true 60): the failure is in choosing λ. (RESULTS.md section 31)""",

# 9 shape fit
"""WHAT IT SHOWS
Left: λ fitted from the shape of the burst-size distribution against the true λ. Right: λ·k0 from the same fit against the true λ·k0. Colours = noise level (SNR 5, 1.5, 0.5); bars = 95% CI; dotted line = perfect recovery.

HOW IT WAS MADE
A different, deliberately favourable fit: the manuscript's model curve exactly as modifiedPoissonFit.m builds it (Poisson at whole numbers of units, zero units removed, joined by straight lines), turned into a probability density and fitted by likelihood to the shape of the distribution only - no bins, no WS score. Left panel: largest 10% of bursts removed; right panel: every burst used. Simulations with known λ (k0 = 60) through the same detection as the real data; CIs from 200 resamplings of the granules.

WHAT IT TELLS US
- Even with the best-behaved fit, λ stays at about 2 for true λ from 0.5 to 5.
- λ·k0 follows the truth within 0.8-1.35x across the range.
- On the real in vivo data the fitted λ is mostly below 1, where the product λ·k0 becomes unstable; there, the plain mean or median burst size is the safer number. (RESULTS.md sections 26-28)""",

# 10 parameter table
"""WHAT IT SHOWS
The proposed replacement for λ and k0: one row per dataset, each value with its 95% CI.

HOW EACH COLUMN WAS MADE
- Granules (movies): granules tracked ≥ 30 min, and the number of movies they come from (movies are the true replicates).
- Brightness: median granule intensity (spot mean, camera offset removed).
- Burst size, granules / noise: the median fraction of the granule lost per detected burst, in the real traces vs in their burst-free copies (slide 4).
- Beyond noise: the difference of those two - the part of the typical burst that noise does not produce.
- Median burst: the typical burst amplitude in manuscript units, all bursts. We use the median rather than the mean because burst sizes are skewed: a few very large (possibly merged) bursts pull the mean up (e.g. in vitro 1:1: median 36, mean 129, largest 5272). The mean is what λ·k0 represents and what adds up to total turnover; the median describes the typical burst and is more robust.
- λ·k0, shape fit: the product from the shape fit of slide 9 (all bursts). Unstable where λ < 1, which is most in vivo rows.
- Turnover beyond noise: % of the granule exchanged per hour in bursts, minus what the burst-free copies give.
- Net change: how much faster (negative) or slower the granule fades than its own surrounding cell background - net loss or gain of material beyond bleaching (same protein, same light).
CIs: resampling movies and then granules within them (in vitro has no movie labels, so granules only).

WHAT IT TELLS US
Red = the CI excludes zero (or, for 4x FUS, the beyond-noise value tested on its own). PP7 shows nothing beyond noise; Qβ shows turnover beyond noise; in vitro shows turnover and net loss that grow with protein:RNA. (RESULTS.md sections 24, 27-29)""",

# 11 FUS 4x
"""WHAT IT SHOWS
4x RNA, PP7 (blue) vs PP7-FUS (green), four panels:
(a) granule brightness of every granule, with the median and 95% CI;
(b) median burst size in the real granules (filled) and in their burst-free copies (open);
(c) the FUS-minus-PP7 gap in burst size split into the part the noise copies also show (grey: set by brightness/noise) and the part beyond noise (green);
(d) turnover beyond noise.

HOW IT WAS MADE
Real traces vs burst-free copies (slide 4); median fraction of the granule shed per burst. p values by shuffling whole movies between the two conditions 10,000 times (movies are the replicates). The hypothesis test: H0 = "noise explains the entire gap", i.e. the real gap equals the gap between the noise copies.

WHAT IT TELLS US
- FUS granules are about half as bright (p = 0.002), so their traces are relatively noisier (22% vs 16% per frame), and the detector finds larger apparent bursts even without biology.
- Burst sizes differ clearly: 17.2% vs 11.0% (p = 0.0006); counting every burst as independent would give p = 8e-11, which overstates it.
- 2.6 of the 6.1-point gap is reproduced by the noise copies; 3.5 points (57%) is beyond noise.
- H0 is rejected one-sided (p = 0.031; one-sided because the direction is the manuscript's own claim), but it is p = 0.061 two-sided and p = 0.088 with a statistic that uses every burst (mean log size) - a trend, not a robust result.
- Turnover beyond noise does not differ. At 30x nothing differs.
- Removing the largest 10% of bursts changes nothing. A firm FUS claim needs ~2-5x more movies. (RESULTS.md sections 20, 25, 26, 28)""",

# 12 robust differences
"""WHAT IT SHOWS
Net change of each granule beyond bleaching (%/h): every grey dot is one granule, red = mean with 95% CI. Zero = the granule fades exactly as fast as its surroundings (bleaching only); below zero = it loses material faster than it bleaches. Left: in vivo, granule spot sum. Right: all datasets, granule intensity above the cell.

HOW IT WAS MADE
For each granule: the rate at which its intensity fades over the movie (robust straight-line fit of the log intensity), minus the rate at which its own cell background fades. Both are the same fluorescent protein under the same light, so they bleach equally. No burst detection is involved, so this is independent of the noise problem.

WHAT IT TELLS US
- In vivo granules are in net balance: every condition fades only by bleaching (within a few %/h); whatever is shed is taken back up.
- In vitro granules lose material steadily, more at higher protein:RNA: about -21, -40 and -71 %/h for 1:1, 10:1 and 100:1 (100:1 vs 1:1: z = -4.5). This holds for any assumed camera offset.
- Qβ shows 9-10 %/h of turnover beyond noise (PP7 ≈ 0; p = 0.005-0.02 vs PP7), and its bursts are larger than its noise copies.
- Caveat: a granule drifting out of focus would also fade faster than the background; for in vitro this cannot be excluded from intensity alone. (RESULTS.md sections 22, 24)""",

# 13 claims
"""Each manuscript claim and what the re-analysis supports.
- Supported: holds with the noise controls and movie-level statistics.
- Reframed: the observation holds, but not in terms of λ and k0.
- Not supported / Replace: does not hold, or the method behind it does not work.

The two claims with the most at stake:
- "FUS raises λ at low valency" becomes "FUS 4x granules are 2x dimmer, with a borderline larger burst size beyond noise (one-sided p = 0.03)".
- The in vitro k0 trend becomes "mean burst size and net loss grow with protein:RNA". The in vitro raw movies are not available, so this rests on the stored traces.""",

# 14 decisions
"""Decisions to make before writing the response letter.
1. Which numbers the paper reports: the recommendation is mean/median burst size (the reliable λ·k0), turnover beyond noise and net change, each with its noise control - not λ and k0 separately.
2. 4x FUS: present it as a trend, or collect more data. About 2-5x more movies per condition (~50) would be needed for 80% power; more granules per movie help less than more movies.
3. A fixed or crosslinked granule movie (same imaging, no exchange) would give a direct, model-free noise control instead of simulated copies.
4. In vitro: the raw movies are missing, so in vitro results cannot be re-extracted or checked with the static-granule control.
5. Extraction: keep the published per-frame thresholds (validated noise controls), or design a stable granule mask (lower noise, but needs its own control).""",
]
