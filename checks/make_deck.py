"""MAKE_DECK  Summary deck of the revision analyses for the PI meeting.
Figures come from analysis_outputs/ (cropped copies in analysis_outputs/deck/assets);
the parameter table is read from analysis_outputs/control_table/control_table.csv.
Usage: python checks/make_deck.py
Writes analysis_outputs/deck/GelMetrics_revision_summary.pptx."""
import os
import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

AO = 'analysis_outputs'; OUT = os.path.join(AO, 'deck'); AS = os.path.join(OUT, 'assets'); os.makedirs(AS, exist_ok=True)
INK = RGBColor(0x1C, 0x24, 0x3A)      # dark slate: title slides, text
RED = RGBColor(0xD9, 0x48, 0x3B)      # mCherry red: key numbers
TEAL = RGBColor(0x1F, 0x7A, 0x74)     # secondary
MUTED = RGBColor(0x5F, 0x67, 0x78)
TINT = RGBColor(0xEE, 0xF2, 0xF5)     # card background
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HEAD, BODY = 'Cambria', 'Calibri'


def crop(src, box, name):
    im = Image.open(os.path.join(AO, src)); w, h = im.size
    im.crop(tuple(int(v) for v in (box[0] * w, box[1] * h, box[2] * w, box[3] * h))).save(os.path.join(AS, name))
    return os.path.join(AS, name)


def pic(slide, path, x, y, w=None, h=None):
    im = Image.open(path); ar = im.size[1] / im.size[0]
    if w is not None and h is not None:          # fit inside the box, keep aspect
        if w * ar > h: w = h / ar
        else: h = w * ar
    elif w is not None: h = w * ar
    else: w = h / ar
    return slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))


def text(slide, x, y, w, h, runs, size=16, color=INK, font=BODY, bold=False, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, spacing=6):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    tf.vertical_anchor = anchor
    for m in ('margin_left', 'margin_right', 'margin_top', 'margin_bottom'):
        setattr(tf, m, 0)
    paras = runs if isinstance(runs, list) else [runs]
    for i, p in enumerate(paras):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); para.alignment = align; para.space_after = Pt(spacing)
        segs = p if isinstance(p, tuple) and isinstance(p[0], tuple) else (p if isinstance(p, tuple) else ((p, {}),))
        if isinstance(segs, tuple) and len(segs) == 2 and isinstance(segs[1], dict):
            segs = (segs,)
        for s, o in segs:
            r = para.add_run(); r.text = s; f = r.font
            f.size = Pt(o.get('size', size)); f.bold = o.get('bold', bold); f.italic = o.get('italic', False)
            f.name = o.get('font', font); f.color.rgb = o.get('color', color)
    return tb


def bullets(slide, x, y, w, h, items, size=16):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for m in ('margin_left', 'margin_right', 'margin_top', 'margin_bottom'):
        setattr(tf, m, 0)
    for i, it in enumerate(items):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); para.space_after = Pt(9)
        pPr = para._p.get_or_add_pPr(); pPr.set('marL', str(Inches(0.25))); pPr.set('indent', str(-Inches(0.25)))
        bu = pPr.makeelement('{http://schemas.openxmlformats.org/drawingml/2006/main}buChar', {'char': '•'}); pPr.append(bu)
        segs = it if isinstance(it, list) else [(it, {})]
        for s, o in segs:
            r = para.add_run(); r.text = s; r.font.size = Pt(o.get('size', size)); r.font.name = BODY
            r.font.bold = o.get('bold', False); r.font.color.rgb = o.get('color', INK)
    return tb


def card(slide, x, y, w, h, fill=TINT):
    s = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    s.adjustments[0] = 0.08; s.fill.solid(); s.fill.fore_color.rgb = fill; s.line.fill.background(); s.shadow.inherit = False
    return s


def stat(slide, x, y, w, big, small, color=RED):
    card(slide, x, y, w, 1.55)
    text(slide, x + 0.25, y + 0.15, w - 0.5, 0.75, big, size=34, color=color, font=HEAD, bold=True)
    text(slide, x + 0.25, y + 0.88, w - 0.5, 0.6, small, size=12, color=MUTED)


def content(prs, title, notes=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    text(s, 0.6, 0.38, 12.1, 0.9, title, size=30, font=HEAD, bold=True, anchor=MSO_ANCHOR.MIDDLE)
    n = len(prs.slides)
    text(s, 12.2, 7.05, 0.6, 0.3, str(n), size=10, color=MUTED, align=PP_ALIGN.RIGHT)
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


def dark(prs, title, sub, notes=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    bg = s.background.fill; bg.solid(); bg.fore_color.rgb = INK
    text(s, 0.9, 2.3, 11.5, 1.6, title, size=42, color=WHITE, font=HEAD, bold=True)
    text(s, 0.9, 4.0, 11.5, 1.2, sub, size=18, color=RGBColor(0xC9, 0xD3, 0xE0))
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


prs = Presentation(); prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)

# 1 title
dark(prs, 'GelMetrics revision: what the re-analysis shows', 'Summary of the code and data checks for the response to the referees\nlab meeting, 20 Oct 2026',
     notes='Purpose: agree on what the revised manuscript can claim, and which parameters it reports.')

# 2 bottom line
s = content(prs, 'Bottom line', notes='Five messages; each has its own slide after this one.')
msgs = [('1', 'A pipeline bug corrupted the PP7 traces', 'Two-channel movies mixed DIC pages into the red signal. Fixed and re-extracted.'),
        ('2', 'Most PP7 "bursts" are noise', 'Burst-free copies of each trace reproduce 88-94% of the detected PP7 turnover.'),
        ('3', 'λ and k0 cannot be separated', 'The manuscript fit returns λ ≈ 2-3 for any true λ; only λ·k0 (≈ mean burst size) is recovered.'),
        ('4', 'FUS at 4x: mostly brightness, partly real', 'FUS granules are 2x dimmer; 57% of the burst-size gap is beyond noise (one-sided p = 0.03).'),
        ('5', 'Robust differences remain', 'Qβ is more dynamic than PP7; in vitro granules lose material faster at higher protein:RNA.')]
for i, (n, h, d) in enumerate(msgs):
    col, row = i % 2, i // 2
    x, y = 0.6 + col * 6.2, 1.55 + row * 1.85
    card(s, x, y, 5.9, 1.6)
    c = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x + 0.25), Inches(y + 0.3), Inches(0.62), Inches(0.62))
    c.fill.solid(); c.fill.fore_color.rgb = RED if i < 4 else TEAL; c.line.fill.background(); c.shadow.inherit = False
    tf = c.text_frame; tf.text = n; p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    p.runs[0].font.size = Pt(18); p.runs[0].font.bold = True; p.runs[0].font.color.rgb = WHITE; p.runs[0].font.name = HEAD
    text(s, x + 1.1, y + 0.22, 4.6, 0.45, h, size=17, bold=True)
    text(s, x + 1.1, y + 0.7, 4.6, 0.85, d, size=13, color=MUTED)

# 3 zigzag
s = content(prs, 'Two-channel movies mixed DIC into the PP7 traces',
            notes='Mosaic walked raw planes (red, DIC, red, ...) of 2-channel movies; the reader used page = frame, so traces alternated red/DIC. '
                  'Confirmed bit-exactly by re-running the original code on movies 005 and 007. Fix: track and read the fluorescence channel only; Mosaic is 0-based (+1); camera offset 200 subtracted.')
pic(s, crop('zigzag_diagnosis/zigzag_mechanism.png', (0, 0.04, 1, 0.49), 'zigzag.png'), 0.6, 1.45, w=6.6)
stat(s, 0.6, 5.25, 3.15, '8 of 9', '4PP7 movies were 2-channel (also all 30PP7 and 30PP7-FUS)', color=MUTED)
stat(s, 4.05, 5.25, 3.15, '59 / 59', 'stored traces reproduced bit-exactly by the original code', color=TEAL)
bullets(s, 7.65, 1.6, 5.05, 5.4, [
    [('What happened: ', {'bold': True}), ('on 2-channel movies the tracker and reader stepped through red and DIC pages alternately.', {})],
    [('Which data: ', {'bold': True}), ('4PP7 (8 of 9 movies), 30PP7, 30PP7-FUS. 4PP7-FUS and in vitro were unaffected.', {})],
    [('Fix: ', {'bold': True}), ('track and read the red channel only, 0-based frame mapping, camera offset removed.', {})],
    [('Re-extracted ', {'bold': True}), ('all PP7 datasets and Qβ from the .nd2 files; Qβ reproduces the stored data.', {})],
    [('Consequence: ', {'bold': True}), ('part of the published 4x PP7 vs PP7-FUS difference came from DIC contamination.', {})]], size=15)

# 4 noise floor
s = content(prs, 'Most PP7 bursts are reproduced by noise alone',
            notes='Each real trace gets a burst-free copy: the dataset\'s typical bleaching curve plus fresh random noise of that trace\'s own frame-to-frame size. '
                  'The same detector finds about as many bursts in the copies as in the real PP7 traces. Lowering SP or changing the smoothing window does not help (SP 3-12 x window 5-21 grid).')
pic(s, crop('noise_background/example_traces.png', (0, 0.04, 1, 0.5), 'traces.png'), 0.6, 1.45, w=7.5)
text(s, 0.6, 4.85, 7.5, 1.9, [(('Left: ', {'bold': True}), ('a real 4PP7-FUS granule. ', {})), (('Right: ', {'bold': True}), ('a simulated copy with the same bleaching and noise level, and no bursts. ', {})), (('Top: ', {'bold': True}), ('raw intensity. ', {})), (('Bottom: ', {'bold': True}), ('what the detector sees, with the bursts it called (red = loss, blue = gain).', {}))], size=14, color=MUTED, spacing=3)
stat(s, 8.5, 1.5, 4.2, '88-94%', 'of the detected PP7 and PP7-FUS turnover is reproduced by burst-free copies')
stat(s, 8.5, 3.3, 4.2, '≈45%', 'for Qβ: a real excess of bursts over noise', color=TEAL)
text(s, 8.5, 5.15, 4.2, 1.6, 'No smoothing window (5-21 frames) or stringency (SP 3-12) separates real bursts from noise in PP7.', size=14, color=MUTED)

# 5 noise source
s = content(prs, 'The noise comes from per-frame thresholding, not the camera',
            notes='Camera + photon noise was predicted from pixels far from any granule. Measured trace noise is 2-6x larger, and the cell background (no bursts) shows the same excess. '
                  'The multithresh thresholds are recomputed every frame, so the granule mask changes by 11-28% every 10 s. A synthetic unchanging granule with only camera noise reproduces the real noise (16.7 vs 16.9%). '
                  'Steady thresholds cut the noise 2-3x but add a slow artefact; not adopted.')
pic(s, os.path.join(AO, 'noise_background/noise_measured_vs_predicted.png'), 0.6, 1.5, w=12.1)
stat(s, 0.6, 5.45, 3.8, '2-6x', 'measured noise vs camera + photon noise')
stat(s, 4.75, 5.45, 3.8, '11-28%', 'frame-to-frame change in the granule mask', color=TEAL)
stat(s, 8.9, 5.45, 3.8, '16.7 vs 16.9%', 'noise of an unchanging synthetic granule vs real 4PP7', color=MUTED)

# 6 degeneracy
s = content(prs, 'λ and k0 cannot be separated from the burst sizes',
            notes='(a) The same 4PP7-FUS histogram is fitted equally well by λ = 0.3/k0 = 85 up to λ = 4/k0 = 23. (b) The likelihood is a long ridge along constant mean burst size. '
                  '(c) Data simulated from different truths look alike. (d) Repeated simulations from λ = 1: the fitted λ lands within 2x of the truth only 44% of the time. The manuscript\'s own Fig. 15b shows the same ridge.')
pic(s, os.path.join(AO, 'lambda_k0_degeneracy/degeneracy.png'), 0.6, 1.45, w=12.1)
bullets(s, 0.6, 4.6, 12.1, 2.4, [
    [('Same histogram, many answers: ', {'bold': True}), ('λ from 0.3 to 4 fits equally well; only λ·k0 (≈ the mean burst size) is fixed.', {})],
    [('Known truth, poor recovery: ', {'bold': True}), ('fitted λ within 2x of the true λ in only 44% of simulated datasets (range 0.05-5.7).', {})],
    [('Integer λ in the manuscript ', {'bold': True}), ('hides this: the grid search picks one point on a flat ridge.', {})]], size=15)

# 6b manuscript approach on simulations
s = content(prs, 'The manuscript’s fit returns λ ≈ 2-3 whatever the truth',
            notes='Simulations as MS section III (1000 signals, k0 = 60, noise 40, SP 15), true λ 0.5-5, 5 datasets each, all negative bursts kept. '
                  'Eq. (1.9) fitted with six scores: MSE, likelihood, KS, Cramer-von Mises, QQ, and the manuscript’s WS. Every score returns λ ≈ 2-3, and k0 drifts from about 30 to 130 (true 60). '
                  'Only λ·k0 tracks the truth. Reasons: the 13-frame smoothing merges close steps and misses small ones, so the detected distribution has the same hump shape for any λ; '
                  'Eq. (1.9) ignores that bursts carry whole units and that 0-unit bursts are invisible; and the answer depends on the score. Near-noise-free simulations fail too.')
pl = crop('ms_sim_test/recovery.png', (0.0, 0.034, 0.333, 0.355), 'ms_rec_lambda.png')
pk = crop('ms_sim_test/recovery.png', (0.0, 0.68, 0.333, 1.0), 'ms_rec_lamk0.png')
pic(s, pl, 0.6, 1.45, w=4.15); pic(s, pk, 4.9, 1.45, w=4.15)
text(s, 0.6, 4.45, 4.0, 0.4, 'Fitted λ (dotted: truth)', size=13, bold=True, color=RED)
text(s, 4.9, 4.45, 4.0, 0.4, 'Fitted λ·k0 (dotted: truth)', size=13, bold=True, color=TEAL)
stat(s, 9.2, 1.5, 3.5, 'λ ≈ 2-3', 'for every true λ from 0.5 to 5, with all six scores')
stat(s, 9.2, 3.25, 3.5, '30 → 130', 'fitted k0 (true value 60) drifts with λ', color=MUTED)
text(s, 0.6, 5.2, 8.4, 1.6, 'Why: smoothing merges and misses bursts, so the detected distribution has the same shape for any λ; and Eq. (1.9) ignores whole units and invisible 0-unit bursts.',
     size=13, color=MUTED)

# 6c Fig. 13 revisited
s = content(prs, 'Fig. 13 revisited: the true λ is not a close runner-up',
            notes='MS Fig. 13 rows (λ 1, 2, 4; λ 2 at SNR 5 and 0.5) plus λ 3 and 5. Candidate λ 1-5, k0 optimised. Granules resampled 100x per dataset (500 fits per row). '
                  'The true λ wins only when it is 2 or 3. At true λ 4 and 5 it is never chosen and the likelihood rejects it (2·ΔNLL 136 and 334). '
                  'The same true λ = 2 is fitted as 1-2 at SNR 5, 2 at SNR 1.5 and 4 at SNR 0.5: λ follows the noise level. '
                  'Success = true λ in the top 3 of 5: 62% overall, vs 60% for a random pick and 71% for always guessing 2, 3 or 4.')
pic(s, os.path.join(AO, 'ms_sim_test/fig13_choice.png'), 0.6, 1.4, w=12.1)
stat(s, 0.6, 4.25, 2.85, '62%', 'true λ in the top 3 of 5 candidates')
stat(s, 3.65, 4.25, 2.85, '60%', 'a random pick of 3 of 5', color=MUTED)
stat(s, 6.7, 4.25, 2.85, '71%', 'always guessing λ = 2, 3 or 4 (ignores the data)', color=MUTED)
stat(s, 9.75, 4.25, 2.95, '1 → 2 → 4', 'fitted λ for the same true λ = 2 as SNR drops from 5 to 0.5', color=TEAL)
text(s, 0.6, 6.05, 12.1, 0.9, 'The method has a built-in answer of λ ≈ 2-3: it is right only when the truth happens to sit there, and the chosen λ follows the noise level, not the number of units.',
     size=14, color=INK)

# 7 shape fit
s = content(prs, 'Shape fitting recovers λ·k0, but not λ',
            notes='The manuscript\'s model curve, fitted to the shape of all bursts by likelihood (no WS, no trimming). In simulations with known λ (k0 = 60), the fitted λ stays near 2 whatever the truth; '
                  'the fitted λ·k0 follows the true λ·k0 within 0.8-1.35x. In vivo λ is mostly below 1, where the product becomes unstable, so the plain mean or median burst size is the safer number.')
pa = crop('lambda_k0_shape_notrim/shape_fit.png', (0.5, 0.0, 0.75, 1), 'shape_c.png')
pb = crop('lambda_k0_shape/shape_fit.png', (0.5, 0.0, 0.75, 1), 'shape_c_lambda.png')
pic(s, pb, 0.6, 1.5, h=4.3); pic(s, pa, 6.7, 1.5, h=4.3)
text(s, 0.6, 6.0, 5.9, 0.9, 'Fitted λ stays ≈ 2 for true λ 0.5-5', size=16, bold=True, color=RED)
text(s, 6.7, 6.0, 5.9, 0.9, 'Fitted λ·k0 tracks the true value (0.8-1.35x)', size=16, bold=True, color=TEAL)

# 8 parameter table
T = pd.read_csv(os.path.join(AO, 'control_table/control_table.csv'))
s = content(prs, 'Proposed parameters, each with its noise control',
            notes='Median vs mean burst size: the distributions are right-skewed (a few very large bursts, possibly merged events). The mean is what λ·k0 represents and what adds up to turnover, '
                  'but it is pulled by rare large bursts (in vitro 1:1: median 36, mean 129, max 5272). The median describes the typical burst and is more robust; use it for comparisons. '
                  'Report the mean when relating to λ·k0 or total turnover. 95% CIs: movie -> granule bootstrap (in vitro: granules).')
cols = ['Dataset', 'Granules\n(movies)', 'Brightness\n(a.u.)', 'Burst size,\ngranules / noise (%)', 'Beyond\nnoise (pts)', 'Median burst\n(a.u.)', 'λ·k0, shape\nfit (a.u.)',
        'Turnover beyond\nnoise (%/h)', 'Net change\n(%/h)']
short = {'PP7 + QCP-5x/PCP-4x': 'PP7, 4x', 'PP7-FUS + QCP-5x/PCP-4x': 'PP7-FUS, 4x', 'PP7 + PCP-30x': 'PP7, 30x', 'PP7-FUS + PCP-30x': 'PP7-FUS, 30x',
         'Qβ + QCP-5x/PCP-4x': 'Qβ, 5x', 'Qβ + QCP-10x': 'Qβ, 10x', 'in vitro, protein:RNA 1:1': 'In vitro 1:1', 'in vitro, protein:RNA 10:1': 'In vitro 10:1',
         'in vitro, protein:RNA 100:1': 'In vitro 100:1'}
ci = lambda r, c, f: f"{r[c]:{f}} [{r[c + '_lo']:{f}}, {r[c + '_hi']:{f}}]"
rows = []
for _, r in T.iterrows():
    rows.append([short[r.label], f"{r.granules}" + (f" ({int(r.movies)})" if r.movies == r.movies else ''), f"{r.brightness:.0f}",
                 f"{r.sizeReal:.1f} / {r.sizeCopy:.1f}", ci(r, 'sizeBeyond', '.1f'), ci(r, 'medianBurst', '.0f'), ci(r, 'lamk0', '.0f'),
                 ci(r, 'turnover', '.1f'), ci(r, 'net', '.1f')])
widths = [1.35, 1.0, 1.0, 1.5, 1.4, 1.45, 1.45, 1.5, 1.45]
tbl = s.shapes.add_table(len(rows) + 1, len(cols), Inches(0.6), Inches(1.45), Inches(sum(widths)), Inches(5.3)).table
for j, w in enumerate(widths):
    tbl.columns[j].width = Inches(w)
for i in range(len(rows) + 1):
    for j in range(len(cols)):
        cell = tbl.cell(i, j); cell.margin_left = cell.margin_right = Inches(0.05); cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        tf = cell.text_frame; tf.text = cols[j] if i == 0 else rows[i - 1][j]
        for p in tf.paragraphs:
            p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            for rr in p.runs:
                rr.font.size = Pt(11 if i == 0 else 11); rr.font.name = BODY; rr.font.bold = (i == 0 or j == 0)
                rr.font.color.rgb = WHITE if i == 0 else INK
        cell.fill.solid(); cell.fill.fore_color.rgb = INK if i == 0 else (TINT if i % 2 == 0 else WHITE)
hl = {('PP7-FUS, 4x', 4), ('Qβ, 5x', 7), ('Qβ, 10x', 7), ('In vitro 10:1', 7), ('In vitro 100:1', 7), ('In vitro 1:1', 8), ('In vitro 10:1', 8), ('In vitro 100:1', 8)}
for i, r in enumerate(rows, 1):
    for j in range(len(cols)):
        if (r[0], j) in hl:
            for rr in tbl.cell(i, j).text_frame.paragraphs[0].runs:
                rr.font.bold = True; rr.font.color.rgb = RED
text(s, 0.6, 6.85, 12.1, 0.4, 'Red: CI excludes 0, or (4x FUS) the beyond-noise gap tested on its own. Median and λ·k0 in spot-mean units; burst size and turnover as % of the granule. 95% CI.',
     size=11, color=MUTED)

# 9 FUS 4x
s = content(prs, '4x FUS: half the gap is brightness, half beyond noise',
            notes='FUS granules are about half as bright, so their traces are relatively noisier and the detector finds larger apparent bursts even without biology. '
                  'Burst-free copies with the same brightness and noise reproduce 2.6 of the 6.1-point gap. H0 "noise explains the entire gap" is rejected one-sided (p = 0.031; the direction is the manuscript\'s claim), '
                  'two-sided p = 0.061, and not with the every-burst statistic (mean log size, p = 0.088). Turnover beyond noise does not differ. At 30x nothing differs.')
pic(s, crop('fus_brightness/fus_brightness.png', (0, 0.045, 1, 0.525), 'fus4x.png'), 0.6, 1.4, w=12.1)
stat(s, 0.6, 4.9, 2.85, '2x', 'dimmer FUS granules (p = 0.002)', color=MUTED)
stat(s, 3.65, 4.9, 2.85, '6.1 pts', 'burst-size gap (p = 0.0006)', color=MUTED)
stat(s, 6.7, 4.9, 2.85, '57%', 'of the gap beyond noise')
stat(s, 9.75, 4.9, 2.95, 'p = 0.031', 'one-sided; 0.061 two-sided; 0.088 every-burst statistic', color=RED)

# 10 robust differences
s = content(prs, 'Robust differences: Qβ vs PP7, and in vitro protein:RNA',
            notes='Net change = granule fade rate minus its own cell-background fade rate (same protein, same light). In vivo granules are in net balance. In vitro granules lose material, more at higher protein:RNA. '
                  'Qβ has turnover beyond noise (9-10 %/h; p = 0.018 and 0.005 vs PP7) and bursts larger than its noise copies. Caveat: focus drift in vitro would also look like net loss.')
pic(s, os.path.join(AO, 'net_loss/net_loss.png'), 0.6, 1.45, w=7.6)
stat(s, 8.6, 1.5, 4.1, '9-10 %/h', 'Qβ turnover beyond noise (PP7: ≈ 0; p = 0.005-0.02)', color=TEAL)
stat(s, 8.6, 3.25, 4.1, '-21 → -71 %/h', 'in vitro net loss from 1:1 to 100:1 (100:1 vs 1:1: z = -4.5)')
text(s, 8.6, 5.0, 4.1, 1.8, 'In vivo, every condition fades only by bleaching: no net loss or gain beyond a few %/h.', size=14, color=MUTED)

# 11 claims
s = content(prs, 'What happens to the manuscript\'s claims',
            notes='Supported = holds with noise control and movie-level statistics. Reframed = the observation holds but not as λ/k0.')
claims = [('PP7 granules: λ = 1, gel-like', 'Not supported', 'λ not identifiable; PP7 dynamics not distinguishable from noise'),
          ('Qβ: λ = 4, brittle', 'Reframed', 'Qβ is clearly more dynamic than PP7 (turnover beyond noise), λ itself not identifiable'),
          ('FUS raises λ at low valency', 'Reframed, weak', 'FUS 4x granules are 2x dimmer; burst size beyond noise is a trend (one-sided p = 0.03)'),
          ('No FUS effect at high valency', 'Supported', 'No difference in any measure at 30x'),
          ('In vitro k0 grows with protein', 'Reframed', 'Mean burst size and net loss grow with protein:RNA (in vitro has no raw movies to re-check)'),
          ('Integer λ from a grid search', 'Replace', 'Report λ·k0 / burst size with CIs and a noise control')]
tb = s.shapes.add_table(len(claims) + 1, 3, Inches(0.6), Inches(1.5), Inches(12.1), Inches(5.0)).table
for j, w in enumerate([4.0, 2.0, 6.1]):
    tb.columns[j].width = Inches(w)
STC = {'Not supported': RED, 'Supported': TEAL, 'Replace': RED}
for i, row in enumerate([('Manuscript claim', 'Status', 'What the re-analysis shows')] + claims):
    for j, v in enumerate(row):
        cell = tb.cell(i, j); cell.text_frame.text = v; cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        cell.fill.solid(); cell.fill.fore_color.rgb = INK if i == 0 else (TINT if i % 2 == 0 else WHITE)
        for rr in cell.text_frame.paragraphs[0].runs:
            rr.font.size = Pt(14); rr.font.name = BODY; rr.font.bold = (i == 0 or j == 1)
            rr.font.color.rgb = WHITE if i == 0 else (STC.get(v, RGBColor(0xB0, 0x6A, 0x10)) if j == 1 else INK)

# 12 decisions
s = dark(prs, 'Decisions for us', '', notes='Items to agree on before writing the response letter.')
items = ['Report λ·k0 / burst size, turnover and net change with noise controls, instead of λ and k0',
         'FUS 4x: present as a trend, or collect ~2-5x more movies (about 50 per condition for 80% power)',
         'A fixed or crosslinked granule movie would give a direct, model-free noise control',
         'In vitro raw movies are missing: in vitro conclusions rest on the stored traces',
         'Extraction: keep per-frame thresholds (published), or design a stable granule mask']
for i, it in enumerate(items):
    y = 3.35 + i * 0.72
    c = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.9), Inches(y), Inches(0.46), Inches(0.46))
    c.fill.solid(); c.fill.fore_color.rgb = RED; c.line.fill.background(); c.shadow.inherit = False
    tf = c.text_frame; tf.text = str(i + 1); p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    p.runs[0].font.size = Pt(14); p.runs[0].font.bold = True; p.runs[0].font.color.rgb = WHITE
    text(s, 1.6, y + 0.02, 11.0, 0.5, it, size=18, color=WHITE)
s.shapes[0].top = Inches(1.3)   # move the title up on this slide
prs.save(os.path.join(OUT, 'GelMetrics_revision_summary.pptx')); print('saved')
