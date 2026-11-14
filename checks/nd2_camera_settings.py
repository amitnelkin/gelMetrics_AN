# Usage (lab PC, needs the raw .nd2): python checks/nd2_camera_settings.py <files.nd2>  -> section 86
# camera and laser settings from the text blocks of an .nd2 (UTF-16), both byte alignments, file head and tail
import re, sys, os
def scan(f):
    sz = os.path.getsize(f); cams, lasers, cfg = [], set(), set()
    with open(f, 'rb') as h:
        for off in (0, max(0, sz - 40_000_000)):
            h.seek(off); b = h.read(40_000_000)
            for a in (0, 1):
                t = b[a:].decode('utf-16-le', 'ignore')
                for m in re.finditer(r'Camera Type: .{0,400}?Temperature', t, re.S):
                    s = m.group(0)
                    g = lambda k: (re.search(k + r':\s*([^\r\n]+)', s) or [None, '?'])[1].strip()
                    cams.append((g('Exposure'), 'EM ' + g('Multiplier'), g('Conversion Gain'), g('Readout Speed')))
                lasers |= set(re.findall(r'Line:\d; ExW:\d+; Power: [\d.]+; Active', t))
    return sorted(set(cams)), sorted(lasers), len(cams)
for f in sys.argv[1:]:
    c, l, n = scan(f)
    print(f.split('/')[-3][-22:], f.split('/')[-1], n, c, l)
