#!/usr/bin/env python3
"""r-profile stats (vr, |B|) from an AMRVAC vtuBCCmpi snapshot.

Usage: python3 rprofile_B.py [file.vtu]   (default: xia_weno30048.vtu)
Prints: piece count + first-piece array layout, then per-cell r computed
from quad centers; stats for outer ring (r>215) and inner region (r<3),
plus global min/med/max of vr and |B|.
"""
import math
import re
import struct
import sys

fn = sys.argv[1] if len(sys.argv) > 1 else 'xia_weno30048.vtu'
raw = open(fn, 'rb').read()
i0 = raw.find(b'<AppendedData')
i1 = raw.find(b'>', i0) + 1
i2 = raw.find(b'</AppendedData>', i1)
app = raw[i1:i2]
app = app[app.find(b'_') + 1:]
xml = raw[:i0].decode('ascii', 'ignore')

pieces = list(re.finditer(r'<Piece\b([^>]*)>(.*?)</Piece>', xml, re.S))
print('pieces:', len(pieces))

W = {'Float32': 4, 'Float64': 8, 'Int32': 4, 'Int64': 8, 'UInt8': 1, 'UInt32': 4}


def rd(off, n, w, ch):
    if w == 8:
        fmt = '<%d%s' % (n, 'd' if ch == 'f' else 'q')
    elif w == 4:
        fmt = '<%d%s' % (n, 'f' if ch == 'f' else 'i')
    else:
        fmt = '<%dB' % n
    return struct.unpack_from(fmt, app, off + 4)


layout_done = False
outer, inner, gvr, gB = [], [], [], []
for pm in pieces:
    attrs, body = pm.group(1), pm.group(2)
    nc = int(re.search(r'NumberOfCells="(\d+)"', attrs).group(1))
    nv = int(re.search(r'NumberOfPoints="(\d+)"', attrs).group(1))
    arrays = {}
    for sm in re.finditer(r'<(Points|Cells|CellData)\b[^>]*>(.*?)</\1>', body, re.S):
        sec = sm.group(1)
        for dm in re.finditer(r'<DataArray\b([^>]*?)/?>', sm.group(2)):
            a = dm.group(1)
            nm = re.search(r'Name="([^"]*)"', a)
            ty = re.search(r'type="([^"]*)"', a)
            nco = re.search(r'NumberOfComponents="(\d+)"', a)
            of = re.search(r'offset="(\d+)"', a)
            if of:
                arrays[(sec, nm.group(1) if nm else '')] = (
                    ty.group(1) if ty else 'Float32',
                    int(nco.group(1)) if nco else 1,
                    int(of.group(1)))
    if not layout_done and len(pieces) > 0:
        for k, v in sorted(arrays.items()):
            print(' layout', k, v)
        layout_done = True
    p = None
    for k, v in arrays.items():
        if k[0] == 'Points':
            p = v
            break
    if p is None:
        continue
    pty, pnco, poff = p
    comp = pnco if pnco else 3
    vals = rd(poff, nv * comp, W.get(pty, 4), 'f')
    X = vals[0::comp]
    Y = vals[1::comp]
    if ('Cells', 'connectivity') not in arrays or ('Cells', 'offsets') not in arrays:
        print('piece missing connectivity/offsets'); continue
    cty, cco, coff = arrays[('Cells', 'connectivity')]
    conn = rd(coff, nc * cco, W.get(cty, 4), 'i')
    oty, oco, ooff = arrays[('Cells', 'offsets')]
    vlen = rd(ooff, nc, W.get(oty, 8), 'i')
    starts = [0] + list(vlen[:-1])
    vr = None
    for name in ('vr',):
        if ('CellData', name) in arrays:
            ty, co, off = arrays[('CellData', name)]
            vr = rd(off, nc, W.get(ty, 4), 'f')
    if ('CellData', 'b1') not in arrays:
        continue
    bs = []
    for name in ('b1', 'b2', 'b3'):
        ty, co, off = arrays[('CellData', name)]
        bs.append(rd(off, nc, W.get(ty, 4), 'f'))
    for i in range(nc):
        s = starts[i]
        e = vlen[i]
        xs = 0.0
        ys = 0.0
        cnt = 0
        for nid in conn[s:e]:
            xs += X[nid]
            ys += Y[nid]
            cnt += 1
        if cnt == 0:
            continue
        r = math.sqrt((xs / cnt) ** 2 + (ys / cnt) ** 2)
        Bm = math.sqrt(bs[0][i] ** 2 + bs[1][i] ** 2 + bs[2][i] ** 2)
        v = vr[i] if vr is not None else float('nan')
        gvr.append(v)
        gB.append(Bm)
        if r > 215.0:
            outer.append((v, Bm, r))
        if r < 3.0:
            inner.append((v, Bm, r))


def st(a):
    a = sorted(a)
    n = len(a)
    if n == 0:
        return '(empty)'
    return 'min=%.5g med=%.5g max=%.5g' % (a[0], a[n // 2], a[-1])


print('global cells:', len(gB))
print('global vr :', st(gvr))
print('global |B|:', st(gB))
if outer:
    print('outer ring n=%d  r=%.3f..%.3f' % (len(outer), min(t[2] for t in outer), max(t[2] for t in outer)))
    print('  vr :', st([t[0] for t in outer]))
    print('  |B|:', st([t[1] for t in outer]))
if inner:
    print('inner r<3 n=%d' % len(inner))
    print('  vr :', st([t[0] for t in inner]))
    print('  |B|:', st([t[1] for t in inner]))
