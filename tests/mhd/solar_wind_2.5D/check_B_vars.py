#!/usr/bin/env python3
"""Analyze B variables & background-field relation in an AMRVAC vtu snapshot.

Usage: python3 check_B_vars.py [file.vtu]

Reports:
  1) array names (first piece) and coordinate ranges; render convention guess
  2) max|b1-br|, max|b2-bth|, max|b3-bphi| over all cells
  3) least-squares fit of (b1-br, b2-bth) to the analytic background dipole
     2.40858*(2*cos(th), sin(th))/r^3 under several orientation conventions
     (Cartesian pole +X/+Y/-X/-Y; raw r,theta). Tiny residual + k~2.41 means
     "state B = total (B0+B1), aux B = perturbation (B1)".
  4) |B| / vr statistics for inner (r=1-3), mid (9-11), outer (213-217)
     rings, plus |diff| = |state - aux| in each ring.
"""
import glob
import math
import re
import struct
import sys

fn = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('xia_weno3*.vtu'))[-1]
raw = open(fn, 'rb').read()
i0 = raw.find(b'<AppendedData')
i1 = raw.find(b'>', i0) + 1
i2 = raw.find(b'</AppendedData>', i1)
app = raw[i1:i2]
app = app[app.find(b'_') + 1:]
xml = raw[:i0].decode('ascii', 'ignore')
pieces = list(re.finditer(r'<Piece\b([^>]*)>(.*?)</Piece>', xml, re.S))
print('file:', fn, '| pieces:', len(pieces))

W = {'Float32': (4, 'f'), 'Float64': (8, 'd'),
     'Int32': (4, 'i'), 'Int64': (8, 'q'),
     'UInt32': (4, 'I'), 'UInt64': (8, 'Q'), 'UInt8': (1, 'B')}


def rdl(key, arr, n):
    ty, co, off = arr[key]
    w, ch = W.get(ty, (4, 'f'))
    return struct.unpack_from('<%d%s' % (n * co, ch), app, off + 4)


def rd_ints(key, arr, n, limit):
    """Read n integer values: try declared type, then int64/int32 variants;
    accept the first whose values are all in [0, limit)."""
    ty, co, off = arr[key]
    for ch in (W.get(ty, (4, 'i'))[1], 'q', 'i'):
        try:
            vals = struct.unpack_from('<%d%s' % (n * (co or 1), ch), app, off + 4)
        except struct.error:
            continue
        if all(isinstance(v, int) and 0 <= v < limit for v in vals):
            return vals
    return None


def parse_body(body):
    arr = {}
    for sm in re.finditer(r'<(Points|Cells|CellData)\b[^>]*>(.*?)</\1>', body, re.S):
        sec = sm.group(1)
        for dm in re.finditer(r'<DataArray\b([^>]*?)/?>', sm.group(2)):
            a = dm.group(1)
            nm = re.search(r'Name="([^"]*)"', a)
            of = re.search(r'offset\s*=\s*"\s*(\d+)"', a)
            ty = re.search(r'type="([^"]*)"', a)
            nco = re.search(r'NumberOfComponents\s*=\s*"\s*(\d+)"', a)
            if of:
                arr[(sec, nm.group(1) if nm else '')] = (
                    ty.group(1) if ty else 'Float32',
                    int(nco.group(1)) if nco else 1,
                    int(of.group(1)))
    return arr


xmin, xmax, ymin, ymax = 1e30, -1e30, 1e30, -1e30
cells = []          # (x, y, d1, d2, d3, |B|state, |B|aux, vr)
maxd = [0.0, 0.0, 0.0]
ncell = 0
badtot = 0

for pi, pm in enumerate(pieces):
    attrs, body = pm.group(1), pm.group(2)
    if pi == 0:
        print('piece0 attrs: %r' % attrs[:200])
    mc = re.search(r'NumberOfCells\s*=\s*"\s*(\d+)"', attrs)
    mp = re.search(r'NumberOfPoints\s*=\s*"\s*(\d+)"', attrs)
    if not mc or not mp:
        if pi < 3:
            print('piece %d attrs unparsed: %r' % (pi, attrs[:220]))
        continue
    nc = int(mc.group(1))
    nv = int(mp.group(1))
    arr = parse_body(body)
    if pi == 0:
        print('first piece: nc=%d nv=%d' % (nc, nv))
        print('cell arrays:', sorted(k[1] for k in arr if k[0] == 'CellData'))
    ptk = [k for k in arr if k[0] == 'Points']
    if (not ptk or ('Cells', 'connectivity') not in arr or
            ('Cells', 'offsets') not in arr or ('CellData', 'b1') not in arr):
        continue
    coP = arr[ptk[0]][1] or 3
    pts = rdl(ptk[0], arr, nv)
    xs = pts[0::coP]
    ys = pts[1::coP] if coP > 1 else [0.0] * nv
    xmin = min(xmin, min(xs))
    xmax = max(xmax, max(xs))
    ymin = min(ymin, min(ys))
    ymax = max(ymax, max(ys))
    vlen = rd_ints(('Cells', 'offsets'), arr, nc, 8 * nc)
    if vlen is None:
        print('piece %d: cannot decode offsets' % pi)
        continue
    total = int(vlen[-1])
    conn = rd_ints(('Cells', 'connectivity'), arr, total, nv)
    if conn is None or len(conn) < total:
        print('piece %d: cannot decode connectivity (total=%d nv=%d)' % (pi, total, nv))
        continue
    if pi == 0:
        print('offsets type=%s (last=%d), conn type=%s'
              % (arr[('Cells', 'offsets')][0], total, arr[('Cells', 'connectivity')][0]))
    b1 = rdl(('CellData', 'b1'), arr, nc)
    b2 = rdl(('CellData', 'b2'), arr, nc)
    b3 = rdl(('CellData', 'b3'), arr, nc)
    have_aux = all(('CellData', k) in arr for k in ('br', 'bth', 'bphi'))
    if have_aux:
        brl = rdl(('CellData', 'br'), arr, nc)
        btl = rdl(('CellData', 'bth'), arr, nc)
        bpl = rdl(('CellData', 'bphi'), arr, nc)
    if ('CellData', 'vr') in arr:
        vrl = rdl(('CellData', 'vr'), arr, nc)
    elif ('CellData', 'm1') in arr and ('CellData', 'rho') in arr:
        m1 = rdl(('CellData', 'm1'), arr, nc)
        rh = rdl(('CellData', 'rho'), arr, nc)
        vrl = [m1[i] / rh[i] if rh[i] != 0 else 0.0 for i in range(nc)]
    else:
        vrl = [0.0] * nc
    starts = [0] + list(vlen[:-1])
    bad = 0
    for i in range(nc):
        s = starts[i]
        e = vlen[i]
        if e <= s or e > len(conn):
            bad += 1
            continue
        try:
            cx = sum(xs[conn[j]] for j in range(s, e)) / (e - s)
            cy = sum(ys[conn[j]] for j in range(s, e)) / (e - s)
        except IndexError:
            bad += 1
            continue
        bm_s = math.sqrt(b1[i] ** 2 + b2[i] ** 2 + b3[i] ** 2)
        if have_aux:
            d1 = b1[i] - brl[i]
            d2 = b2[i] - btl[i]
            d3 = b3[i] - bpl[i]
            bm_a = math.sqrt(brl[i] ** 2 + btl[i] ** 2 + bpl[i] ** 2)
        else:
            d1 = d2 = d3 = 0.0
            bm_a = 0.0
        maxd[0] = max(maxd[0], abs(d1))
        maxd[1] = max(maxd[1], abs(d2))
        maxd[2] = max(maxd[2], abs(d3))
        cells.append((cx, cy, d1, d2, d3, bm_s, bm_a, vrl[i]))
        ncell += 1
    if bad:
        print('piece %d: skipped %d bad cells' % (pi, bad))
        badtot += bad

print('coord ranges: X %.6g..%.6g   Y %.6g..%.6g' % (xmin, xmax, ymin, ymax))
is_cart = ymax > 3.2
print('render: %s' % ('Cartesian' if is_cart else 'raw (r,theta)'))
print('cells scanned: %d (bad %d)' % (ncell, badtot))
print('max|b1-br| = %.6g   max|b2-bth| = %.6g   max|b3-bphi| = %.6g' % tuple(maxd))
if max(maxd) == 0.0:
    print('NOTE: state B and aux B are numerically identical in this file.')
allb = sorted(c[5] for c in cells)
print('global |B| (state): min=%.6g med=%.6g max=%.6g'
      % (allb[0], allb[len(allb) // 2], allb[-1]))


def fit(name, gen):
    num = 0.0
    den = 0.0
    pairs = []
    for c in cells:
        g = gen(c[0], c[1])
        if g is None:
            continue
        r, th = g
        if r < 1e-9 or r > 10.0:
            continue
        ar = 2.0 * math.cos(th) / r ** 3
        at = math.sin(th) / r ** 3
        num += c[2] * ar + c[3] * at
        den += ar * ar + at * at
        pairs.append((c[2], c[3], ar, at))
    if not pairs or den == 0.0:
        print('fit %-18s: no data' % name)
        return
    k = num / den
    R = max(math.hypot(d1 - k * ar, d2 - k * at) for d1, d2, ar, at in pairs)
    print('fit %-18s: k=%9.5f  max resid=%.3g  (n=%d)' % (name, k, R, len(pairs)))


print('dipole-relation fits on (b1-br, b2-bth), r<10:')
fit('cart pole +X', lambda x, y: (math.hypot(x, y), math.atan2(y, x)))
fit('cart pole +Y', lambda x, y: (math.hypot(x, y), math.pi / 2 - math.atan2(y, x)))
fit('cart pole -X', lambda x, y: (math.hypot(x, y), math.pi - math.atan2(y, x)))
fit('cart pole -Y', lambda x, y: (math.hypot(x, y), math.pi / 2 + math.atan2(y, x)))
fit('raw (r,theta)', lambda x, y: (x, y))


def ring(tag, rmin, rmax):
    sel = []
    for c in cells:
        r = math.hypot(c[0], c[1]) if is_cart else c[0]
        if rmin <= r <= rmax:
            sel.append(c)
    if not sel:
        print('  %-18s none' % tag)
        return
    n = len(sel)
    bs = sorted(c[5] for c in sel)
    ba = sorted(c[6] for c in sel)
    dv = sorted(c[7] for c in sel)
    dd = sorted(max(abs(c[2]), abs(c[3]), abs(c[4])) for c in sel)
    print('  %-18s n=%-6d |B|s med=%.5g max=%.5g | |B|a med=%.5g | vr med=%.5g max=%.5g | |diff| med=%.3g max=%.3g'
          % (tag, n, bs[n // 2], bs[-1], ba[n // 2], dv[n // 2], dv[-1], dd[n // 2], dd[-1]))


print('rings:')
ring('inner 1-3', 0.9, 3.0)
ring('mid 9-11', 9.0, 11.0)
ring('outer 213-217', 213.0, 217.0)
ring('outer 215-217', 215.0, 217.0)
