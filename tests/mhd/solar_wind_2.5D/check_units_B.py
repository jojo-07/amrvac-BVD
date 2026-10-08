#!/usr/bin/env python3
"""Extract |B| statistics from AMRVAC appended-binary .vtu snapshots.

Usage: python3 check_units_B.py [snapshot.vtu]   (default: latest xia_weno3*.vtu)
"""
import glob
import math
import re
import struct
import sys

fn = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('xia_weno3*.vtu'))[-1]
raw = open(fn, 'rb').read()
m = re.search(rb'<AppendedData[^>]*>(.*)</AppendedData>', raw, re.S)
app = m.group(1)
app = app[app.find(b'_') + 1:]
xml = raw[:m.start()]
names = sorted({n.decode() for n in re.findall(rb'Name="([^"]+)"', xml)})
print('file:', fn)
print('arrays:', names)

B = []
for pm in re.finditer(rb'<Piece[^>]*NumberOfCells="\s*(\d+)"[^>]*>(.*?)</Piece>', xml, re.S):
    nc = int(pm.group(1))
    body = pm.group(2)
    vals = []
    for nm in (b'b1', b'b2', b'b3'):
        om = re.search(rb'Name="' + nm + rb'"[^>]*offset="(\d+)"', body)
        if om:
            vals.append(struct.unpack_from('<%df' % nc, app, int(om.group(1)) + 4))
    if len(vals) == 3:
        B.extend(math.sqrt(x * x + y * y + z * z) for x, y, z in zip(*vals))

B.sort()
n = len(B)
print('cells:', n)
if n:
    print('|B| min / med / max: %.6g / %.6g / %.6g' % (B[0], B[n // 2], B[-1]))
    print('|B| p99: %.6g' % B[int(n * 0.99)])
