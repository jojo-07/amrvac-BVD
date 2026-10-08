#!/usr/bin/env python3
"""对比 vac logo 各限制器最终快照的 rho 统计。"""
import struct, re, glob

def stats(path):
    raw = open(path, 'rb').read()
    m = re.search(rb'<AppendedData[^>]*>(.*)</AppendedData>', raw, re.S)
    app = m.group(1)[m.group(1).find(b'_')+1:]
    xml = raw[:m.start()]
    allv = []
    for pm in re.finditer(rb'<Piece NumberOfPoints="\s*(\d+)" NumberOfCells="\s*(\d+)">(.*?)</Piece>', xml, re.S):
        ncl = int(pm.group(2))
        body = pm.group(3)
        off = int(re.search(rb'Name="rho"[^>]*offset="(\d+)"', body).group(1))
        allv += struct.unpack_from('<%df' % ncl, app, off+4)
    nan = sum(1 for v in allv if v != v)
    return len(allv), min(allv), max(allv), nan

for name in ['bvd', 'weno3']:
    path = sorted(glob.glob('vaclogo_%s*.vtu' % name))[-1]
    n, lo, hi, nan = stats(path)
    print('%-6s 快照=%s  单元=%d  min=%.4f  max=%.4f  NaN=%d'
          % (name, path, n, lo, hi, nan))
