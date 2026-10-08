#!/usr/bin/env python3
"""打印 vac logo 快照的 rho 统计（min/max/NaN 数）。

用法: python3 snap_stats.py [glob1 glob2 ...]
      无参数时默认统计两级 BVD 的全部快照。
"""
import struct
import re
import sys
import glob


def stats(path):
    raw = open(path, 'rb').read()
    m = re.search(rb'<AppendedData[^>]*>(.*)</AppendedData>', raw, re.S)
    app = m.group(1)[m.group(1).find(b'_') + 1:]
    xml = raw[:m.start()]
    vals = []
    for pm in re.finditer(rb'<Piece NumberOfPoints="\s*(\d+)" NumberOfCells="\s*(\d+)">(.*?)</Piece>', xml, re.S):
        ncl = int(pm.group(2))
        body = pm.group(3)
        off = int(re.search(rb'Name="rho"[^>]*offset="(\d+)"', body).group(1))
        vals += struct.unpack_from('<%df' % ncl, app, off + 4)
    nan = sum(1 for v in vals if v != v)
    return len(vals), min(vals), max(vals), nan


if __name__ == '__main__':
    pats = sys.argv[1:] or ['vaclogo_bvd_ts*.vtu']
    for pat in pats:
        files = sorted(glob.glob(pat))
        if not files:
            print('%s: 无匹配文件' % pat)
            continue
        print('== %s （%d 个文件）==' % (pat, len(files)))
        for f in files:
            n, lo, hi, nan = stats(f)
            print('%-28s min=%.4f max=%.4f NaN=%d' % (f, lo, hi, nan))
