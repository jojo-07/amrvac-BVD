#!/usr/bin/env python3
"""检查 vac BVD 运行最终快照的 rho 数据是否有 NaN/越界。"""
import struct, re, glob

path = sorted(glob.glob('vaclogo_bvd*.vtu'))[-1]
raw = open(path, 'rb').read()
m = re.search(rb'<AppendedData[^>]*>(.*)</AppendedData>', raw, re.S)
app = m.group(1)[m.group(1).find(b'_')+1:]
xml = raw[:m.start()]
allvals = []
for pm in re.finditer(rb'<Piece NumberOfPoints="\s*(\d+)" NumberOfCells="\s*(\d+)">(.*?)</Piece>', xml, re.S):
    ncl = int(pm.group(2))
    body = pm.group(3)
    off = int(re.search(rb'Name="rho"[^>]*offset="(\d+)"', body).group(1))
    (sz,) = struct.unpack_from('<I', app, off)
    allvals += struct.unpack_from('<%df' % ncl, app, off+4)
nan = sum(1 for v in allvals if v != v)
print('快照: %s' % path)
print('rho: 总单元=%d  min=%.4f  max=%.4f  NaN数=%d' % (len(allvals), min(allvals), max(allvals), nan))
