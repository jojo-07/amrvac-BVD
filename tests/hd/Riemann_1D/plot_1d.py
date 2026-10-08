#!/usr/bin/env python3
"""plot_1d.py —— 纯 Python 解析 AMRVAC 1D vtu（appended 格式）并生成 SVG 对比图。

用法:
    plot_1d.py [输出.svg] 变量名... 标签=vtu文件 标签=vtu文件 ...

例（在 tests/hd/Riemann_1D 下）:
    ./plot_1d.py sod_compare.svg rho v1 p \
        bvd=rm_1d_3step_hllc_bvdviz0000.vtu \
        ppm=rm_1d_3step_hllc_ppm0000.vtu

无 matplotlib/gnuplot 依赖，生成的 SVG 可直接用浏览器/图片查看器打开。
"""
import struct, re, sys, math

# ---------- vtu 解析（appended raw 格式） ----------

def read_vtu(path):
    raw = open(path, 'rb').read()
    m = re.search(rb'<AppendedData[^>]*>(.*)</AppendedData>', raw, re.S)
    blob = m.group(1)
    app = blob[blob.find(b'_')+1:]
    return raw[:m.start()], app

def grab(app, off, n, fmt):
    (sz,) = struct.unpack_from('<I', app, off)
    assert sz == n * struct.calcsize(fmt), f"size mismatch {sz}"
    return struct.unpack_from('<%d%s' % (n, fmt), app, off+4)

def load_1d(path, varnames):
    """返回 {变量名: [值...]}, x: [坐标...]"""
    xml, app = read_vtu(path)
    xs, data = [], {v: [] for v in varnames}
    for pm in re.finditer(
            rb'<Piece NumberOfPoints="\s*(\d+)" NumberOfCells="\s*(\d+)">(.*?)</Piece>',
            xml, re.S):
        npt, ncl = int(pm.group(1)), int(pm.group(2))
        body = pm.group(3)
        offs = {}
        for v in varnames:
            mm = re.search(rb'Name="%s"[^>]*offset="(\d+)"' % v.encode(), body)
            if mm:
                offs[v] = int(mm.group(1))
        pts_off = int(re.search(rb'NumberOfComponents="3"[^>]*offset="(\d+)"',
                                body).group(1))
        pts = grab(app, pts_off, 3*npt, 'f')
        xs += [(pts[3*i] + pts[3*i+3])/2 for i in range(npt-1)]
        for v in varnames:
            if v in offs:
                data[v] += list(grab(app, offs[v], ncl, 'f'))
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    return [xs[i] for i in order], {v: [data[v][i] for i in order]
                                    for v in varnames if v in data}

# ---------- SVG 绘图 ----------

PALETTE = ['#d62728', '#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd', '#8c564b']
W, H = 800, 520
ML, MR, MT, MB = 70, 16, 10, 42   # 边距

def panel(svg, x, ys, labels, title, ymin, ymax, xmin, xmax, px, py, pw, ph):
    def X(v):  return px + (v - xmin) / (xmax - xmin) * pw
    def Y(v):  return py + ph - (v - ymin) / (ymax - ymin) * ph
    # 轴与刻度
    svg.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" '
               f'fill="none" stroke="#999" stroke-width="1"/>')
    for i in range(6):
        xt = xmin + (xmax - xmin) * i / 5
        yt = ymin + (ymax - ymin) * i / 5
        svg.append(f'<line x1="{X(xt)}" y1="{py}" x2="{X(xt)}" y2="{py+4}" '
                   f'stroke="#666"/>')
        svg.append(f'<line x1="{px}" y1="{Y(yt)}" x2="{px+4}" y2="{Y(yt)}" '
                   f'stroke="#666"/>')
        svg.append(f'<text x="{X(xt)}" y="{py+ph+14}" font-size="11" '
                   f'text-anchor="middle" fill="#333">{xt:.2f}</text>')
        svg.append(f'<text x="{px-6}" y="{Y(yt)+4}" font-size="11" '
                   f'text-anchor="end" fill="#333">{yt:.3f}</text>')
    svg.append(f'<text x="{px+pw/2}" y="{py-4}" font-size="13" '
               f'text-anchor="middle" font-weight="bold">{title}</text>')
    # 曲线
    for k, (lab, arr) in enumerate(zip(labels, ys)):
        col = PALETTE[k % len(PALETTE)]
        pts = ' '.join(f'{X(x[i])},{Y(arr[i])}' for i in range(len(x)))
        svg.append(f'<polyline points="{pts}" fill="none" stroke="{col}" '
                   f'stroke-width="1.6"/>')
    return labels

def main():
    args = sys.argv[1:]
    if len(args) < 2 or '=' not in args[-1]:
        print(__doc__)
        sys.exit(1)
    out = args[0]
    varnames = args[1].split(',')
    series = []          # (label, path)
    for a in args[2:]:
        lab, path = a.split('=', 1)
        series.append((lab, path))

    x0, d0 = load_1d(series[0][1], varnames)
    xmin, xmax = min(x0), max(x0)
    labels = [l for l, _ in series]

    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" '
           f'height="{H}" viewBox="0 0 {W} {H}">',
           '<rect width="100%" height="100%" fill="white"/>']
    npan = len(varnames)
    gap = 12
    ph = (H - MT - MB - gap * (npan - 1)) / npan
    for j, v in enumerate(varnames):
        ys, ymin, ymax = [], None, None
        for _, path in series:
            x, d = load_1d(path, [v])
            arr = d[v]
            ys.append(arr)
            lo, hi = min(arr), max(arr)
            ymin = lo if ymin is None else min(ymin, lo)
            ymax = hi if ymax is None else max(ymax, hi)
        pad = (ymax - ymin) * 0.06 or 0.1
        py = MT + j * (ph + gap)
        panel(svg, x0, ys, labels, v, ymin-pad, ymax+pad,
              xmin, xmax, ML, py, W-ML-MR, ph)
    # 图例
    lx = ML + 14
    ly = H - 12
    for k, lab in enumerate(labels):
        svg.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+22}" y2="{ly}" '
                   f'stroke="{PALETTE[k % len(PALETTE)]}" stroke-width="2.4"/>')
        svg.append(f'<text x="{lx+28}" y="{ly+4}" font-size="12" '
                   f'fill="#333">{lab}</text>')
        lx += 60 + len(lab) * 7
    svg.append('</svg>')
    with open(out, 'w') as f:
        f.write('\n'.join(svg))
    print(f'写入 {out}（面板: {varnames}，系列: {labels}）')

if __name__ == '__main__':
    main()
