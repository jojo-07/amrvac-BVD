#!/usr/bin/env python3
"""Compare two AMRVAC .log files at selected times (max relative difference).

Usage: python3 cmp_logs.py A.log B.log [t1 t2 ...]
Columns: 0=it 1=time 2=dt 3..10=rho m1 m2 m3 e b1 b2 b3 ...
"""
import sys


def load(fn):
    rows = {}
    for line in open(fn):
        p = line.split()
        if len(p) < 12:
            continue
        try:
            t = float(p[1])
        except ValueError:
            continue
        rows[t] = p
    return rows


def at(rows, t, tol):
    if not rows:
        return None
    best = min(rows, key=lambda x: abs(x - t))
    return rows[best] if abs(best - t) < tol else None


def main(fa, fb, ts):
    A = load(fa)
    B = load(fb)
    print(f"A={fa}  tmax={max(A):.2f}")
    print(f"B={fb}  tmax={max(B):.2f}")
    for t in ts:
        tol = 1e-12 if t == 0.0 else 0.01
        ra, rb = at(A, t, tol), at(B, t, tol)
        if ra and rb:
            mx, who = 0.0, -1
            for i in range(3, 11):
                a, b = float(ra[i]), float(rb[i])
                if max(abs(a), abs(b)) > 1e-30:
                    r = abs(a - b) / max(abs(b), 1e-30)
                    if r > mx:
                        mx, who = r, i
            print(f"t={t:<6}: max rel diff = {mx:.3e} (col {who})")
        else:
            print(f"t={t:<6}: row missing (A={bool(ra)}, B={bool(rb)})")


if __name__ == "__main__":
    fa, fb = sys.argv[1], sys.argv[2]
    ts = [float(x) for x in sys.argv[3:]] if len(sys.argv) > 3 else [0.0, 1.0, 5.0, 9.0, 10.0]
    main(fa, fb, ts)
