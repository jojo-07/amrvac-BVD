#!/usr/bin/env python3
"""收敛阶计算：sin 平流 t=1.0 的 L1 误差。"""
import math

bvd = [2.06441e-2, 2.65390e-3, 3.34049e-4, 4.18319e-5, 5.23474e-6]
ppm = [6.67946e-2, 1.57549e-2, 3.64743e-3, 7.27401e-4, 1.29144e-4]
mm  = [1.51404e-1, 5.94000e-2, 1.96601e-2, 5.70702e-3, 1.56342e-3]

for name, errs in [('bvd', bvd), ('ppm', ppm), ('minmod', mm)]:
    orders = [math.log(errs[i]/errs[i+1], 2) for i in range(len(errs)-1)]
    avg = sum(orders)/len(orders)
    print('%-7s 各段阶数: %s  平均=%.2f' % (name, [round(o,2) for o in orders], avg))

print()
print('N=256 L1 误差比  bvd:ppm:minmod = 1 : %.0f : %.0f' % (ppm[-1]/bvd[-1], mm[-1]/bvd[-1]))
