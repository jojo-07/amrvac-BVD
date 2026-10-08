#!/bin/bash
# BVD 收敛阶测试：sin 平流，5 档分辨率 × 3 种限制器
# 用法: ./conv_test.sh  （在 tests/rho/analytic_test_1d 目录下）
set -e
export AMRVAC_DIR=/home/amrvac_bvd/amrvac-BVD
cd /home/amrvac_bvd/amrvac-BVD/tests/rho/analytic_test_1d

# 一次性设置（幂等）
if [ ! -f makefile ]; then
  $AMRVAC_DIR/setup.pl -d=1 -arch=default
fi
make -j4 2>&1 | tail -1

for lim in bvd ppm minmod; do
  for N in 16 32 64 128 256; do
    cat > conv_${lim}_${N}.par << EOF
 &methodlist
        time_stepper='threestep'
        flux_scheme='hll'
        limiter='${lim}'
 /
 &meshlist
        refine_max_level=1
        domain_nx1=${N}
        iprob=1
 /
 &stoplist
        time_max=1.0d0
 /
EOF
    echo "=== ${lim} N=${N} ==="
    mpirun -np 1 ./amrvac -i amrvac.par conv_${lim}_${N}.par 2>&1 | grep "CONV_TEST" | tail -1
  done
done
