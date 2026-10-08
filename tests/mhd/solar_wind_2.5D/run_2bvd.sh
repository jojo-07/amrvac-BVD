#!/bin/bash
# Two-stage BVD production run on the 2.5D solar-wind case.
# Data files use base name "xia_2bvd" to stay distinct from:
#   xia_bvd   (single-stage BVD)   xia_weno3 (WENO3)   xia (koren reference)
cd /home/amrvac_bvd/amrvac-BVD/tests/mhd/solar_wind_2.5D || exit 1
rm -f xia_2bvd.log
mpirun -np 8 ./amrvac -i xia2018_2bvd.par > 2bvd_prod.out 2>&1 &
MPID=$!
sleep 120
echo "== start check $(date '+%F %H:%M') =="
tail -2 xia_2bvd.log 2>/dev/null | cut -c1-100
echo "NaN=$(grep -c NaN xia_2bvd.log 2>/dev/null)"
k=0
while kill -0 $MPID 2>/dev/null; do
  sleep 600
  k=$((k+1))
  N=$(grep -c NaN xia_2bvd.log 2>/dev/null || true)
  if [ -n "$N" ] && [ "$N" != "0" ]; then
    echo "NAN_BVD2 rows=$N $(date '+%F %H:%M')"
    kill $(ps -eo pid,args | awk '/xia2018_2bvd[.]par/ && !/awk/ {print $1}') 2>/dev/null
    sleep 15
    break
  fi
  if [ $((k % 6)) -eq 0 ]; then
    echo "progress $(date '+%F %H:%M') t=$(tail -1 xia_2bvd.log | awk '{print $2}')"
  fi
done
wait $MPID
echo "BVD2_END $(date '+%F %H:%M')"
echo "rows=$(wc -l < xia_2bvd.log) NaN=$(grep -c NaN xia_2bvd.log)"
tail -1 xia_2bvd.log | cut -c1-110
