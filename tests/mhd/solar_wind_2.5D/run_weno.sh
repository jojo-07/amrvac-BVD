#!/bin/bash
# WENO5 production run on the 2.5D solar-wind case (xia2018 configuration).
# Data files use base name "xia_weno".
cd /home/amrvac_bvd/amrvac-BVD/tests/mhd/solar_wind_2.5D || exit 1
rm -f xia_weno.log
mpirun -np 8 ./amrvac -i xia2018_weno.par > weno_prod.out 2>&1 &
MPID=$!
sleep 120
echo "== start check $(date '+%F %H:%M') =="
tail -2 xia_weno.log 2>/dev/null | cut -c1-100
echo "NaN=$(grep -c NaN xia_weno.log 2>/dev/null)"
k=0
while kill -0 $MPID 2>/dev/null; do
  sleep 600
  k=$((k+1))
  N=$(grep -c NaN xia_weno.log 2>/dev/null || true)
  if [ -n "$N" ] && [ "$N" != "0" ]; then
    echo "NAN_WENO rows=$N $(date '+%F %H:%M')"
    kill $(ps -eo pid,args | awk '/xia2018_weno[.]par/ && !/awk/ {print $1}') 2>/dev/null
    sleep 15
    break
  fi
  if [ $((k % 6)) -eq 0 ]; then
    echo "progress $(date '+%F %H:%M') t=$(tail -1 xia_weno.log | awk '{print $2}')"
  fi
done
wait $MPID
echo "WENO_END $(date '+%F %H:%M')"
echo "rows=$(wc -l < xia_weno.log) NaN=$(grep -c NaN xia_weno.log)"
tail -1 xia_weno.log | cut -c1-110
