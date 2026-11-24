#!/bin/sh
# F14 analyses at 110, 115 and 120 terms on the order-120 wall series (run in parallel).
cd "$(dirname "$0")/.."
PY=${PY:-analysis/.venv/Scripts/python}
S=results/series/wall120.json; K=bond_wall_120; A=results/analysis; L=results/logs
"$PY" analysis/wall_analysis.py $S $K --terms 110 115 120 --names T S Rt Rx --K 2 3 \
  --Lmax 9 --Lstep 3 --mstep 2 --xc 0.64470018 0.644700185 0.64470019 \
  -o $A/wall120_biased.json > $L/wall120_biased.log 2>&1 &
"$PY" analysis/wall_analysis.py $S $K --terms 110 115 120 --names T S --K 1 2 3 4 \
  --Lmax 12 --Lstep 2 --xc 0.644700185 \
  -o $A/wall120_biased_robust.json > $L/wall120_robust.log 2>&1 &
"$PY" analysis/wall_analysis.py $S $K --terms 110 115 120 --names T S Rt Rx --K 2 3 \
  -o $A/wall120_unbiased.json > $L/wall120_unbiased.log 2>&1 &
wait
echo "all done $(date)"
