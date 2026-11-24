#!/bin/bash
# Wall series to order 130 with 31-bit residues and u128 keys (u32 + wide
# build), spilling layers above 256 MiB to disk. Coefficients reach about
# 126 bits, so five primes carry them and a sixth verifies. 6 threads.
#
# Pause: create $SPILL/PAUSE; the run stops at the next block boundary.
# Resume, after a pause, kill or reboot: run this script again. Finished
# primes are skipped and the current prime continues from its checkpoint.
cd "$(dirname "$0")/.."
SPILL=${SPILL:-C:/Users/timja/dpspill}
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
log "start wall130 (u32+wide, 6 primes, 6 threads, spill $SPILL)"
results/dpwall_v3_u32_wide.exe -n 130 --threads 6 --block 6 --shards 256 --primes 6 \
    --spill "$SPILL" -v -o results/runs/wall130.jsonl >> results/logs/wall130.log 2>&1
code=$?
if [ $code -eq 3 ]; then log "wall130 paused"; else log "wall130 done (exit $code)"; fi
