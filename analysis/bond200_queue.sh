#!/bin/bash
# Bond series to order 200, all four moments, 31-bit residues (u32 build),
# out of core. The largest coefficients (mu02) reach about 193 bits, so seven
# primes carry them and an eighth verifies. 768 shards keep each shard in
# flight near 90 MiB at the projected 32 GB peak layer. 6 threads.
#
# Pause: create $SPILL/PAUSE; the run stops at the next phase boundary.
# Resume, after a pause, kill or reboot: run this script again. Finished
# primes are skipped and the current prime continues from its checkpoint.
cd "$(dirname "$0")/.."
SPILL=${SPILL:-C:/Users/timja/dpspill_bond}
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
log "start bond200 (u32, all moments, 8 primes, 6 threads, spill $SPILL)"
results/dpenum_v15_u32.exe -n 200 --bond --threads 6 --shards 768 --primes 8 \
    --spill "$SPILL" -v -o results/runs/bond200.jsonl >> results/logs/bond200.log 2>&1
code=$?
if [ $code -eq 3 ]; then log "bond200 paused"; else log "bond200 done (exit $code)"; fi
