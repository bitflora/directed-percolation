#!/bin/bash
# Wall series (T, S, x and t moments) to order 110: three 62-bit primes (the
# third verifies), 6 threads. Each prime is appended to the output when done;
# resume with: analysis/wall110_queue.sh <first prime index>.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
FIRST=${1:-0}
COUNT=$((3 - FIRST))
log "start wall110 (u64, primes $FIRST..2, 6 threads)"
results/dpwall_v1.exe -n 110 --threads 6 --block 6 --first-prime "$FIRST" --primes "$COUNT" -v \
    -o results/runs/wall110.jsonl >> results/logs/wall110.log 2>&1
log "wall110 done (exit $?)"
