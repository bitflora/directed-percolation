#!/bin/bash
# Site series to order 158 (Jensen 1999's site length), all moments:
# fourteen 15-bit primes (about 177-bit coefficients expected, so 13 primes
# reconstruct and the 14th verifies), u16 build, one process with 6 threads.
# Resume with: analysis/site158_queue.sh <first-prime-index>
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
FIRST=${1:-0}
COUNT=$((14 - FIRST))
log "start site158 (all moments, u16, primes $FIRST..13, 6 threads)"
results/dpenum_v14u16.exe -n 158 --site --threads 6 --block 6 --first-prime "$FIRST" --primes "$COUNT" -v \
    -o results/runs/site158.jsonl >> results/logs/site158.log 2>&1
log "site158 done (exit $?)"
