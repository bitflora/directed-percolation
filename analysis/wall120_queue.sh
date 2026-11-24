#!/bin/bash
# Wall series to order 120 with 31-bit residues (u32 build, compact u64 keys):
# coefficients reach about 116 bits, so four primes carry them and a fifth
# verifies. 6 threads. Each prime is appended to the output when done; resume
# with: analysis/wall120_queue.sh <first prime index>.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
FIRST=${1:-0}
COUNT=$((5 - FIRST))
log "start wall120 (u32, primes $FIRST..4, 6 threads)"
results/dpwall_v2_u32.exe -n 120 --threads 6 --block 6 --first-prime "$FIRST" --primes "$COUNT" -v \
    -o results/runs/wall120.jsonl >> results/logs/wall120.log 2>&1
log "wall120 done (exit $?)"
