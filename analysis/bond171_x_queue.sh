#!/bin/bash
# Bond order 171 with all moments (adds mu20 to the verified S, mu01, mu02):
# twelve 15-bit primes (u16 build, about 2 GB heap per prime with the
# split-layer parallel steps), one process with 6 threads. Each prime is
# appended to the output as soon as it is done, so the queue can be resumed
# with --first-prime k after an interruption.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
FIRST=${1:-0}
COUNT=$((12 - FIRST))
log "start bond171_x (all moments, u16, primes $FIRST..11, 6 threads)"
results/dpenum_v14u16.exe -n 171 --bond --threads 6 --block 6 --first-prime "$FIRST" --primes "$COUNT" -v \
    -o results/runs/bond171_x.jsonl >> results/logs/bond171_x.log 2>&1
log "bond171_x done (exit $?)"
