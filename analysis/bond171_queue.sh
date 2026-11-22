#!/bin/bash
# After site146 completes: bond S, mu01, mu02 to order 171 (Jensen 1999's
# length), 7 31-bit primes, one process with 6 threads (block-parallel).
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
until [ "$(grep -c '"prime"' results/runs/site146.jsonl 2>/dev/null)" -ge 7 ]; do sleep 20; done
log "start bond171 (--no-x, u32, 7 primes, 6 threads)"
results/dpenum_v12u32.exe -n 171 --bond --no-x --threads 6 --block 6 --primes 7 -v -o results/runs/bond171.jsonl > results/logs/bond171.log 2>&1
log "bond171 done"
