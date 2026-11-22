#!/bin/bash
# S, mu01, mu02 to order 166 (no mu20), 7 31-bit primes in two processes.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
E=results/dpenum_v11u32.exe
log "start bond166 P (primes 0-3, --no-x)"
$E -n 166 --bond --no-x --threads 1 --primes 4 --first-prime 0 -v -o results/runs/bond166.jsonl > results/logs/bond166_p.log 2>&1 &
until [ "$(grep -c '"prime"' results/runs/bond156_v.jsonl 2>/dev/null)" -ge 1 ]; do sleep 20; done
log "start bond166 Q (primes 4-6, --no-x)"
$E -n 166 --bond --no-x --threads 1 --primes 3 --first-prime 4 -v -o results/runs/bond166.jsonl > results/logs/bond166_q.log 2>&1 &
wait
log "bond166 done"
