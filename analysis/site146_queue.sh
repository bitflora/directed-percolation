#!/bin/bash
# Site series to order 146 (all moments, 7 31-bit primes in three processes),
# started as the order-166 bond processes finish.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
E=results/dpenum_v12u32.exe
cnt() { grep -c '^prime' "$1" 2>/dev/null || echo 0; }
until [ "$(cnt results/logs/bond166_q.log)" -ge 3 ] || [ "$(cnt results/logs/bond166_p.log)" -ge 4 ]; do sleep 20; done
log "start site146 S1 (primes 0-2)"
$E -n 146 --site --threads 1 --primes 3 --first-prime 0 -v -o results/runs/site146.jsonl > results/logs/site146_s1.log 2>&1 &
until [ "$(cnt results/logs/bond166_q.log)" -ge 3 ] && [ "$(cnt results/logs/bond166_p.log)" -ge 4 ]; do sleep 20; done
log "start site146 S2 (primes 3-4) and S3 (primes 5-6)"
$E -n 146 --site --threads 1 --primes 2 --first-prime 3 -v -o results/runs/site146.jsonl > results/logs/site146_s2.log 2>&1 &
$E -n 146 --site --threads 1 --primes 2 --first-prime 5 -v -o results/runs/site146.jsonl > results/logs/site146_s3.log 2>&1 &
wait
log "site146 done"
