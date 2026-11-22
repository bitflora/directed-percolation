#!/bin/bash
# Overnight queue after the order-146 bond runs: bond to 156 with 31-bit primes
# (7 primes in two processes) and the remaining two site-136 primes.
cd "$(dirname "$0")/.."
U=results/dpenum_v8u32.exe; V=results/dpenum_v8.exe
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
waitfor() { # file pattern count
  until [ "$(grep -c '"prime"' "$1" 2>/dev/null)" -ge "$2" ]; do sleep 20; done
}
# process X starts when bond146 has 3 primes (process b done), Y when it has 4
( waitfor results/runs/bond146.jsonl 3; log "start bond156 X (primes 0-3)";
  $U -n 156 --bond --threads 1 --primes 4 --first-prime 0 -v -o results/runs/bond156.jsonl > results/logs/bond156_x.log 2>&1;
  log "bond156 X done" ) &
( waitfor results/runs/bond146.jsonl 4; log "start bond156 Y (primes 4-6)";
  $U -n 156 --bond --threads 1 --primes 3 --first-prime 4 -v -o results/runs/bond156.jsonl > results/logs/bond156_y.log 2>&1;
  log "bond156 Y done" ) &
( waitfor results/runs/site136.jsonl 2; log "start site136 primes 2-3";
  $V -n 136 --site --threads 1 --primes 2 --first-prime 2 -v -o results/runs/site136.jsonl > results/logs/site136_b.log 2>&1;
  log "site136 b done" ) &
wait
log "queue done"
