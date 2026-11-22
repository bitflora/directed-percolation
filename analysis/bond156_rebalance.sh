#!/bin/bash
# Rebalance the remaining bond156 primes (indices 2 and 6) onto fresh v11u32
# processes; stop the original processes once their current prime is written.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
E=results/dpenum_v11u32.exe
killpid() { powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='dpenum_v8u32.exe'\" | Where-Object { \$_.CommandLine -like '*$1*' } | ForEach-Object { Stop-Process -Id \$_.ProcessId -Force; \$_.ProcessId }"; }
log "start bond156 W (prime 2, v11u32)"
$E -n 156 --bond --threads 1 --primes 1 --first-prime 2 -v -o results/runs/bond156_w.jsonl > results/logs/bond156_w.log 2>&1 &
until [ "$(grep -c '^prime' results/logs/bond156_x.log)" -ge 2 ]; do sleep 20; done
log "stopping X after prime 1: $(killpid 'first-prime 0')"
log "start bond156 V (prime 6, v11u32)"
$E -n 156 --bond --threads 1 --primes 1 --first-prime 6 -v -o results/runs/bond156_v.jsonl > results/logs/bond156_v.log 2>&1 &
until [ "$(grep -c '^prime' results/logs/bond156_y.log)" -ge 2 ]; do sleep 20; done
log "stopping Y after prime 5: $(killpid 'first-prime 4')"
wait
log "bond156 rebalance done"
