#!/bin/bash
# Run bond156 prime index 3 in its own process once site136 is complete, and stop
# the original 4-prime process (X) after it has written its third prime.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
until [ "$(grep -c '"prime"' results/runs/site136.jsonl 2>/dev/null)" -ge 4 ]; do sleep 30; done
log "start bond156 Z (prime 3, v10u32)"
results/dpenum_v10u32.exe -n 156 --bond --threads 1 --primes 1 --first-prime 3 -v -o results/runs/bond156_z.jsonl > results/logs/bond156_z.log 2>&1 &
# X writes primes 0,1,2 then starts 3; stop it once the third X prime is written
until [ "$(grep -c '^prime' results/logs/bond156_x.log 2>/dev/null)" -ge 3 ]; do sleep 30; done
pid=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='dpenum_v8u32.exe'\" | Where-Object { \$_.CommandLine -like '*--primes 4*' }).ProcessId")
log "stopping bond156 X (pid $pid) after its third prime"
[ -n "$pid" ] && powershell -NoProfile -Command "Stop-Process -Id $pid -Force"
wait
log "bond156 Z done"
