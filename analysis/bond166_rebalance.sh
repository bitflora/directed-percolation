#!/bin/bash
# Finish bond166 (S, mu01, mu02) sooner: one 15-bit verification prime (u16
# build) now, then u32 prime index 2 in its own process; stop P after its
# second prime (indices 0 and 1). Q keeps indices 4-6. Index 3 is not needed:
# six 31-bit primes (186 bits) reconstruct, the 15-bit prime verifies.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
cnt() { grep -c '^prime' "$1" 2>/dev/null || echo 0; }
log "start bond166 R1 (u16 prime 0, --no-x)"
results/dpenum_v12u16.exe -n 166 --bond --no-x --threads 1 --primes 1 --first-prime 0 -v -o results/runs/bond166_r.jsonl > results/logs/bond166_r1.log 2>&1
log "start bond166 R2 (u32 prime 2, --no-x)"
results/dpenum_v12u32.exe -n 166 --bond --no-x --threads 1 --primes 1 --first-prime 2 -v -o results/runs/bond166_r.jsonl > results/logs/bond166_r2.log 2>&1 &
until [ "$(cnt results/logs/bond166_p.log)" -ge 2 ]; do sleep 20; done
pid=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='dpenum_v11u32.exe'\" | Where-Object { \$_.CommandLine -like '*--primes 4 --first-prime 0*' }).ProcessId")
log "stopping bond166 P after primes 0-1 (pid $pid)"
[ -n "$pid" ] && powershell -NoProfile -Command "Stop-Process -Id $pid -Force"
wait
log "bond166 rebalance done"
