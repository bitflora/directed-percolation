#!/bin/bash
# Wait for the running wall130 job to exit, then resume bond200, so the two
# never share the machine. bond200 starts only if wall130 actually finished
# (all 6 primes in results/runs/wall130.jsonl); after a pause or crash this
# logs that and exits without starting anything.
#
# Polls every 5 minutes. Not reboot-proof: after a reboot, rerun
# wall130_queue.sh and then this script.
cd "$(dirname "$0")/.."
log() { echo "[$(date '+%F %T')] $*" >> results/logs/queue.log; }
log "after_wall130: waiting for wall130 to exit, then bond200"
while tasklist //FI "IMAGENAME eq dpwall_v4_u32_wide.exe" 2>/dev/null | grep -q dpwall_v4_u32_wide; do
    sleep 300
done
lines=$(grep -c . results/runs/wall130.jsonl 2>/dev/null || echo 0)
if [ "$lines" -ge 6 ]; then
    log "after_wall130: wall130 has $lines primes; resuming bond200"
    exec bash analysis/bond200_queue.sh
else
    log "after_wall130: wall130 exited with $lines of 6 primes; bond200 not started"
fi
