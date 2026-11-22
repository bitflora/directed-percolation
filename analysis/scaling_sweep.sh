#!/bin/bash
# Engine scaling sweep: peak states, peak MiB and time against order n (one prime, all moments).
E=${E:-enum/target/release/dpenum.exe}
out=results/scaling.csv
echo "model,n,peak_states,peak_mib,seconds,bound_calls,state_steps" > $out
for model in bond site; do
  for n in 60 70 80 90 100 110 120; do
    line=$($E -n $n --$model --threads 1 -v 2>&1 >/dev/null | grep "^done:")
    ps=$(echo "$line" | sed -E 's/.*peak states ([0-9]+).*/\1/')
    mib=$(echo "$line" | sed -E 's/.*peak MiB ([0-9.]+).*/\1/')
    sec=$(echo "$line" | sed -E 's/.*MiB [0-9.]+ ([0-9.]+)s.*/\1/')
    bc=$(echo "$line" | sed -E 's/.*bound calls ([0-9]+).*/\1/')
    ss=$(echo "$line" | sed -E 's/.*state-steps ([0-9]+).*/\1/')
    echo "$model,$n,$ps,$mib,$sec,$bc,$ss" >> $out
  done
done
