#!/bin/bash

#TODO: need two modes. dry and calculation mode
# we work only on the calculation mode

start_time=$(grep -m1 -E 'Start calculation|Start dry-Run' cc4s.log | cut -d: -f1)

ranks=$(grep -m1 -E 'Start calculation|Start dry-Run' cc4s.log | awk '
  /employing [0-9]+ ranks/ {
    match($0, /employing ([0-9]+) ranks/, arr)
    print arr[1]
    exit
  }
')


(
  if grep -q "Start calculation" cc4s.log; then
    sed -n '/Start calculation/,$p' cc4s.log
  else
    sed -n '/Start dry-Run/,$p' cc4s.log
  fi
) | grep "Contraction:" | awk  -F'Contraction:' '{print $2}' |
  awk -v ranks="$ranks" -F'\\|' '
function _sqrt(x) { return x<0 ? 0 : sqrt(x) }
{
  key = $1
  perf = $2


  if (match(perf, /([0-9.]+) GF, ([0-9.]+) s/, arr)) {
    gf = arr[1] + 0
    t = arr[2] + 0
    has_time = 1
  } else if (match(perf, /([0-9.]+) GF/, arr)) {
    gf = arr[1] + 0
    t = 0
    has_time = 0
  } else {
    gf = 0; t = 0; has_time = 0
  }

  count[key]++
  if (!(key in gf_val)) gf_val[key] = gf

  if (has_time) {
    sum_t[key] += t
    sum_t2[key] += t*t
  }
  if (length(key) > maxlen) maxlen = length(key)  # Track longest key
  has_time_arr[key] = has_time
}
END {
  for (k in count) {
    n = count[k]

    if (has_time_arr[k]) {
      avg_t = sum_t[k]/n
      sd_t = _sqrt(sum_t2[k]/n - avg_t*avg_t)
      printf "%f %-*s|  %.6f GF, %.6f ± %.6f s (n=%d) -> Total %.6f GF, %.6f s, %.6f GF/s/rank\n",
             gf_val[k]*n, maxlen+2, k, gf_val[k], avg_t, sd_t, n, gf_val[k]*n, avg_t *n, gf_val[k]/avg_t/ranks
    } else {
      # no timing info for dry-run lines
      printf "%f %-*s|  %.6f GF (n=%d) -> Total %.6f GF\n",
              gf_val[k]*n, maxlen+2, k, gf_val[k], n, gf_val[k] *n
    }
  }
}' | sort -nk1 | cut -d' ' -f2-

