#!/bin/bash


SCANDRY=0
MODE=summary
STEP=
ALGO=
logfile=cc4s.log

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --dry             Scan the dry run's log.
  --algo <name>     Select algorithm to investigate
  --mode=<mode>     Set mode (default: summary)
  --step=<step>     Specify step to execute
  --file=<file>     Log file name (default: cc4s.log)
  -h, --help        Show this help message and exit

Examples:
  $(basename "$0") --mode=full --step=2 --file=run.log
  $(basename "$0") --dry --algo myalgo
EOF
}

# Parse command line arguments
for arg in "$@"; do
  case $arg in
    --dry)
      SCANDRY=true
      shift
      ;;
    --algo)
      ALGO=algo
      shift
      ;;
    --mode=*)
      MODE="${arg#--mode=}"
      ;;
    --step=*)
      STEP="${arg#--step=}"
      ;;
    --file=*)
      logfile="${arg#--file=}"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $arg"
      usage
      exit 1
      ;;
  esac
done

start_time=$(grep -m1 -E 'Start calculation|Start dry-Run' $logfile | cut -d: -f1)

ranks=$(grep -m1 -E 'Start calculation|Start dry-Run' $logfile | awk '
  /employing [0-9]+ ranks/ {
    match($0, /employing ([0-9]+) ranks/, arr)
    print arr[1]
    exit
  }
')



log_extract=$(
  if [ "$SCANDRY" = true ]; then
    sed -n '1,/Start calculation/p' $logfile
  else
    if grep -q "Start calculation" $logfile; then
      sed -n '/Start calculation/,$p' $logfile
    else
      sed -n '/Start dry-Run/,$p' $logfile
    fi
  fi
)

numberContractions=$(echo "$log_extract" | grep "Contraction:" | wc -l)
printf "Total number of contractions: %d\n\n" $numberContractions

if [ "$ALGO" = "algo" ]; then

  echo "$log_extract" | grep "end step:" | awk -v dry=$SCANDRY '
{
  match($0, /end step: *[0-9]+/, m1)
  match($0, /name: *[^,]+/, m2)
  match($0, /operations: *[0-9.eE+-]+/, m3)

  step = (m1[0] ~ /[0-9]+/ ? gensub(/.* ([0-9]+)$/, "\\1", "g", m1[0]) : "?")
  name = (m2[0] ~ /name: */ ? gensub(/name: */, "", "g", m2[0]) : "?")
  ops = (m3[0] ~ /operations: */ ? gensub(/operations: */, "", "g", m3[0]) : "?")

  printf("step: %s, name: %s, Operations: %s GFLOP", step, name, ops)
  if (!dry) {
    match($0, /realtime: *[0-9.eE+-]+/, m4)
    time = (m4[0] ~ /realtime: */ ? gensub(/realtime: */, "", "g", m4[0]) : "?")

    printf(", Realtime: %s", time)
  }
  printf("\n")
}'


  log_base="$log_extract"
  [ -z "$STEP" ] && read -p "Enter step number to analyze (e.g., 1; 0 or return to choose all algos): " STEP

  if [[ ! -z "$STEP" && "$STEP" != 0 ]]; then
    STEP_BEGIN="begin step: $STEP"
    STEP_END="end step: $STEP"
    log_extract=$(echo "$log_base" | awk -v step_begin="$STEP_BEGIN" -v step_end="$STEP_END" '
      $0 ~ step_begin { printing=1 }
      printing { print }
      $0 ~ step_end { printing=0 }
    ')
    echo "Analyzing step $STEP:"
  fi
fi





echo "$log_extract" | grep "Contraction:" | awk  -F'Contraction:' '{print $2}' |
  awk -v ranks="$ranks" -v mode="$MODE" -F'\\|' '
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

  if (mode != "summary") {
    print NR, $0
    total_gf+=gf
    if (has_time) total_sec+=t
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
  if (mode == "summary") {
    for (k in count) {
      n = count[k]
      total_flop += gf_val[k] * n
      if (has_time_arr[k]) {
        total_time += sum_t[k]
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

    # dirty hack: we sort afterwards...so we have to ensure that this is the last appearance

    if (total_time > 0) {
      printf "%f ---\n%f \n%f \n%fTotal Flops: %.6f GF, Total time: %.6f s, Overall: %.6f GF/s/ranks\n",
             1e12, 2e12, 3e12, 4e12, total_flop, total_time, total_flop / total_time / ranks
    }
  } else {
    if (total_sec > 0) {
      printf "%f Total Flops: %.6f GF, Total time: %.6f s, Overall: %.6f GF/s/ranks\no\n",
              5e12, total_gf, total_sec, total_gf / total_sec / ranks
    }
  }

}' | sort -nk1 | cut -d' ' -f2-

