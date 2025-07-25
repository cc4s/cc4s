#!/bin/bash

file="cc4s.log"
thresh=0.05  # displaying only memory allocation/frees for tensors which are significant - compared to the largest tensor in the calculation
algo="all"

# Function to display usage information
usage() {
  echo "Usage: $0 [-f file] [-t threshold] [-h]"
  echo "  -f file      Specify the log file to process (default: cc4s.log)"
  echo "  -t threshold Specify the threshold for significant memory allocations (default: 0.05)"
  echo "  -a algo-id   Specify which algo you want to have detailed output (default: all)"
  echo "  -h           Display this help message"
  exit 0
}


# Parse command-line arguments
while getopts "f:t:a:h" opt; do
  case $opt in
    f) file="$OPTARG" ;;
    t) thresh="$OPTARG" ;;
    a) algo="$OPTARG" ;;
    h) usage ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

if [[ "$algo" =~ ^[0-9]+$ || "$algo" == "all" ]]; then
  :
else
  echo "Error: algo must be an integer or 'all'"
  exit
fi


if [ ! -f "$file" ]; then
  echo "Error: File '$file' not found."
  exit 1
fi


largestTensor=0.0



tmpfile=$(mktemp -p .)
# Only pick relevant lines
grep -E "(Allocate tensor|Free tensor|Allocate object|Free object|begin step|end step|Start calculation|Start dry|Finished calculation)" "$file" | \
#grep -E "(Allocate tensor|Free tensor|begin step|end step|Start calculation|Start dry|Finished calculation)" "$file" | \
awk -F ':' \
' \
BEGIN {
  current = 0.0;
  hwm = 0.0;
  largest=0.0;
  largestName="xxx"
  running = 0
}
{

  time = $1;
  location = $2;
  event = $3;

  // watch out!! the line reads "begin step:..." which introduces a further seperator
  if (event == "begin step") {
    print "\n=== Step " $4 $5" ==="
    next;
  }

  if (event == "end step" || event == "Finished calculation") {
    printf("\n=== High water mark: %8.2f MB\n", hwm);
    next;
  }

  split($3, msg, " ");

  # if we are already in a dry run and hit Start calculaton we exit
  if (msg[1] == "Start" && (msg[2] == "calculation")  && running) {
    exit;
  }

  # Determine the number of processes.
  if (msg[1] == "Start" && (msg[2] ~ /^(calculation|dry-Run)$/)) {
    running = 1
    n_ranks = msg[5]
    next
  }
  if (msg[1] == "Allocate" || msg[1] == "Free") {
    tensor_name = msg[3];
    elements = msg[5]
    bytes = msg[7]
    if ( n_ranks < 1 ) {
      print "n_ranks not detected. Exit!"
      print n_ranks
      exit 1
    }
    size_mb = elements / (1024.0 * 1024.0 * n_ranks) * msg[7]
    if (size_mb > largest) { largest = size_mb; largestName = tensor_name }
    if (msg[1] == "Allocate") {
      action="+"
      current += size_mb;
      if (current > hwm) hwm = current;
    } else if (msg[1] == "Free") {
      action="-"
      current -= size_mb;
    }
    printf("%10.6f  %30s  %s%8.2f MB  (Total/HWM: %8.2f / %8.2f MB)\n", time, tensor_name, action, size_mb, current, hwm);
  }
}
END {
  printf("\n--\n=== Largest tensor in this calculation is %s of size %8.2f MB\n", largestName, largest)
}
' > $tmpfile


echo "==============="
echo "Summary"
echo "==============="

awk '{

  if ($0 ~ /^=/ || $0 ~ /^$/ || $0 ~ /^-/) {
    print $0;
  }
}' $tmpfile

nsteps=$(grep "=== Step" $tmpfile | wc -l)

#TODO:
# We want to let the user decide what to do now...
# either show all the data or only the results of one particular algorithm
#read -p "Enter your choice [1-3]: " choice

largest=$(tail -n1 $tmpfile | awk '{print $11}')
minimum=$(python3 -c "print($largest*$thresh)")


while true; do

prompt="

Declare if you want to see further details of the allocation.

Just hit return to procede.
Type e if you want to leave.
Type a if you want to change the algorithm(s) you want to see details from (current value: $algo)
Type m if you want to change the size of the smallest tensor shown (current value: $minimum)
"


read -p "$prompt" choice

#proceed
[ -z "$choice" ] && break
#exit
[[ "$choice" == e* ]] && rm $tmpfile && exit
#change algo
if [[ "$choice" == a* ]]; then
  read -p $'Provide new value for algo\n' nalgo
  [[ "$nalgo" =~ ^[0-9]+$ || "$nalgo" == "all" ]] && algo=$nalgo
  if (( algo > nsteps || algo < 1 )); then
  echo "==========================================="
  echo "===== Warning: algo should be chosen between 1 and $nsteps"
  echo "==========================================="
  fi
  continue;
fi
#change minimum
if [[ "$choice" == m* ]]; then
  read -p $'Provide new value for minimum\n' nmin
  [[ "$nmin" =~ ^([0-9]*[.])?[0-9]+$ ]] && minimum=$nmin
  continue;
fi

done


awk -v m="$minimum" -v algo="$algo" '
 BEGIN {
    show_all = (algo == "all")
  }

  $0 ~ /^=== Step[ \t]+[0-9]+/ {
    step = $3 + 0
    in_block = show_all || (step == algo)
  }
  in_block {
    if ($0 ~ /^=/ || $0 ~ /^$/ || $0 ~ /^-/) {
      print $0
    } else if ($4 > m) {
      print $0
    }
  }
' "$tmpfile"




rm $tmpfile
