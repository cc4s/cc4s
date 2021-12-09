#!/usr/bin/env bash
set -eu

CONFIGS=( `echo ${@%.mk} | xargs -n1 basename` )

for config in ${CONFIGS[@]}; do
  env="etc/env/${ENVIRONMENT}/${config}.sh"
  cat <<EOF
            =======
        CONFIG = $config
EOF
  [[ -f "$env" ]] ||
    { echo -e "\n\n\tWARNING: $env does not exist\n\n"; continue; }

  module purge
  source "$env"
  make clean-all CONFIG=$config
  make clean CONFIG=$config
  make -sj CONFIG=$config extern
  make -sj CONFIG=$config
  {
    # you should have python and so on
    make -C test run CONFIG=$config CC4S_RUN='mpirun -np 48 ${CC4S_PATH}'
    # TODO: add check phase too
  }

done


