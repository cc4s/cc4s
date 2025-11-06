#!/usr/bin/env bash

set -e

if [[ -n "${1}" ]]; then
  compilers=(
    "${1}"
  )
else
  compilers=(
    clang_12 clang_13 clang_14 clang_15 clang_16 clang_17 clang_18 clang_19 clang_20
    gcc9 gcc10 gcc11 gcc12 gcc13 gcc14
  )
fi

for compiler in ${compilers[@]}; do
  build_dir="nix-build/${compiler}"
  if [[ -f ${build_dir}/Makefile ]]; then
    nix-shell etc/env/nix/shell.nix \
      --argstr compiler $compiler \
      --pure \
      --run "cd ${build_dir}; make -j1"
  else
    nix-shell etc/env/nix/shell.nix \
      --argstr compiler $compiler \
      --pure \
      --run "./build.sh nix-build/${compiler} Debug yes"
  fi
done
