#!/usr/bin/env bash
set -eu

default_build_directory=build
build_dir=${1:-${default_build_directory}}

if [[ ${build_dir} = -h || ${build_dir} == --help ]]; then
  cat <<EOF
usage: ./build.sh
usage: ./build.sh other-build-directory
usage: ./build.sh other-build-directory Debug
usage: ./build.sh other-build-directory Debug notest

- The first argument is a build directory, by default ${default_build_directory}
- The second argument is the build type, by default Release.
  It must be one of { Release, Debug }.
- When the third argument is different from 'no', it will run the
  tests in addition.

EOF
  exit 0
fi

# this could be Release or Debug
build_type=${2:-Release}

if [[ "${build_type}" != Release && "${build_type}" != Debug ]]; then
  echo "Build type must be one of {Release, Debug}"
  exit 1
fi

notest=${3:-no}

mkdir -p "${build_dir}"
root_dir="$(realpath "$0" | xargs dirname)"
cd "${build_dir}"

if command -v mpiicpc >/dev/null 2>&1; then
  CMAKE_CXX_COMPILER=$(command -v mpiicpc)
else
  CMAKE_CXX_COMPILER=$(command -v mpicxx)
fi

cat <<EOF


  BUILDING CC4S in   '${build_dir}'
  CMAKE_BUILD_TYPE = '${build_type}'
  Running tests    = '${notest}'


EOF

cmake "${root_dir}" \
  -DCMAKE_BUILD_TYPE="${build_type}" \
  -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
make -j "$(nproc --all)"

# run the tests
if [[ "${notest}" = no ]]; then
  ctest -V
fi
