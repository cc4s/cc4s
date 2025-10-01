#!/usr/bin/env bash
set -eu

# usage: ./build.sh
# usage: ./build.sh other-build-directory
# usage: ./build.sh other-build-directory Debug
# usage: ./build.sh other-build-directory Debug notest

build_dir=${1:-build}
# this could be Release or Debug
build_type=${2:-Release}
notest=${3:-no}

mkdir -p ${build_dir}
cd ${build_dir}

if $(which mpiicpc 2>&1 >/dev/null); then
  CMAKE_CXX_COMPILER=$(which mpiicpc)
else
  CMAKE_CXX_COMPILER=$(which mpicxx)
fi

cat <<EOF


  BUILDING CC4S in '${build_dir}'


EOF

cmake .. \
  -DCMAKE_BUILD_TYPE=${build_type} \
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
make -j $(nproc --all)

# run the tests
if [[ ! ${notest} = no ]]; then
  ctest -V
fi
