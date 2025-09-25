cmake_minimum_required(VERSION 3.15)
project(Cc4s LANGUAGES CXX)

include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Add external dependencies
find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Using Intel C++ (classic icpc)")
else()
    message(STATUS "Using C++ ${CMAKE_CXX_COMPILER_ID}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Debug)
endif()


FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://gitlab.cc4s.org/cc4s/yaml-cpp
  GIT_TAG        c9460110e072df84b7dee3eb651f2ec5df75fb18
)
FetchContent_MakeAvailable(yaml-cpp)

FetchContent_Declare(
  ctf
  GIT_REPOSITORY https://gitlab.cc4s.org/cc4s/ctf
  GIT_TAG        01ea778f8883f872fc009dee5872069a3ffebe50
)
FetchContent_MakeAvailable(ctf)


find_package(LAPACK REQUIRED)

################################################################################
# MKL and ScaLAPACK

# If MKL_DIR is not set by the user, but MKLROOT exists, use it
if(NOT MKL_DIR AND DEFINED ENV{MKLROOT})
    set(MKL_DIR "$ENV{MKLROOT}/lib/cmake/mkl")
endif()

if(NOT DEFINED ENV{MKLROOT})
  find_library(SCALAPACK_LIB
      NAMES scalapack
      HINTS /usr/local/lib/
      REQUIRED)
  set(CC4S_SCALAPACK ${SCALAPACK_LIB})
else()
  set(CC4S_SCALAPACK)
endif()

# set(MKL_DIR MKLROOT)
# find_package(MKL CONFIG)
# set(CC4S_SCALAPACK MKL::ScaLAPACK)

# if (MKL_FOUND)
#   set(CC4S_SCALAPACK MKL::ScaLAPACK)
# else()
#   find_library(SCALAPACK_LIB
#       NAMES scalapack
#       HINTS /usr/local/lib/
#       REQUIRED)
#   set(CC4S_SCALAPACK ${SCALAPACK_LIB})
# endif()

################################################################################
# Executable

# Main executable
add_executable(Cc4s
  src/main/Cc4s.cxx
  src/main/GC.cxx
  src/main/Options.cxx
  src/main/Writer.cxx
  src/main/Reader.cxx
  src/main/Setting.cxx
  src/main/Log.cxx
  src/main/Timer.cxx
  src/main/TensorIo.cxx
  src/main/TensorSet.cxx
  src/main/tcc/Tcc.cxx
  src/main/engines/DryTensor.cxx

  # Mixers (don't change the order)
  src/main/mixers/Mixer.cxx
  src/main/mixers/LinearMixer.cxx
  src/main/mixers/DiisMixer.cxx

  #
  # Put Algorithm.cxx as the first file of all algorithms
  #
  src/main/algorithms/Algorithm.cxx
  src/main/algorithms/BasisSetCorrection.cxx
  src/main/algorithms/CcsdChannels.cxx
  src/main/algorithms/CcsdtNonit.cxx
  src/main/algorithms/CcsdTTT.cxx
  # Warning:
  # Always put CoupledClusterMethod.cxx before all other
  # files that depend on it to declare statically the methods!
  # Otherwise the intel linker will not find it.
  src/main/algorithms/coupledcluster/CoupledClusterMethod.cxx
  src/main/algorithms/coupledcluster/Ccsd.cxx
  src/main/algorithms/coupledcluster/CcsdReference.cxx
  src/main/algorithms/coupledcluster/Ccsdt.cxx
  src/main/algorithms/CoupledCluster.cxx
  src/main/algorithms/coupledcluster/Drccd.cxx
  src/main/algorithms/coupledcluster/LambdaCcsdReference.cxx

  src/main/algorithms/DefineHolesAndParticles.cxx
  src/main/algorithms/DimensionProperty.cxx
  src/main/algorithms/FiniteSizeCorrection.cxx
  src/main/algorithms/MrccReader.cxx
  src/main/algorithms/NonZeroCondition.cxx
  src/main/algorithms/PerturbativeCT.cxx
  src/main/algorithms/PerturbativeTriples.cxx
  src/main/algorithms/PerturbativeTriplesReference.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Atrip.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Blas.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Complex.cxx
  src/main/algorithms/perturbative-triples/src/atrip/DatabaseCommunicator.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Equations.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Mpi.cxx
  src/main/algorithms/perturbative-triples/src/atrip/RankMap.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Reader.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Slice.cxx
  src/main/algorithms/perturbative-triples/src/atrip/SliceUnion.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Tuples.cxx
  src/main/algorithms/SecondOrderPerturbationTheory.cxx
  src/main/algorithms/SliceOperator.cxx
  src/main/algorithms/TensorSetToTensor.cxx
  src/main/algorithms/UegVertexGenerator.cxx
  src/main/algorithms/VertexCoulombIntegrals.cxx

  src/main/algorithms/Read.cxx
  src/main/algorithms/Write.cxx
)


# Include directories within the source code
target_include_directories(Cc4s PRIVATE ${CMAKE_SOURCE_DIR}/src/main)
target_include_directories(Cc4s PRIVATE ${CMAKE_SOURCE_DIR}/src/main/algorithms/perturbative-triples/src/)

################################################################################
# Several versions and compilation date

string(TIMESTAMP CC4S_DATE)
execute_process(COMMAND bash -c "mpirun --version | head -n1"
  OUTPUT_VARIABLE CC4S_MPI_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND bash -c "git describe --all --dirty --long"
  OUTPUT_VARIABLE CC4S_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND bash -c "${CMAKE_CXX_COMPILER} --version | head -n1"
  OUTPUT_VARIABLE CC4S_COMPILER_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE)

target_compile_definitions(Cc4s
  PRIVATE
  COMPILER_VERSION="${CC4S_COMPILER_VERSION}"
  CC4S_MPI_VERSION="${CC4S_MPI_VERSION}"
  CC4S_DATE="${CC4S_DATE}"
  CC4S_VERSION="${CC4S_VERSION}"

  _POSIX_C_SOURCE=200112L
  __STDC_LIMIT_MACROS
  FTN_UNDERSCORE=1
)

target_compile_options(Cc4s
  PRIVATE
  -g
  $<$<CONFIG:Debug>:-O0>
  $<$<CONFIG:Release>:-O3>
  -Wall
  -pedantic
  --all-warnings
  -fmax-errors=1
)

# if (CMAKE_BUILD_TYPE EQUAL Debug)
#   target_compile_options(Cc4s
#     PRIVATE
#     -O0
#   )
# else()
#   target_compile_options(Cc4s
#     PRIVATE
#     -O3
#   )
# endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(Cc4s PRIVATE
    -Wno-vla -Wno-int-in-bool-context)
endif()

# Link against external dependencies
target_link_libraries(Cc4s
  PRIVATE
  yaml-cpp ctf
  OpenMP::OpenMP_CXX
  ${CC4S_SCALAPACK}
  ${LAPACK_LIBRARIES}
  quadmath
)


message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("CMAKE_BUILD_TYPE   = ${CMAKE_BUILD_TYPE}")
message("yaml-cpp           = ${yaml-cpp}")
message("ctf                = ${ctf}")
message("CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("SCALAPACK_LIB      = ${CC4S_SCALAPACK}")
message("LAPACK_LIBRARIES   = ${LAPACK_LIBRARIES}")
