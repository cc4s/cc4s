cmake_minimum_required(VERSION 3.15)
project(Cc4s LANGUAGES CXX)

include(FetchContent)

# Set C++ standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
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
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Building ${CMAKE_BUILD_TYPE} profile")

if(EXISTS "${CMAKE_SOURCE_DIR}/thirdparty/yaml-cpp")
  message(STATUS "Using local fallback for yaml-cpp")
  FetchContent_Declare(
    yaml-cpp
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/thirdparty/yaml-cpp"
  )
else()
  message(STATUS "Cloning yaml-cpp from the internet")
  set(YAML_CPP_GIT_TAG "c9460110e072df84b7dee3eb651f2ec5df75fb18")
  FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://gitlab.cc4s.org/cc4s/yaml-cpp
    GIT_TAG        ${YAML_CPP_GIT_TAG}
  )
endif()
FetchContent_MakeAvailable(yaml-cpp)

if(EXISTS "${CMAKE_SOURCE_DIR}/thirdparty/ctf")
  message(STATUS "Using local fallback for CTF")
  FetchContent_Declare(
    ctf
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/thirdparty/ctf"
  )
else()
  message(STATUS "Cloning CTF from the internet")
  set(CTF_GIT_TAG "01ea778f8883f872fc009dee5872069a3ffebe50")
  FetchContent_Declare(
    ctf
    GIT_REPOSITORY https://gitlab.cc4s.org/cc4s/ctf
    GIT_TAG        ${CTF_GIT_TAG}
  )
endif()
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
  src/main/StepCheck.cxx
  src/main/Options.cxx

  # Set writer and reader in front of all files that use its
  # functionality (TensorIO)
  src/main/Writer.cxx
  src/main/Reader.cxx
  src/main/TensorIo.cxx
  src/main/TensorSet.cxx

  src/main/Setting.cxx
  src/main/Log.cxx
  src/main/Timer.cxx
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
  src/main/algorithms/coupledcluster/LambdaCcsd.cxx

  src/main/algorithms/DefineHolesAndParticles.cxx
  src/main/algorithms/DimensionProperty.cxx
  src/main/algorithms/FiniteSizeCorrection.cxx
  src/main/algorithms/MrccReader.cxx
  src/main/algorithms/NonZeroCondition.cxx
  src/main/algorithms/PerturbativeCT.cxx
  src/main/algorithms/PerturbativeTriples.cxx
  src/main/algorithms/PerturbativeTriplesReference.cxx
  src/main/algorithms/SecondOrderPerturbationTheory.cxx
  src/main/algorithms/SliceOperator.cxx
  src/main/algorithms/TensorSetToTensor.cxx
  src/main/algorithms/UegVertexGenerator.cxx
  src/main/algorithms/VertexCoulombIntegrals.cxx

  src/main/algorithms/Read.cxx
  src/main/algorithms/Write.cxx

  # atrip source
  src/main/algorithms/perturbative-triples/src/atrip/Atrip.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Blas.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Complex.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Checkpoint.cxx
  src/main/algorithms/perturbative-triples/src/atrip/DatabaseCommunicator.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Equations.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Mpi.cxx
  src/main/algorithms/perturbative-triples/src/atrip/RankMap.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Reader.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Slice.cxx
  src/main/algorithms/perturbative-triples/src/atrip/SliceUnion.cxx
  src/main/algorithms/perturbative-triples/src/atrip/Tuples.cxx

)

set_source_files_properties(
    src/main/algorithms/coupledcluster/LambdaCcsdReference.cxx
    src/main/algorithms/coupledcluster/LambdaCcsd.cxx
    PROPERTIES
    COMPILE_FLAGS "-O0"
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
  CTF_SWITCH_REDISTRIBUTION
  CTF_LOWMEM_CONTRACTION
  CTF_MEMORY_LOGGER
  ATRIP_NO_OUTPUT
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

get_target_property(CC4S_FLAGS Cc4s COMPILE_OPTIONS)
get_target_property(CC4S_DEFS Cc4s COMPILE_DEFINITIONS)

# Enable testing
option(ENABLE_TESTING "Enable building/running tests" ON)
set(CTEST_OUTPUT_ON_FAILURE ON)
if(ENABLE_TESTING)
    # User can optionally override number of processes
    set(NP 4 CACHE STRING "Number of MPI processes")
    set(CC4S_PATH "${CMAKE_BINARY_DIR}/Cc4s")
    set(FULL_RUN_CMD "mpirun -np ${NP} ${CC4S_PATH}"
      CACHE STRING "Command to run CC4S")
    set(TEST_DIR ${CMAKE_BINARY_DIR}/test)
    message(STATUS "CC4S for integration tests = ${FULL_RUN_CMD}")

    enable_testing()

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    message(STATUS "Using Python3              = ${Python3_EXECUTABLE}")
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import yaml"
        RESULT_VARIABLE PY_YAML_FOUND
        OUTPUT_QUIET
        ERROR_QUIET
    )

    set(CC4S_PYTHONPATH "${CMAKE_SOURCE_DIR}/test/:${CMAKE_SOURCE_DIR}/thirdparty/pyyaml/lib/")
    if(NOT PY_YAML_FOUND EQUAL 0)
      if(NOT EXISTS "${CMAKE_SOURCE_DIR}/thirdparty/pyyaml")
        message(STATUS
          "PyYAML is not installed for Python3 (${Python3_EXECUTABLE}). "
          "CC4S will build the pyyaml from source")

        execute_process(
          COMMAND bash -c
          "
          mkdir -p ${CMAKE_SOURCE_DIR}/thirdparty/ &&
          git clone https://github.com/yaml/pyyaml ${CMAKE_SOURCE_DIR}/thirdparty/pyyaml
          "
          RESULT_VARIABLE PY_YAML_BUILT)
      else()
        message(STATUS "Locak version of pyyaml will be used")
      endif()
    endif()

    set(TEST_BUILD_DIR ${CMAKE_BINARY_DIR}/test/tests)
    file(MAKE_DIRECTORY ${TEST_BUILD_DIR})

    # copy testis script and global sources
    #file(COPY "${CMAKE_SOURCE_DIR}/test/testis.py" DESTINATION "${CMAKE_BINARY_DIR}/test")
    file(CREATE_LINK "${CMAKE_SOURCE_DIR}/test/testis.py"
      "${CMAKE_BINARY_DIR}/test/testis.py")
    file(CREATE_LINK
      "${CMAKE_SOURCE_DIR}/test/sources.json"
      "${CMAKE_BINARY_DIR}/test/sources.json")
    file(CREATE_LINK
      "${CMAKE_SOURCE_DIR}/test/testis-store"
      "${CMAKE_BINARY_DIR}/test/testis-store"
      SYMBOLIC)

    file(GLOB_RECURSE TEST_JSON_FILES ${CMAKE_SOURCE_DIR}/test/tests/*/test.json)

    foreach(TEST_FILE IN LISTS TEST_JSON_FILES)
        # Get source directory of the test
        get_filename_component(TEST_SRC_DIR ${TEST_FILE} DIRECTORY)
        get_filename_component(TEST_NAME ${TEST_SRC_DIR} NAME)

        # Destination folder in the build tree
        set(TEST_DST_DIR ${TEST_BUILD_DIR}/${TEST_NAME})
        file(MAKE_DIRECTORY ${TEST_DST_DIR})

        # Copy only files in the source directory (non-recursive)
        file(GLOB TEST_FILES "${TEST_SRC_DIR}/*")
        foreach(f IN LISTS TEST_FILES)
            if(EXISTS "${f}" AND NOT IS_DIRECTORY "${f}" AND NOT IS_SYMLINK "${f}")
                get_filename_component(fname ${f} NAME)
                file(CREATE_LINK "${f}" "${TEST_DST_DIR}/${fname}" SYMBOLIC)
            endif()
        endforeach()
        set(TESTIS_SRC "${CMAKE_BINARY_DIR}/test/testis.py")
        set(TESTIS_LINK "${TEST_DST_DIR}/testis.py")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${TESTIS_SRC} ${TESTIS_LINK}
        )

    endforeach()
    # Register the test
    # message("what now ${CMAKE_BINARY_DIR}")
    string(MD5 CC4S_HASH "${CC4S_PATH}")
    set(OUTPUT_NAME "${CMAKE_BUILD_TYPE}-${CC4S_HASH}")

    add_test(
      NAME integration-tests
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test/run-integration-tests 4)

    set_tests_properties(integration-tests PROPERTIES
        ENVIRONMENT "CC4S_RUN=${FULL_RUN_CMD}")

      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/test/run-integration-tests.in
        ${CMAKE_CURRENT_BINARY_DIR}/test/run-integration-tests
        @ONLY)

      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/test/run-integration-tests-check.in
        ${CMAKE_CURRENT_BINARY_DIR}/test/run-integration-tests-check
        @ONLY)

      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/test/run-integration-tests-run.in
        ${CMAKE_CURRENT_BINARY_DIR}/test/run-integration-tests-run
        @ONLY)

      configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/test/run.in
        ${CMAKE_CURRENT_BINARY_DIR}/test/run
        @ONLY)

    message(STATUS "Run tests with `make test` or `ctest -V` for verbose output")

endif()

message(STATUS "CXX FLAGS                 = ${CC4S_FLAGS}")
message(STATUS "CXX DEFS                  = ${CC4S_DEFS}")
message(STATUS "CMAKE_BUILD_TYPE          = ${CMAKE_BUILD_TYPE}")
message(STATUS "Using CTF at GIT_TAG      = ${CTF_GIT_TAG}")
message(STATUS "Using yaml-cpp at GIT_TAG = ${YAML_CPP_GIT_TAG}")
message(STATUS "CMAKE_CXX_COMPILER        = ${CMAKE_CXX_COMPILER}")
message(STATUS "SCALAPACK_LIB             = ${CC4S_SCALAPACK}")
message(STATUS "LAPACK_LIBRARIES          = ${LAPACK_LIBRARIES}")



