# Currently Loaded Modulefiles:
# 1) intel/19.1.2   2) intel-mpi/2019.8   3) intel-mkl/2020.2 4) cmake/3.15.1-intel-19.1.1.217-nkzxutj 

include Cc4s.mk
include etc/make/ctf.mk
include etc/make/yaml.mk

# compiler and linker
CXX = mpiicc

# general and language options (for preprocessing, compiling and linking)
CC4S_OPTIONS = \
-mkl -lpthread  -std=c++11 \
-Wall -pedantic -fmax-errors=3 \
-qopenmp \
-qoverride-limits -DINTEL_COMPILER \
-Qoption,cpp,--extended_float_types

# optimization options (only for compiling and linking)
OPTIMIZE = -O3

# libraries provided by the enviornment
# BLAS and LAPACK library contained in mkl
# ScaLAPACK libarary, expects mkl and intelmpi to be loaded
MKL_LIB += -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

CTF_CONFIG_FLAGS = CXX=mpicxx \
                   CXXFLAGS="-O3" \
                   --no-dynamic

LINK_LIBS = \
-mkl \
${MKL_LIB} \
${YAML_LIB} \
${CTF_LIB} \

INCLUDE_FLAGS = \
${YAML_INCLUDE} \
${CTF_INCLUDE}
