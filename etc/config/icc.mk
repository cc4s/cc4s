include Libraries.mk
include etc/make/ctf.mk
include etc/make/yaml.mk

# main target
TARGET = Cc4s

# compiler and linker
CXX = mpicxx

# general and language options (for preprocessing, compiling and linking)
OPTIONS = -mkl -lpthread  -std=c++11 \
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

# destination path for installation
INSTALL=~/bin/cc4s/${CONFIG}

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

# vim: ft=make
