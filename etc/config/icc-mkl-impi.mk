include Cc4s.mk
include Extern.mk

# compiler and linker
CXX = mpiicc

# general and language options (for preprocessing, compiling and linking)
CXXFLAGS += \
-mkl -lpthread  -std=c++11 \
-Wall -pedantic -fmax-errors=3 \
-qopenmp \
-qoverride-limits -DINTEL_COMPILER \
-Qoption,cpp,--extended_float_types

# optimization options (only for compiling and linking)
CXXFLAGS += -O3

# libraries provided by the enviornment
# BLAS and LAPACK library contained in mkl
MKL_LIB += -lmkl_blacs_intelmpi_lp64

CTF_CONFIG_FLAGS = CXX=mpiicc \
                   CXXFLAGS="-O3 -std=c++11" \
                   --no-dynamic

LDFLAGS += \
-lpthread -qopenmp \
-mkl \
${MKL_LIB} \
${STATIC_LIBS} \
