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
-Qoption,cpp,--extended_float_types \
-DDEBUG

# optimization options (only for compiling and linking)
CXXFLAGS += -O0 -g

# libraries provided by the enviornment
# BLAS and LAPACK library contained in mkl
# ScaLAPACK libarary, expects mkl and intelmpi to be loaded
MKL_LIB += -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

CTF_CONFIG_FLAGS = CXX=mpiicc \
                   CXXFLAGS="-O0 -g" \
                   --no-dynamic

LDFLAGS += \
-lpthread -qopenmp \
-mkl \
${MKL_LIB} \
${STATIC_LIBS} \
