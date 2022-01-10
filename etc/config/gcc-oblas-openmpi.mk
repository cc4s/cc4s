include Cc4s.mk
include Extern.mk
include etc/make/blas.mk

# compiler and linker
CXX = mpicxx

# general and language options (for preprocessing, compiling and linking)
CXXFLAGS += \
-fopenmp -std=c++11 \
-Wall -pedantic --all-warnings -fmax-errors=3 \
-Wno-vla \
-Wno-int-in-bool-context

# optimization options (only for compiling and linking)
CXXFLAGS += -Ofast -march=native -fno-lto

CTF_CONFIG_FLAGS = CXX=${CXX} \
                   AR=gcc-ar \
                   CXXFLAGS="-Ofast -march=native -fno-lto" \
                   LIBS="-L$(BLAS_PATH)" \
                   --no-dynamic

LDFLAGS += \
-Wl,-Bstatic \
${STATIC_LIBS} \
-lgfortran -lquadmath \
-Wl,-Bdynamic \
${DYNAMIC_LIBS} \
