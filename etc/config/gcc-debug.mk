include Cc4s.mk
include Extern.mk
include etc/make/scalapack.mk
include etc/make/blas.mk

# compiler and linker
CXX = mpicxx

# general and language options (for preprocessing, compiling and linking)
CXXFLAGS += \
-fopenmp -std=c++11 \
-Wall -pedantic --all-warnings -fmax-errors=1 \
-Wno-vla \
-Wno-int-in-bool-context \
-DDEBUG

# optimization options (only for compiling and linking)
CXXFLAGS += -O0 -g -fno-lto

CTF_CONFIG_FLAGS = CXX=$(CXX) \
                   AR=gcc-ar \
                   CXXFLAGS="-Ofast -march=native -fno-lto -fopenmp -DPROFILE" \
                   LIBS="-L$(BLAS_PATH)/lib" \
                   --no-dynamic

LDFLAGS += \
-Wl,-Bstatic \
${STATIC_LIBS} \
-lgfortran -lquadmath \
-Wl,-Bdynamic \
${DYNAMIC_LIBS} \
