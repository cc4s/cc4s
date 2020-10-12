#THIS CONFIG SHOULD WORK ON CQC

include Cc4s.mk
include etc/make/ctf.mk
include etc/make/yaml.mk

# compiler and linker
CXX = mpicxx

# general and language options (for preprocessing, compiling and linking)
CC4S_OPTIONS = \
-fopenmp -std=c++11 \
-Wall -pedantic --all-warnings -fmax-errors=1 \
-Wno-vla \
-Wno-int-in-bool-context \
-DDEBUG

# optimization options (only for compiling and linking)
OPTIMIZE = -O0 -g -fno-lto

CTF_CONFIG_FLAGS = CXX=$(CXX) \
                   AR=gcc-ar \
                   CXXFLAGS="-Ofast -march=native -fno-lto -fopenmp -DPROFILE" \
                   LIBS="-L$(BLAS_PATH)/lib" \
                   --no-dynamic

LINK_LIBS = \
-Wl,-Bstatic \
${YAML_LIB} \
${CTF_LIB} \
${BLAS_LIB} \
${SCALAPACK_LIB} \
-lgfortran -lquadmath \
-Wl,-Bdynamic

INCLUDE_FLAGS = \
${YAML_INCLUDE} \
${BLAS_INCLUDE} \
${CTF_INCLUDE}
