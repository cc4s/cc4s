#Currently Loaded Modulefiles:
# 1) gcc/9.1.0-gcc-4.8.5-mj7s6dg        4) openmpi/4.0.3-gcc-9.1.0-it4lnqb           
# 2) cmake/3.15.1-gcc-9.1.0-4pcxi7s     5) netlib-scalapack/2.1.0-gcc-9.1.0-ja5stbc  
# 3) openblas/0.3.9-gcc-9.1.0-eqtjopu  

include Cc4s.mk
include etc/make/ctf.mk
include etc/make/yaml.mk

# compiler and linker
CXX = mpicxx

# general and language options (for preprocessing, compiling and linking)
CC4S_OPTIONS = \
-fopenmp -std=c++11 \
-Wall -pedantic --all-warnings -fmax-errors=3

# optimization options (only for compiling and linking)
OPTIMIZE = -Ofast -march=native -fno-lto

CTF_CONFIG_FLAGS = CXX=${CXX} \
                   AR=gcc-ar \
                   CXXFLAGS="-Ofast -march=native -fno-lto" \
                   LIBS="-L$(BLAS_PATH)" \
                   --no-dynamic

LINK_LIBS = \
-Wl,-Bstatic \
${YAML_LIB} \
${CTF_LIB} \
${BLAS_LIB} \
-lgfortran -lquadmath \
-Wl,-Bdynamic \
${SCALAPACK_LIB} \



INCLUDE_FLAGS = \
${YAML_INCLUDE} \
${BLAS_INCLUDE} \
${CTF_INCLUDE}
