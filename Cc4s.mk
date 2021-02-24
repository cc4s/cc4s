# Note: Do not edit this file, if you want to override a variable
# do it in your configuration file before you include this file

# Cyclops Tensor Framework ====================================================
CTF_COMMIT     ?= 968f8f9eb6aab1d6b67d2fcc1a70c9fc3b98adfa
CTF_BUILD_PATH ?= $(abspath lib/build/${CONFIG}/ctf/$(CTF_COMMIT))
CTF_LIB        ?= -L${CTF_BUILD_PATH}/lib -lctf
CTF_INCLUDE    ?= -I${CTF_BUILD_PATH}/include

# yaml-cpp ====================================================================
YAML_COMMIT     ?= c9460110e072df84b7dee3eb651f2ec5df75fb18
YAML_BUILD_PATH ?= $(abspath lib/build/${CONFIG}/yaml-cpp/$(YAML_COMMIT))
YAML_LIB        ?= -L${YAML_BUILD_PATH} -lyaml-cpp
YAML_INCLUDE    ?= -I${YAML_BUILD_PATH}/include

# BLAS ========================================================================
BLAS_INCLUDE ?= -I${BLAS_PATH}/include
BLAS_LIB     ?= -L${BLAS_PATH}/lib -lopenblas

# ScaLAPACK ===================================================================
SCALAPACK_LIB ?= -L${SCALAPACK_PATH}/lib -lscalapack

# General settings ============================================================
# destination path for installation
CC4S_INSTALL = ~/bin/cc4s/$(CONFIG)
# main target
CC4S_TARGET = Cc4s
