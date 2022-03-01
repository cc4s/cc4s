{ compiler, pkgs ? import <nixpkgs> {} }:

let
  make-config = {CXX, BLAS_PATH}:
    pkgs.writeText "config.mk"
    ''
    include Cc4s.mk
    include Extern.mk

    # setup openmpi cxx
    OMPI_CXX = ${CXX}
    export OMPI_CXX

    # compiler and linker
    CXX = mpicxx

    # blas path
    BLAS_PATH = ${BLAS_PATH}

    # general and language options (for preprocessing, compiling and linking)
    CXXFLAGS += \
    -fopenmp -std=c++11 \
    -Wall -pedantic --all-warnings -fmax-errors=3 \
    -Wno-vla \
    -Wno-int-in-bool-context

    # optimization options (only for compiling and linking)
    CXXFLAGS += -O0 -g -march=native -fno-lto

    CTF_CONFIG_FLAGS = CXX=$(CXX) \
                       AR=gcc-ar \
                       CXXFLAGS="-O0 -g --std=c++11 -march=native -fno-lto" \
                       LIBS="-L$(BLAS_PATH)" \
                       --no-dynamic \
                       --without-scalapack

    LDFLAGS += \
    -Wl,-Bstatic \
    $(STATIC_LIBS) \
    -lquadmath \
    -Wl,-Bdynamic \
    -L$(BLAS_PATH) \
    -lopenblas \
    $(DYNAMIC_LIBS) \
    '';
in


pkgs.mkShell rec {

  buildInputs = with pkgs; [
    python3
    gdb
    openssh
    cmake
    git
    coreutils
    binutils
    openblas
    scalapack
    gfortran
    openmpi
  ];

  compiler-pkg
    = if compiler    == "gcc11" then pkgs.gcc11
    else if compiler == "gcc10" then pkgs.gcc10
    else if compiler == "gcc9" then pkgs.gcc9
    else if compiler == "gcc8" then pkgs.gcc8
    else if compiler == "gcc7" then pkgs.gcc7
    else if compiler == "gcc6" then pkgs.gcc6
    else if compiler == "clang9" then pkgs.clang_9
    else pkgs.gcc;


  BLAS_PATH = "${pkgs.openblas}";
  SCALAPACK_PATH = "${pkgs.scalapack}";
  GFORTRAN_LIB_PATH = "${pkgs.gfortran.cc.lib}";

  CXX = "${compiler-pkg}/bin/c++";
  CC = "${compiler-pkg}/bin/cc";
  LD = "${compiler-pkg}/bin/ld";

  config = make-config { inherit CXX; inherit BLAS_PATH; };

  shellHook = ''
    export LD_LIBRARY_PATH=${pkgs.gfortran.cc.lib}/lib:$LD_LIBRARY_PATH

    export OMPI_CXX=${CXX}
    export OMPI_CC=${CC}
    CXX=${CXX}
    CC=${CC}
    LD=${LD}

    type -a mpicxx
    mpicxx --version

    export CONFIG=nix-${compiler}
    echo creating etc/config/nix-${compiler}.mk
    ln -fs ${config} etc/config/nix-${compiler}.mk
  '';


}
