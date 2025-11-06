{ compiler, pkgs ? import <nixpkgs> { } }:

pkgs.mkShell rec {

  buildInputs = with pkgs; [
    python3
    gdb
    openssh
    openssl
    cacert
    curl
    cmake
    git
    coreutils
    binutils
    openblas
    scalapack
    gfortran
    openmpi
    llvmPackages.openmp
  ];
  compiler-pkg = pkgs.${compiler};

  CXX = "${compiler-pkg}/bin/c++";
  CC = "${compiler-pkg}/bin/cc";
  LD = "${compiler-pkg}/bin/ld";

  shellHook = ''
    export LD_LIBRARY_PATH=${pkgs.gfortran.cc.lib}/lib:$LD_LIBRARY_PATH
    export OPENSSL_MODULES=${pkgs.openssl.out}/lib/ossl-modules
    export SSL_CERT_FILE=${pkgs.cacert}/etc/ssl/certs/ca-bundle.crt


    openssl version -a

    export OMPI_CXX=${CXX}
    export OMPI_CC=${CC}
    CXX=${CXX}
    CC=${CC}
    LD=${LD}

    type -a mpicxx
    mpicxx --version
  '';

}
