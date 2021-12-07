{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell rec {

  buildInputs = with pkgs; [
    openblas
    scalapack
    gfortran
    openmpi
    gcc
  ];

  BLAS_PATH = "${pkgs.openblas}";
  SCALAPACK_PATH = "${pkgs.scalapack}";
  GFORTRAN_LIB_PATH = "${pkgs.gfortran.cc.lib}";
  shellHook = ''
    export LD_LIBRARY_PATH=${pkgs.gfortran.cc.lib}/lib:$LD_LIBRARY_PATH
  '';

}
