#include <IterativePseudoInverter.hpp>

IterativePseudoInverter::IterativePseudoInverter(
  Matrix<complex> const &matrix_
):
  matrix_(matrix),
  inverse(matrix_.lens[1], matrix_.lens[0], *matrix.wrld),
  square(matrix_.lens[0], matrix_.lens[0], *matrix.wrld)
{
  Matrix<complex> conjugate(matrix_.lens[1], matrix_.lens[0], *matrix_.wrld),
  Univar_Function<> fconj(&MathFunctions::conj);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fconj);
  square["ij"] = matrix_["ik"] * conjugate["kj"];
}

void invert(CTF::Matrix<complex> &pseudoInverse, double accuracy);


