#ifndef COMPLEX_TENSOR_DEFINED
#define COMPLEX_TENSOR_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Decomposes the tensor of complex elements into
   * two tensors containing the real and imaginary parts,
   * respectively.
   */
  void fromComplexTensor(
    CTF::Tensor<complex> const &c,
    CTF::Tensor<double> &r, CTF::Tensor<double> &i
  );

  /**
   * \brief Discards the real part of a complex tensor.
   */
  void fromComplexTensor(
    CTF::Tensor<complex> const &c,
    CTF::Tensor<double> &r
  );

  /**
   * \brief Composes a tensor of complex elements
   * containing of the given tensors of real and imaginary parts.
   * Note that in this overload the imaginary part may be redistributed
   * during reading.
   */
  void toComplexTensor(
    CTF::Tensor<double> const &r, CTF::Tensor<double> &i,
    CTF::Tensor<complex> &c
  );

  /**
   * \brief Composes a tensor of complex elements
   * containing of the given tensors of real and imaginary parts.
   * Note this may redistribute the tensor c during writing
   * if the distributions of r and i differ.
   */
  void toComplexTensor(
    CTF::Tensor<double> const &r, CTF::Tensor<double> const &i,
    CTF::Tensor<complex> &c
  );
}

#endif

