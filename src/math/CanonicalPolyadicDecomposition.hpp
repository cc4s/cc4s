#ifndef CANONICAL_POLYADIC_DECOMPOSITION_DEFINED
#define CANONICAL_POLYADIC_DECOMPOSITION_DEFINED

#include <tcc/DryTensor.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Calculates \f$T_{ijk}=A_{iR}B_{jR}C_{kR}\f$ with minimal memory
   * footprint.
   */
  template <typename F=double>
  void composeCanonicalPolyadicDecompositionTensors(
    CTF::Tensor<F> &A, CTF::Tensor<F> &B, CTF::Tensor<F> &C,
    CTF::Tensor<F> &T
  );
  /**
   * \brief Performs a dry run of the calculation
   * \f$T_{ijk}=A_{iR}B_{jR}C_{kR}\f$ with minimal memory footprint.
   */
  template <typename F=double>
  void dryComposeCanonicalPolyadicDecompositionTensors(
    DryTensor<F> &A, DryTensor<F> &B, DryTensor<F> &C,
    DryTensor<F> &T
  );

  /**
   * \brief Calculates \f$A_{iR} = T_{ijk}B^{jR}C^{kR}\f$
   * using a contraction order with minimal memory footprint.
   * \param[in] T the tensor \f$T_{ijk}\f$
   * \param[in] indicesT the indices to compose \f$T\f$, e.g. "ijk"
   * \param[in] B the tensor \f$B_{jR}\f$
   * \param[in] idxB the index letter \f$B\f$ has in common with \f$T\f$
   * \param[in] C the tensor \f$C_{kR}\f$
   * \param[in] idxC the index letter \f$C\f$ has in common with \f$T\f$
   * \param[out] A the tensor \f$A_{iR}\f$
   * \param[in] idxA the index letter \f$A\f$ has in common with \f$T\f$
   */
  template <typename F=double>
  void contractWithCanonicalPolyadicDecompositionTensors(
    CTF::Tensor<F> &T, char const *indicesT,
    CTF::Tensor<F> &B, char const idxB,
    CTF::Tensor<F> &C, char const idxC,
    CTF::Tensor<F> &A, char const idxA
  );
  /**
   * \brief Performs a dry run of the calculation
   * \f$A_{iR} = T_{ijk}B^{jR}C^{kR}\f$
   * using a contraction order with minimal memory footprint.
   * \param[in] T the tensor \f$T_{ijk}\f$
   * \param[in] indicesT the indices to compose \f$T\f$, e.g. "ijk"
   * \param[in] B the tensor \f$B_{jR}\f$
   * \param[in] idxB the index letter \f$B\f$ has in common with \f$T\f$
   * \param[in] C the tensor \f$C_{kR}\f$
   * \param[in] idxC the index letter \f$C\f$ has in common with \f$T\f$
   * \param[out] A the tensor \f$A_{iR}\f$
   * \param[in] idxA the index letter \f$A\f$ has in common with \f$T\f$
   */
  template <typename F=double>
  void dryContractWithCanonicalPolyadicDecompositionTensors(
    DryTensor<F> &T, char const *indicesT,
    DryTensor<F> &B, char const idxB,
    DryTensor<F> &C, char const idxC,
    DryTensor<F> &A, char const idxA
  );
}

#endif

