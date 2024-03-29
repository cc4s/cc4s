#ifndef COMPLETE_TRIPLES
#define COMPLETE_TRIPLES

#include <map>
#include <TensorSet.hpp>
#include <tcc/Tcc.hpp>

namespace cc4s {

namespace ct {

// F = Field Value, double, complex etc.
// TE = TensorEngine (CTF)
/*
 * These equations are taken from
 *
 *     Efficient computer implementation of the renormalized coupled-cluster methods:
 *     The R-CCSD[T], R-CCSD(T), CR-CCSD[T], and CR-CCSD(T) approaches
 *
 *     https://doi.org/10.1016/S0010-4655(02)00598-2
 * 
 */
template <typename F, typename TE>
std::shared_ptr<TensorSet<F, TE>>
 getCompleteTriples(
  std::shared_ptr<TensorSet<F, TE>> coulombIntegrals,
  std::shared_ptr<TensorSet<F, TE>> amplitudes
); 

template <typename F, typename TE>
double getDenominator(
 std::shared_ptr<TensorSet<F, TE>> amplitudes,
 Ptr<Tensor<F,TE>> tabcijk,
 Ptr<Tensor<F,TE>> zabcijk
);

}

}

#endif
