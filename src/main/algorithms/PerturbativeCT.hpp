#ifndef PERTURBATIVE_CT
#define PERTURBATIVE_CT

#include <map>
#include <TensorSet.hpp>
#include <tcc/Tcc.hpp>

namespace cc4s {

/**
 * \brief Caclulates the interemediates for the complete Triples method
 */
class PerturbativeCT: public Algorithm {
 public:
  ALGORITHM_REGISTRAR_DECLARATION(PerturbativeCT)

  Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
 protected:
  template <typename F, typename TE>
  Ptr<MapNode> run(const Ptr<MapNode> &arguments);
  template <typename F, typename TE>
  Ptr<MapNode> run_tcc(const Ptr<MapNode> &arguments);



};


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
 getPerturbativeCTReference(
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
