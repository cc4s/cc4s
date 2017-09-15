/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel, all rights reserved.*/
#ifndef THERMAL_EOM_DEFINED
#define THERMAL_EOM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <vector>

namespace cc4s {
  template <typename F=complex>
  class ThermalHamiltonian {
  public:
    ThermalHamiltonian(
      CTF::Tensor<F> *E0_,
      CTF::Tensor<F> *epsi_,
      CTF::Tensor<F> *epsa_,
      CTF::Tensor<F> *Vabij_
    );

    FockVector<F> rightApply(FockVector<F> &v);
    FockVector<F> leftApply(FockVector<F> &v);

  protected:
    CTF::Tensor<F> *E0;
    CTF::Tensor<F> *epsi, *epsa;
    CTF::Tensor<F> *Vabij;
  };

  /**
   * \brief Implements the diagonal preconditionar for the davidson method
   * \tparam F It is the field variable to be used, in general it will be
   * complex
   */
  template <typename F=complex>
  class ThermalHamiltonianPreConditioner {
  public:
    typedef FockVector<F> V;

    /**
     * \brief Constructor for the preconditioner.
     */
    ThermalHamiltonianPreConditioner (
      CTF::Tensor<F> &E0,
      CTF::Tensor<F> &epsi,
      CTF::Tensor<F> &epsa,
      CTF::Tensor<F> &Vabij
    );

    std::vector<V> getInitialBasis(int eigenVectorsCount);

    V getCorrection(const complex eigenValue, V &residuum);

    V getDiagonalH() const {
      return diagonalH;
    }

  protected:
    V diagonalH;
  };

  class ThermalEquationOfMotion: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalEquationOfMotion);
    ThermalEquationOfMotion(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalEquationOfMotion();

    virtual void run();
  };
}

#endif

