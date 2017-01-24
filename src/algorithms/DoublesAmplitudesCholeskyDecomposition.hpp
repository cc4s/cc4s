/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DOUBLES_AMPLITUDES_CHOLESKY_DECOMPOSITION_DEFINED
#define DOUBLES_AMPLITUDES_CHOLESKY_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <memory>

namespace cc4s {
  class DoublesAmplitudesCholeskyDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DoublesAmplitudesCholeskyDecomposition);
    DoublesAmplitudesCholeskyDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~DoublesAmplitudesCholeskyDecomposition();
    /**
     * \brief Calculates left singular vectors of the particle-hole Coulomb vertex
     * \f$\tilde\Gamma^a_{iG}\f$.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the left singular vectors of
     * \f$\Gamma^a_{iG}\f$.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_REDUCTION = 2.0;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;

  protected:
    void diagonlizeAmplitudes();
    void sliceLargestEigenValues();

    int Nv, No, NvNo, NF, lower, upper;
    double *lambdas;
    complex *sqrtLambdas;
    int64_t lambdasCount, *lambdaIndices;

    std::shared_ptr<CTF::Tensor<>> Taibj, UaiF;
    std::shared_ptr<CTF::Tensor<complex>> LFai, sqrtLambdaF;
    CTF::Tensor<> *LambdaF;
  };
}

#endif

