/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REDUCE_ENERGY_MATRIX_DEFINED
#define REDUCE_ENERGY_MATRIX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ReduceEnergyMatrix: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ReduceEnergyMatrix);
    ReduceEnergyMatrix(std::vector<Argument> const &argumentList);
    virtual ~ReduceEnergyMatrix();
    /**
     * \brief calculates the eigenvalue decomposition of the given
     * energy matrix \f$E_G^{G'}=U_G^F \lambda_F {U^\ast}_F^{G'}\f$,
     * truncates \f$F\f$ to a minimal set of indices to reproduce the
     * the energy \f${\rm Tr}\{E\}\approx\sum_F\lambda_F\f$ within
     * the required accuracy specified.
     * Subsequently, the particle hole Coulomb vertex is reduced by
     * \f$\Gamma^{qF}_r = \Gamma^{qG}_r U_G^F\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of reducing the Coulomb vertex.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_REDUCTION = 0.5;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;

  protected:
    void preconditionEnergyMatrix();
    void readEnergyMatrix();
    void diagonalizeEnergyMatrix();
    void truncateUnitaryTransform();
    void writeUnitaryTransform();
    void writeEnergySpectrum();

    CTF::Tensor<complex> *EGH;
    int NG;
    int64_t elementsCount;
    int64_t *indices;
    complex *elements;
    double *eigenValues;
    double shift;

    int NF;
    complex *transformElements;
  };
}

#endif

