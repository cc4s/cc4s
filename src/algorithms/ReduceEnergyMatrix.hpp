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
     * energy matrix \f$E_G^{G'}=U_G^H \lambda_H {U^\ast}_H^{G'}\f$,
     * truncates \f$H\f$ to a minimal set of indices to reproduce the
     * the energy \f${\rm Tr}\{E\}\approx\sum_H\lambda_H\f$ within
     * the required accuracy specified.
     * Subsequently, the particle hole Coulomb vertex is reduced by
     * \f$\Gamma^{qH}_r = \Gamma^{qG}_r U_G^H\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of reducing the Coulomb vertex.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_REDUCTION = 0.5;

  protected:
    void readEnergyMatrix();
    void diagonalizeEnergyMatrix();
    void truncateUnitaryTransform();
    void writeUnitaryTransform();
    void writeEnergySpectrum();

    CTF::Tensor<complex> *EGH;
    int nG;
    int64_t elementsCount;
    int64_t *indices;
    complex *elements;
    double *eigenValues;

    int ng;
    complex *transformElements;
  };
}

#endif

