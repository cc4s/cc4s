/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DIIS_MIXER_DEFINED
#define DIIS_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class DiisMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(DiisMixer)
    DiisMixer(const Ptr<MapNode> &arguments);
    virtual ~DiisMixer();

    virtual void append(
      const Ptr<TensorUnion<F,TE>> &A, const Ptr<TensorUnion<F,TE>> &R
    );
    virtual Ptr<const TensorUnion<F,TE>> get();
    virtual Real<64> getResiduumNorm();

    Ptr<TensorUnion<F,TE>> next;
    Ptr<TensorUnion<F,TE>> nextResiduum; // AI: dont know what this is good for

    /**
     * \brief The amplitudes associated to each residuum in the overlap matrix B
     **/
    std::vector<Ptr<TensorUnion<F,TE>>> amplitudes;
    /**
     * \brief The residua contained in the overlap matrix B
     **/
    std::vector<Ptr<TensorUnion<F,TE>>> residua;
    /**
     * \brief Overlap matrix \f$B(i,j) = 2*Re(\langleR(i)|R(j)\rangle),
       B(i,N) = 1, B(N,j) = 1, B(N,N) = 0\f$ where \f$i,j < N\f$.
     * The coefficients \f$c_j\f$ and the Lagrangian multiplyer are then given
     * by \f$c_j = B^+(j,N), lambda = B^+(N,N)\f$, where $\fB^+\f$ denotes
     * the Moore--Pensore pseudo inverse of \f$B\f$.
     **/
    std::vector<F> B;


    size_t N;
    size_t nextIndex;
    size_t count;
    Real<64> residuumNorm;
    std::vector<Real<64>> inverse(std::vector<Real<64>> matrix, size_t N);
    std::vector<Complex<64>> inverse(std::vector<Complex<64>> matrix, size_t N);

  };
}

#endif

