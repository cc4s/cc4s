/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DIIS_MIXER_DEFINED
#define DIIS_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class DiisMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(DiisMixer);
    DiisMixer(const Ptr<MapNode> &arguments);
    virtual ~DiisMixer();

    virtual void append(
      const Ptr<FockVector<F,TE>> &A, const Ptr<FockVector<F,TE>> &R
    );
    virtual Ptr<const FockVector<F,TE>> get();
    virtual Ptr<const FockVector<F,TE>> getResiduum();

    Ptr<FockVector<F,TE>> next;
    Ptr<FockVector<F,TE>> nextResiduum; // AI: dont know what this is good for

    /**
     * \brief The amplitudes associated to each residuum in the overlap matrix B
     **/
    std::vector<Ptr<FockVector<F,TE>>> amplitudes;
    /**
     * \brief The residua contained in the overlap matrix B
     **/
    std::vector<Ptr<FockVector<F,TE>>> residua;
    /**
     * \brief Overlap matrix \f$B(i,j) = 2*Re(\langleR(i)|R(j)\rangle),
       B(i,N) = 1, B(N,j) = 1, B(N,N) = 0\f$ where \f$i,j < N\f$.
     * The coefficients \f$c_j\f$ and the Lagrangian multiplyer are then given
     * by \f$c_j = B^+(j,N), lambda = B^+(N,N)\f$, where $\fB^+\f$ denotes
     * the Moore--Pensore pseudo inverse of \f$B\f$.
     **/
    std::vector<F> B;


    int N;
    int nextIndex;
    int count;

    std::vector<double> inverse(std::vector<double> matrix, int N);
    std::vector<complex<double>> inverse(std::vector<complex<double>> matrix, int N);

  };
}

#endif

