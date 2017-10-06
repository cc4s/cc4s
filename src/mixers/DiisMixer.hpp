/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DIIS_MIXER_DEFINED 
#define DIIS_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

#include <ctf.hpp>

namespace cc4s {
  template <typename F>
  class DiisMixer: public Mixer<F> {
  public:
    MIXER_REGISTRAR_DECLARATION(DiisMixer);
    DiisMixer(Algorithm *algorithm);
    virtual ~DiisMixer();

    virtual void append(
      const PTR(FockVector<F>) &A, const PTR(FockVector<F>) &R
    );
    virtual PTR(FockVector<F>) get();
    virtual PTR(FockVector<F>) getResiduum();

    /**
     * \brief The best estimate for the next amplitudes as returned by get().
     **/
    PTR(FockVector<F>) next;

    /**
     * \brief The best estimate for the next residuum as returned by
     * getResiduum().
     **/
    PTR(FockVector<F>) nextResiduum;

    /**
     * \brief The amplitudes associated to each residuum in the overlap matrix B
     **/
    std::vector<PTR(FockVector<F>)> amplitudes;
    /**
     * \brief The residua contained in the overlap matrix B
     **/
    std::vector<PTR(FockVector<F>)> residua;

    /**
     * \brief Overlap matrix \f$B(i,j) = 2*Re(\langleR(i)|R(j)\rangle),
       B(i,N) = 1, B(N,j) = 1, B(N,N) = 0\f$ where \f$i,j < N\f$.
     * The coefficients \f$c_j\f$ and the Lagrangian multiplyer are then given
     * by \f$c_j = B^+(j,N), lambda = B^+(N,N)\f$, where $\fB^+\f$ denotes
     * the Moore--Pensore pseudo inverse of \f$B\f$.
     **/
    PTR(CTF::Tensor<F>) B;

    /**
     * \brief The index of the residuum in the residua vector and in the
     * overlap matrix to be replaced with the next one.
     **/
    int nextIndex;

    /**
     * \brief The number of non-zero amplitudes in the vector.
     **/
    int count;
  };
}

#endif

