/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LINEAR_MIXER_DEFINED 
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F>
  class LinearMixer: public Mixer<F> {
  public:
    MIXER_REGISTRAR_DECLARATION(LinearMixer);
    LinearMixer(Algorithm *algorithm);
    virtual ~LinearMixer();

    virtual void append(
      const PTR(FockVector<F>) &A, const PTR(FockVector<F>) &R
    );
    virtual PTR(FockVector<F>) get();
    virtual PTR(FockVector<F>) getResiduum();

    PTR(FockVector<F>) last;
    PTR(FockVector<F>) lastResiduum;

    double ratio;
  };
}

#endif

