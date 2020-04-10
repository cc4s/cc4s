/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LINEAR_MIXER_DEFINED 
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class LinearMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(LinearMixer);
    LinearMixer(Algorithm *algorithm);
    virtual ~LinearMixer();

    virtual void append(
      const PTR(ESC(FockVector<F,TE>)) &A, const PTR(ESC(FockVector<F,TE>)) &R
    );
    virtual PTR(ESC(const FockVector<F,TE>)) get();
    virtual PTR(ESC(const FockVector<F,TE>)) getResiduum();

    PTR(ESC(FockVector<F,TE>)) last;
    PTR(ESC(FockVector<F,TE>)) lastResiduum;

    F ratio;
  };
}

#endif

