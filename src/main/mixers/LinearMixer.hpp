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
    LinearMixer(const Ptr<MapNode> &arguments);

    virtual void append(
      const Ptr<FockVector<F,TE>> &A, const Ptr<FockVector<F,TE>> &R
    );
    virtual Ptr<const FockVector<F,TE>> get();
    virtual Real<> getResiduumNorm();

    Ptr<FockVector<F,TE>> last;
    Ptr<FockVector<F,TE>> lastResiduum;

    F ratio;
    Real<> residuumNorm;
  };
}

#endif

