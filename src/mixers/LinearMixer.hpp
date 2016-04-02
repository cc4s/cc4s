/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LINEAR_MIXER_DEFINED 
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

namespace cc4s {
  template <typename F>
  class LinearMixer: public Mixer<F> {
  public:
    MIXER_REGISTRAR_DECLARATION(LinearMixer);
    LinearMixer(Algorithm *algorithm);
    virtual ~LinearMixer();

    virtual void append(CTF::Tensor<F> &A);
    virtual CTF::Tensor<F> &getNext();

    CTF::Tensor<F> *last;
    double ratio;
  };
}

#endif

