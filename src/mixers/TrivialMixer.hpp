/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TRIVIAL_MIXER_DEFINED 
#define TRIVIAL_MIXER_DEFINED

#include <mixers/Mixer.hpp>

namespace cc4s {
  template <typename F>
  class TrivialMixer: public Mixer<F> {
  public:
    MIXER_REGISTRAR_DECLARATION(TrivialMixer);
    TrivialMixer(Algorithm *algorithm);
    virtual ~TrivialMixer();

    virtual void append(CTF::Tensor<F> &A);
    virtual CTF::Tensor<F> &getNext();

    CTF::Tensor<F> *last;
    double ratio;
  };
}

#endif

