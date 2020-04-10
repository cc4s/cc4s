/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <mixers/Mixer.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>

using namespace cc4s;

template <typename F, typename TE>
Mixer<F,TE>::Mixer(Algorithm *algorithm_): algorithm(algorithm_) {
}

template <typename F, typename TE>
Mixer<F,TE>::~Mixer() {
}

// instantiate
template class Mixer<cc4s::Real<64>,DryEngine>;
template class Mixer<cc4s::Complex<64>,DryEngine>;
template class Mixer<cc4s::Real<64>,CtfEngine>;
template class Mixer<cc4s::Complex<64>,CtfEngine>;


template <typename F, typename TE>
std::map<
  std::string,
  std::function<PTR(ESC(Mixer<F,TE>)) (Algorithm *algorithm)>
> *MixerFactory<F,TE>::mixerMap;

// instantiate
template class MixerFactory<cc4s::Real<64>,DryEngine>;
template class MixerFactory<cc4s::Complex<64>,DryEngine>;
template class MixerFactory<cc4s::Real<64>,CtfEngine>;
template class MixerFactory<cc4s::Complex<64>,CtfEngine>;

