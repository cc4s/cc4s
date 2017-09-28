/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <mixers/Mixer.hpp>

using namespace cc4s;

template <typename F>
Mixer<F>::Mixer(Algorithm *algorithm_): algorithm(algorithm_) {
}

template <typename F>
Mixer<F>::~Mixer() {
}

// instantiate
template class Mixer<double>;
template class Mixer<complex>;


template <typename F>
std::map<
  std::string,
  std::function<Mixer<F> *(Algorithm *algorithm)>
> *MixerFactory<F>::mixerMap;

// instantiate
template class MixerFactory<double>;
template class MixerFactory<complex>;

