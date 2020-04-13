/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <mixers/Mixer.hpp>
#include <Data.hpp>

using namespace cc4s;

template <typename F, typename TE>
Mixer<F,TE>::Mixer(const Ptr<Algorithm> &algorithm_): algorithm(algorithm_) {
}

template <typename F, typename TE>
Mixer<F,TE>::~Mixer() {
}

// instantiate
template class Mixer<Real<64>,DryTensorEngine>;
template class Mixer<Complex<64>,DryTensorEngine>;
template class Mixer<Real<64>,DefaultTensorEngine>;
template class Mixer<Complex<64>,DefaultTensorEngine>;


template <typename F, typename TE>
Ptr<typename MixerFactory<F,TE>::MixerMap> MixerFactory<F,TE>::mixerMap;

// instantiate
template class MixerFactory<Real<64>,DryTensorEngine>;
template class MixerFactory<Complex<64>,DryTensorEngine>;
template class MixerFactory<Real<64>,DefaultTensorEngine>;
template class MixerFactory<Complex<64>,DefaultTensorEngine>;

