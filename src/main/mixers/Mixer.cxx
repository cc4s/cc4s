
#include <mixers/Mixer.hpp>
#include <algorithms/Algorithm.hpp>

using namespace cc4s;

template <typename F, typename TE>
Mixer<F,TE>::Mixer(const Ptr<MapNode> &) {
}

template <typename F, typename TE>
Mixer<F,TE>::~Mixer() {
}

// instantiate
template class cc4s::Mixer<Real<64>,DefaultDryTensorEngine>;
template class cc4s::Mixer<Complex<64>,DefaultDryTensorEngine>;
template class cc4s::Mixer<Real<64>,DefaultTensorEngine>;
template class cc4s::Mixer<Complex<64>,DefaultTensorEngine>;


template <typename F, typename TE>
Ptr<typename MixerFactory<F,TE>::MixerMap> MixerFactory<F,TE>::mixerMap;

// instantiate
template class cc4s::MixerFactory<Real<64>,DefaultDryTensorEngine>;
template class cc4s::MixerFactory<Complex<64>,DefaultDryTensorEngine>;
template class cc4s::MixerFactory<Real<64>,DefaultTensorEngine>;
template class cc4s::MixerFactory<Complex<64>,DefaultTensorEngine>;

