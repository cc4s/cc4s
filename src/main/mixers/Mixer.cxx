/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

