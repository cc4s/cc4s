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

#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>
#include <algorithms/Algorithm.hpp>

using namespace cc4s;

template <typename F, typename TE>
CoupledClusterMethod<F,TE>::CoupledClusterMethod(
  const Ptr<MapNode> &arguments_
): arguments(arguments_) {
}

template <typename F, typename TE>
CoupledClusterMethod<F,TE>::~CoupledClusterMethod() {
}

// instantiate
template class cc4s::CoupledClusterMethod<Real<64>,DefaultDryTensorEngine>;
template class cc4s::CoupledClusterMethod<Complex<64>,DefaultDryTensorEngine>;
template class cc4s::CoupledClusterMethod<Real<64>,DefaultTensorEngine>;
template class cc4s::CoupledClusterMethod<Complex<64>,DefaultTensorEngine>;


template <typename F, typename TE>
Ptr<
  typename CoupledClusterMethodFactory<F,TE>::CoupledClusterMethodMap
> CoupledClusterMethodFactory<F,TE>::methodMap;

// instantiate
template class cc4s::CoupledClusterMethodFactory<Real<64>,DefaultDryTensorEngine>;
template class cc4s::CoupledClusterMethodFactory<Complex<64>,DefaultDryTensorEngine>;
template class cc4s::CoupledClusterMethodFactory<Real<64>,DefaultTensorEngine>;
template class cc4s::CoupledClusterMethodFactory<Complex<64>,DefaultTensorEngine>;

