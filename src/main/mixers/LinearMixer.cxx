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

#include <mixers/LinearMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <algorithms/Algorithm.hpp>
#include <Node.hpp>

using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(LinearMixer)

template <typename F, typename TE>
LinearMixer<F,TE>::LinearMixer(
  const Ptr<MapNode> &arguments
):
  Mixer<F,TE>(arguments), last(nullptr), lastResiduum(nullptr)
{
  ratio = arguments->getValue<Real<>>("ratio", 1.0);
  LOG() << "ratio=" << real(ratio) << std::endl;
}

template <typename F, typename TE>
void LinearMixer<F,TE>::append(
  const Ptr<TensorUnion<F,TE>> &next,
  const Ptr<TensorUnion<F,TE>> &nextResiduum
) {
  if (last) {
    // mix accordingly
    *last *= F(1)-ratio;
    *next *= ratio;
    *next += *last;

    *lastResiduum *= F(1)-ratio;
    *nextResiduum *= ratio;
    *nextResiduum += *lastResiduum;
  }
  last = next;
  lastResiduum = nextResiduum;
  residuumNorm = sqrt(real(lastResiduum->dot(*lastResiduum)));
}

template <typename F, typename TE>
Ptr<TensorUnion<F,TE>> LinearMixer<F,TE>::get() {
  return last;
}

template <typename F, typename TE>
double LinearMixer<F,TE>::getResiduumNorm() {
  return residuumNorm;
}

// instantiate
template class cc4s::LinearMixer<Real<64>, DefaultDryTensorEngine>;
template class cc4s::LinearMixer<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::LinearMixer<Real<64>, DefaultTensorEngine>;
template class cc4s::LinearMixer<Complex<64>, DefaultTensorEngine>;

