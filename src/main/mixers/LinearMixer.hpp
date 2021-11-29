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

#ifndef LINEAR_MIXER_DEFINED 
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class LinearMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(LinearMixer)
    LinearMixer(const Ptr<MapNode> &arguments);

    std::string describeOptions() override;

    void append(
      const Ptr<TensorUnion<F,TE>> &A, const Ptr<TensorUnion<F,TE>> &R
    ) override ;
    Ptr<TensorUnion<F,TE>> get() override;
    Real<> getResiduumNorm() override;

    Ptr<TensorUnion<F,TE>> last;
    Ptr<TensorUnion<F,TE>> lastResiduum;

    F ratio;
    Real<> residuumNorm;
  };
}

#endif

