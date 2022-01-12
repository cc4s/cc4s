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

#ifndef DIIS_MIXER_DEFINED
#define DIIS_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class DiisMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(DiisMixer)
    DiisMixer(const Ptr<MapNode> &arguments);
    virtual ~DiisMixer();

    std::string describeOptions() override;

    void append(
      const Ptr<TensorSet<F,TE>> &A, const Ptr<TensorSet<F,TE>> &R
    ) override ;
    Ptr<TensorSet<F,TE>> get() override;
    Real<> getResiduumNorm() override;

    Ptr<TensorSet<F,TE>> next;
    Ptr<TensorSet<F,TE>> nextResiduum;

    /**
     * \brief The amplitudes associated to each residuum in the overlap matrix B
     **/
    std::vector<Ptr<TensorSet<F,TE>>> amplitudes;
    /**
     * \brief The residua contained in the overlap matrix B
     **/
    std::vector<Ptr<TensorSet<F,TE>>> residua;
    /**
     * \brief Overlap matrix \f$B(i,j) = 2*Re(\langleR(i)|R(j)\rangle),
       B(i,N) = 1, B(N,j) = 1, B(N,N) = 0\f$ where \f$i,j < N\f$.
     * The coefficients \f$c_j\f$ and the Lagrangian multiplyer are then given
     * by \f$c_j = B^+(j,N), lambda = B^+(N,N)\f$, where $\fB^+\f$ denotes
     * the Moore--Pensore pseudo inverse of \f$B\f$.
     **/
    std::vector<F> B;


    size_t N;
    size_t nextIndex;
    size_t count;
    Real<64> residuumNorm;
    std::vector<Real<64>> inverse(std::vector<Real<64>> matrix, size_t N);
    std::vector<Complex<64>> inverse(std::vector<Complex<64>> matrix, size_t N);

  };
}

#endif

