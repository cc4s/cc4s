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

#ifndef COULOMB_INTEGRALS_FROM_VERTEX_DEFINED
#define COULOMB_INTEGRALS_FROM_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates the Coulomb Integrals \f$V_{ij}^{ab}, V_{bj}^{ai},
   * V_{kl}^{ij}, V_{cd}^{ab}, V_{ka}^{ij}, V_{ci}^{ab}\f$ (if given) from the
   * Coulomb Vertex \f$\Gamma_{rG}^q\f$ and stores them in CTF Tensors Vabij,
   * Vaibj, Vijkl, Vabcd, Vijka, and Vabci respectively. The arguments of the
   * integrals are PPPP, PPHH, HHHH, PHPH, HHHP, and PPPHCoulombIntegrals.
   */
  class CoulombIntegralsFromVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombIntegralsFromVertex)
    /**
     * \brief Calculates Coulomb integrals Vabcd,Vabij,Vaibj,Vabci,Vijka,Vijkl
     * from GammaGai,GammaGab,GammaGij Coulomb Vertices.
     * Arguments can be any combination of the sort
     * PPPP, PHPH, PPHH, HHHH, HHHP, PPPHCoulombIntegrals.
    */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename TE>
    Ptr<MapNode> calculateRealIntegrals(const Ptr<MapNode> &arguments);
    template <typename TE>
    Ptr<MapNode> calculateComplexIntegrals(const Ptr<MapNode> &arguments);

/*
    CTF::Tensor<complex> *GammaGab;
    CTF::Tensor<complex> *GammaGij;
    std::array<int,4> syms, vvvv, vovo, vvoo, voov, oovv, oooo, ooov, vooo,
                            vvvo, ovoo, ovov, ovvv, vvov, ovvo, oovo, vovv;
*/
  };
}

#endif

