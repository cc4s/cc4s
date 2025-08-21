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

#ifndef CCSDTNONIT_DEFINED
#define CCSDTNONIT_DEFINED

#include <SharedPointer.hpp>
#include <algorithms/Algorithm.hpp>
#include <tcc/Tcc.hpp>


namespace cc4s {
  /**
   * \brief Implements the non iterative ccsdt methods
   */
  class CcsdtNonit: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdtNonit)

    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

    template <typename F, typename TE>
    std::vector<Real<>> getEnergy(std::map<std::string, Ptr<TensorExpression<F,TE>>>,
                                  Ptr<Tensor<F,TE>> Rppphhh,
                                  std::string name);

    template <typename F, typename TE>
    void getAmplitudes(std::map<std::string, Ptr<TensorExpression<F,TE>>>,
                       Ptr<Tensor<F,TE>> Rppphhh,
                       Ptr<Tensor<F,TE>> Tppphhh);

    template <typename F, typename TE>
    void r3_VT3T1_a(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                    Ptr<Tensor<F,TE>> Tph,
                    Ptr<Tensor<F,TE>> Tpphh,
                    Ptr<Tensor<F,TE>> Tppphhh,
                    Ptr<Tensor<F,TE>> Rppphhh,
                    std::string &t3t3mode);

    template <typename F, typename TE>
    void r3_VT3T2_b(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                    Ptr<Tensor<F,TE>> Tph,
                    Ptr<Tensor<F,TE>> Tpphh,
                    Ptr<Tensor<F,TE>> Tppphhh,
                    Ptr<Tensor<F,TE>> Rppphhh,
                    std::string &t3mode);

    template <typename F, typename TE>
    void r3_VT2T1(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                  Ptr<Tensor<F,TE>> Tph,
                  Ptr<Tensor<F,TE>> Tpphh,
                  Ptr<Tensor<F,TE>> Rppphhh,
                  std::string &t3mode);

    template <typename F, typename TE>
    void r3_VT2T1_2nd(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                  Ptr<Tensor<F,TE>> Tph,
                  Ptr<Tensor<F,TE>> Tpphh,
                  Ptr<Tensor<F,TE>> Rppphhh,
                  std::string &t3mode);


  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);

  };
}

#endif

