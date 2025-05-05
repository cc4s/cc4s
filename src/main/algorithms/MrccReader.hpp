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

#ifndef MRCC_READER_DEFINED
#define MRCC_READER_DEFINED

#include <string>
#include <tcc/Tcc.hpp>
#include <algorithms/Algorithm.hpp>
#include <tcc/TensorExpression.hpp>
#include <DryRun.hpp>

namespace cc4s {
  class MrccReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(MrccReader)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:

    // Tensor wrapper to be able to write a real tensor
    // irrespective of the TensorEngine easier
    // TODO: maybe we want to adapt such a concept in the whole code
    template <typename TE>
    using ___tensor_wrapper = Ptr<Tensor<Real<>, TE>>;

    template <typename TE, enum DryMode dry_mode>
    std::vector<Real<>> get_gamma_ab(bool fnoApprox,
                                     ___tensor_wrapper<TE> Gamma_igm,
                                     ___tensor_wrapper<TE> Epsp,
                                     ___tensor_wrapper<TE> Epsh,
                                     ___tensor_wrapper<TE> Mos_vir);

    template <typename TE, enum DryMode dry_mode>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);

    template <typename T>
    std::vector<T> fileRead(const std::string& , size_t nint = 0);

    template <typename TE>
    std::vector<Real<>> eval_mp2(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                 Ptr<Tensor<Real<>,TE>> Epsp,
                                 Ptr<Tensor<Real<>,TE>> Epsh,
                                 Ptr<Tensor<Real<>,TE>> Mos_vir,
                                 bool eval_gamma_ab = false);

    template <typename TE>
    std::vector<Real<>> eval_mp2_approx_ref(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                            Ptr<Tensor<Real<>,TE>> Epsp,
                                            Ptr<Tensor<Real<>,TE>> Epsh,
                                            Ptr<Tensor<Real<>,TE>> Mos_vir);

    template <typename TE>
    std::vector<Real<>> eval_mp2_approx(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                        Ptr<Tensor<Real<>,TE>> Epsp,
                                        Ptr<Tensor<Real<>,TE>> Epsh,
                                        Ptr<Tensor<Real<>,TE>> Mos_vir);

    template <typename TE>
    std::vector<Real<>> eval_mp2_approx_dry(Ptr<Tensor<Real<>,TE>> Gamma_igm,
                                            Ptr<Tensor<Real<>,TE>> Epsp,
                                            Ptr<Tensor<Real<>,TE>> Epsh,
                                            Ptr<Tensor<Real<>,TE>> Mos_vir);


    #pragma pack(1)
    struct Data1 {
      Real<> value;
      int index;
    };

    struct Data2 {
      Real<> value;
      int index1;
      int index2;
    };

    struct Data3 {
      int64_t index1;
      int64_t index2;
      int64_t index3;
    };
    #pragma pack()

    struct SystemDimensions {
        int conventionalBasis;
        int scf_fitting;
        int correlation_fitting;
        int electrons;
        int core_electrons;
    };


    std::vector<Real<>> halfTransform(
      const std::vector<Real<>> &mos,
      const std::vector<Data3>  &intpos,
      const size_t norb,
      const size_t nbf,
      const size_t ndfbasis
    );


    SystemDimensions extractSystemDimensions(const std::string filename);


  };
}

#endif
