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

#ifndef _PERTURBATIVE_TRIPLES_HEADER_FILE
#define _PERTURBATIVE_TRIPLES_HEADER_FILE

#include <TensorSet.hpp>
#include <algorithms/Algorithm.hpp>
namespace cc4s {
struct PerturbativeTriples : public Algorithm {
  ALGORITHM_REGISTRAR_DECLARATION(PerturbativeTriples)
  /**
   * \brief run routine as always
   */
  Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
};

template <typename F, typename TE>
CTF::Tensor<F> &get_ctf_tensor(Ptr<TensorSet<F, TE>> tensorSet,
                               const std::string &key) {
  // Get the tensor expression from the TensorSet by key
  auto tensorExpression = tensorSet->get(key);

  if (!tensorExpression) {
    throw std::runtime_error("Tensor with key " + key + " not found.");
  }
  // Evaluate the tensor expression (this will return a Tensor)
  auto evaluatedTensor = tensorExpression->evaluate();
  return evaluatedTensor->getMachineTensor()->tensor;
}

} // namespace cc4s

#endif
