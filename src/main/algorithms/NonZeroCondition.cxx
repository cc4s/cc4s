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

#include <algorithms/NonZeroCondition.hpp>

#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <Parser.hpp>
#include <Scanner.hpp>
#include <Real.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>

#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(NonZeroCondition)

/**
 * \brief Extracts a Delta tensor from the given tensor
 * representing a given non-zero condition. The
 * Delta has indices according to the properties used
 * by the given non-zero condition. Its entries are one
 * for those combinations where the non-zero condition of
 * the given tensor is met and zero otherwise.
 * E.g the non-zero condition "spin" of a
 * spin-conserving one-body Hamiltonian with two state indices is a
 * 2x2 diagonal matrix, indicating that the spin component of the
 * Hamiltonian's indices must match.
 */
Ptr<MapNode> NonZeroCondition::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  Ptr<MapNode> result;
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  } else {
    using TE = DefaultTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  }
  ASSERT_LOCATION(
    result, "expecting operator to be a tensor", arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> NonZeroCondition::run(
  const Ptr<MapNode> &arguments
) {
  typedef TensorExpression<F,TE> T;
  auto tensorExpression(arguments->getPtr<T>("operator"));
  if (!tensorExpression) return nullptr;
  auto tensor(tensorExpression->inspect());

  auto conditionName(arguments->getValue<std::string>("conditionName"));
  for (auto condition: tensor->nonZeroConditions->all) {
    if (condition->name == conditionName) {
      // find dimensions of Delta tensor
      std::vector<Natural<>> conditionLens;
      std::vector<std::string> conditionDimensions;
      for (auto dimProperty: condition->dimensionPropertyReferences) {
        conditionLens.push_back(dimProperty.property->indicesOfProperty.size());
        conditionDimensions.push_back(dimProperty.property->name);
      }
      // create Delta tensor for non-zero condition
      auto Delta(
        Tcc<TE>::template tensor<F>(
          conditionLens, tensor->getName() + "." + conditionName
        )
      );

      // enter one at the coordinates found in the non-zero tuples list
      Natural<> indicesCount(
        Cc4s::world->getRank() == 0 ? condition->tuples.size() : 0
      );
      std::vector<Natural<>> indices(indicesCount);
      std::vector<F> values(indicesCount);
      for (Natural<> t(0); t < indicesCount; ++t) {
        Natural<> index(0);
        for (Natural<> d(conditionLens.size()-1); d > 0; --d) {
          index += condition->tuples[t][d];
          index *= conditionLens[d-1];
        }
        index += condition->tuples[t][0];
        indices[t] = index;
        values[t] = F(1);
      }
      Delta->write(indicesCount, indices.data(), values.data());
      Delta->getUnit() = 1.0;
      // generate result
      auto result(New<MapNode>(SOURCE_LOCATION));
      result->setPtr("nonZeroCondition", Delta);
      return result;
    }
  }

  THROW_LOCATION(
    "Tensor "
      + tensor->getName()
      + " does not have a non-zero condition named "
      + conditionName,
    arguments->sourceLocation
  );
}

