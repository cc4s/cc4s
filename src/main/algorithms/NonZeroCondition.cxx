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

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <util/Scanner.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
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
  auto op(arguments->getMap("operator"));
  auto scalarType(op->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {  
      return run<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {  
      return run<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>(arguments);
    }
  }
  ASSERT_LOCATION(
    false, "scalar type '" + scalarType + "' not supported",
    op->get("scalarType")->sourceLocation
  );
}

template <typename F, typename TE>
Ptr<MapNode> NonZeroCondition::run(
  const Ptr<MapNode> &arguments
) {
  typedef TensorExpression<F,TE> T;
  auto op(arguments->getMap("operator"));
  auto tensorExpression(op->getPtr<T>("data"));
  ASSERT_LOCATION(
    tensorExpression, "expecting operator to be a tensor", op->sourceLocation
  );
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

      // build meta data of Delta tensor
      auto nonZeroConditionNode( New<MapNode>(op->sourceLocation) );
      nonZeroConditionNode->setValue<>("version", std::string("v1.0"));
      nonZeroConditionNode->setValue<>("scalarType", TypeTraits<F>::getName());
      auto dimensionsNode( New<MapNode>(SOURCE_LOCATION) );
      for (Natural<> d(0); d < conditionLens.size(); ++d) {
        auto dimensionNode( New<MapNode>(SOURCE_LOCATION) );
        dimensionNode->setValue<>("length", conditionLens[d]);
        dimensionNode->setValue<>("type", conditionDimensions[d]);
        dimensionsNode->get(d) = dimensionNode;
      }
      nonZeroConditionNode->get("dimensions") = dimensionsNode;
      nonZeroConditionNode->setValue("data", Delta);

      auto result(New<MapNode>(SOURCE_LOCATION));
      result->get("nonZeroCondition") = nonZeroConditionNode;
      return result;
    }
  }

  THROW_LOCATION(
    "Tensor "
      + tensor->getName()
      + " does not have a non-zero condition named "
      + conditionName,
    op->sourceLocation
  );
}

