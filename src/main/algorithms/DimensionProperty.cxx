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

#include <algorithms/DimensionProperty.hpp>

#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <Parser.hpp>
#include <Scanner.hpp>
#include <Real.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>

#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DimensionProperty)

/**
 * \brief Extracts a property matrix for a property of a
 * selected dimension of the given tensor.
 * The first dimension of the matrix corresponds
 * to the selected dimension of the given tensor,
 * the second dimension corresponds to each possible property value.
 * The matrix contains ones for all combinations where
 * the property of the index in the first dimension matches the
 * property value of the second dimension.
 * The matrix is zero otherwise.
 */
Ptr<MapNode> DimensionProperty::run(const Ptr<MapNode> &arguments) {
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
Ptr<MapNode> DimensionProperty::run(
  const Ptr<MapNode> &arguments
) {
  auto tensorExpression(arguments->getPtr<TensorExpression<F,TE>>("operator"));
  if (!tensorExpression) return nullptr;

  auto tensor(tensorExpression->inspect());
  auto dimensionIndex(arguments->getValue<Natural<>>("dimension"));
  auto propertyName(arguments->getValue<std::string>("property"));

  ASSERT_LOCATION(
    dimensionIndex < tensor->getLens().size(),
    "Tensor "
      + tensor->getName()
      + " has no dimension #"
      + std::to_string(dimensionIndex),
    arguments->sourceLocation
  );

  auto propertyIt(
    tensor->dimensions[dimensionIndex]->properties.find(propertyName)
  );
  ASSERT_LOCATION(
    propertyIt != tensor->dimensions[dimensionIndex]->properties.end(),
    "Tensor "
      + tensor->getName()
      + " dimension "
      + std::to_string(dimensionIndex)
      + " does not have a property "
      + propertyName,
    arguments->sourceLocation
  );
  auto property(propertyIt->second);

  // create property matrix
  std::vector<Natural<>> lens(
    {tensor->getLens()[dimensionIndex], property->indicesOfProperty.size()}
  );
  std::vector<std::string> types(
    {tensor->dimensions[dimensionIndex]->name, propertyName}
  );
  auto P( Tcc<TE>::template tensor<F>(lens, propertyName) );

  // enter one at the respective positions
  std::vector<Natural<>> indices(0);
  std::vector<F> values(0);
  for (Natural<> p(0); p < property->indicesOfProperty.size(); ++p) {
    for (Natural<> index: property->indicesOfProperty[p]) {
      if (Cc4s::world->getRank() == 0) {
        indices.push_back(index + lens[0] * p);
        values.push_back(F(1));
      }
    }
  }
  P->write(indices.size(), indices.data(), values.data());

  // build result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("dimensionProperty", P);
  return result;
}

