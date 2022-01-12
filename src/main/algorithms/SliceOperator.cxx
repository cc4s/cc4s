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

#include <algorithms/SliceOperator.hpp>

#include <tcc/Tcc.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(SliceOperator)

Ptr<MapNode> SliceOperator::run(const Ptr<MapNode> &arguments) {
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
Ptr<MapNode> SliceOperator::run(
  const Ptr<MapNode> &arguments
) {
  typedef TensorExpression<F,TE> T;
  auto op(arguments->getMap("operator"));
  auto tensorExpression(op->getPtr<T>("data"));
  ASSERT_LOCATION(
    tensorExpression,
    "expecting operator to be a tensor expression", op->sourceLocation
  );

  // read dimensions from eigen energies meta data
  auto slicedEigenEnergies(arguments->getMap("slicedEigenEnergies"));
  No = slicedEigenEnergies->getValue<size_t>("holesCount");
  Nv = slicedEigenEnergies->getValue<size_t>("particlesCount");

  auto dimensions(op->getMap("dimensions"));
  size_t d(0);
  for (auto key: dimensions->getKeys()) {
    auto dimension(dimensions->getMap(key));
    if (dimension->getValue<std::string>("type") == "orbital") {
      // only slice dimensions of type 'orbital'
      dims.push_back(d);
    }
    ++d;
  }

  slices = New<MapNode>(op->sourceLocation);
  OUT() <<
    "Slicing " << tensorExpression->inspect()->getName() <<
    " into holes and particles." << std::endl;
  slice(tensorExpression, "");

  // create result
  auto slicedOperator(New<MapNode>(op->sourceLocation));
  // copy all meta tensorExpression from original operator
  for (auto key: op->getKeys()) {
    if (key != "data") {
      slicedOperator->get(key) = op->get(key);
    }
  }
  // enter slices
  slicedOperator->get("slices") = slices;
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("slicedOperator") = slicedOperator;
  return result;
}

template <typename F, typename TE>
void SliceOperator::slice(
  const Ptr<TensorExpression<F,TE>> &tensorExpression, const std::string &parts
) {
  if (parts.length() < dims.size()) {
    slice(tensorExpression, parts + "h");
    slice(tensorExpression, parts + "p");
  } else {
    auto tensor(tensorExpression->inspect());
    std::vector<size_t> begins(tensor->getLens().size());
    std::vector<size_t> ends(tensor->getLens());
    std::string index("");
    for (unsigned int i(0); i < tensor->getLens().size(); ++i) {
      index += ('a' + i);
    }
    for (unsigned int i(0); i < dims.size(); ++i) {
      if (parts[i] == 'h') {
        ends[dims[i]] = No;
      } else {
        begins[dims[i]] = tensor->getLens()[dims[i]] - Nv;
      }
    }
    auto result(Tcc<TE>::template tensor<F>(tensor->getName() + parts ));
    slices->setPtr(
      parts,
      COMPILE_RECIPE(result,
        (*result)[index] <<= (*(*tensorExpression)(begins,ends))[index]
      )
    );
  }
}

