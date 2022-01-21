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
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    (result = run<Real<>,TE>(arguments)) ||
      (result = run<Complex<>,TE>(arguments));
  } else {
    using TE = DefaultTensorEngine;
    (result = run<Real<>,TE>(arguments)) ||
      (result = run<Complex<>,TE>(arguments));
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> SliceOperator::run(
  const Ptr<MapNode> &arguments
) {
  auto tensorExpression(
    arguments->getPtr<TensorExpression<F,TE>>("operator")
  );
  if (!tensorExpression) return nullptr;
  auto tensor(tensorExpression->inspect());
  auto slicedEigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  // TODO: determine units in tcc
  ASSERT_LOCATION(
    slicedEigenEnergies, "expecting TensorSet of Real64 as 'slicedEigenEnergies'",
    arguments->sourceLocation
  );

  No = slicedEigenEnergies->get("h")->inspect()->getLen(0);
  Nv = slicedEigenEnergies->get("p")->inspect()->getLen(0);

  for (Natural<> d(0); d < tensor->dimensions.size(); ++d) {
    auto dimension(tensor->dimensions[d]);
    if (dimension->name == "State") {
      // only slice dimensions of type 'State'
      dims.push_back(d);
    }
  }

  OUT() <<
    "Slicing " << tensorExpression->inspect()->getName() <<
    " into holes and particles." << std::endl;
  auto slices( New<TensorSet<F,TE>>() );
  slice(tensorExpression, "", slices);

  // create result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("slicedOperator", slices);
  return result;
}

template <typename F, typename TE>
void SliceOperator::slice(
  const Ptr<TensorExpression<F,TE>> &tensorExpression,
  const std::string &parts,
  const Ptr<TensorSet<F,TE>> &slices
) {
  if (parts.length() < dims.size()) {
    slice(tensorExpression, parts + "h", slices);
    slice(tensorExpression, parts + "p", slices);
  } else {
    auto tensor(tensorExpression->inspect());
    std::vector<Natural<>> begins(tensor->getLens().size());
    std::vector<Natural<>> ends(tensor->getLens());
    std::string index("");
    for (Natural<> i(0); i < tensor->getLens().size(); ++i) {
      index += ('a' + i);
    }
    for (Natural<> i(0); i < dims.size(); ++i) {
      if (parts[i] == 'h') {
        ends[dims[i]] = No;
      } else {
        begins[dims[i]] = tensor->getLens()[dims[i]] - Nv;
      }
    }
    auto result(Tcc<TE>::template tensor<F>(tensor->getName() + parts ));
    slices->get(parts) =
      COMPILE_RECIPE(result,
        (*result)[index] <<= (*(*tensorExpression)(begins,ends))[index]
      );
  }
}

