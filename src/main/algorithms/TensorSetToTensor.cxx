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

#include <algorithms/TensorSetToTensor.hpp>

#include <tcc/Tcc.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorSetToTensor)

Ptr<MapNode> TensorSetToTensor::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
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
    result, "algorithm not successful",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> TensorSetToTensor::run(
  const Ptr<MapNode> &arguments
) {


  auto tensorSet(arguments->template getPtr<TensorSet<F,TE>>("tensorSet"));

  if (!tensorSet) return nullptr;

  auto component(arguments->template getValue<std::string>("component"));
  auto outputTensor(tensorSet->get(component));
  ASSERT_LOCATION(
   outputTensor, "given tensor Set does not have a component named " + component,
   SOURCE_LOCATION
  )

  OUT() << "Taking out the component " << component << " from tensorSet\n";

  // create result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("tensor", outputTensor);
  return result;
}

