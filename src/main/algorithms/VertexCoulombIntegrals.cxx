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

#include <algorithms/VertexCoulombIntegrals.hpp>

#include <tcc/Tcc.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(VertexCoulombIntegrals)

template<typename F, typename TE>
using TensorMap = std::map<std::string, Ptr<TensorExpression<F, TE>>>;


Ptr<MapNode> VertexCoulombIntegrals::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  Ptr<MapNode> result;
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (result = run<Real<>,TE>(arguments)) || (result = run<Complex<>, TE>(arguments));
  } else {
    using TE = DefaultTensorEngine;
    (result = run<Real<>,TE>(arguments)) || (result = run<Complex<>, TE>(arguments));
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type for CoulombVertex",
    arguments->sourceLocation
  );
  return result;

}

template <typename F, typename TE>
Ptr<MapNode> VertexCoulombIntegrals::run(const Ptr<MapNode>& arguments) {
  auto slicedCoulombVertex = arguments->getPtr<TensorSet<F, TE>>("slicedCoulombVertex");
  if (!slicedCoulombVertex) return nullptr;
  TensorMap<F,TE> gammaG;
  TensorMap<F,TE> conjTG;
  auto coulombIntegrals = New<TensorSet<F,TE>>();

  std::vector<std::string> keys = {"pp","ph","hp","hh"};
  for (auto& key : keys) {
    gammaG[key] = slicedCoulombVertex->get(key);
  }

  //TODO: dont conj when we are use real integrals
  for (auto& key : keys) {
    auto t = Tcc<TE>::template tensor<F>(std::string("conjTGammaG") + key);
    std::string rKey{key[1], key[0]};  // reversed order
    EXECUTE(
        (*t)["Gqr"] <<= map<F>(conj<F>, (*gammaG.at(rKey))["Grq"])
    );
    conjTG[key] = t;
  }

  auto NF(conjTG["pp"]->inspect()->lens[0]);
  context["NF"] = NF;
  const char ckeys[] = {'p', 'h'};
  for (char LO : ckeys) {
    for (char RO : ckeys) {
      for (char LI : ckeys) {
        for (char RI : ckeys) {
          std::string leftKey = {LO, LI};
          std::string rightKey = {RO, RI};
          std::string sliceName = {LO, RO, LI, RI};
          LOG() << leftKey << " " << rightKey << " " << sliceName << std::endl;
          auto result = Tcc<TE>::template tensor<F>("V" + sliceName);
          result->getUnit() = pow(gammaG["hh"]->inspect()->getUnit(), 2.0);

          auto recipe = COMPILE_RECIPE(result, (
            (*result)["pqsr"] <<= (*conjTG.at(leftKey))["Gps"] * (*gammaG.at(rightKey))["Gqr"]
          ));
          coulombIntegrals->get(sliceName) = recipe;
        }
      }
    }
  }

  // create result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("coulombIntegrals", coulombIntegrals);
  return result;

}


