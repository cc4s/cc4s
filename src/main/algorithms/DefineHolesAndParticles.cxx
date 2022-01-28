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

#include <algorithms/DefineHolesAndParticles.hpp>
#include <tcc/Tcc.hpp>
#include <TensorSet.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DefineHolesAndParticles)

Ptr<MapNode> DefineHolesAndParticles::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    return run<DefaultDryTensorEngine>(arguments);
  } else {
    return run<DefaultTensorEngine>(arguments);
  }
}

template <typename TE>
Ptr<MapNode> DefineHolesAndParticles::run(
  const Ptr<MapNode> &arguments
) {
  auto eps(arguments->getPtr<TensorExpression<Real<>,TE>>("eigenEnergies"));
  Ptr<MapNode> metaData(eps->inspect()->getMetaData());
  auto energies(metaData->getMap("energies"));
  auto Np(energies->getSize());

  // find fermi energy to determine No and Nv
  auto fermiEnergy(metaData->getValue<Real<>>("fermiEnergy"));
  size_t No(0);
  while (No < Np && energies->getValue<Real<>>(No) < fermiEnergy) { ++No; }
  ASSERT_LOCATION(
    0 < No, "Fermi energy below all eigen energies.",
    metaData->sourceLocation
  );
  ASSERT_LOCATION(
    No < Np, "Fermi energy above all eigen energies.",
    metaData->sourceLocation
  );

  auto Nv(Np-No);
  OUT() << "number of holes     No: " << No << std::endl;
  OUT() << "number of particles Nv: " << Nv << std::endl;
  OUT() << "number of states    Np: " << Np << std::endl;

  auto epsh(Tcc<TE>::template tensor<Real<>>("epsh"));
  auto epshRecipe(
    COMPILE_RECIPE(epsh,
      (*epsh)["i"] <<= (*(*eps)({0},{No}))["i"]
    )
  );

  auto epsp(Tcc<TE>::template tensor<Real<>>("epsp"));
  auto epspRecipe(
    COMPILE_RECIPE(epsp,
      (*epsp)["a"] <<= (*(*eps)({No},{Np}))["a"]
    )
  );

  // create result
  auto slicedEigenEnergies(
    New<TensorSet<Real<>,TE>>(
      std::map<std::string,Ptr<TensorExpression<Real<>,TE>>>(
        {{"h",epshRecipe}, {"p",epspRecipe}}
      )
    )
  );
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("slicedEigenEnergies", slicedEigenEnergies);
  return result;
}

