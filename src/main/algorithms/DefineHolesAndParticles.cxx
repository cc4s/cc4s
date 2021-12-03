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

#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
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
  auto eigenEnergies(arguments->getMap("eigenEnergies"));
  auto energies(eigenEnergies->getMap("energies"));
  auto eps(eigenEnergies->getPtr<TensorExpression<Real<>,TE>>("data"));
  ASSERT_LOCATION(
    eps, "expecting list of eigenEnergies",
    eigenEnergies->sourceLocation
  );

  auto Np(energies->size());

  // find fermi energy to determine No and Nv
  auto fermiEnergy(eigenEnergies->getValue<Real<>>("fermiEnergy"));
  size_t No(0);
  while (No < Np && energies->getValue<Real<>>(No) < fermiEnergy) { ++No; }
  ASSERT_LOCATION(
    0 < No, "Fermi energy below all eigen energies.",
    eigenEnergies->sourceLocation
  );
  ASSERT_LOCATION(
    No < Np, "Fermi energy above all eigen energies.",
    eigenEnergies->sourceLocation
  );

  auto Nv(Np-No);
  OUT() << "number of holes     No: " << No << std::endl;
  OUT() << "number of particles Nv: " << Nv << std::endl;
  OUT() << "number of states    Np: " << Np << std::endl;

  auto slices(New<MapNode>(eigenEnergies->sourceLocation));
  {
    auto epsi(Tcc<TE>::template tensor<Real<>>("epsi"));
    slices->setPtr(
      "h",
      COMPILE_RECIPE(epsi,
        (*epsi)["i"] <<= (*(*eps)({0},{No}))["i"]
      )
    );
  }
  {
    auto epsa(Tcc<TE>::template tensor<Real<>>("epsa"));
    slices->setPtr(
      "p",
      COMPILE_RECIPE(epsa,
        (*epsa)["a"] <<= (*(*eps)({No},{Np}))["a"]
      )
    );
  }

  // create result
  auto slicedEigenEnergies(New<MapNode>(eigenEnergies->sourceLocation));
  slicedEigenEnergies->get("scalarType") = eigenEnergies->get("scalarType");
  slicedEigenEnergies->get("indices") = eigenEnergies->get("indices");
  slicedEigenEnergies->get("dimensions") = eigenEnergies->get("dimensions");
  slicedEigenEnergies->get("unit") = eigenEnergies->get("unit");
  slicedEigenEnergies->setValue("holesCount", No);
  slicedEigenEnergies->setValue("particlesCount", Nv);
  slicedEigenEnergies->get("slices") = slices;
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("slicedEigenEnergies") = slicedEigenEnergies;
  return result;
}

