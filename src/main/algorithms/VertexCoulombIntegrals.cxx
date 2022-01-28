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

Ptr<MapNode> VertexCoulombIntegrals::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    return run<DefaultDryTensorEngine>(arguments);
  } else {
    return run<DefaultTensorEngine>(arguments);
  }
}

template <typename TE>
Ptr<MapNode> VertexCoulombIntegrals::run(const Ptr<MapNode> &arguments) {
  auto slicedCoulombVertex(
    arguments->getPtr<TensorSet<Complex<>,TE>>("slicedCoulombVertex")
  );
  Ptr<MapNode> metaData(
    slicedCoulombVertex->get("hh")->inspect()->getMetaData()
  );
  auto halfGrid(metaData->getValue<bool>("halfGrid", 0));
  if (halfGrid) {
    OUT() << "using real coulomb vertex" << std::endl;
    return calculateRealIntegrals<TE>(slicedCoulombVertex);
  } else {
    OUT() << "using complex coulomb vertex" << std::endl;
    return calculateComplexIntegrals<TE>(slicedCoulombVertex);
  }
}

template <typename TE>
Ptr<MapNode> VertexCoulombIntegrals::calculateRealIntegrals(
  const Ptr<TensorSet<Complex<>,TE>> &slicedCoulombVertex
) {
  // get input recipes
  auto GammaGhh(slicedCoulombVertex->get("hh"));
  auto GammaGhp(slicedCoulombVertex->get("hp"));
  auto GammaGph(slicedCoulombVertex->get("ph"));
  auto GammaGpp(slicedCoulombVertex->get("pp"));

  auto NF(GammaGhh->inspect()->lens[0]);
  OUT() << "number of field variables NF: " << NF << std::endl;

#define DEFINE_VERTEX_PART(PART, SLICE) \
  Ptr<TensorRecipe<Real<>,TE>> PART##GammaG##SLICE; \
  { \
    auto result(Tcc<TE>::template tensor<Real<>>( \
      std::string(#PART) + "GammaG" + #SLICE) \
    ); \
    PART##GammaG##SLICE = COMPILE_RECIPE(result, (\
      (*result)["Gqr"] <<= \
        map<Real<>>(PART<Complex<>>, (*GammaG##SLICE)["Gqr"]) \
      ) \
    ); \
  }
  // define intermediate recipes
  DEFINE_VERTEX_PART(real, pp)
  DEFINE_VERTEX_PART(real, ph)
  DEFINE_VERTEX_PART(real, hp)
  DEFINE_VERTEX_PART(real, hh)
  DEFINE_VERTEX_PART(imag, pp)
  DEFINE_VERTEX_PART(imag, ph)
  DEFINE_VERTEX_PART(imag, hp)
  DEFINE_VERTEX_PART(imag, hh)
#undef DEFINE_VERTEX_PART

#define DEFINE_REAL_INTEGRALS_SLICE(LO,RO,LI,RI) \
  { \
    auto sliceName(std::string(#LO) + #RO + #LI + #RI); \
    auto result( \
      Tcc<TE>::template tensor<Real<>>(std::string("V") + sliceName) \
    ); \
    result->getUnit() = pow(realGammaG##RO##RI->inspect()->getUnit(),2.0); \
    coulombIntegrals->get(sliceName) = \
      COMPILE_RECIPE(result, ( \
        (*result)["pqsr"] <<= \
          (*realGammaG##LO##LI)["Gps"] * (*realGammaG##RO##RI)["Gqr"], \
        (*result)["pqsr"] += \
          (*imagGammaG##LO##LI)["Gps"] * (*imagGammaG##RO##RI)["Gqr"] \
        ) \
      ); \
  }
  // define recipes for basic integral slices
  auto coulombIntegrals(New<TensorSet<Real<>,TE>>());
  // TODO: symmetry considerations
  DEFINE_REAL_INTEGRALS_SLICE(h,h,h,h);
  DEFINE_REAL_INTEGRALS_SLICE(p,h,h,h);
  DEFINE_REAL_INTEGRALS_SLICE(h,p,h,h);
  DEFINE_REAL_INTEGRALS_SLICE(p,p,h,h);
  DEFINE_REAL_INTEGRALS_SLICE(h,h,p,h);
  DEFINE_REAL_INTEGRALS_SLICE(p,h,p,h);
  DEFINE_REAL_INTEGRALS_SLICE(h,p,p,h);
  DEFINE_REAL_INTEGRALS_SLICE(p,p,p,h);
  DEFINE_REAL_INTEGRALS_SLICE(h,h,h,p);
  DEFINE_REAL_INTEGRALS_SLICE(p,h,h,p);
  DEFINE_REAL_INTEGRALS_SLICE(h,p,h,p);
  DEFINE_REAL_INTEGRALS_SLICE(p,p,h,p);
  DEFINE_REAL_INTEGRALS_SLICE(h,h,p,p);
  DEFINE_REAL_INTEGRALS_SLICE(p,h,p,p);
  DEFINE_REAL_INTEGRALS_SLICE(h,p,p,p);
  DEFINE_REAL_INTEGRALS_SLICE(p,p,p,p);
#undef DEFINE_REAL_INTEGRALS_SLICE

  // create result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("coulombIntegrals", coulombIntegrals);
  return result;
}

template <typename TE>
Ptr<MapNode> VertexCoulombIntegrals::calculateComplexIntegrals(
  const Ptr<TensorSet<Complex<>,TE>> &slicedCoulombVertex
) {
  // get input recipes
  auto GammaGhh(slicedCoulombVertex->get("hh"));
  auto GammaGhp(slicedCoulombVertex->get("hp"));
  auto GammaGph(slicedCoulombVertex->get("ph"));
  auto GammaGpp(slicedCoulombVertex->get("pp"));

  auto NF(GammaGhh->inspect()->lens[0]);
  OUT() << "number of field variables NF: " << NF << std::endl;

#define DEFINE_VERTEX_CONJT(O,I) \
  Ptr<TensorRecipe<Complex<>,TE>> conjTGammaG##O##I; \
  { \
    auto result(Tcc<TE>::template tensor<Complex<>>( \
      std::string("conjTGammaG") + #O + #I) \
    ); \
    conjTGammaG##O##I = COMPILE_RECIPE(result, \
      ( \
        (*result)["Gqr"] <<= \
          map<Complex<>>(conj<Complex<>>, (*GammaG##I##O)["Grq"]) \
      ) \
    ); \
  }
  // define intermediate recipes
  DEFINE_VERTEX_CONJT(p,p)
  DEFINE_VERTEX_CONJT(p,h)
  DEFINE_VERTEX_CONJT(h,p)
  DEFINE_VERTEX_CONJT(h,h)
#undef DEFINE_VERTEX_CONJT

#define DEFINE_COMPLEX_INTEGRALS_SLICE(LO,RO,LI,RI) \
  { \
    auto sliceName(std::string(#LO) + #RO + #LI + #RI); \
    auto result( \
      Tcc<TE>::template tensor<Complex<>>(std::string("V") + sliceName) \
    ); \
    result->getUnit() = pow(GammaG##RO##RI->inspect()->getUnit(),2.0); \
    coulombIntegrals->get(sliceName) = \
      COMPILE_RECIPE(result, (\
        (*result)["pqsr"] <<= \
          (*conjTGammaG##LO##LI)["Gps"] * (*GammaG##RO##RI)["Gqr"] \
        ) \
      ); \
  }
  // define recipes for integral slices
  auto coulombIntegrals(New<TensorSet<Complex<>,TE>>());
  // TODO: symmetry considerations
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,h,h,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,h,h,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,p,h,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,p,h,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,h,p,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,h,p,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,p,p,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,p,p,h);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,h,h,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,h,h,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,p,h,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,p,h,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,h,p,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,h,p,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(h,p,p,p);
  DEFINE_COMPLEX_INTEGRALS_SLICE(p,p,p,p);
#undef DEFINE_COMPLEX_INTEGRALS_SLICE

  // create result
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("coulombIntegrals", coulombIntegrals);
  return result;
}

