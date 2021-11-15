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

#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(VertexCoulombIntegrals)

Ptr<MapNode> VertexCoulombIntegrals::run(const Ptr<MapNode> &arguments) {
  auto slicedCoulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto momentumType(
    slicedCoulombVertex->getMap(
      "indices"
    )->getMap("momentum")->getValue<std::string>("type")
  );

  // multiplex calls to template methods
  if (momentumType == "halfGrid") {
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      return calculateRealIntegrals<TE>(arguments);
    } else {
      using TE = DefaultTensorEngine;
      return calculateRealIntegrals<TE>(arguments);
    }
  } else if (momentumType == "fullGrid") {
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      return calculateComplexIntegrals<TE>(arguments);
    } else {
      using TE = DefaultTensorEngine;
      return calculateComplexIntegrals<TE>(arguments);
    }
  }
  ASSERT_LOCATION(
    false, "momentum 'type' must specify either 'halfGrid' or 'fullGrid'",
    slicedCoulombVertex->getMap(
      "indices"
    )->getMap("momentum")->get("type")->sourceLocation
  );
}

template <typename TE>
Ptr<MapNode> VertexCoulombIntegrals::calculateRealIntegrals(
  const Ptr<MapNode> &arguments
) {
  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaGpp(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("pp"));
  auto GammaGph(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("ph"));
  auto GammaGhp(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("hp"));
  auto GammaGhh(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("hh"));

  auto NF(GammaGhh->inspect()->lens[0]);
  OUT() << "number of field variables NF= " << NF << std::endl;

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
    integralSlices->setValue(sliceName, \
      COMPILE_RECIPE(result, ( \
        (*result)["pqsr"] <<= \
          (*realGammaG##LO##LI)["Gps"] * (*realGammaG##RO##RI)["Gqr"], \
        (*result)["pqsr"] += \
          (*imagGammaG##LO##LI)["Gps"] * (*imagGammaG##RO##RI)["Gqr"] \
        ) \
      ) \
    ); \
  }
  // define recipes for basic integral slices
  auto integralSlices(New<MapNode>(SOURCE_LOCATION));
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
  auto coulombIntegrals(New<MapNode>(SOURCE_LOCATION));
  coulombIntegrals->get("slices") = integralSlices;
  coulombIntegrals->setValue<std::string>("scalarType", "real64");
  coulombIntegrals->setValue<Real<>>(
    "unit", pow(coulombVertex->getValue<Real<>>("unit"),2.0)
  );
  // create indices entry
  auto indices(New<MapNode>(SOURCE_LOCATION));
  indices->get("orbital") = coulombVertex->getMap("indices")->get("orbital");
  coulombIntegrals->get("indices") = indices;
  // create dimensions entry
  auto dimensions(New<MapNode>(SOURCE_LOCATION));
  dimensions->get(0) = coulombVertex->getMap("dimensions")->get(1);
  dimensions->get(1) = coulombVertex->getMap("dimensions")->get(1);
  dimensions->get(2) = coulombVertex->getMap("dimensions")->get(2);
  dimensions->get(3) = coulombVertex->getMap("dimensions")->get(2);
  coulombIntegrals->get("dimensions") = dimensions;
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("coulombIntegrals") = coulombIntegrals;
  return result;
}

template <typename TE>
Ptr<MapNode> VertexCoulombIntegrals::calculateComplexIntegrals(
  const Ptr<MapNode> &arguments
) {
  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaGpp(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("pp"));
  auto GammaGph(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("ph"));
  auto GammaGhp(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("hp"));
  auto GammaGhh(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("hh"));

  auto NF(GammaGhh->inspect()->lens[0]);
  OUT() << "number of field variables NF= " << NF << std::endl;

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

#define DEFINE_COMPLEX_INTEGRALS_SLICE(LO,RO,LI,RI) \
  { \
    auto sliceName(std::string(#LO) + #RO + #LI + #RI); \
    auto result( \
      Tcc<TE>::template tensor<Complex<>>(std::string("V") + sliceName) \
    ); \
    integralSlices->setValue(sliceName, \
      COMPILE_RECIPE(result, (\
        (*result)["pqsr"] <<= \
          (*conjTGammaG##LO##LI)["Gps"] * (*GammaG##RO##RI)["Gqr"] \
        ) \
      ) \
    ); \
  }
  // define recipes for integral slices
  auto integralSlices(New<MapNode>(SOURCE_LOCATION));
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

  // create result
  auto coulombIntegrals(New<MapNode>(SOURCE_LOCATION));
  coulombIntegrals->get("slices") = integralSlices;
  coulombIntegrals->setValue<Real<>>(
    "unit", pow(coulombVertex->getValue<Real<>>("unit"),2.0)
  );
  coulombIntegrals->setValue<std::string>("scalarType", "complex64");
  // create indices entry
  auto indices(New<MapNode>(SOURCE_LOCATION));
  indices->get("orbital") = coulombVertex->getMap("indices")->get("orbital");
  coulombIntegrals->get("indices") = indices;
  // create dimensions entry
  auto dimensions(New<MapNode>(SOURCE_LOCATION));
  dimensions->get(0) = coulombVertex->getMap("dimensions")->get(1);
  dimensions->get(1) = coulombVertex->getMap("dimensions")->get(1);
  dimensions->get(2) = coulombVertex->getMap("dimensions")->get(2);
  dimensions->get(3) = coulombVertex->getMap("dimensions")->get(2);
  coulombIntegrals->get("dimensions") = dimensions;
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("coulombIntegrals") = coulombIntegrals;
  return result;
}

