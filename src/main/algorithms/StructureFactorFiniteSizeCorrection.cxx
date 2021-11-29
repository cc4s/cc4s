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

#include <algorithms/StructureFactorFiniteSizeCorrection.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <gte/TricubicInterpolation.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(StructureFactorFiniteSizeCorrection)

Ptr<MapNode> StructureFactorFiniteSizeCorrection::run(
  const Ptr<MapNode> &arguments
) {
  auto result(New<MapNode>(SOURCE_LOCATION));

  auto slicedCoulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto momentumType(
    slicedCoulombVertex->getMap(
      "indices"
    )->getMap("momentum")->getValue<std::string>("type")
  );


  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (momentumType == "halfGrid") {
//   i   return calculateStructureFactor<Real<>, TE>(arguments);
      calculateStructureFactor<Real<>, TE>(arguments, result);
    } else if (momentumType == "fullGrid") {
//      return calculateStructureFactor<Complex<>,TE>(arguments);
      calculateStructureFactor<Complex<>, TE>(arguments, result);
    }
    interpolation<TE>(arguments, result);
  }
  else {
    using TE = DefaultTensorEngine;
    if (momentumType == "halfGrid") {
//   i   return calculateStructureFactor<Real<>, TE>(arguments);
      calculateStructureFactor<Real<>, TE>(arguments, result);
    } else if (momentumType == "fullGrid") {
//      return calculateStructureFactor<Complex<>,TE>(arguments);
      calculateStructureFactor<Complex<>, TE>(arguments, result);
    }

    interpolation<TE>(arguments, result);
  } 
  return result;
}

// This algorithm works as follows:
// - undo the SVD of the CoulombVertex and bring it back to reciprocal mesh
// - divide te CoulombVertex by the Coulomb potential to obtain the codensities

template <typename F, typename TE>
void StructureFactorFiniteSizeCorrection::calculateStructureFactor(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {
  using Tc = TensorExpression<Complex<>, TE>;
  using Tr = TensorExpression<Real<>, TE>;

  // READ input to compute structure factors

  auto coulombVertexSingularVectors(
    arguments->getMap("coulombVertexSingularVectors")
  );
  auto singularVectors(coulombVertexSingularVectors->getPtr<Tc>("data"));

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input tensor expression
  auto GammaFph(slices->getPtr<Tc>("ph"));
  auto GammaFhp(slices->getPtr<Tc>("hp"));
  auto GammaFhh(slices->getPtr<Tc>("hh"));
  auto GammaGph( Tcc<TE>::template tensor<Complex<>>("Gph"));
  auto GammaGhp( Tcc<TE>::template tensor<Complex<>>("Ghp"));
  auto GammaGhh( Tcc<TE>::template tensor<Complex<>>("Ghh"));
  COMPILE(
    (*GammaGph)["Gai"] <<= (*GammaFph)["Fai"] * (*singularVectors)["GF"],
    (*GammaGhp)["Gia"] <<= (*GammaFhp)["Fia"] * (*singularVectors)["GF"],
    (*GammaGhh)["Gij"] <<= (*GammaFhh)["Fij"] * (*singularVectors)["GF"]
  )->execute();


  //We have to calculate the overlap coefficients Cpq(G) from Γpq(G) by dividing
  //by the reciprocal Coulomb kernel
  //Finally the StructureFactor reads: S(G)=Cai(G) (Cjb(G))^(*) ( Tabij + Tia Tjb )
  auto CGph   = ( Tcc<TE>::template tensor<Complex<>>("CGph"));
  auto cTCGph = ( Tcc<TE>::template tensor<Complex<>>("cTCGph"));
 
  //units of Coulomb potential [Energy*Volume]
  auto CoulombPotential(arguments->getMap("coulombPotential"));
  auto VofG(CoulombPotential->getPtr<Tr>("data"));
  auto invSqrtCoulombPotential(
    Tcc<TE>::template tensor<Complex<>>("invSqrtCoulombPotential")
  );

  COMPILE(
    (*invSqrtCoulombPotential)["G"] <<=
      map<Complex<>>(inverseSqrt<Complex<>>, (*VofG)["G"]),
    // PH codensities
    (*CGph)["Gai"]    <<= (*GammaGph)["Gai"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGph)["Gai"]  <<= map<Complex<>>(conj<Complex<>>, (*GammaGhp)["Gia"]),
    (*cTCGph)["Gai"]  <<= (*cTCGph)["Gai"] * (*invSqrtCoulombPotential)["G"]
  )->execute();

  auto StructureFactor( Tcc<TE>::template tensor<Real<>>("StructureFactor"));
  //prepare T amplitudes


  //THESE TWO LINES ARE SEGFAULTING
  auto amplitudes(arguments->getPtr<TensorUnion<F,TE>>("amplitudes"));

  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );
  auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));
  auto Tai( Tcc<TE>::template tensor<Complex<>>("Tai"));


  COMPILE(
//    (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
    (*Tabij)["abij"] <<= map<Complex<>>(toComplex<F>, (*Tpphh)["abij"]),
    (*Tai)["ai"] <<= map<Complex<>>(toComplex<F>, (*Tph)["ai"]),
    (*Tabij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"] 
  )->execute();

  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*cTCGph)["Gai"] * (*CGph)["Gbj"] * (*Tabij)["abij"],
    (*SofG)["G"]  += (-1.0) * (*cTCGph)["Gaj"] * (*CGph)["Gbi"] * (*Tabij)["abij"],
    (*StructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

//  auto result(New<MapNode>(SOURCE_LOCATION));

  result->setPtr<TensorExpression<Real<>, TE>>(
    "structureFactor", StructureFactor
  );
//  return result;
}

template <typename TE>
void StructureFactorFiniteSizeCorrection::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = TensorExpression<Real<>, TE>;

  Real<> sum3D(0.0), inter3D(0.0);
  Natural<> countNO(0), countNOg(0);
  // resolution for fine grid used for interpolating the transition structure factor
  auto N(arguments->getValue<size_t>("interpolationGridSize", 20));

  auto gridVectors(arguments->getMap("gridVectors"));


  // READ THE MOMENTUM GRID
  auto grid(gridVectors->getPtr<T>("data")->evaluate());
  ASSERT_LOCATION(
    grid, "expecting the reciprocal Grid",
    gridVectors->sourceLocation
  );

  Natural<> NG(grid->lens[1]);
  std::vector<Vector<>> cartesianGrid(NG);
  std::vector<Real<>> output(NG*3);
  output = grid->readAll();
  for (Natural<> i(0); i < NG; i++){
    cartesianGrid[i][0] = output[3*i+0];
    cartesianGrid[i][1] = output[3*i+1];
    cartesianGrid[i][2] = output[3*i+2];
  }


  // reciprocal lattice vectors ( 2pi/a )
  std::vector<Vector<>> B(3);

  auto Gi(gridVectors->getMap("Gi"));
  auto Gj(gridVectors->getMap("Gj"));
  auto Gk(gridVectors->getMap("Gk"));

  auto Gix(Gi->getValue<Real<>>(0));
  auto Giy(Gi->getValue<Real<>>(1));
  auto Giz(Gi->getValue<Real<>>(2));

  auto Gjx(Gj->getValue<Real<>>(0));
  auto Gjy(Gj->getValue<Real<>>(1));
  auto Gjz(Gj->getValue<Real<>>(2));

  auto Gkx(Gk->getValue<Real<>>(0));
  auto Gky(Gk->getValue<Real<>>(1));
  auto Gkz(Gk->getValue<Real<>>(2));

  B[0][0] = Gix; B[0][1] = Giy; B[0][2] = Giz ;
  B[1][0] = Gjx; B[1][1] = Gjy; B[1][2] = Gjz ;
  B[2][0] = Gkx; B[2][1] = Gky; B[2][2] = Gkz ;


  // READ THE Structure Factor
  std::vector<Real<>> SofG;
  auto structureFactor(result->getPtr<T>("structureFactor")->evaluate());
  SofG = structureFactor->readAll();

  // the tricubic interpolation requires a rectangular grid
  // ---> transformation

  //construct the transformation matrix, which is the real cell
  std::vector<Vector<>> A(3);
  Real<> Omega((B[0].cross(B[1])).dot(B[2]));
  A[0] = B[1].cross(B[2])/Omega;
  A[1] = B[2].cross(B[0])/Omega;
  A[2] = B[0].cross(B[1])/Omega;

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto CoulombPotential(arguments->getMap("coulombPotential"));
// convert all computed emergies to units used for CoulombVertex
  auto toCoulombVertexUnits = pow(coulombVertex->getValue<Real<>>("unit"),2.0) /
    pow(gridVectors->getValue<Real<>>("unit"),3.0) / CoulombPotential->getValue<Real<>>("unit");
  // 4.583594607547605 = 4*PI*2*RYTOEV*AUTOA/4/pi**2 =EDEPS/4/pi**2
  //Real<> factor(4.583594607547605*Omega/2.0/M_PI);
  Real<> factor(gridVectors->getValue<Real<>>("unit")/pow(coulombVertex->getValue<Real<>>("unit"),2.0)*Omega/2.0/M_PI/M_PI);
  


  // determine bounding box in direct coordinates (in reciprocal space)
  Vector<> directMin, directMax;
  for (Natural<> g(0); g < NG; ++g)
  for (Natural<> d(0); d < 3; ++d) {
    Real<> directComponent(A[d].dot(cartesianGrid[g]));
    directMin[d] = std::min(directMin[d], directComponent);
    directMax[d] = std::max(directMax[d], directComponent);
  }
  // build grid for the entire bounding box
  Vector<> boxDimensions, boxOrigin;
  Natural<> boxSize(1);
  for (Natural<> d(0); d < 3; ++d) {
    boxSize *= boxDimensions[d] = std::floor(directMax[d] - directMin[d] + 1.5);
    boxOrigin[d] = std::floor(directMin[d] + 0.5);
  }

  // allocate and initialize direct-grid StructureFactor
  std::vector<Real<>> directSofG(boxSize, 0.0);

  for (Natural<> g(0); g < NG; ++g) {
    Natural<> index(0);
    Vector<> directG;
    for (Natural<> d(0); d < 3; ++d) {
      directG[d] = A[d].dot(cartesianGrid[g]);
      index *= boxDimensions[d];
      index += std::floor(directG[d] + 0.5) - boxOrigin[d];
    }
    directSofG[index] = SofG[g];
  }

  // create trilinear or tricubic interpolator
  gte::IntpTricubic3<Real<>> interpolatedSofG(
    (int)boxDimensions[0], (int)boxDimensions[1], (int)boxDimensions[2],
    boxOrigin[0], 1.0, boxOrigin[1], 1.0, boxOrigin[2], 1.0,
    directSofG.data(),
    true
  );

  // Real<> check: calculate the structure Factor on the regular grid
  for (Natural<> i(0); i < NG; ++i) {
    if (cartesianGrid[i].length() < 1e-8) continue;
    sum3D += toCoulombVertexUnits * factor/cartesianGrid[i].sqrLength()*SofG[i];
  }


  MpiCommunicator communicator(
    Cc4s::world->getRank(), Cc4s::world->getProcesses(), Cc4s::world->getComm()
  );
  communicator.barrier();
  for (Natural<> a(0); a < N; ++a)
  for (Natural<> b(0); b < N; ++b)
  for (Natural<> c(0); c < N; ++c) {
    if ((a+b+c) % Cc4s::world->getProcesses() != Cc4s::world->getRank()) continue;
    Vector<Real<>> ga(((B[0]/Real<>(N))*Real<>(a)));
    Vector<Real<>> gb(((B[1]/Real<>(N))*Real<>(b)));
    Vector<Real<>> gc(((B[2]/Real<>(N))*Real<>(c)));
    Vector<Real<>> g(ga+gb+gc);
    countNOg++;
    // loop over coarse grid. i.e. shift fine grid onto every coarse grid-point
    for (Natural<> i(0); i < NG; i++){
      g  = ga+gb+gc;
      g += cartesianGrid[i];
      Vector<Real<>> directG;
      for (Natural<> d(0); d < 3; d++) directG[d] = A[d].dot(g);
      countNO++;
      if (g.length() < 1e-8) continue;
      Real<> interpol(interpolatedSofG(directG[0], directG[1], directG[2]));
      inter3D += toCoulombVertexUnits * interpol * factor/g.sqrLength();
    }
  }
  Real<> totalInter3D(0.0);
  communicator.allReduce(inter3D, totalInter3D);


  OUT() << "Uncorrected correlation energy: " << std::setprecision(10) << sum3D << "\n";
  OUT() << "Finite-size energy correction:   " << std::setprecision(10) << totalInter3D/N/N/N-sum3D << "\n";


  result->setValue("corrected", inter3D);
  result->setValue("uncorrected", sum3D);
}
