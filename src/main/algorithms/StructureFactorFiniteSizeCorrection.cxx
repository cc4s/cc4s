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
  using TRc = TensorRecipe<Complex<>, TE>;
  using Tc = Tensor<Complex<>, TE>;
  using Tr = Tensor<Real<>, TE>;

  // READ input to compute structure factors

  auto coulombVertexSingularVectors = arguments->getMap("coulombVertexSingularVectors");
  auto singularVectors(coulombVertexSingularVectors->getValue<Ptr<Tc>>("data"));

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaFph(slices->getValue<Ptr<TRc>>("ph"));
  auto GammaFhp(slices->getValue<Ptr<TRc>>("hp"));
  auto GammaFhh(slices->getValue<Ptr<TRc>>("hh"));
  auto GammaGph( Tcc<TE>::template tensor<Complex<>>("Gph"));
  auto GammaGhp( Tcc<TE>::template tensor<Complex<>>("Ghp"));
  auto GammaGhh( Tcc<TE>::template tensor<Complex<>>("Ghh"));
  COMPILE(
    (*GammaGph)["Gai"] <<= (*GammaFph)["Fai"] * (*singularVectors)["GF"],
    (*GammaGhp)["Gia"] <<= (*GammaFhp)["Fia"] * (*singularVectors)["GF"],
    (*GammaGhh)["Gij"] <<= (*GammaFhh)["Fij"] * (*singularVectors)["GF"]
  )->execute();


  //We have to take out calculate the overlap coefficients Cpq(G) from Γpq(G) by taking
  //out the reciprocal Coulomb kernel
  //Finally the StructureFactor reads: S(G)=Cai(G)*Cbj*(G)*Tabij
  auto CGph   = ( Tcc<TE>::template tensor<Complex<>>("CGph"));
  auto cTCGph = ( Tcc<TE>::template tensor<Complex<>>("cTCGph"));

  auto CoulombPotential(arguments->getMap("coulombPotential"));
  auto VofG(CoulombPotential->getValue<Ptr<Tr>>("data"));
  auto invSqrtCoulombPotential( Tcc<TE>::template tensor<Complex<>>
    ("invSqrtCoulombPotential"));

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
  auto amplitudesNode(
    arguments->get("amplitudes")->toAtom<Ptr<const TensorUnion<F,TE>>>()
  );
  auto amplitudes(amplitudesNode->value);

  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );
  auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));


  COMPILE(
//    (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
    (*Tabij)["abij"] <<= map<Complex<>>(toComplex<F>, (*Tpphh)["abij"])
  )->execute();

  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*cTCGph)["Gai"] * (*CGph)["Gbj"] * (*Tabij)["abij"],
    (*SofG)["G"]  += (-1.0) * (*cTCGph)["Gaj"] * (*CGph)["Gbi"] * (*Tabij)["abij"],
    (*StructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

//  auto result(New<MapNode>(SOURCE_LOCATION));

  result->setValue<Ptr<Tensor<Real<>, TE>>>("structureFactor", StructureFactor);
//  return result;
}

template <typename TE>
void StructureFactorFiniteSizeCorrection::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = Tensor<Real<>, TE>;

  Real<> sum3D(0.0), inter3D(0.0);
  Natural<> countNO(0), countNOg(0);
  // hard coded resolution of the fine grid
  Natural<> N(20);
  // this is hard coded for vasp:
  // the coulomb Energy in eV using the lattice dimensions of vasp

  auto gridVectors(arguments->getMap("gridVectors"));
  auto volume(gridVectors->getValue<Real<>>("volume"));

// constant factor is still unclear to me
//  Real<> factor(4.5835494674469/volume);
  Real<> factor(0.00020880076563/volume);
  // READ THE MOMENTUM GRID
  auto grid(gridVectors->getValue<Ptr<T>>("data"));
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


  // READ THE RECIPROCAL CELL
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

/*
  auto Gix(gridVectors->getValue<Real<>>("Gix"));
  auto Giy(gridVectors->getValue<Real<>>("Giy"));
  auto Giz(gridVectors->getValue<Real<>>("Giz"));
  auto Gjx(gridVectors->getValue<Real<>>("Gjx"));
  auto Gjy(gridVectors->getValue<Real<>>("Gjy"));
  auto Gjz(gridVectors->getValue<Real<>>("Gjz"));
  auto Gkx(gridVectors->getValue<Real<>>("Gkx"));
  auto Gky(gridVectors->getValue<Real<>>("Gky"));
  auto Gkz(gridVectors->getValue<Real<>>("Gkz"));
*/

  B[0][0] = Gix; B[0][1] = Giy; B[0][2] = Giz ;
  B[1][0] = Gjx; B[1][1] = Gjy; B[1][2] = Gjz ;
  B[2][0] = Gkx; B[2][1] = Gky; B[2][2] = Gkz ;


  // READ THE Structure Factor
  std::vector<Real<>> SofG;
  auto structureFactor(result->getValue<Ptr<T>>("structureFactor"));
//  ASSERT_LOCATION(
//    structureFactor, "expecting the structureFactor",
//    structureFactor->sourceLocation
//  );
  SofG = structureFactor->readAll();

  // the tricubic interpolation requires a rectangular grid
  // ---> transformation

  //construct the transformation matrix, which is the real cell
  std::vector<Vector<>> A(3);
  Real<> Omega((B[0].cross(B[1])).dot(B[2]));
  A[0] = B[1].cross(B[2])/Omega;
  A[1] = B[2].cross(B[0])/Omega;
  A[2] = B[0].cross(B[1])/Omega;


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

//  OUT() << boxSize << std::endl;
  // enter known SG values
  for (Natural<> g(0); g < NG; ++g) {
    Natural<> index(0);
    Vector<> directG;
    for (Natural<> d(2); d-- > 0; ) {
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
    sum3D += factor/cartesianGrid[i].sqrLength()*SofG[i];
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
      inter3D += interpol * factor/g.sqrLength();
    }
  }
  Real<> totalInter3D(0.0);
  communicator.allReduce(inter3D, totalInter3D);


  OUT() << "Uncorrected correlation energy: " << sum3D << "\n";
  OUT() << "Basis-set energy correction:   " << totalInter3D/N/N/N-sum3D << "\n";


  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("corrected", inter3D);
  energy->setValue<Real<>>("uncorrected", sum3D);
  result->get("energy") = energy;
}
