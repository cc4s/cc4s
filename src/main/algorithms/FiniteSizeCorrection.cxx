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

#include <algorithms/FiniteSizeCorrection.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <TensorSet.hpp>
#include <gte/TricubicInterpolation.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(FiniteSizeCorrection)

Ptr<MapNode> FiniteSizeCorrection::run(
  const Ptr<MapNode> &arguments
) {
  auto result(New<MapNode>(SOURCE_LOCATION));

  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    calculateTransitionStructureFactor<Real<>, TE>(arguments, result)
      || calculateTransitionStructureFactor<Complex<>, TE>(arguments, result);
    interpolation<TE>(arguments, result);
  }
  else {
    using TE = DefaultTensorEngine;
    calculateTransitionStructureFactor<Real<>, TE>(arguments, result)
      || calculateTransitionStructureFactor<Complex<>, TE>(arguments, result);
    interpolation<TE>(arguments, result);
  }
  return result;
}

// This algorithm works as follows:
// - undo the SVD of the CoulombVertex and bring it back to reciprocal mesh
// - divide te CoulombVertex by the Coulomb potential to obtain the codensities

template <typename F, typename TE>
bool FiniteSizeCorrection::calculateTransitionStructureFactor(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {
  using Tc = TensorExpression<Complex<>, TE>;
  using Tr = TensorExpression<Real<>, TE>;
  using TSc = TensorSet<Complex<>, TE>;

  // READ input to compute structure factors
  auto amplitudes(arguments->getPtr<TensorSet<F,TE>>("amplitudes"));
  if (!amplitudes) return false;

  auto singularVectors(
    arguments->getPtr<Tc>("coulombVertexSingularVectors")
  );

  auto coulombVertex(arguments->getPtr<TSc>("slicedCoulombVertex"));
  // get input tensor expression
  auto GammaFph(coulombVertex->get("ph"));
  auto GammaFhp(coulombVertex->get("hp"));
  auto GammaFhh(coulombVertex->get("hh"));
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
  //Finally the TransitionStructureFactor reads: S(G)=Cai(G) (Cjb(G))^(*) ( Tabij + Tia Tjb )
  auto CGhp   = ( Tcc<TE>::template tensor<Complex<>>("CGhp"));
  auto cTCGhp = ( Tcc<TE>::template tensor<Complex<>>("cTCGhp"));

  //units of Coulomb potential [Energy*Volume]
  auto VofG(arguments->getPtr<Tr>("coulombPotential"));
  auto invSqrtCoulombPotential(
    Tcc<TE>::template tensor<Complex<>>("invSqrtCoulombPotential")
  );

  auto inverseSqrt( [](const Real<> x) { return 1.0 / Complex<>(sqrt(x)); } );
  COMPILE(
    (*invSqrtCoulombPotential)["G"] <<=
      map<Complex<>>(inverseSqrt, (*VofG)["G"]),
    // PH codensities
    (*CGhp)["Gia"]    <<= (*GammaGhp)["Gia"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGhp)["Gia"]  <<= map<Complex<>>(conj<Complex<>>, (*GammaGph)["Gai"]),
    (*cTCGhp)["Gia"]  <<= (*cTCGhp)["Gia"] * (*invSqrtCoulombPotential)["G"]
  )->execute();

  auto TransitionStructureFactor( Tcc<TE>::template tensor<Real<>>("TransitionStructureFactor"));
  //prepare T amplitudes

  auto Tph( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );
  auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));
  auto Tai( Tcc<TE>::template tensor<Complex<>>("Tai"));


  auto toComplex( [](F x) { return Complex<>(x); } );
  COMPILE(
//    (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
    (*Tabij)["abij"] <<= map<Complex<>>(toComplex, (*Tpphh)["abij"]),
    (*Tai)["ai"] <<= map<Complex<>>(toComplex, (*Tph)["ai"]),
    (*Tabij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"]
  )->execute();

  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*cTCGhp)["Gia"] * (*CGhp)["Gjb"] * (*Tabij)["abij"],
    (*SofG)["G"]  += (-1.0) * (*cTCGhp)["Gja"] * (*CGhp)["Gib"] * (*Tabij)["abij"],
    (*TransitionStructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

  // TODO: determinde unit in tcc
  TransitionStructureFactor->getUnit() = 1.0;
//  auto result(New<MapNode>(SOURCE_LOCATION));

  result->setPtr("transitionStructureFactor", TransitionStructureFactor);
  return true;
}

template <typename TE>
void FiniteSizeCorrection::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = TensorExpression<Real<>, TE>;

  Real<> sum3D(0.0), inter3D(0.0);
  Natural<> countNO(0), countNOg(0);
  // resolution for fine grid used for interpolating the transition structure factor
  auto N(arguments->getValue<size_t>("interpolationGridSize", 20));

  // READ THE MOMENTUM GRID
  auto grid(arguments->getPtr<T>("gridVectors"));



  Natural<> NG(grid->inspect()->lens[1]);
  std::vector<Vector<>> cartesianGrid(NG);
  std::vector<Real<>> output(NG*3);
  output = grid->inspect()->readAll();
  for (Natural<> i(0); i < NG; i++){
    cartesianGrid[i][0] = output[3*i+0];
    cartesianGrid[i][1] = output[3*i+1];
    cartesianGrid[i][2] = output[3*i+2];
  }


  // reciprocal lattice vectors ( 2pi/a )
  std::vector<Vector<>> B(3);

  Ptr<MapNode> metaData(grid->inspect()->getMetaData());
  auto Gi(metaData->getMap("Gi"));
  auto Gj(metaData->getMap("Gj"));
  auto Gk(metaData->getMap("Gk"));

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


  // READ THE Transition Structure Factor
  std::vector<Real<>> SofG;
  auto transitionStructureFactor(result->getPtr<T>("transitionStructureFactor")->evaluate());
  SofG = transitionStructureFactor->readAll();

  // the tricubic interpolation requires a rectangular grid
  // ---> transformation

  //construct the transformation matrix, which is the real cell
  std::vector<Vector<>> A(3);
  Real<> Omega(abs((B[0].cross(B[1])).dot(B[2])));

  A[0] = B[1].cross(B[2])/Omega;
  A[1] = B[2].cross(B[0])/Omega;
  A[2] = B[0].cross(B[1])/Omega;

  auto coulombVertex(
    arguments->getPtr<TensorSet<Complex<>,TE>>("slicedCoulombVertex")
  );
  auto GammaFhh(coulombVertex->get("hh"));
  auto coulombPotential(arguments->getPtr<T>("coulombPotential"));
// convert all computed emergies to units used for CoulombVertex
  auto toCoulombVertexUnits(
    pow(GammaFhh->inspect()->getUnit(),2.0)
      / pow(grid->inspect()->getUnit(),3.0)
      / coulombPotential->inspect()->getUnit()
  );
  auto factor(
      grid->inspect()->getUnit()
      / pow(GammaFhh->inspect()->getUnit(),2.0)
      * Omega/2.0/Pi<>()/Pi<>()
  );



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

  // allocate and initialize direct-grid TransitionStructureFactor
  std::vector<Real<>> directSofG(boxSize, 0.0);

  for (Natural<> g(0); g < NG; ++g) {
    Natural<> index(0);
    Vector<> directG;
    for (Natural<> d(0); d < 3; ++d) {
      directG[2-d] = A[2-d].dot(cartesianGrid[g]);
      index *= boxDimensions[2-d];
      index += std::floor(directG[2-d] + 0.5) - boxOrigin[2-d];
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
  for (int a(-int(N)); a < int(N); ++a)
  for (int b(-int(N)); b < int(N); ++b)
  for (int c(-int(N)); c < int(N); ++c) {
    if ((a+b+c) % Cc4s::world->getProcesses() != Cc4s::world->getRank()) continue;
    Vector<Real<>> ga(((B[0]/Real<>(N*2))*Real<>(a)));
    Vector<Real<>> gb(((B[1]/Real<>(N*2))*Real<>(b)));
    Vector<Real<>> gc(((B[2]/Real<>(N*2))*Real<>(c)));
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

//  OUT() << "Uncorrected correlation energy:   " << std::setprecision(10) << sum3D << "\n";
  OUT() << "Finite-size energy correction:    " << std::setprecision(10) << totalInter3D/N/N/N/8.0-sum3D << "\n";


  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("correction", totalInter3D/N/N/N/8.0-sum3D);
  energy->setValue("uncorrected", sum3D);
  energy->setValue("corrected", totalInter3D/N/N/N/8.0);
  // TODO: should match units of eigenenergies
  energy->setValue("unit",
    coulombPotential->inspect()->getUnit()
      * pow(grid->inspect()->getUnit(),3.0)
  );
  result->get("energy") = energy;

}
