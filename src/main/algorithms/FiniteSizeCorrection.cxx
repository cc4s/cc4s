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

  context["method"] = "Finite-size correction";

  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    calculateTransitionStructureFactor<Real<>, TE>(arguments, result)
      || calculateTransitionStructureFactor<Complex<>, TE>(arguments, result);
  }
  else {
    using TE = DefaultTensorEngine;
    calculateTransitionStructureFactor<Real<>, TE>(arguments, result)
      || calculateTransitionStructureFactor<Complex<>, TE>(arguments, result);

    interpolation<Real<>, TE>(arguments, result)
      || interpolation<Complex<>, TE>(arguments, result);
  }
  return result;
}

// This algorithm works as follows:
// - undo the SVD of the CoulombVertex and bring it back to reciprocal mesh
// - divide te CoulombVertex by the Coulomb potential to obtain the codensities


template<typename F, typename TE>
typename std::enable_if<std::is_same<Complex<>,F>::value>::type
rotateVertices(
  Ptr<Tensor<Complex<>,TE>> GammaGph,
  Ptr<TensorExpression<F,TE>> GammaFph,
  Ptr<TensorExpression<Complex<>,TE>> singularVectors
) {
  EXECUTE((*GammaGph)["Gai"] <<= (*GammaFph)["Fai"] * (*singularVectors)["GF"]);
}

// helper for non-Complex case (do nothing)
template<typename F, typename TE>
typename std::enable_if<!std::is_same<Complex<>,F>::value>::type
rotateVertices(
  Ptr<Tensor<Complex<>,TE>> /*GammaGph*/,
  Ptr<TensorExpression<F,TE>> /*GammaFph*/,
  Ptr<TensorExpression<Complex<>,TE>> /*singularVectors*/
) {
  // nothing to do
}



template <typename F, typename TE>
bool FiniteSizeCorrection::calculateTransitionStructureFactor(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {

  // READ input to compute structure factors
  auto amplitudes(arguments->getPtr<TensorSet<F,TE>>("amplitudes"));
  if (!amplitudes) return false;

  auto singularVectors(
    arguments->getPtr<TensorExpression<Complex<>,TE>>("coulombVertexSingularVectors")
  );
  ASSERT_LOCATION(singularVectors, "singular Vectors are of wrong type", arguments->sourceLocation);

  auto coulombVertex(arguments->getPtr<TensorSet<F,TE>>("slicedCoulombVertex"));
  // get input tensor expression: in auxiliary basis
  auto GammaFph(coulombVertex->get("ph"));
  auto GammaFhp(coulombVertex->get("hp"));
  //These are the rotated vertices (in reciprocal mesh basis)

  auto no = GammaFhp->evaluate()->lens[1];
  auto nv = GammaFhp->evaluate()->lens[2];
  auto ng = singularVectors->evaluate()->lens[0];
  auto nf = singularVectors->evaluate()->lens[1];

  context["occ"] = no;
  context["vir"] = nv;
  context["NF"]  = GammaFhp->evaluate()->lens[0];
  context["type"] = TypeTraits<F>::getName();

  //we have to provide the length of these objects for the Real<> version
  auto GammaGph( Tcc<TE>::template tensor<Complex<>>({ng, nv, no}, "Gph"));
  auto GammaGhp( Tcc<TE>::template tensor<Complex<>>({ng, no, nv}, "Ghp"));


  // If F is complex, everything is straighforward and we can just rotate the vertex

  if (std::is_same<Complex<>,F>::value) {
    // as long as we dont have constexpr we have to do it the hard way
    rotateVertices<F,TE>(GammaGph, GammaFph, singularVectors);
    rotateVertices<F,TE>(GammaGhp, GammaFhp, singularVectors);
  } else {

    // Now, the CoulombVertex is of same type as the amplitudes. This means that we
    // have to transform it into complex in order to perform the unitary rotation
    // with the (complex) singular Vectors

    auto np = Cc4s::getProcessesCount();
    auto rank = Cc4s::world->getRank();

    GammaGhp->getMachineTensor();
    auto cGammaFhp( Tcc<TE>::template tensor<Complex<>>({nf, no, nv}, "cGammaFhp"));

    size_t els = nf * no;
    std::vector<size_t> idx(2*els), c_idx(els);
    std::vector<F> data(els*2);
    std::vector<Complex<>> c_data(els);

    size_t rounds = (nv + np - 1) / np; // this is a nice ceil()
    for (size_t r(0); r < rounds; r++) {
      size_t a = r * np + rank;
      if (a >= nv) els = 0;
      std::iota(idx.begin(), idx.end(), a*els*2);
      GammaFhp->evaluate()->read(els*2, idx.data(), data.data());
      // Intentional bitwise copy: interpret pairs of doubles as std::complex<double>
      // Warning -Wclass-memaccess is safe to ignore here.
      std::memcpy(c_data.data(), data.data(), c_data.size() * sizeof(Complex<>));
      std::iota(c_idx.begin(), c_idx.end(), a * els);
      cGammaFhp->write(els, c_idx.data(), c_data.data());
    }
    EXECUTE((*GammaGhp)["Gia"] <<= (*cGammaFhp)["Fia"] * (*singularVectors)["GF"]);

    // i think we need another conj to undo to conj, performed later
    // however, we are only taking the real value of the strucutre factor.....
    //EXECUTE((*GammaGph)["Gai"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGhp)["Gia"]));
    EXECUTE((*GammaGph)["Gai"] <<= (*GammaGhp)["Gia"]);

// The following shows an optimized version which does not need read/write
/*
#if !defined(CTF_SWITCH_REDISTRIBUTION)
#error "CTF_SWITCH_REDISTRIBUTION must be defined to compile this code"
#endif

    // store the partition of the orginal vertex
    auto partFph = GammaFph->evaluate()->getPartition();
    auto partFhp = GammaFhp->evaluate()->getPartition();

    // We just look for a distribution where we distribute the vertex NOT
    // with respect to the auxiliary field variable. --> all aux field entries
    // are consecutive on each rank
    std::vector<int> part(4, 0);
    auto np = Cc4s::getProcessesCount();
    int i, j=1;
    for (i = 1; i < np; i++) {
      if (np % i == 0) {
        j = np / i;
        if ( j < i ) break;
      }
    }
    part[1] = 1; part[2] = i; part[3] = j;
    LOG() << "Custom Partition: " << part[1] <<  " " << part[2] << " " << part[3] << '\n';

    GammaFph->evaluate()->switchDistribution(part);
    GammaFhp->evaluate()->switchDistribution(part);

    // now we can work with the ctf data pointer - because the raw data is corretly aligned
*/
  }

  //We have to calculate the overlap coefficients Cpq(G) from Γpq(G) by dividing
  //by the reciprocal Coulomb kernel
  //Finally the TransitionStructureFactor reads: S(G)=Cai(G) (Cjb(G))^(*) ( Tabij + Tia Tjb )

  auto VofG(arguments->getPtr<Tensor<Real<>, TE>>("coulombPotential"));
  auto iVofG( Tcc<TE>::template tensor<Complex<>>("inverseCoulombPotential"));
  auto inverse( [](const Real<> x) { return 1.0 / Complex<>(x); } );

  COMPILE(
    (*iVofG)["G"] <<= map<Complex<>>(inverse, (*VofG)["G"]),
    // PH codensities
    (*GammaGhp)["Gia"]  <<= (*GammaGhp)["Gia"] * (*iVofG)["G"],
    (*GammaGph)["Gai"]  <<= map<Complex<>>(conj<Complex<>>, (*GammaGph)["Gai"])
  )->execute();

  //prepare T amplitudes

  auto Tph( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );
  auto Xpphh( Tcc<TE>::template tensor<Complex<>>("Xpphh"));

  auto toComplex( [](F x) { return Complex<>(x); } );
  COMPILE(
    (*Xpphh)["abij"] <<= map<Complex<>>(toComplex, (*Tpphh)["abij"]),
    (*Xpphh)["abij"]  += map<Complex<>>(toComplex, (*Tph)["ai"] * (*Tph)["bj"])
  )->execute();

  auto TransitionStructureFactor( Tcc<TE>::template tensor<Real<>>("TransitionStructureFactor"));
  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*GammaGph)["Gai"] * (*GammaGhp)["Gjb"] * (*Xpphh)["abij"],
    (*SofG)["G"]  += (-1.0) * (*GammaGph)["Gaj"] * (*GammaGhp)["Gib"] * (*Xpphh)["abij"],
    (*TransitionStructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

  TransitionStructureFactor->getUnit() = 1.0;
  result->setPtr("transitionStructureFactor", TransitionStructureFactor);
  return true;
}

template <typename F, typename TE>
bool FiniteSizeCorrection::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = TensorExpression<Real<>, TE>;
  auto coulombVertex(
    arguments->getPtr<TensorSet<F,TE>>("slicedCoulombVertex")
  );
  if (!coulombVertex) return false;

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


  Cc4s::world->barrier();
  for (int a(-int(N)); a < int(N); ++a)
  for (int b(-int(N)); b < int(N); ++b)
  for (int c(-int(N)); c < int(N); ++c) {
    if ((a+b+c+3*int(N)) % Cc4s::world->getProcesses() != Cc4s::world->getRank()) continue;
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
  Cc4s::world->allReduce(inter3D, totalInter3D);

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

  return true;
}
