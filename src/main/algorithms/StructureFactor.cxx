#include <algorithms/StructureFactor.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <gte/TricubicInterpolation.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(StructureFactor)

Ptr<MapNode> StructureFactor::run(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  // multiplex calls to template methods
  bool success(false);
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    success =
      StructureFactorCalculator<Real<>,TE>::run(arguments, result) ||
      StructureFactorCalculator<Complex<>,TE>::run(arguments, result);
  } else {
    using TE = DefaultTensorEngine;
    success =
      StructureFactorCalculator<Real<>,TE>::run(arguments, result) ||
      StructureFactorCalculator<Complex<>,TE>::run(arguments, result);
  }
  ASSERT(
    success, "unsupported orbitals type in amplitudes"
  );
  return result;
}

// code for real and hope that complex works
template <typename F, typename TE>
bool StructureFactorCalculator<F,TE>::run(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
){
  auto amplitudesNode(
    arguments->get("amplitudes")->toAtom<Ptr<const TensorUnion<F,TE>>>()
  );
  // amplitudesNode is nullptr if the node was of different type
  if (!amplitudesNode) return false;

  // create calculator object
  StructureFactorCalculator<F,TE> calculator(amplitudesNode->value);

  calculator.calculate(arguments, result);
  return true;
}

template <typename F, typename TE>
void StructureFactorCalculator<F,TE>::calculate(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {
  using TRc = TensorRecipe<Complex<>, TE>;
  using Tc = Tensor<Complex<>, TE>;
  using Tr = Tensor<Real<>, TE>;
  auto amplitudesNode(arguments->get("amplitudes")->toAtom<Ptr<const TensorUnion<F,TE>>>());
  //amplitudesNode is nullptr if the node was of different type
  if (!amplitudesNode) return 0;
  auto amplitudes(amplitudesNode->value);
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );
  //Note: the field variable F is the SVD reduced version of the reciprocal grid G.
  // We work with the full reciprocal grid G in this algorithm

  auto coulombVertexSingularVectors = arguments->getMap("coulombVertexSingularVectors");
  auto singularVectors(coulombVertexSingularVectors->getValue<Ptr<Tc>>("data"));

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaFph(slices->getValue<Ptr<TRc>>("ph"));
  auto GammaFhp(slices->getValue<Ptr<TRc>>("hp"));
  auto GammaGph( Tcc<TE>::template tensor<Complex<>>("Gph"));
  auto GammaGhp( Tcc<TE>::template tensor<Complex<>>("Ghp"));
  COMPILE(
    (*GammaGph)["Gai"] <<= (*GammaFph)["Fai"] * (*singularVectors)["GF"],
    (*GammaGhp)["Gia"] <<= (*GammaFhp)["Fia"] * (*singularVectors)["GF"]
  )->execute();


  //We have to take out calculate the overlap coefficients Cpq(G) from Γpq(G) by taking
  //out the reciprocal Coulomb kernel
  //Finally the StructureFactor reads: S(G)=Cai(G)*Cbj*(G)*Tabij
  auto CGph = ( Tcc<TE>::template tensor<Complex<>>("CGph"));
  auto cTCGhp = ( Tcc<TE>::template tensor<Complex<>>("CGhp"));
  auto CoulombPotential(arguments->getMap("coulombPotential"));
  auto VofG(CoulombPotential->getValue<Ptr<Tr>>("data"));
  auto invSqrtCoulombPotential( Tcc<TE>::template tensor<Complex<>>
    ("invSqrtCoulombPotential"));

  COMPILE(
    (*invSqrtCoulombPotential)["G"] <<=
      map<Complex<>>(inverseSqrt<Complex<>>, (*VofG)["G"]),
    (*CGph)["Gai"]   <<= (*GammaGph)["Gai"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGhp)["Gia"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGph)["Gai"]),
    (*cTCGhp)["Gia"] <<= (*cTCGhp)["Gia"] * (*invSqrtCoulombPotential)["G"],
    (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"]
  )->execute();

  auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));
  COMPILE(
    (*Tabij)["abij"] <<= map<Complex<>>(toComplex<F>, (*Tpphh)["abij"])
  )->execute();

  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*cTCGhp)["Gia"] * (*CGph)["Gbj"] * (*Tabij)["abij"],
    (*SofG)["G"]  += (-1.0) * (*cTCGhp)["Gja"] * (*CGph)["Gbi"] * (*Tabij)["abij"]
  )->execute();


  auto StructureFactor( Tcc<TE>::template tensor<Real<>>("StructureFactor"));
  COMPILE(
    (*StructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

  auto structureFactor(New<MapNode>(SOURCE_LOCATION));
  result->setValue<Ptr<Tensor<Real<>,TE>>>("structureFactor", StructureFactor);
  return 1;
}

template <typename TE>
void StructureFactor::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = Tensor<Real<>, TE>;

  double sum3D(0.0), inter3D(0.0);
  int64_t countNO(0), countNOg(0);
  // hard coded resolution of the fine grid
  int Na(51), Nb(51), Nc(51);
  // this is hard coded for vasp:
  // the coulomb Energy in eV using the lattice dimensions of vasp

  auto gridVectors(arguments->getMap("gridVectors"));
  auto volume(gridVectors->getValue<Real<>>("volume"));

  double factor(4.5835494674469/volume);
  OUT() << factor << std::endl;
  // READ THE MOMENTUM GRID
  auto grid(gridVectors->getValue<Ptr<T>>("data"));
  ASSERT_LOCATION(
    grid, "expecting the reciprocal Grid",
    gridVectors->sourceLocation
  );

  int64_t NG(grid->lens[1]);
  std::vector<Vector<>> cartesianGrid(NG);
  std::vector<Real<>> output(NG*3);
  output = grid->readAll();
  for (int64_t i(0); i < NG; i++){
    cartesianGrid[i][0] = output[3*i+0];
    cartesianGrid[i][1] = output[3*i+1];
    cartesianGrid[i][2] = output[3*i+2];
  }


  // READ THE RECIPROCAL CELL
  std::vector<Vector<>> B(3);
  auto Gix(gridVectors->getValue<Real<>>("Gix"));
  auto Giy(gridVectors->getValue<Real<>>("Giy"));
  auto Giz(gridVectors->getValue<Real<>>("Giz"));
  auto Gjx(gridVectors->getValue<Real<>>("Gjx"));
  auto Gjy(gridVectors->getValue<Real<>>("Gjy"));
  auto Gjz(gridVectors->getValue<Real<>>("Gjz"));
  auto Gkx(gridVectors->getValue<Real<>>("Gkx"));
  auto Gky(gridVectors->getValue<Real<>>("Gky"));
  auto Gkz(gridVectors->getValue<Real<>>("Gkz"));


  B[0][0] = Gix; B[0][1] = Giy; B[0][2] = Giz ;
  B[1][0] = Gjx; B[1][1] = Gjy; B[1][2] = Gjz ;
  B[2][0] = Gkx; B[2][1] = Gky; B[2][2] = Gkz ;
  OUT() << "here you can die\n";

  // READ THE Structure Factor
  std::vector<Real<>> SofG;
  auto structureFactor(arguments->getMap("structureFactor"));
  auto sfactor(structureFactor->getValue<Ptr<T>>("data"));
  ASSERT_LOCATION(
    sfactor, "expecting the reciprocal Grid",
    structureFactor->sourceLocation
  );
  SofG = sfactor->readAll();

  OUT() << "here you can diee\n";
  // the tricubic interpolation requires a rectangular grid
  // ---> transformation

  //construct the transformation matrix, which is the real cell
  std::vector<Vector<>> A(3);
  double Omega((B[0].cross(B[1])).dot(B[2]));
  A[0] = B[1].cross(B[2])/Omega;
  A[1] = B[2].cross(B[0])/Omega;
  A[2] = B[0].cross(B[1])/Omega;


  // determine bounding box in direct coordinates (in reciprocal space)
  Vector<> directMin, directMax;
  for (int g(0); g < NG; ++g)
  for (int d(0); d < 3; ++d) {
    double directComponent(A[d].dot(cartesianGrid[g]));
    directMin[d] = std::min(directMin[d], directComponent);
    directMax[d] = std::max(directMax[d], directComponent);
  }
  // build grid for the entire bounding box
  Vector<> boxDimensions, boxOrigin;
  int64_t boxSize(1);
  for (int d(0); d < 3; ++d) {
    boxSize *= boxDimensions[d] = std::floor(directMax[d] - directMin[d] + 1.5);
    boxOrigin[d] = std::floor(directMin[d] + 0.5);
  }
  OUT() << "here you can dieee\n";

  // allocate and initialize direct-grid StructureFactor
  std::vector<Real<>> directSofG(boxSize, 0.0);

  // enter known SG values
  for (int g(0); g < NG; ++g) {
    int64_t index(0);
    Vector<> directG;
    for (int d(2); d >= 0; --d) {
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
/*

  // dont fully understand the concept of smallBZ
  std::vector<Vector<>> smallBZ;
  smallBZ.push_back(cartesianGrid[1]);
  for (int t(2); t<NG; t++){
    if (IsInSmallBZ(cartesianGrid[t], 1., smallBZ))
      smallBZ.push_back(cartesianGrid[t]);
  }
*/

  // double check: calculate the structure Factor on the regular grid
  for (int64_t i(0); i < NG; ++i) {
    if (cartesianGrid[i].length() < 1e-8) continue;
    sum3D += factor/cartesianGrid[i].sqrLength()*SofG[i];
  }


  MpiCommunicator communicator(
    Cc4s::world->getRank(), Cc4s::world->getProcesses(), Cc4s::world->getComm()
  );
  communicator.barrier();
  for (int t0(-Na); t0 <= Na; ++t0)
  for (int t1(-Nb); t1 <= Nb; ++t1)
  for (int t2(-Nc); t2 <= Nc; ++t2) {
    if ((std::abs(t0+t1+t2)) % Cc4s::world->getProcesses() != Cc4s::world->getRank()) continue;
    Vector<double> ga(((B[0]/double(Na))*double(t0)));
    Vector<double> gb(((B[1]/double(Nb))*double(t1)));
    Vector<double> gc(((B[2]/double(Nc))*double(t2)));
    Vector<double> g(ga+gb+gc);
    // check if g is within smallBZ
//    if (!IsInSmallBZ(g, 2, smallBZ)) continue;
    countNOg++;
    // loop over coarse grid. i.e. shift fine grid onto every coarse grid-point
    for ( size_t i(0); i < NG; i++){
      g  = ga+gb+gc;
      g += cartesianGrid[i];
      Vector<double> directG;
      for (size_t d(0); d < 3; d++) directG[d] = A[d].dot(g);
      countNO++;
      if (g.length() < 1e-8) continue;
      double interpol(interpolatedSofG(directG[0], directG[1], directG[2]));
      inter3D += interpol * factor/g.sqrLength();
    }
  }

  OUT() << "uncorrected: " << sum3D << "\n";
  OUT() << "corrected:   " << inter3D << "\n";

//  double inter3D, sum3D;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("corrected", inter3D);
  energy->setValue<Real<>>("uncorrected", sum3D);
  result->get("energy") = energy;
}


//scale=1 is used to search for the vectors which
//define smallBZ; scale = 2 is used to tell if a vector is
//within the smallBZ or not. For a vector that is within the smallBZ,
//its projection on any vectors which define smallBZ must be less than 1/2
/*bool FiniteSizeCorrection::IsInSmallBZ(
  Vector<> point, double scale, std::vector<Vector<>> smallBZ
)
{
  double epsilon(smallBZ[0].length() * 1e-10);
  size_t countVector(0);
  for (auto &bz: smallBZ) {
    if (
      abs(abs(bz.dot(point))/bz.sqrLength()*scale -1.0) < epsilon
       || abs(bz.dot(point))/bz.sqrLength()*scale -1.0 < -epsilon
    ) {
      countVector++;
    }
    else {
      break;
    }
  }

  if (countVector == smallBZ.size()){
    return true;
  }
  return false;
}
*/

