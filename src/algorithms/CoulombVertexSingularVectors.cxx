#include <algorithms/CoulombVertexSingularVectors.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CoulombVertexSingularVectors);

CoulombVertexSingularVectors::CoulombVertexSingularVectors(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CoulombVertexSingularVectors::~CoulombVertexSingularVectors() {
}

void CoulombVertexSingularVectors::run() {
  // read the Coulomb vertex GammaGqr
  // Its singular value decomposition is U.Sigma.W*
  // where W* is a matrix with the compound orbital index (q,r)
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("FullCoulombVertex"));

  // construct the conjugate of the Coulomb vertex GammaGqr
  Tensor<complex> conjGammaGqr(*GammaGqr);
  conjugate(conjGammaGqr);

  // contract away the orbitals away, leaving U.Sigma^2.U*
  int NG(GammaGqr->lens[0]);
  Matrix<complex> USSUT(NG, NG, *GammaGqr->wrld, "USSUT");
  LOG(1, "CoulombVertexSingularVectors")
    << "Contracting over orbitals of " << GammaGqr->get_name() 
    << " to get U.Sigma^2.U*, with NG=" << NG << std::endl;
  // FIXME: shouldn't it be Grq? doesn't matter with real orbitals
  USSUT["GH"] = conjGammaGqr["Gqr"] * (*GammaGqr)["Hqr"];

  // use ScaLapack routines to diagonalise the USSUT matrix, i.e. find U
  BlacsWorld world(USSUT.wrld->rank, USSUT.wrld->np);
  ScaLapackMatrix<complex> scaUSSUT(USSUT, &world);
  ScaLapackMatrix<complex> scaU(scaUSSUT);
  ScaLapackHermitianEigenSystemDc<complex> eigenSystem(&scaUSSUT, &scaU);
  double *SS(new double[NG]);
  eigenSystem.solve(SS);

  // get number of field variables
  int NF(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  // write singular vectors back to CTF
  Matrix<complex> U(USSUT);
  scaU.write(U);
  // slice singular vectors U corresponding to NF largest singular values S
  int start[] = {0, NG-NF}, end[] = {NG, NG};
  allocatedTensorArgument<complex>(
    "CoulombVertexSingularVectors", new Tensor<complex>(U.slice(start, end))
  );

  // TODO: also write out the singular values
  delete[] SS;
}

void CoulombVertexSingularVectors::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex, DryTensor<complex>>("FullCoulombVertex")
  );

  DryTensor<complex> conjGammaGqr(*GammaGqr, SOURCE_LOCATION);

  int NG(GammaGqr->lens[0]);
  // get number of field variables
  int NF(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(NG * reduction + 0.5);
  }

  allocatedTensorArgument<complex, DryTensor<complex>>(
    "CoulombVertexSingularVectors",
    new DryMatrix<complex>(NG, NF, NS, SOURCE_LOCATION)
  );
}

