#include <algorithms/ProjectorEnergyMatrixFromCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ProjectorEnergyMatrixFromCoulombVertex);

ProjectorEnergyMatrixFromCoulombVertex::ProjectorEnergyMatrixFromCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ProjectorEnergyMatrixFromCoulombVertex::~ProjectorEnergyMatrixFromCoulombVertex() {
}

void ProjectorEnergyMatrixFromCoulombVertex::run() {
  // read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("CoulombVertex"));

  // construct the conjugate of the Coulomb vertex GammaGqr
  Tensor<complex> conjGammaGqr(*GammaGqr);
  conjugate(conjGammaGqr);

  // compute the energy matrix
  Matrix<complex> *energyMatrix = new Matrix<complex>(
    GammaGqr->lens[0], GammaGqr->lens[0], *GammaGqr->wrld
  );
  allocatedTensorArgument<complex>("EnergyMatrix", energyMatrix);

  LOG(1, "EnergyMatrix") << "Computing projector energy matrix from Coulomb vertex" << GammaGqr->get_name() 
			 << ", with NG=" << GammaGqr->lens[0] << std::endl;

  (*energyMatrix)["GH"] = conjGammaGqr["Gqr"] * (*GammaGqr)["Hqr"];
}

void ProjectorEnergyMatrixFromCoulombVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(getTensorArgument<complex, 
			       DryTensor<complex>>("CoulombVertex"));

  DryTensor<complex> conjGammaGqr(*GammaGqr);

  LOG(1, "EnergyMatrix") << "Computing projector energy matrix from Coulomb vertex with NG=" << GammaGqr->lens[0] << std::endl;

  DryMatrix<complex> *energyMatrix = new DryMatrix<complex>(
    GammaGqr->lens[0], GammaGqr->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "EnergyMatrix", energyMatrix
  );
}

