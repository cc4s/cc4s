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
  // read the Coulomb vertex GammaGpq
  Tensor<complex> *GammaGpq( getTensorArgument<complex>("CoulombVertex"));

  // construct the conjugate of the Coulomb vertex GammaGpq
  Tensor<complex> conjGammaGpq(*GammaGpq);
  conjugate(conjGammaGpq);

  // compute the energy matrix
  Matrix<complex> *energyMatrix = new Matrix<complex>(
    GammaGpq->lens[0], GammaGpq->lens[0], *GammaGpq->wrld
  );
  allocatedTensorArgument<complex>("EnergyMatrix", energyMatrix);

  LOG(1, "EnergyMatrix") << "Computing energy matrix from Coulomb vertex" << GammaGpq->get_name() 
			 << ", with NG=" << GammaGpq->lens[0] << std::endl;

  (*energyMatrix)["GH"] = conjGammaGpq["Gqr"] * (*GammaGpq)["Hqr"];
}

void ProjectorEnergyMatrixFromCoulombVertex::dryRun() {
  // Read the Coulomb vertex GammaGpq
  DryTensor<complex> *GammaGpq(getTensorArgument<complex, 
			       DryTensor<complex>>("CoulombVertex"));

  DryTensor<complex> conjGammaGpq(*GammaGpq);

  LOG(1, "EnergyMatrix") << "Computing energy matrix from Coulomb vertex with NG=" << GammaGpq->lens[0] << std::endl;

  DryMatrix<complex> *energyMatrix = new DryMatrix<complex>(
    GammaGpq->lens[0], GammaGpq->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "EnergyMatrix", energyMatrix
  );
}

