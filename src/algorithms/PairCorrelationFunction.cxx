#include <algorithms/PairCorrelationFunction.hpp>

#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(PairCorrelationFunction);

PairCorrelationFunction::PairCorrelationFunction(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PairCorrelationFunction::~PairCorrelationFunction() {
}

void PairCorrelationFunction::run() {

  GtoRFourier();

}



void PairCorrelationFunction::GtoRFourier(){
  Tensor<> *ctfMomenta(getTensorArgument<>("GGrid"));
  Tensor<> *ctfStructureFactor(getTensorArgument<>("StructureFactor"));
  Tensor<> *ctfReciprocalLattice(getTensorArgument<>("ReciprocalLattice"));
  Tensor<> *ctfRealLattice(getTensorArgument<>("RealLattice"));

  std::vector<Vector<>> momenta;
  std::vector<double> structureFactor;
  std::vector<Vector<>> reciprocalLattice;
  std::vector<Vector<>> realLattice;

  int NG(ctfMomenta->lens[1]);
  momenta.resize(NG);
  structureFactor.resize(NG);
  reciprocalLattice.resize(3);
  realLattice.resize(3);

  ctfMomenta->read_all(momenta.data()->coordinate);
  ctfStructureFactor->read_all(structureFactor.data());
  ctfReciprocalLattice->read_all(reciprocalLattice.data()->coordinate);
  ctfRealLattice->read_all(realLattice.data()->coordinate);

  cc4s::Vector<> check_grid;
  for (int g(0); g<NG; ++g){
    check_grid+=momenta[g];
  }
  if (check_grid.length() > 1e-10){
    throw new EXCEPTION("Only full Grid implemented");
  }

  Vector<> directMin, directMax;
  for (int g(0); g < NG; ++g) {
    double directComponentx(momenta[g].dot(realLattice[0]));
    double directComponenty(momenta[g].dot(realLattice[1]));
    double directComponentz(momenta[g].dot(realLattice[2]));

    directMin[0] = std::min(directMin[0], directComponentx);
    directMax[0] = std::max(directMax[0], directComponentx);
    directMin[1] = std::min(directMin[1], directComponenty);
    directMax[1] = std::max(directMax[1], directComponenty);
    directMin[2] = std::min(directMin[2], directComponentz);
    directMax[2] = std::max(directMax[2], directComponentz);
  }

  Vector<int> boxDimension;
  for (int d(0);d < 3; ++d){
    boxDimension[d] = std::floor(directMax[d]-directMin[d] + 1.5);
  }
  int64_t boxSize(boxDimension[0]*boxDimension[1]*boxDimension[2]);

  int minx((int) (directMin[0]-0.5));   int maxx((int) (directMax[0]+0.5));
  int miny((int) (directMin[1]-0.5));   int maxy((int) (directMax[1]+0.5));
  int minz((int) (directMin[2]-0.5));   int maxz((int) (directMax[2]+0.5));
  std::vector<double> realStructureFactor;
  std::vector<Vector<double>> realSpaceMesh;
  realStructureFactor.resize(boxSize);
  realSpaceMesh.resize(boxSize);
  int u(0);
  // set up real space mesh
  for ( int i(minx) ; i <= maxx ; ++i ){
    for ( int j(miny) ; j <= maxy ; ++j ){
      for ( int k(minz) ; k <= maxz ; ++k){
	realSpaceMesh[u][0]  = static_cast<double>(i)/static_cast<double>(boxDimension[0])*realLattice[0][0];
	realSpaceMesh[u][0] += static_cast<double>(j)/static_cast<double>(boxDimension[1])*realLattice[1][0];
	realSpaceMesh[u][0] += static_cast<double>(k)/static_cast<double>(boxDimension[2])*realLattice[2][0];
	realSpaceMesh[u][1]  = static_cast<double>(i)/static_cast<double>(boxDimension[0])*realLattice[0][1];
	realSpaceMesh[u][1] += static_cast<double>(j)/static_cast<double>(boxDimension[1])*realLattice[1][1];
	realSpaceMesh[u][1] += static_cast<double>(k)/static_cast<double>(boxDimension[2])*realLattice[2][1];
	realSpaceMesh[u][2]  = static_cast<double>(i)/static_cast<double>(boxDimension[0])*realLattice[0][2];
	realSpaceMesh[u][2] += static_cast<double>(j)/static_cast<double>(boxDimension[1])*realLattice[1][2];
	realSpaceMesh[u][2] += static_cast<double>(k)/static_cast<double>(boxDimension[2])*realLattice[2][2];
	u++;
      }						
    }
  }

  LOG(0,"Dimensions:") << "x: " << minx << " " << maxx << std::endl;
  LOG(0,"Dimensions:") << "y: " << miny << " " << maxz << std::endl;
  LOG(0,"Dimensions:") << "z: " << minz << " " << maxz << std::endl;
  LOG(0,"Dimensions:") << "#X*#Y*#Z = " << maxx-minx+1 << " " << maxy-miny+1 << " " << maxz-minz +1 << std::endl;
  LOG(0,"Dimensions:") << "Real space mesh:" << std::endl;
  LOG(0,"Dimensions:") << realLattice[0][0]/static_cast<double>(boxDimension[0]) << " "
		       << realLattice[0][1]/static_cast<double>(boxDimension[0]) << " "
		       << realLattice[0][2]/static_cast<double>(boxDimension[0]) << std::endl;
  LOG(0,"Dimensions:") << realLattice[1][0]/static_cast<double>(boxDimension[1]) << " "
		       << realLattice[1][1]/static_cast<double>(boxDimension[1]) << " "
		       << realLattice[1][2]/static_cast<double>(boxDimension[1]) << std::endl;
  LOG(0,"Dimensions:") << realLattice[2][0]/static_cast<double>(boxDimension[2]) << " "
		       << realLattice[2][1]/static_cast<double>(boxDimension[2]) << " "
		       << realLattice[2][2]/static_cast<double>(boxDimension[2]) << std::endl;

  for ( int i(0); i<boxSize; ++i ) {
    for ( int j(0); j<NG; ++j) {
      double phase(2.*M_PI*realSpaceMesh[i].dot(momenta[j]));
      realStructureFactor[i] += cos(phase)*structureFactor[j];
    }
    LOG(2,"Structure") << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << realStructureFactor[i] << std::endl;
  }
}

