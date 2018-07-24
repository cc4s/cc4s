#include <algorithms/BasisSetExtrapolation.hpp>

#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <tcc/DryTensor.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>
#include <util/MpiCommunicator.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(BasisSetExtrapolation);

BasisSetExtrapolation::BasisSetExtrapolation(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

BasisSetExtrapolation::~BasisSetExtrapolation() {
}

void BasisSetExtrapolation::run() {

  bool lHaveReciprocalMesh(false);
  bool lHaveRealSpaceMesh(false);

  int fFourierSF(getIntegerArgument("FourierStructureFactor",0));
  if ( fFourierSF > 0){
    if (!lHaveReciprocalMesh) {
      readReciprocalGridFromFile();
      lHaveReciprocalMesh = true;
    }
    if (!lHaveRealSpaceMesh) {
      constructRealSpaceMesh(1,realSpaceMesh);
      lHaveRealSpaceMesh = true;
    }
    LOG(0,"run FourierStructureFactor") << std::endl;
    GtoRFourier();
  }

  int fFourierPCF(getIntegerArgument("FourierPairCorrelationFunction",0));
  if ( fFourierPCF > 0){
    LOG(0,"run FourierPairCorrelationFunction") << std::endl;
    if (!lHaveReciprocalMesh) {
      readReciprocalGridFromFile();
      lHaveReciprocalMesh = true;
    }
    if (!lHaveRealSpaceMesh) {
      constructRealSpaceMesh(1,realSpaceMesh);
      lHaveRealSpaceMesh = true;
    }
    RtoGFourier();
  }

  int fFourierCompleteness(getIntegerArgument("FourierCompleteness",0));
  if ( fFourierCompleteness > 0){
    if (!lHaveReciprocalMesh) {
      readReciprocalGridFromFile();
      lHaveReciprocalMesh = true;
    }
    if (!lHaveRealSpaceMesh) {
      constructRealSpaceMesh(1,realSpaceMesh);
      lHaveRealSpaceMesh = true;
    }
    fourierCompleteness();
  }

  int fQGG(getIntegerArgument("calculateQGG",0));
  if(fQGG == 1){
    int iStart(getIntegerArgument("iStart",-1));
    int iEnd(getIntegerArgument("iEnd",-1));
    QGG(iStart,iEnd);
  }

  int fFitF12(getIntegerArgument("fitF12",-1));
  if (fFitF12 >= 0){
    real minG(getRealArgument("minG",-1));
    real maxG(getRealArgument("maxG",-1));
    if ( minG > maxG || minG < 0. ) throw new EXCEPTION("need fitting range:minG and maxG");
    fitF12(1,fFitF12,minG,maxG);
  }

  int fFitF12Slater(getIntegerArgument("fitF12Slater",-1));
  if (fFitF12Slater >= 0){
    real minG(getRealArgument("minG",-1));
    real maxG(getRealArgument("maxG",-1));
    if ( minG > maxG || minG < 0. ) throw new EXCEPTION("need fitting range:minG and maxG");
    fitF12(2,fFitF12Slater,minG,maxG);
  }

  int CGi(getIntegerArgument("CGi",0));
  //not really important right now. but maybe we will continue on it!
  if (CGi > 0){
    realCGi();
  }
}

void BasisSetExtrapolation::dryRun(){
 int fQGG(getIntegerArgument("calculateQGG",0));
 if (fQGG > 0){
   int iStart(getIntegerArgument("iStart",-1));
   int iEnd(getIntegerArgument("iEnd",-1));
   dryQGG(iStart,iEnd);
 }
}

void BasisSetExtrapolation::GtoRFourier(){

  NR = realSpaceMesh.size();
  std::vector<real> localPairCorrelationFunction(NR);
  std::vector<real> localdrPairCorrelationFunction(NR);
  std::vector<real> localddrPairCorrelationFunction(NR);
  std::vector<real> localdddrPairCorrelationFunction(NR);
  std::vector<real> localddddrPairCorrelationFunction(NR);
  std::vector<real> localdddddrPairCorrelationFunction(NR);
  MpiCommunicator communicator(
    Cc4s::world->rank, Cc4s::world->np, Cc4s::world->comm
  );
  for ( int j(communicator.getRank()); j<NG; j+=communicator.getProcesses() ) {
    for ( int i(0); i < NR; ++i) {
      real phase(realSpaceMesh[i].dot(reciprocalGrid[j]));
      localPairCorrelationFunction[i] += cos(phase)*structureFactor[j];
      real deri(phase/realSpaceMesh[i].length());
      localdrPairCorrelationFunction[i] -= sin(phase)*structureFactor[j]*deri;
      localddrPairCorrelationFunction[i] -= cos(phase)*structureFactor[j]*deri*deri;
      localdddrPairCorrelationFunction[i] += sin(phase)*structureFactor[j]*deri*deri*deri;
      localddddrPairCorrelationFunction[i] += cos(phase)*structureFactor[j]*deri*deri*deri*deri;
      localdddddrPairCorrelationFunction[i] -= sin(phase)*structureFactor[j]*deri*deri*deri*deri*deri;
    }
  }

  pairCorrelationFunction.resize(NR);
  std::vector<real> drPairCorrelationFunction(NR);
  std::vector<real> ddrPairCorrelationFunction(NR);
  std::vector<real> dddrPairCorrelationFunction(NR);
  std::vector<real> ddddrPairCorrelationFunction(NR);
  std::vector<real> dddddrPairCorrelationFunction(NR);
  for ( int i(0); i < NR; ++i) {
    communicator.allReduce(localPairCorrelationFunction[i], pairCorrelationFunction[i]);
    communicator.allReduce(localdrPairCorrelationFunction[i], drPairCorrelationFunction[i]);
    communicator.allReduce(localddrPairCorrelationFunction[i], ddrPairCorrelationFunction[i]);
    communicator.allReduce(localdddrPairCorrelationFunction[i], dddrPairCorrelationFunction[i]);
    communicator.allReduce(localddddrPairCorrelationFunction[i], ddddrPairCorrelationFunction[i]);
    communicator.allReduce(localdddddrPairCorrelationFunction[i], dddddrPairCorrelationFunction[i]);
  }

  for ( int i(0); i<NR; ++i ) {
    LOG(2,"Structure") << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << pairCorrelationFunction[i] << std::endl;
    LOG(2,"drS")       << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << drPairCorrelationFunction[i] << std::endl;
    LOG(2,"d2rS")      << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << ddrPairCorrelationFunction[i] << std::endl;
    LOG(2,"d3rS")      << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << dddrPairCorrelationFunction[i] << std::endl;
    LOG(2,"d4rS")      << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << ddddrPairCorrelationFunction[i] << std::endl;
    LOG(2,"d5rS")      << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << dddddrPairCorrelationFunction[i] << std::endl;
  }

  auto ctfRealSpaceMesh(
    new CTF::Tensor<>(2, std::vector<int>( {3,NR} ).data() )
  );
  auto ctfPairCorrelationFunction(new CTF::Vector<>(NR));
  std::vector<int64_t> indices2(ctfRealSpaceMesh->wrld->rank == 0 ? 3*NR : 0);
  std::vector<int64_t> indices(ctfPairCorrelationFunction->wrld->rank == 0 ? NR : 0);
  for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; }
  for (size_t i(0); i < indices2.size(); ++i) { indices2[i] = i; }
  ctfRealSpaceMesh->write(
    indices2.size(), indices2.data(), realSpaceMesh.data()->coordinate
  );
  ctfPairCorrelationFunction->write(
    indices.size(), indices.data(), pairCorrelationFunction.data()
  );
  allocatedTensorArgument<>("PairCorrelationFunction",ctfPairCorrelationFunction);
  allocatedTensorArgument<>("RealSpaceMesh",ctfRealSpaceMesh);
}

void BasisSetExtrapolation::RtoGFourier(){
  NR = realSpaceMesh.size();
  std::vector<real> localffStructureFactor(NG);

  MpiCommunicator communicator(
    Cc4s::world->rank, Cc4s::world->np, Cc4s::world->comm
  );
  for ( int j(communicator.getRank()); j<NR; j+=communicator.getProcesses() ) {
    //    if ( realSpaceMesh[i].length() < 3.){
    for ( int i(0); i < NG ; ++i) {
      real phase(realSpaceMesh[j].dot(reciprocalGrid[i]));
      localffStructureFactor[i] += cos(phase)*pairCorrelationFunction[j];
    }
  }

  std::vector<real> ffStructureFactor(NG);
  for (int i(0); i < NG; ++i) {
    communicator.allReduce(localffStructureFactor[i], ffStructureFactor[i]);
  }

  for ( int i(0); i<NG; ++i ) {
    LOG(2,"ffstructure") << reciprocalGrid[i][0] << " " << reciprocalGrid[i][1] << " "
		         << reciprocalGrid[i][2] << " "
			 << ffStructureFactor[i]/static_cast<real>(NR) << std::endl;
  } 

}

void BasisSetExtrapolation::constructRealSpaceMesh(int augmentationFactor,std::vector<Vector<>> &realSpaceMesh){
  Tensor<> *ctfReciprocalLattice(getTensorArgument<>("ReciprocalLattice"));
  Tensor<> *ctfRealLattice(getTensorArgument<>("RealLattice"));

  std::vector<Vector<>> reciprocalLattice;
  std::vector<Vector<>> realLattice;

  reciprocalLattice.resize(3);
  realLattice.resize(3);

  ctfReciprocalLattice->read_all(reciprocalLattice.data()->coordinate);
  ctfRealLattice->read_all(realLattice.data()->coordinate);

  Vector<> directMin, directMax;
  for (int g(0); g < NG; ++g) {
    real directComponentx(reciprocalGrid[g].dot(realLattice[0])/2./Pi<real>());
    real directComponenty(reciprocalGrid[g].dot(realLattice[1])/2./Pi<real>());
    real directComponentz(reciprocalGrid[g].dot(realLattice[2])/2./Pi<real>());
    
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

  //int64_t boxSize(boxDimension[0]*boxDimension[1]*boxDimension[2]);

  int minx((int) (directMin[0]-0.5));   int maxx((int) (directMax[0]+0.5));
  int miny((int) (directMin[1]-0.5));   int maxy((int) (directMax[1]+0.5));
  int minz((int) (directMin[2]-0.5));   int maxz((int) (directMax[2]+0.5));

  minx *= augmentationFactor;
  maxx *= augmentationFactor;
  miny *= augmentationFactor;
  maxy *= augmentationFactor;
  minz *= augmentationFactor;
  maxz *= augmentationFactor;
  boxDimension[0] *= augmentationFactor*augmentationFactor;
  boxDimension[1] *= augmentationFactor*augmentationFactor;
  boxDimension[2] *= augmentationFactor*augmentationFactor;

  NR =(maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1);
  realSpaceMesh.resize(NR);

  int u(0);
  // set up real space mesh
  for ( int i(minx) ; i <= maxx ; ++i ){
    for ( int j(miny) ; j <= maxy ; ++j ){
      for ( int k(minz) ; k <= maxz ; ++k){
	for ( int xyz(0) ; xyz < 3; ++xyz){
	  realSpaceMesh[u][xyz]  = static_cast<real>(i)/
	    static_cast<real>(boxDimension[0])*realLattice[0][xyz];
	  realSpaceMesh[u][xyz] += static_cast<real>(j)/
	    static_cast<real>(boxDimension[1])*realLattice[1][xyz];
	  realSpaceMesh[u][xyz] += static_cast<real>(k)/
	    static_cast<real>(boxDimension[2])*realLattice[2][xyz];
	}
	u++;
      }						
    }
  }

  LOG(0,"Dimensions:") << "x: " << minx << " " << maxx << std::endl;
  LOG(0,"Dimensions:") << "y: " << miny << " " << maxz << std::endl;
  LOG(0,"Dimensions:") << "z: " << minz << " " << maxz << std::endl;
  LOG(0,"Dimensions:") << "#X*#Y*#Z = " << maxx-minx+1 << " * " << maxy-miny+1
		       << " * " << maxz-minz +1 << " = " << NR << std::endl;
  LOG(0,"Dimensions:") << "Real space mesh:" << std::endl;
  LOG(0,"Dimensions:") << realLattice[0][0]/static_cast<real>(boxDimension[0]) << " "
		       << realLattice[0][1]/static_cast<real>(boxDimension[0]) << " "
		       << realLattice[0][2]/static_cast<real>(boxDimension[0]) << std::endl;
  LOG(0,"Dimensions:") << realLattice[1][0]/static_cast<real>(boxDimension[1]) << " "
		       << realLattice[1][1]/static_cast<real>(boxDimension[1]) << " "
		       << realLattice[1][2]/static_cast<real>(boxDimension[1]) << std::endl;
  LOG(0,"Dimensions:") << realLattice[2][0]/static_cast<real>(boxDimension[2]) << " "
		       << realLattice[2][1]/static_cast<real>(boxDimension[2]) << " "
		       << realLattice[2][2]/static_cast<real>(boxDimension[2]) << std::endl;
}

void BasisSetExtrapolation::readReciprocalGridFromFile(){
  Tensor<> *ctfReciprocalGrid(getTensorArgument<>("Momenta"));
  Tensor<> *ctfStructureFactor(getTensorArgument<>("StructureFactor"));

  NG = ctfStructureFactor->lens[0];
  reciprocalGrid.resize(NG);
  structureFactor.resize(NG);
  ctfReciprocalGrid->read_all(reciprocalGrid.data()->coordinate);
  ctfStructureFactor->read_all(structureFactor.data());

  cc4s::Vector<> check_grid;
  for (int g(0); g<NG; ++g){
    check_grid+=reciprocalGrid[g];
  }
  // for convenience always work with the full grid 
  if (check_grid.length() > 1e-10){ 
    for (int g(1); g<NG; ++g){
      structureFactor.push_back(structureFactor[g]);
      Vector<> minusG;
      minusG -= reciprocalGrid[g];
      reciprocalGrid.push_back(minusG);
    }
    NG = NG*2-1;
    for (int g(0); g<NG; ++g){
      structureFactor[g] *= 0.5;
    }
  }
  // work with the real G-vectors
  for (int g(0); g<NG; ++g){
    reciprocalGrid[g] *= 2.*Pi<real>();
  }

}

void BasisSetExtrapolation::fourierCompleteness(){

  Tensor<> *ctfBasisSetCompletenessNumber(getTensorArgument<>("BasisSetCompletenessNumber"));
  std::vector<real> basisSetCompletenessNumber(NG);
  ctfBasisSetCompletenessNumber->read_all(basisSetCompletenessNumber.data());
  int NGcheck(ctfBasisSetCompletenessNumber->lens[0]);

  basisSetCompletenessNumber[0] = 1.;

  if ( NG == NGcheck*2-1){
    // blow up half grid // 0th element G=0
    for (int g(1); g<NGcheck; ++g){
      basisSetCompletenessNumber[NGcheck+g-1] = basisSetCompletenessNumber[g];
    }
    for (int g(1); g<NG; ++g){
      basisSetCompletenessNumber[g] /= 2.;
    }
  }
  else if ( NGcheck != NG){
    throw new EXCEPTION("Dimension Problems!");
  }

  std::vector<real> localRealSpaceBasisSetCompleteness(NR);

  MpiCommunicator communicator(
    Cc4s::world->rank, Cc4s::world->np, Cc4s::world->comm
  );
  for ( int j(communicator.getRank()); j<NG; j+=communicator.getProcesses() ) {
    for ( int i(0); i < NR ; ++i) {
      real phase(realSpaceMesh[i].dot(reciprocalGrid[j]));
      localRealSpaceBasisSetCompleteness[i] += cos(phase)*basisSetCompletenessNumber[j];
    }
  }

  realSpaceBasisSetCompleteness.resize(NR);

  for ( int i(0); i < NR; ++i) {
    communicator.allReduce(localRealSpaceBasisSetCompleteness[i], realSpaceBasisSetCompleteness[i]);
  }

  for (int i(0); i<NR; ++i){
    LOG(2,"completeR") << realSpaceMesh[i][0] << " " << realSpaceMesh[i][1] << " "
		       << realSpaceMesh[i][2] << " " << realSpaceBasisSetCompleteness[i] << std::endl;
  }

}

void BasisSetExtrapolation::inverseConvolution(){
  Tensor<> *ctfBasisSetCompletenessNumber(getTensorArgument<>("BasisSetCompletenessNumber"));
  Tensor<> *ctfCorrectedStructureFactor(getTensorArgument<>("StructureFactor"));

  if ( ctfBasisSetCompletenessNumber->lens[0] !=
       ctfCorrectedStructureFactor->lens[0]) throw new EXCEPTION("Dimension problems!");

  Transform<real, real>(
    std::function<void(real, real &)>(
      [](real cN, real &s) {
	cN == 0 ? s = 0. : s /= cN;  }
    )
  )(
    (*ctfBasisSetCompletenessNumber)["G"],
    (*ctfCorrectedStructureFactor)["G"]
  );

  allocatedTensorArgument<>("CorrectedStructureFactor",ctfCorrectedStructureFactor);  
}

void BasisSetExtrapolation::QGG(int iStart, int iEnd){

  PTR(Tensor<complex>) GammaGai;
  GammaGai = NEW(Tensor<complex>,getTensorArgument<complex>("ParticleHoleCoulombVertex"));

  int NF(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);
 
  Tensor<> *ctfCoulombKernel(getTensorArgument<>("CoulombKernel"));

  int NFF[] = {NF};
  Tensor<complex> invSqrtVG(1, NFF);

  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
    [](real vG, complex &invVG){
      if ( std::isinf(vG) ){
        invVG = 0.;
      }
      else{
        invVG = std::sqrt(1./vG);
      }
    }
   )
  )(
    (*ctfCoulombKernel)["G"],invSqrtVG["G"]
  );

  (*GammaGai)["Gai"] *= invSqrtVG["G"];
  auto conjCGai(new Tensor<complex>(*GammaGai));
  Univar_Function<complex> fConj(conj<complex>);

  conjCGai->sum(1.0, *GammaGai, "Gai", 0.0, "Gai", fConj);

  if (iStart < 0 || iStart > No){
   iStart = 0;
  }
  if (iEnd <= iStart || iEnd > No){
   iEnd = No;
  }
  LOG(0,"QGG Sliced") << "iStart: " << iStart << " iEnd: " << iEnd << std::endl;
  int CGaiStart[] = {0, 0, iStart};
  int CGaiEnd[] =   {NF, Nv, iEnd};
  auto CGaiSliced(new Tensor<complex>(GammaGai->slice(CGaiStart,CGaiEnd)));
  auto conjCGaiSliced(new Tensor<complex>(conjCGai->slice(CGaiStart,CGaiEnd)));
  CGaiSliced->set_name("CGaiSliced");
  conjCGaiSliced->set_name("conjCGaiSliced");

  int NGG[] = {NF, NF};
  int NFG[] = {NF, NF, iEnd-iStart, iEnd-iStart}; 

  auto QGG( new CTF::Tensor<complex>(2,NGG));

  PTR(Tensor<complex>) FGij;
  PTR(Tensor<complex>) FGji;
  PTR(Tensor<complex>) FGone;
  PTR(Tensor<complex>) FGtwo;

  FGij = NEW(Tensor<complex>,4,NFG);
  FGji = NEW(Tensor<complex>,4,NFG);
  FGone = NEW(Tensor<complex>,2,NGG);
  FGtwo = NEW(Tensor<complex>,2,NGG);

  (*FGone)["GF"] =  (*conjCGaiSliced)["Gai"] * (*CGaiSliced)["Fai"];
  (*FGtwo)["GF"] =  (*conjCGaiSliced)["Fbj"] * (*CGaiSliced)["Gbj"];
  (*QGG)["GF"]   = (2.0) * (*FGone)["GF"] * (*FGtwo)["GF"];

  (*FGij)["GFij"] = (*conjCGaiSliced)["Gai"] * (*CGaiSliced)["Faj"];
  (*FGji)["GFij"] = (*conjCGaiSliced)["Fbi"] * (*CGaiSliced)["Gbj"];
  (*QGG)["GF"] += (-1.0) * (*FGij)["GFij"] * (*FGji)["GFij"];

  auto realQGG(new CTF::Tensor<>(2,NGG));
  auto absQGG(new CTF::Tensor<>(2,NGG));

  fromComplexTensor(*QGG,*realQGG);

  CTF::Transform<complex, double>(
    std::function<void(complex, double &)>(
      [](complex qgg, double &absqgg){
       absqgg = static_cast<double>(std::abs(qgg));
      }
    )
  )(
   (*QGG)["GF"], (*absQGG)["GF"]
  );

  allocatedTensorArgument<complex>("QGG",QGG);
  allocatedTensorArgument<>("AbsQGG", absQGG);
  allocatedTensorArgument<>("RealQGG", realQGG);

}

void BasisSetExtrapolation::calculateNewSF(
  int type, real gamma, Tensor<> *coulombKernel, Tensor<> *newSF, Tensor<> *resNewSF){

  CTF::Tensor<complex> *QGG(getTensorArgument<complex>("QGG"));
  int NG(QGG->lens[0]);
  int NFF[] = {NG};
  
  auto cK(new Tensor<>(1,NFF));
  (*cK) = (*coulombKernel);

  real volume(getRealArgument("volume",-1.));

  auto reciprocalYC(new Tensor<complex>(1,NFF));

  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [&type, &volume, &gamma](real cK, complex &YC){
        if ( cK == 0 ){
          YC = 0.;
        }
        else{
          if (type ==1){
            YC = 1./cK/cK; //*volume/4.5835494674469; //Ha(eV)*a0(A)*4*pi/(2*pi)**2
  //        rYC = rYC/150.4121/(gamma*gamma+150.4121/rYC); hbar*hbar/2*me*10**20meter*(2pi)**2 in eV
  //aG      rYC = constantFactor/180.956*4*Pi<>()/(2*Pi<>())/(2*Pi<>())*2*rYC/(gamma*gamma+1.0/rYC);
  //aI      rYC = rYC*4.5835494674469/volume/(gamma*gamma+2.*150.4121/rYC);   // GAMMA[Energy]
            YC *= (-0.015237)/volume/(gamma*gamma+1./YC);               // GAMMA[1/A]
          }
          else if (type ==2){
            YC = (-0.015237)/volume/(gamma*gamma+cK*cK)/(gamma*gamma+cK*cK);
          }
       }
      }
    )
  )(
   (*cK)["G"],(*reciprocalYC)["G"]
  );

  auto dReciprocalYC(new Tensor<complex>(1,NFF));

  CTF::Transform<real, complex, complex >(
    std::function<void(real, complex, complex &)>(
      [&type, &volume, &gamma](real cK, complex YC, complex &dYC){
        if ( cK == 0){
          dYC = 0.; 
        }
        else{
          if(type==1){
            dYC = (-2.)*YC*gamma/(gamma*gamma+cK*cK);
          }
          else if(type==2){
            dYC = (-4.)*YC*gamma/(gamma*gamma+cK*cK);
          }
        }
       }
     ) 
  )(
   (*cK)["G"],(*reciprocalYC)["G"],(*dReciprocalYC)["G"]
  );
  auto complexNewSF(new Tensor<complex>(1,NFF));
  auto complexResNewSF(new Tensor<complex>(1,NFF));

  (*complexNewSF)["G"] = (*QGG)["GF"] * (*reciprocalYC)["F"];
  (*complexResNewSF)["G"] = (*QGG)["GF"] * (*dReciprocalYC)["F"];

  fromComplexTensor(*complexNewSF,*newSF);
  fromComplexTensor(*complexResNewSF,*resNewSF);

}

void BasisSetExtrapolation::fitF12(int type,int iter, real minG, real maxG){

  real volume(getRealArgument("volume",-1));
  if (volume < 0. ) throw new EXCEPTION("Set volume");

  CTF::Tensor<> *structureFactor(getTensorArgument<>("StructureFactor"));
  CTF::Tensor<> *coulombKernel(getTensorArgument<>("CoulombKernel"));

  int NG(coulombKernel->lens[0]);
  int NFF[] = {NG};

  auto resNewSF(new Tensor<>(1,NFF));
  auto newSF(new Tensor<>(1,NFF));

  CTF::Transform<real>(
    std::function<void(real &)>(
      [&volume](real &cK){
        if( std::isinf(cK)){
          cK = 0.;
        }
        else{
          cK *= volume/4.5835494674469;
          cK = 1./std::sqrt(cK);
        }
      }
    )
  )(
   (*coulombKernel)["G"]
  );

  real gamma(getRealArgument("gamma",1));

  for (int i(0); i<=iter; ++i){  

    calculateNewSF(type,gamma,coulombKernel,newSF,resNewSF);

    auto dummy(new Tensor<>(1,NFF));
    int counter(0);
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&counter, &maxG, &minG](real cK, real res, real &dummy){
          if ( cK > maxG || cK < minG ){
            dummy = 0.;
          }
          else{
            dummy = res*res;
            counter++;
          }
        }
     )
    )(
     (*coulombKernel)["G"],(*resNewSF)["G"],(*dummy)["G"]
    );

    CTF::Scalar<double> denom;
    denom[""] = (*dummy)["G"];
    
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [](real newsf, real sf, real &dummy){
          dummy = newsf - sf;
        }
      )
    )(
     (*newSF)["G"],(*structureFactor)["G"],(*dummy)["G"]
    );

    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&maxG, &minG](real cK, real resnewsf, real &dummy){
          if ( cK > maxG || cK < minG ){
            dummy = 0.;
          }
          else{
            dummy *= resnewsf;
          }
        }
      )
    )(
     (*coulombKernel)["G"],(*resNewSF)["G"],(*dummy)["G"]
    );

    CTF::Scalar<double> numer;
    numer[""] = (*dummy)["G"];
    gamma -= numer.get_val()/denom.get_val();

    (*resNewSF) = (*structureFactor);
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&maxG, &minG](real cK, real newsf, real &resnewsf){
          if ( cK > maxG || cK < minG){
            resnewsf = 0.;
          }
          else{
            resnewsf -= newsf;
          }
         }
      )
    )(
     (*coulombKernel)["G"],(*newSF)["G"],(*resNewSF)["G"]
    );

    LOG(0,"norm") << gamma << " " << resNewSF->norm2() << std::endl;

  }
  allocatedTensorArgument<>("ResNewSF",resNewSF);
  allocatedTensorArgument<>("NewSF",newSF);
}

void BasisSetExtrapolation::dryQGG(int iStart, int iEnd){

  DryTensor<> *epsi(getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));

  int No(epsi->lens[0]);
  if (iStart >= 0 && iEnd <= No && iStart < iEnd){
    No = iEnd - iStart;
  }
  DryTensor<complex> *GammaGai(getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex"));

  int NF(GammaGai->lens[0]);
  int NFF[] = {NF, NF};
  int NGG[] = {NF, NF, No};
  int NFG[] = {NF, NF, No, No};
  int syms0[] = {NS, NS};
  int syms1[] = {NS, NS, NS};
  int syms2[] = {NS, NS, NS, NS}; 

  DryTensor<complex> FGijd(4,NFG,syms2);
  DryTensor<complex> FGij(4,NFG,syms2);
  DryTensor<complex> FGji(4,NFG,syms2);
  DryTensor<complex> FGi(3,NGG,syms1);

  allocatedTensorArgument("QGG", new DryTensor<>(2, NFF, syms0, SOURCE_LOCATION));

}

void BasisSetExtrapolation::realCGi(){

  Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));

  // Get the Particle Hole Coulomb Vertex
  No = GammaFqr->lens[2];
  Np = GammaFqr->lens[1];
  NF = GammaFqr->lens[0];

  LOG(0,"dimensions:") << No << " <-No " << Np << " <- Np " << NF << " <-NF" << std::endl;
  int FzeropStart[] = {0, 0, 0};
  int FzeropEnd[]   = {NF, Np, 1};  

  Tensor<complex> GammaGzerop(GammaFqr->slice(FzeropStart, FzeropEnd));

  int NNF[] = {NF};
  Tensor<> *ctfCoulombKernel(getTensorArgument<>("CoulombKernel"));
  auto invSqrtVG(new Tensor<complex>(1,NNF));

  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
    [](real vG, complex &invVG){
      if ( std::isinf(vG) ){
        invVG = 0.;
      }
      else{
        invVG = std::sqrt(1./vG);
      }
    }
   )
  )(
    (*ctfCoulombKernel)["G"],(*invSqrtVG)["G"]
  ); 

  Tensor<complex> CGzerop(false,GammaGzerop);

  CGzerop["Gai"] = GammaGzerop["Gai"] * (*invSqrtVG)["G"];

  auto realCGzerop(new CTF::Tensor<> (3,FzeropEnd));
  auto imagCGzerop(new CTF::Tensor<> (3,FzeropEnd));

  fromComplexTensor(CGzerop,*realCGzerop,*imagCGzerop);

  allocatedTensorArgument<>("RealCGp",realCGzerop);
  allocatedTensorArgument<>("ImagCGp",imagCGzerop);

}
