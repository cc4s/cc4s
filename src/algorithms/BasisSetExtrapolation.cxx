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
  
  int fbasisSetExtrapolation(getIntegerArgument("basisSetExtrapolation",0));
  if (fbasisSetExtrapolation > 0) {
    LOG(0,"run BasisSetExtrapolation") << std::endl;
    if (!lHaveReciprocalMesh) {
      readReciprocalGridFromFile();
      lHaveReciprocalMesh = true;
    }
    real minG(getRealArgument("minG",-2.));
    real maxG(getRealArgument("maxG",-1.));
    if ( (minG > maxG) || (minG <=  0.) ) fbasisSetExtrapolation = 0;
    extrapolation(minG,maxG,fbasisSetExtrapolation);
  }
  
  int fbasisSetCompleteness(getIntegerArgument("basisSetCompleteness",0));
  if (fbasisSetCompleteness == 1) {
    LOG(0,"run basisSetCompleteness full FTODDUMP") << std::endl;
    basisSetCompletenessFullFTODDUMP();
  }
  else if ( fbasisSetCompleteness == 2){
    LOG(0,"run basiSetCompleteness FTODDUMPIA") << std::endl;
    basisSetCompletenessFTODDUMPIA();
  }
 
  //  int fInverseConvolution(getIntegerArgument("InverseConvolution",0));
  //  if ( fInverseConvolution > 0){
  //    inverseConvolution();
  //  }


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

  int fQGG(getIntegerArgument("fQGG",0));
  if (fQGG == 1){
    QGG();
  }
  else if(fQGG == 2){
    int iStart(getIntegerArgument("iStart",-1));
    int iEnd(getIntegerArgument("iEnd",-1));
    QGGSliced(iStart,iEnd);
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
  if (CGi > 0){
    realCGi();
  }
}

void BasisSetExtrapolation::dryRun(){
 int fQGG(getIntegerArgument("QGG",0));
 if (fQGG > 0){
   dryQGG();
 }
}

void BasisSetExtrapolation::extrapolation(
  real Gmin, real Gmax, int itype){

  real constantFactor(getRealArgument("constantFactor",-1.));
  if ( constantFactor < 0.) throw new EXCEPTION("Set constant Factor");
  //  real augmentationFactor(getRealArgument("augmentationFactor",1.0));


  real gridEnergy(0.);
  
  for (int d(1); d < NG; ++d) {
    real invG(1/reciprocalGrid[d].length());
    gridEnergy  += invG*invG*structureFactor[d];
  }
  gridEnergy *=constantFactor;
  LOG(0,"ReferenceEnergy:") << gridEnergy << std::endl;

  /*

  real coeff(0.);
  std::vector<real> fitabsG; std::vector<real> fitSF;

  if ( itype > 0){
    for (int d(1); d < NG; ++d) {
      real fac(1.);
      fac = simplestWindow(Gmin,Gmax,cartesianGrid[d].l);
      longRangeEnergyExplicit += fac*cartesianGrid[d].vg * cartesianGrid[d].s;
      if ( (fac > 0.) && ( fac < 1.)){
	fitabsG.push_back(cartesianGrid[d].l);
	fitSF.push_back(cartesianGrid[d].s);
      }
    }
    LOG(2,"LongRangeEnergyExplicit:") << longRangeEnergyExplicit << std::endl;
    if ( fitSF.size() == 0) throw new EXCEPTION("No |G|-values in fitting range");
    //  for (unsigned int d(0); d< fitSF.size(); ++d){
    //  LOG(2,"PrepFit") << fitabsG[d] << " " << fitSF[d] << std::endl;
    //}

    coeff = leastSquareFit(fitabsG,fitSF);

    LOG(2,"Fitting") << "coeff: " << coeff << std::endl;
    real residuum(0.);
    for (unsigned int d(0); d < fitabsG.size(); ++d){
      real res;
      res = coeff/pow(fitabsG[d],4) - fitSF[d];
      residuum += res*res;
    }
    LOG(2,"Fitting") << "Residuum: " << residuum << std::endl;
  }
  
  // Scheme 1: Fill the full reciprocal mesh with perfect 1/G**4 values

  std::vector<real> extrapolatedSF;
  std::vector<Vector<>> extrapolatedGrid;


  if ( itype == 0) {
    LOG(0,"IO") << "Writing unmodified Grid" << std::endl;
    for (int d(0); d < NG; ++d){
      extrapolatedSF.push_back(cartesianGrid[d].s);
      extrapolatedGrid.push_back(cartesianGrid[d].v);
    }
    auto ctfExtrapolatedSF(
      new CTF::Tensor<>(1, std::vector<int>( {NG}).data())
    );
    std::vector<int64_t> indices(ctfExtrapolatedSF->wrld->rank == 0 ? NG : 0);
    for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; }

    ctfExtrapolatedSF->write(
      indices.size(), indices.data(), extrapolatedSF.data()
    );

    allocatedTensorArgument<>("ExtrapolatedSF", ctfExtrapolatedSF);

    auto ctfExtrapolatedGrid(
      new CTF::Tensor<>(2, std::vector<int>( {3,NG} ).data())
    );
    std::vector<int64_t> indices3(ctfExtrapolatedGrid->wrld->rank == 0 ? 3*NG : 0);
    for (size_t i(0); i < indices3.size(); ++i) { indices3[i] = i; }

    ctfExtrapolatedGrid->write(
      indices3.size(), indices3.data(), extrapolatedGrid.data()->coordinate
    );

    allocatedTensorArgument<>("ExtrapolatedGrid", ctfExtrapolatedGrid);
  }

  if ( itype > 0) {
    for (int d(0); d < NG; ++d){
      if ( cartesianGrid[d].l < Gmax){
	extrapolatedSF.push_back(cartesianGrid[d].s);
	extrapolatedGrid.push_back(cartesianGrid[d].v);
      }
      else{
	extrapolatedSF.push_back(coeff/pow(cartesianGrid[d].l,4));
	extrapolatedGrid.push_back(cartesianGrid[d].v);
      }
    }
    real basisSetExtrapolatedEnergy(0.);


    for (int d(1); d < NG; ++d){
      basisSetExtrapolatedEnergy += cartesianGrid[d].vg*extrapolatedSF[d];
    }
    LOG(0,"ExtrapolatedEnergy:") << basisSetExtrapolatedEnergy << std::endl;
  }

  if ( itype == 1) {
    LOG(0,"IO") << "Writing Extrapolated Grid" << std::endl;
    auto ctfExtrapolatedSF(
      new CTF::Tensor<>(1, std::vector<int>( {NG}).data())
    );
    std::vector<int64_t> indices(ctfExtrapolatedSF->wrld->rank == 0 ? NG : 0);
    for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; }

    ctfExtrapolatedSF->write(
      indices.size(), indices.data(), extrapolatedSF.data()
    );

    allocatedTensorArgument<>("ExtrapolatedSF", ctfExtrapolatedSF);

    auto ctfExtrapolatedGrid(
      new CTF::Tensor<>(2, std::vector<int>( {3,NG} ).data())
    );
    std::vector<int64_t> indices3(ctfExtrapolatedGrid->wrld->rank == 0 ? 3*NG : 0);
    for (size_t i(0); i < indices3.size(); ++i) { indices3[i] = i; }

    ctfExtrapolatedGrid->write(
      indices3.size(), indices3.data(), extrapolatedGrid.data()->coordinate
    );

    allocatedTensorArgument<>("ExtrapolatedGrid", ctfExtrapolatedGrid);
  }
  



  // Scheme 2: Augment the reciprocal lattice and calculate the extrapolated
  //           Structure Factor

  Tensor<> *ctfReciprocalLattice(getTensorArgument<>("ReciprocalLattice"));


  std::vector<Vector<>> reciprocalLattice(3);
  reciprocalLattice.resize(3);
  ctfReciprocalLattice->read_all(reciprocalLattice.data()->coordinate);

  reciprocalLattice[0] /= M_PI*2.;
  reciprocalLattice[1] /= M_PI*2.;
  reciprocalLattice[2] /= M_PI*2.;

  real maxabsG(cartesianGrid[0].l);
  for ( int d(1); d < NG ; ++d){
    maxabsG = std::max(maxabsG,cartesianGrid[d].l);
  }
  //  maxabsG = sqrt(2.)*maxabsG;  // real the plane-wave cutoff

  real maxAugmentedAbsG(0.);

  maxAugmentedAbsG = sqrt(augmentationFactor)*maxabsG;

  LOG(0,"MaxAbsG") << maxabsG << " " << maxAugmentedAbsG<< std::endl;


  // Purpose: construct augmented reciprocal lattice with larger G-vectors
  //          than the used lattice.
  // This is high quality Fortran influenced coding
  // enlarge the triples until no further elements are included

  int m(1);
  std::vector<Vector<>> augmentedGrid;
  std::vector<real> augmentedSF;

  int current(1);
  while(current > 0){
    current = 0;
    for ( int i(-m); i<=m; ++i){
      for ( int j(-m); j<=m; ++j){
	for ( int k(-m); k<=m; ++k){
	  if ( ( std::abs(i) < m ) && ( std::abs(j) < m )
	       &&  ( std::abs(k) < m ) ) continue;
	  Vector<> GridPoint(real(i)*reciprocalLattice[0]+
			     real(j)*reciprocalLattice[1]+real(k)*reciprocalLattice[2]);
	  if ( GridPoint.length() < maxAugmentedAbsG ) {
	    current++;
	    if ( GridPoint.length() > maxabsG ) {
	      augmentedGrid.push_back(GridPoint);
	      real locSF(coeff/pow(GridPoint.length(),4));
	      augmentedSF.push_back(locSF);
	    }
	  }
	}
      }
    }
    m++;
  }

  real ExtrapolatedAugmentedEnergy(0.);
  for ( unsigned int d(0); d< augmentedGrid.size(); ++d){
    real vg(constantFactor);
    vg /= augmentedGrid[d].length()*augmentedGrid[d].length();
    ExtrapolatedAugmentedEnergy += vg*augmentedSF[d];
  }

  LOG(0,"ExtrapolatedEnergy:") << " Grid Points: " <<  NG
			       << " Augmented Grid points: " << NG + augmentedGrid.size()
			       << " Energy: " << ExtrapolatedAugmentedEnergy << std::endl;



  if ( itype == 2) {
    int realdNG(extrapolatedGrid.size()+augmentedGrid.size());

    std::vector<real> totalSF;
    std::vector<Vector<>> totalGrid;

    totalSF.reserve(realdNG);
    totalGrid.reserve(realdNG);

    totalSF.insert(totalSF.end(),extrapolatedSF.begin(),extrapolatedSF.end());
    totalSF.insert(totalSF.end(),augmentedSF.begin(),augmentedSF.end());

    totalGrid.insert(totalGrid.end(),extrapolatedGrid.begin(),extrapolatedGrid.end());
    totalGrid.insert(totalGrid.end(),augmentedGrid.begin(),augmentedGrid.end());

    LOG(0,"IO") << "Writing Augmented Extrapolated Grid" << std::endl;

    auto ctfExtrapolatedSF(
      new CTF::Tensor<>(1, std::vector<int>( {realdNG}).data())
    );
    std::vector<int64_t> indices(ctfExtrapolatedSF->wrld->rank == 0 ? realdNG : 0);
    for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; }

    ctfExtrapolatedSF->write(
      indices.size(), indices.data(), totalSF.data()
    );

    allocatedTensorArgument<>("ExtrapolatedSF", ctfExtrapolatedSF);

    auto ctfExtrapolatedGrid(
      new CTF::Tensor<>(2, std::vector<int>( {3,realdNG} ).data())
    );
    std::vector<int64_t> indices3(ctfExtrapolatedGrid->wrld->rank == 0 ? 3*realdNG : 0);
    for (size_t i(0); i < indices3.size(); ++i) { indices3[i] = i; }

    ctfExtrapolatedGrid->write(
      indices3.size(), indices3.data(), totalGrid.data()->coordinate
    );
    allocatedTensorArgument<>("ExtrapolatedGrid", ctfExtrapolatedGrid);
  }

  */
}

real BasisSetExtrapolation::leastSquareFit(
  std::vector<real> fitabsG, std::vector<real> fitSF){
  // simple fit of data points with the curve f(x) = a/x**4
  real num(0.); real denum(0.);
  for (unsigned int d(0); d < fitabsG.size(); ++d){
    num += fitSF[d]/pow(fitabsG[d],4);
    denum += pow(fitabsG[d],-8);
  }
  return num/denum;
}

real BasisSetExtrapolation::simplestWindow(
  real Gmin, real Gmax, real G)
{
  //  real x(G-Gmin/(Gmax-Gmin));
  //return 6.*x**5-15.*x**4+10.*x**3;
  real output(1.);
  if ( G > Gmin ) output = (Gmax-G)/(Gmax-Gmin);
  if ( G > Gmax ) output = 0.;
  return output;
}

real BasisSetExtrapolation::integrateSimplestWindow(
  real Gmin, real Gmax)
{
  return log(Gmax/Gmin)/(Gmax-Gmin);
}

void BasisSetExtrapolation::basisSetCompletenessFTODDUMPIA(){

  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);

  PTR(Tensor<complex>) GammaGai;
  GammaGai = NEW(Tensor<complex>,getTensorArgument<complex>("ParticleHoleCoulombVertex"));
  
  int NF(GammaGai->lens[0]);
  
  int Nocc[] = {NF, No};



  
  Tensor<> *realInfVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realVG(new Tensor<>(false, *realInfVG));
  //Define take out inf funciton                                                                                 
  class TakeOutInf {
  public:
    real operator ()(real x){
      return std::isinf(x) ? 0.0 : x;
    }
  };
  //Take out the inf from realVG.                                                                                
  TakeOutInf takeOutInf;
  Univar_Function<> fTakeOutInf(takeOutInf);
  realVG->sum(1.0, *realInfVG, "G", 0.0, "G", fTakeOutInf);
  realVG->set_name("realVG");
  Tensor<complex> VG(
    1, realVG->lens, realVG->sym, *realVG->wrld, "VG"
  );
  toComplexTensor(*realVG, VG);
  Tensor<> realInvSqrtVG(false, *realVG);
  Tensor<complex> invSqrtVG(
    1, realInvSqrtVG.lens, realInvSqrtVG.sym,
    *realInvSqrtVG.wrld, "invSqrtVG"
  );

  //Starting a new space whose memory will be erased after operation                                             
  //Define operation inverse square root                                                                         
  class InvSqrt {
  public:
    real operator ()(real x){
    return std::sqrt(1.0 / x);
    }
  };

  //Get the inverted square root of VG                                                                           
  InvSqrt invSqrt;
  Univar_Function<> fInvSqrt(invSqrt);
  realInvSqrtVG.sum(1.0, *realInfVG, "G", 0.0, "G", fInvSqrt);
  toComplexTensor(realInvSqrtVG, invSqrtVG);

  (*GammaGai)["Gai"] *= invSqrtVG["G"];

  Tensor<complex> conjCGai(false, *GammaGai);
  Univar_Function<complex> fConj(conj<complex>);
  conjCGai.sum(1.0, *GammaGai, "Gai", 0.0, "Gai", fConj);

  auto sumUnoccupied(new Tensor<complex>(2,Nocc));
  (*sumUnoccupied)["Gi"] += (*GammaGai)["Gai"]*conjCGai["Gai"];
  auto realSumUnoccupied(new Tensor<real>(2,Nocc));
  fromComplexTensor(*sumUnoccupied,*realSumUnoccupied);

  allocatedTensorArgument<>("SumUnoccupied", realSumUnoccupied);

  int ctfNG[] = { NF };
  
  auto unoccupiedCompletenessNumber(new CTF::Vector<real>(NF));
  
  (*unoccupiedCompletenessNumber)["G"] += (*realSumUnoccupied)["Gi"];
  // TODO: WATCH OUT..:ONLY RESTRICTED!!!!!!                   
  real invOcc(1./static_cast<real>(No));

  (*unoccupiedCompletenessNumber)["G"] *= invOcc;

  allocatedTensorArgument<>("UnoccupiedCompletenessNumber", unoccupiedCompletenessNumber);

  
}

void BasisSetExtrapolation::basisSetCompletenessFullFTODDUMP(){

  // Read the Particle/Hole Eigenenergies + CoulombVertex
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));

  // Get the Particle Hole Coulomb Vertex
  No = epsi->lens[0];
  Nv = epsa->lens[0];
  Np = GammaFqr->lens[1];
  NF = GammaFqr->lens[0];

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int FaiStart[] = {0, aStart,iStart};
  int FaiEnd[]   = {NF,aEnd,  iEnd};
  int FijStart[] = {0, iStart, iStart};
  int FijEnd[]   = {NF,iEnd, iEnd};

  PTR(Tensor<complex>) GammaGij;
  // Definition of ParticleHole Coulomb Vertex
  PTR(Tensor<complex>) GammaGai;
  {
    Tensor<complex> GammaFai(GammaFqr->slice(FaiStart, FaiEnd));
    Tensor<complex> GammaFij(GammaFqr->slice(FijStart, FijEnd));

    if (isArgumentGiven("CoulombVertexSingularVectors")) {
      Tensor<complex> *UGF(
	getTensorArgument<complex>("CoulombVertexSingularVectors")
      );
      int lens[]= {UGF->lens[0], Nv, No};
      GammaGai = NEW(Tensor<complex>,
        3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
      );
      (*GammaGai)["Gai"] = GammaFai["Fai"] * (*UGF)["GF"];

      int lensocc[] = {UGF->lens[0], No, No};
      GammaGij = NEW(Tensor<complex>,
        3, lensocc, GammaFqr->sym, *GammaFqr-> wrld, "GammaGij"
      );
      (*GammaGij)["Gij"] = GammaFij["Fij"] * (*UGF)["GF"];

    } else {
      int lens[]= {NF, Nv, No};
      GammaGai = NEW(Tensor<complex>,
        3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
      );
      (*GammaGai) = GammaFai;

      int lensocc[]= {NF, No, No};
      GammaGij = NEW(Tensor<complex>,
        3, lensocc, GammaFqr->sym, *GammaFqr->wrld, "GammaGij"
      );
      (*GammaGij) = GammaFij;
    }
  }
  
  Tensor<> *realInfVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realVG(new Tensor<>(false, *realInfVG));
  //Define take out inf funciton
  class TakeOutInf {
  public:
    real operator ()(real x){
      return std::isinf(x) ? 0.0 : x;
    }
  };
  //Take out the inf from realVG.
  TakeOutInf takeOutInf;
  Univar_Function<> fTakeOutInf(takeOutInf);
  realVG->sum(1.0, *realInfVG, "G", 0.0, "G", fTakeOutInf);
  realVG->set_name("realVG");
  Tensor<complex> VG(
    1, realVG->lens, realVG->sym, *realVG->wrld, "VG"
  );
  toComplexTensor(*realVG, VG);
  Tensor<> realInvSqrtVG(false, *realVG);
  Tensor<complex> invSqrtVG(
    1, realInvSqrtVG.lens, realInvSqrtVG.sym,
    *realInvSqrtVG.wrld, "invSqrtVG"
  );

  //Starting a new space whose memory will be erased after operation
  //Define operation inverse square root
  class InvSqrt {
  public:
    real operator ()(real x){
      return std::sqrt(1.0 / x);
    }
  };

  //Get the inverted square root of VG
  InvSqrt invSqrt;
  Univar_Function<> fInvSqrt(invSqrt);
  realInvSqrtVG.sum(1.0, *realInfVG, "G", 0.0, "G", fInvSqrt);
  toComplexTensor(realInvSqrtVG, invSqrtVG);

  //Define CGai
  Tensor<complex> CGai(*GammaGai);
  CGai["Gai"] *= invSqrtVG["G"];

  //Conjugate of CGai
  Tensor<complex> conjCGai(false, CGai);
  Univar_Function<complex> fConj(conj<complex>);
  conjCGai.sum(1.0, CGai, "Gai", 0.0, "Gai", fConj);

  Tensor<complex> CGij(*GammaGij);
  CGij["Gij"] *= invSqrtVG["G"];



  Tensor<complex> conjCGij(false, CGij);
  conjCGij.sum(1.0, CGij, "Gij", 0.0, "Gij", fConj);

  int Nocc[] = {NF, No};
  int ctfNF[] = {NF};
  auto spectrumGi(new Tensor<complex>(2,Nocc));
  auto sumOccupied(new CTF::Vector<complex>(NF));
  auto sumUnoccupied(new CTF::Vector<complex>(NF));
  (*spectrumGi)["Gi"] += CGij["Gii"]*conjCGij["Gii"];
  (*sumOccupied)["G"] += CGij["Gij"]*conjCGij["Gij"];
  (*sumUnoccupied)["G"] += CGai["Gai"]*conjCGai["Gai"];

  auto realSpectrumGi(new Tensor<real>(2,Nocc));
  auto realSumOccupied(new CTF::Vector<real>(NF));
  auto realSumUnoccupied(new CTF::Vector<real>(NF));
                                 
  fromComplexTensor(*spectrumGi,*realSpectrumGi);
  fromComplexTensor(*sumOccupied,*realSumOccupied);
  fromComplexTensor(*sumUnoccupied,*realSumUnoccupied);


  auto basisSetCompletenessNumber(new CTF::Vector<real>(NF));
  
  
  (*basisSetCompletenessNumber)["G"] += (*realSumOccupied)["G"] + (*realSumUnoccupied)["G"];
  // TODO: WATCH OUT..:ONLY RESTRICTED!!!!!!
  real invOcc(1./static_cast<real>(No));
  
  (*basisSetCompletenessNumber)["G"] *= invOcc;
  (*realSumOccupied)["G"] *= invOcc;
  (*realSumUnoccupied)["G"] *= invOcc;
  allocatedTensorArgument<>("SpectrumGi", realSpectrumGi);
  allocatedTensorArgument<>("SumOccupied", realSumOccupied);
  allocatedTensorArgument<>("SumUnoccupied", realSumUnoccupied);
  allocatedTensorArgument<>("BasisSetCompletenessNumber", basisSetCompletenessNumber);
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
    //    if ( realSpaceMesh[i].length() < 3.){
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
  //std::vector<real> pairCorrelationFunction(boxSize);
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
    //    if ( realSpaceMesh[i].length() < 3.){
    for ( int i(0); i < NR ; ++i) {
      real phase(realSpaceMesh[i].dot(reciprocalGrid[j]));
      localRealSpaceBasisSetCompleteness[i] += cos(phase)*basisSetCompletenessNumber[j];
    }
  }
  //std::vector<real> pairCorrelationFunction(boxSize);

  realSpaceBasisSetCompleteness.resize(NR);
  
  //  real normConvoluter(0.);
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
void BasisSetExtrapolation::QGGSliced(int iStart, int iEnd){

  PTR(Tensor<complex>) GammaGai;
  GammaGai = NEW(Tensor<complex>,getTensorArgument<complex>("ParticleHoleCoulombVertex"));

  int NF(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);
  int Np(No+Nv);
  Tensor<> *ctfCoulombKernel(getTensorArgument<>("CoulombKernel"));

//  LOG(0,"dims:") << NF << " " << Nv << " " << No << std::endl;

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
  
//  LOG(0,"dims2:") << CGaiSliced->lens[0] << " " << CGaiSliced->lens[1] << " " << CGaiSliced->lens[2] << std::endl;

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

  int NGR((NF+1)/2);
  int halfs[] = {0, 0};
  int halfe[] = {NGR, NGR};
  auto QGGhalf(new Tensor<complex>(QGG->slice(halfs,halfe)));
  auto absQGGhalf(new Tensor<>(absQGG->slice(halfs,halfe)));
  auto realQGGhalf(new Tensor<>(realQGG->slice(halfs,halfe)));

  allocatedTensorArgument<complex>("QGG",QGGhalf);
  allocatedTensorArgument<>("AbsQGG", absQGGhalf);
  allocatedTensorArgument<>("RealQGG", realQGGhalf);



}


void BasisSetExtrapolation::QGG(){

  PTR(Tensor<complex>) GammaGai;
  GammaGai = NEW(Tensor<complex>,getTensorArgument<complex>("ParticleHoleCoulombVertex"));

  int NF(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]); 
  int Np(No+Nv);
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
  Tensor<complex> conjCGai(false, *GammaGai);
  Univar_Function<complex> fConj(conj<complex>);
  
  conjCGai.sum(1.0, *GammaGai, "Gai", 0.0, "Gai", fConj);
  

  int NGG[] = {NF, NF};
  int NFG[] = {NF, NF, No, No};
  int NFGi[] = {NF, NF, No};
  auto QGG( new CTF::Tensor<complex>(2,NGG));
//ij  PTR(Tensor<complex>) FGijd;
  PTR(Tensor<complex>) FGij;
  PTR(Tensor<complex>) FGji;
  PTR(Tensor<complex>) FGone;
  PTR(Tensor<complex>) FGtwo;
//ij  FGijd = NEW(Tensor<complex>,4,NFG);
  FGij = NEW(Tensor<complex>,4,NFG);
  FGji = NEW(Tensor<complex>,4,NFG);
  FGone = NEW(Tensor<complex>,2,NGG);
  FGtwo = NEW(Tensor<complex>,2,NGG);

//  PTR(Tensor<complex>) FGi;
//  FGi = NEW(Tensor<complex>,3,NFGi); 

// aI
  
//ij    (*FGi)["GFi"] = (*GammaGai)["Gai"] * conjCGai["Fai"];
//ij    (*FGijd)["GFij"] = (2.0) * (*FGi)["GFi"] * conjCGai["Gbj"] * (*GammaGai)["Fbj"];
    
//ij    (*FGij)["GFij"] = conjCGai["Gbj"] * conjCGai["Fbi"]; 
//ij    (*FGji)["GFij"] = (*GammaGai)["Gai"] * (*GammaGai)["Faj"];
//ij    (*FGijd)["GFij"] += (-1.0) * (*FGij)["GFij"] * (*FGji)["GFij"];

//aI  (*FGone)["GF"] = (*GammaGai)["Gai"] * conjCGai["Fai"];
//aI  (*FGtwo)["GF"] = conjCGai["Gbj"] * (*GammaGai)["Fbj"];
//aI  (*QGG)["GF"] = (2.0) * (*FGone)["GF"] * (*FGtwo)["GF"];
 
//aI  (*FGij)["GFij"] = (*GammaGai)["Gai"] * (*GammaGai)["Faj"];
//aI  (*FGji)["GFij"] = conjCGai["Gbj"] * conjCGai["Fbi"];
//aI  (*QGG)["GF"] += (-1.0) * (*FGij)["GFij"] * (*FGji)["GFij"];


// aG


  (*FGone)["GF"] =  conjCGai["Gai"] * (*GammaGai)["Fai"];
  (*FGtwo)["GF"] =  conjCGai["Fbj"] * (*GammaGai)["Gbj"];
  (*QGG)["GF"]   = (2.0) * (*FGone)["GF"] * (*FGtwo)["GF"];

  
  (*FGij)["GFij"] = conjCGai["Gai"] * (*GammaGai)["Faj"];
  (*FGji)["GFij"] = conjCGai["Fbi"] * (*GammaGai)["Gbj"];
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


  int NGR((NF+1)/2);
  int halfs[] = {0, 0};
  int halfe[] = {NGR, NGR};
  auto QGGhalf(new Tensor<complex>(QGG->slice(halfs,halfe)));
  auto absQGGhalf(new Tensor<>(absQGG->slice(halfs,halfe)));
  auto realQGGhalf(new Tensor<>(realQGG->slice(halfs,halfe)));

  allocatedTensorArgument<complex>("QGG",QGGhalf);
  allocatedTensorArgument<>("AbsQGG", absQGGhalf);
  allocatedTensorArgument<>("RealQGG", realQGGhalf);

}

void BasisSetExtrapolation::calculateNewSF(
  int type, real gamma, Tensor<> *coulombKernel, Tensor<> *newSF, Tensor<> *resNewSF){

  CTF::Tensor<complex> *QGG(getTensorArgument<complex>("QGG"));
  int NG(QGG->lens[0]);
  int NFF[] = {NG};
// Multiply QGG' f12_G' for an approximative structure factor
  auto cK(new Tensor<>(1,NFF));
  if ( NG != coulombKernel->lens[0]){
   int NGR(NG*2-1);
   if ( NGR  != coulombKernel->lens[0]) throw new EXCEPTION("dimension missmatch");
   int gstart[] = {0};
   int gend[]  = {NG};
   (*cK) = coulombKernel->slice(gstart,gend);
  } 
  else{
   (*cK) = (*coulombKernel);
  }
 
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

//  allocatedTensorArgument<complex>("cK",complexNewSF);

 
  if (NG == coulombKernel->lens[0]){
    fromComplexTensor(*complexNewSF,*newSF);
    fromComplexTensor(*complexResNewSF,*resNewSF);
  }
  else if(NG*2-1 == coulombKernel->lens[0]){
    auto realNewSF(new Tensor<>(1,NFF));
    auto realResNewSF(new Tensor<>(1,NFF));
    fromComplexTensor(*complexNewSF,*realNewSF);
    fromComplexTensor(*complexResNewSF,*realResNewSF);
    LOG(0,"JAAAA") << std::endl; 
    allocatedTensorArgument<>("cK",realNewSF);
    //write to vector...double it...write to tensor...
    //or ask felix
    std::vector<real> vectorSF;
    std::vector<real> vectorResSF;
    std::vector<int64_t> indices(realNewSF->wrld->rank == 0 ? NG : 0);
    for (size_t i(0); i < indices.size(); ++i) { indices[i] = i;}
    vectorSF.resize(NG);
    vectorResSF.resize(NG);
    realNewSF->write(
      indices.size(), indices.data(), vectorSF.data()
    );
    realResNewSF->write(
      indices.size(), indices.data(), vectorResSF.data()
    );

    for (size_t i(1); i< indices.size(); ++i){
//      LOG(0,"sf") << vectorSF[i] << std::endl;
      vectorSF.push_back(vectorSF[i]);
      vectorResSF.push_back(vectorResSF[i]);
    }     
    newSF->read_all(vectorSF.data());
    resNewSF->read_all(vectorResSF.data());
    
  }
  else{
    throw new EXCEPTION("dimension missmatch");
  }



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

//  gamma = 0.001;
  for (int i(0); i<=iter; ++i){  
//    gamma = 0.001 + static_cast<real>(i)*0.001;
   // LOG(0,"gamma") << gamma << std::endl;
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
  //  LOG(0,"counter") << counter << std::endl;
    CTF::Scalar<double> denom;
    denom[""] = (*dummy)["G"];
  //  LOG(0,"denom") << denom.get_val() << std::endl;
 
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
//    LOG(0,"numer") << numer.get_val() << std::endl;
    gamma -= numer.get_val()/denom.get_val();

//    LOG(0,"gamma") << gamma << std::endl;
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
//    LOG(0,"-----") << std::endl; 

  }
  allocatedTensorArgument<>("ResNewSF",resNewSF);
  allocatedTensorArgument<>("NewSF",newSF);
}

void BasisSetExtrapolation::dryQGG(){

  DryTensor<> *epsi(getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));

  DryTensor<> *epsa(getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));
  
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);

  DryTensor<complex> *GammaGai(getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex"));
//  GammaGai = NEW(DryTensor<complex>,getTensorArgument<complex, DryTensor<double>>("ParticleHoleCoulombVertex"));

  int NF(GammaGai->lens[0]);

 // DryTensor<> *realInfVG(getTensorArgument<double, DryTensor<double>>("CoulombKernel"));


  int NGG[] = {NF, NF, No};
  int NFG[] = {NF, NF, No, No};
  int syms1[] = {NS, NS, NS};
  int syms2[] = {NS, NS, NS, NS}; 
//  PTR(DryTensor<complex>) QGG;
//  DryTensor<> FGbij(5,NFG,syms2);
 

  DryTensor<complex> FGijd(4,NFG,syms2);
  DryTensor<complex> FGij(4,NFG,syms2);
  DryTensor<complex> FGji(4,NFG,syms2);

  DryTensor<complex> FGi(3,NGG,syms1);


//  QGG = NEW(DryTensor<complex>,2,NGG);

//  PTR(DryTensor<complex>) FGij;
//  PTR(DryTensor<complex>) FGji;
//  FGij = NEW(DryTensor<complex>,4,NFG);
//  FGji = NEW(DryTensor<complex>,4,NFG);


}
//void BasisSetExtrapolation::convolutionInR(){
//}
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
