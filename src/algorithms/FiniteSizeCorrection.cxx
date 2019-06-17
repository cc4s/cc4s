#include <algorithms/FiniteSizeCorrection.hpp>

#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Vector.hpp>
#include <math/Interpolation.hpp>
#include <gte/TricubicInterpolation.hpp>
#include <gte/TrilinearInterpolation.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/SharedPointer.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <iostream>
// FIXME: use common way for math constants
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <sstream>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(FiniteSizeCorrection);

FiniteSizeCorrection::FiniteSizeCorrection(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

FiniteSizeCorrection::~FiniteSizeCorrection() {
}

void FiniteSizeCorrection::run() {
  int fReadFromFile(getIntegerArgument("fReadFromFile", 0));
  if (fReadFromFile == 1) {
    LOG(0,"FiniteSize") << "Reading structure factor from file" << std::endl;
    readFromFile();
  } else {
    Data *Tabij(getArgumentData("DoublesAmplitudes"));
    TensorData<double> *realTabij(dynamic_cast<TensorData<double>*>(Tabij));
    if (realTabij) {
      LOG(0,"FiniteSize") << "Calculating real structure factor" << std::endl;
      calculateRealStructureFactor();
    } else {
      LOG(0,"FiniteSize") << "Calculating complex structure factor" << std::endl;
      calculateComplexStructureFactor();
    }
  }

  int fFiniteSizeCorrection(getIntegerArgument("FiniteSize",0));
  if (fFiniteSizeCorrection == 0){
    LOG(0,"FiniteSize") << "Interpolating and integrating" << std::endl;
    interpolation3D();
    LOG(0,"FiniteSize") << "Caclulating finite size correction" << std::endl;
    calculateFiniteSizeCorrection();
  }
}

void FiniteSizeCorrection::dryRun() {
  int fReadFromFile(getIntegerArgument("fReadFromFile", 0));
  if (fReadFromFile == 1) {
    LOG(0,"FiniteSize") << "Reading structure factor from file" << std::endl;
  } else {
    LOG(0,"FiniteSize") << "Calculating structure factor" << std::endl;
    dryCalculateStructureFactor();
  }
  //constructFibonacciGrid();
  LOG(0,"FiniteSize") << "Interpolating and integrating" << std::endl;
  dryInterpolation3D();
  LOG(0,"FiniteSize") << "Caclulating finite size correction" << std::endl;
  dryCalculateFiniteSizeCorrection();
}


class FiniteSizeCorrection::Momentum {
  public:
  cc4s::Vector<> v;
  double s;
  double l;
  double vg;
  Momentum(): s(0.0), l(0.0), vg(0.) {
  }
  Momentum(
    cc4s::Vector<> v_, double s_=0., double vg_=0.
  ): v(v_), s(s_), l(v_.length()), vg(vg_) {
  }
  double locate(Momentum *m, int const n) {
    cc4s::Vector<> u(v);
    //if (v[3] < 0.) u= v*(-1.);
    for (int d(0); d < n; ++d) {
      if (u.approximately(m[d].v)) {
	return m[d].s;
      }
    }
    return 0;
  }

  static bool sortByLength (Momentum const &n, Momentum const &m) {
    return n.l < m.l;
  }
  static bool sortByVector (Momentum const &n, Momentum const &m) {
    return n.v < m.v;
  }
};

void FiniteSizeCorrection::readFromFile(){
  LOG(0,"readFromFile") << "reading " << std::endl;
  Tensor<> *realVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realSG(getTensorArgument<>("StructureFactor"));
  LOG(1,"readFromFile") << "success\n Loading into Vectors " << std::endl;
  NG=realVG->lens[0];
  VofG.resize(NG);
  realVG->read_all(VofG.data());
  LOG(1,"readFromFile") << "VofG Finished" << std::endl;
  structureFactors.resize(NG);
  realSG->read_all(structureFactors.data());
  LOG(1,"readFromFile") << "Finished" << std::endl;
}

void FiniteSizeCorrection::calculateRealStructureFactor() {
  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);
  
  int orbitalPairStart(getIntegerArgument("orbitalPairStart",-1));
  int orbitalPairEnd(getIntegerArgument("orbitalPairEnd",-1));


  if (orbitalPairStart < 0 || orbitalPairStart >= orbitalPairEnd){
    orbitalPairStart = 0;
  }
  if (orbitalPairEnd > No || orbitalPairEnd <= orbitalPairStart){
    orbitalPairEnd = No;
  }
  
  int numberOrbitalPairs(orbitalPairEnd - orbitalPairStart);
  bool orbitalPairs(numberOrbitalPairs != No);
  if ( orbitalPairs ){
    LOG(0,"Orbital Pair Analysis") << "Treating only pairs " << orbitalPairStart
                                  << " to " << orbitalPairEnd << std::endl;
  }
  
  PTR(Tensor<complex>) GammaGai;
  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  if (orbitalPairs){
    iStart = orbitalPairStart;
    iEnd = orbitalPairEnd;
  }

  // Read the Coulomb vertex GammaGqr
  if ( isArgumentGiven("CoulombVertex")){
    Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));
    // Get the Particle Hole Coulomb Vertex
    int NF(GammaFqr->lens[0]);
    
    int FaiStart[] = {0, aStart,iStart};
    int FaiEnd[]   = {NF,aEnd,  iEnd};

    Tensor<complex> GammaFai(GammaFqr->slice(FaiStart, FaiEnd));

    if (isArgumentGiven("CoulombVertexSingularVectors")) {
      Tensor<complex> *UGF(
        getTensorArgument<complex>("CoulombVertexSingularVectors")
      );
      int lens[]= {UGF->lens[0], Nv, No};
      GammaGai = NEW(Tensor<complex>,
        3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
      );
      (*GammaGai)["Gai"] = GammaFai["Fai"] * (*UGF)["GF"];
    }
    else {
      int lens[]= {NF, Nv, No};
      if (orbitalPairs){
        lens[0] = NF; lens[1] = Nv; lens[2] = numberOrbitalPairs;
      }
      GammaGai = NEW(Tensor<complex>,
        3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
      );
      (*GammaGai) = GammaFai;
    }
  }
  else if (isArgumentGiven("ParticleHoleCoulombVertex")){
    if (orbitalPairs){
      Tensor<complex> *GGai(getTensorArgument<complex>("ParticleHoleCoulombVertex"));
      int NF(GGai->lens[0]);
      int GGaiStart[] = {0, aStart,iStart};
      int GGaiEnd[]   = {NF,aEnd, iEnd };
      int sGGaiStart[] = {0, aStart, 0};
      int sGGaiEnd[] = {NF,aEnd,numberOrbitalPairs};
      int lens[] = {NF,Nv,numberOrbitalPairs};
      GammaGai = NEW(Tensor<complex>,3,lens, GGai->sym, *GGai->wrld, "GammaGai");
      GammaGai->slice(sGGaiStart,sGGaiEnd,0.0,*GGai,GGaiStart,GGaiEnd,1.0);
    }
    else{
      GammaGai = NEW(Tensor<complex>,getTensorArgument<complex>("ParticleHoleCoulombVertex"));
    }
  }
  else {
    throw new EXCEPTION("Need Appropriate Coulomb Vertex");    
  }

  Tensor<> *realInfVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realVG(new Tensor<>(false, *realInfVG));
  //Define take out inf funciton
  class TakeOutInf {
  public:
    double operator ()(double x){
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
    double operator ()(double x){
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

  //Split CGai and conjCGai into real and imag parts
  Tensor<> realCGai(3, GammaGai->lens, GammaGai->sym,
                        *GammaGai->wrld, "RealCGai");
  Tensor<> imagCGai(3, GammaGai->lens, GammaGai->sym,
                        *GammaGai->wrld, "ImagCGai");
  fromComplexTensor(CGai, realCGai, imagCGai);

  Tensor<> realConjCGai(3, GammaGai->lens, GammaGai->sym,
                        *GammaGai->wrld, "RealConjCGai");
  Tensor<> imagConjCGai(3, GammaGai->lens, GammaGai->sym,
                        *GammaGai->wrld, "ImagConjCGai");
  fromComplexTensor(conjCGai, realConjCGai, imagConjCGai);

  Tensor<> *realTabij;
  if (orbitalPairs){
    Tensor<> *Tabij(getTensorArgument("DoublesAmplitudes"));
    if (isArgumentGiven("SinglesAmplitudes") ) {
      //Get Tai
      Tensor<> *realTai(getTensorArgument("SinglesAmplitudes"));
      (*Tabij)["abij"] += (*realTai)["ai"] * (*realTai)["bj"];
    }
    int aiSliceStart[] = { 0,  0, orbitalPairStart, orbitalPairStart};
    int aiSliceEnd[]   = {Nv, Nv, orbitalPairEnd,   orbitalPairEnd  };
    int lens[] = {Nv, Nv, numberOrbitalPairs, numberOrbitalPairs};
    realTabij = new Tensor<>(3, lens, Tabij->sym, *Tabij->wrld, "realTabij");
    (*realTabij) = Tabij->slice(aiSliceStart,aiSliceEnd); 
  }
  else{
    //Get Tabij
    realTabij = getTensorArgument("DoublesAmplitudes");

    if (isArgumentGiven("SinglesAmplitudes") ) {
      //Get Tai
      Tensor<> *realTai(getTensorArgument("SinglesAmplitudes"));
      (*realTabij)["abij"] += (*realTai)["ai"] * (*realTai)["bj"];
    }
  }
  
 //Construct SG
  NG = GammaGai->lens[0];
  CTF::Vector<> *realSG(new CTF::Vector<>(NG, *GammaGai->wrld, "realSG"));
  CTF::Vector<> *realSGd(new CTF::Vector<>(NG, *GammaGai->wrld, "realSGd"));
  CTF::Vector<> *realSGx(new CTF::Vector<>(NG, *GammaGai->wrld, "realSGx"));
  (*realSGd)["G"]  = ( 1.0) * realConjCGai["Gai"] * realCGai["Gbj"] * (*realTabij)["abij"];
  (*realSGd)["G"] += (-1.0) * imagConjCGai["Gai"] * imagCGai["Gbj"] * (*realTabij)["abij"];

  (*realSGx)["G"] += ( 1.0) * realConjCGai["Gaj"] * realCGai["Gbi"] * (*realTabij)["abij"];
  (*realSGx)["G"] += (-1.0) * imagConjCGai["Gaj"] * imagCGai["Gbi"] * (*realTabij)["abij"];

  (*realSG)["G"] = (2.0) * (*realSGd)["G"] + (-1.0) * (*realSGx)["G"];
  allocatedTensorArgument<>("StructureFactor", realSG);
  if(isArgumentGiven("StructureFactors")){
    CTF::Vector<> *realSGs(new CTF::Vector<>(NG, *GammaGai->wrld, "realSGs"));
    (*realSGs)["G"] = (0.5) * (*realSGd)["G"] + (0.5) * (*realSGx)["G"];
    allocatedTensorArgument<>("StructureFactors",realSGs);
    Scalar <> senergy(*Cc4s::world);
    senergy[""] = (*realVG)["G"] * (*realSGs)["G"];
    double _senergy(senergy.get_val());
    LOG(0,"Singlet energy:") << _senergy << std::endl;

  }

  if(isArgumentGiven("StructureFactort")){
    CTF::Vector<> *realSGt(new CTF::Vector<>(NG, *GammaGai->wrld, "realSGt"));
    (*realSGt)["G"] = (1.5) * (*realSGd)["G"] + (-1.5) * (*realSGx)["G"];
    allocatedTensorArgument<>("StructureFactort",realSGt);
    Scalar <> tenergy(*Cc4s::world);
    tenergy[""] = (*realVG)["G"] * (*realSGt)["G"];
    double _tenergy(tenergy.get_val());
    LOG(0,"Triplet energy:") << _tenergy << std::endl;
  }

  VofG.resize(NG);
  realVG->read_all(VofG.data());
  structureFactors.resize(NG);
  realSG->read_all(structureFactors.data());
}


void FiniteSizeCorrection::calculateComplexStructureFactor() {
  Tensor<> *realInfVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realVG(new Tensor<>(false, *realInfVG));
  // Define take out inf funciton
  class TakeOutInf {
  public:
    double operator ()(double x){
      return std::isinf(x) ? 0.0 : x;
    }
  };
  // Take out the inf from realVG.
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
    double operator ()(double x){
      return std::sqrt(1.0 / x);
    }
  };

  //Get the inverted square root of VG
  InvSqrt invSqrt;
  Univar_Function<> fInvSqrt(invSqrt);
  realInvSqrtVG.sum(1.0, *realInfVG, "G", 0.0, "G", fInvSqrt);
  toComplexTensor(realInvSqrtVG, invSqrtVG);

  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaFqr( getTensorArgument<complex>("CoulombVertex"));
  Tensor<complex> *GammaGqr;

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  if (isArgumentGiven("CoulombVertexSingularVectors")) {
    Tensor<complex> *UGF(
      getTensorArgument<complex>("CoulombVertexSingularVectors")
    );
    int lens[]= {UGF->lens[0], GammaFqr->lens[1], GammaFqr->lens[2]};
    GammaGqr = new Tensor<complex>(
      3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
    );
    (*GammaGqr)["Gqr"] = (*GammaFqr)["Fqr"] * (*UGF)["GF"];
  } else {
    int lens[]= {GammaFqr->lens[0], GammaFqr->lens[1], GammaFqr->lens[2]};
    GammaGqr = new Tensor<complex>(
      3, lens, GammaFqr->sym, *GammaFqr->wrld, "GammaGqr"
    );
    (*GammaGqr) = (*GammaFqr);
  }

  // Compute the No,Nv,NG,Np
  NG=GammaGqr->lens[0];
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(GammaGqr->lens[1]);

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int GiaStart[] = {0, iStart,aStart};
  int GiaEnd[]   = {NG,iEnd,  aEnd};
  int GaiStart[] = {0, aStart,iStart};
  int GaiEnd[]   = {NG,aEnd,  iEnd};
  //  GammaGia = new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd));
  //  GammaGai = new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd));
  Tensor<complex> GammaGia(GammaGqr->slice(GiaStart, GiaEnd));
  Tensor<complex> GammaGai(GammaGqr->slice(GaiStart, GaiEnd));

  delete GammaGqr;

  //Define CGia
  Tensor<complex> CGia(GammaGia);
  CGia["Gia"] *= invSqrtVG["G"];

  Tensor<complex> conjTransposeCGia(false, GammaGia);
  Univar_Function<complex> fConj(conj<complex>);
  conjTransposeCGia.sum(1.0,GammaGai,"Gai", 0.0,"Gia", fConj);
  conjTransposeCGia["Gia"] *= invSqrtVG["G"];

  Tensor<complex> conjTransposeGammaGia(false, GammaGia);
  conjTransposeGammaGia.sum(1.0,GammaGai,"Gai", 0.0,"Gia", fConj);

  /*
  Tensor<complex> conjCGai(false, GammaGai);
  Univar_Function<complex> fConj(conj<complex>);
  conjCGai.sum(1.0,GammaGai,"Gai", 0.0,"Gai", fConj);
  conjCGai["Gai"] *= invSqrtVG["G"];
  */

  //Get Tabij
  Tensor<complex> *Tabij(getTensorArgument<complex>("DoublesAmplitudes"));

  if (isArgumentGiven("SinglesAmplitudes") ) {
    //Get Tai
    Tensor<complex> *Tai(getTensorArgument<complex>("SinglesAmplitudes"));
    (*Tabij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  }

  //construct SG

  CTF::Vector<> *realSG(new CTF::Vector<>(NG, *CGia.wrld, "realSG"));
  CTF::Vector<complex> *SG(new CTF::Vector<complex>(NG, *CGia.wrld, "SG"));
  (*SG)["G"]  = ( 2.0) * conjTransposeCGia["Gia"] * CGia["Gjb"] * (*Tabij)["abij"];
  (*SG)["G"] += (-1.0) * conjTransposeCGia["Gja"] * CGia["Gib"] * (*Tabij)["abij"];
  fromComplexTensor(*SG, *realSG);
  allocatedTensorArgument<>("StructureFactor", realSG);

  VofG.resize(NG);
  realVG->read_all(VofG.data());
  structureFactors.resize(NG);
  realSG->read_all(structureFactors.data());
}


void FiniteSizeCorrection::dryCalculateStructureFactor() {
  DryTensor<complex> *GammaGai(nullptr);
  //Definition of the variables
  if ( isArgumentGiven("CoulombVertex")){
    DryTensor<complex> *GammaFai(
      getTensorArgument<complex, DryTensor<complex>>("CoulombVertex")
    );
    // don't slice in dry run
    GammaGai = GammaFai;
  }
  else if( isArgumentGiven("ParticleHoleCoulombVertex")){
    DryTensor<complex> *GammaFai(
      getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
    );
    GammaGai = GammaFai;
  }

  int symo[] = { NS, NS };

  int NG=GammaGai->lens[0];
  int len[]={NG,1};

  DryTensor<> realInfVG(2,len,symo);
  DryTensor<> realVG(realInfVG);
  DryTensor<> VG(realVG);
  DryTensor<> realInvSqrtVG(realVG);
  DryTensor<> invSqrtVG(realInfVG);

  //Define CGai
  DryTensor<complex> CGai(*GammaGai);

  //Conjugate of CGai
  DryTensor<complex> conjCGai(CGai);


  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms4[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms4);

  DryTensor<> rTabij(Tabij);
  DryTensor<> iTabij(Tabij);
  int syms[] = {NS, NS};
  allocatedTensorArgument("StructureFactor", new DryTensor<>(2, len, syms, SOURCE_LOCATION));
}


void FiniteSizeCorrection::constructFibonacciGrid(double R, int N) {
  //This function construct a Fibonacci grid on a sphere with a certain radius.
  //Returns a vector of vectors: {x,y,z}
  //The N should be fixed and R should be a vector which is selected by another
  //function which determines the R's
  //N = 128; N is the number of points on the sphere, defined in .cxx file
  double inc = M_PI * (3 - std::sqrt(5));
  fibonacciGrid = new Momentum[N];

  for (int k(0); k < N; ++k) {
    double z((2.0*k+1)/N - 1.0);
    double r(R * std::sqrt(1.0 - z*z));
    double phi(k * inc);
    fibonacciGrid[k].v[0] = r * std::cos(phi);
    fibonacciGrid[k].v[1] = r * std::sin(phi);
    fibonacciGrid[k].v[2] = R * z;
  }
}

void FiniteSizeCorrection::interpolation3D() {
  Tensor<> *momenta(getTensorArgument<>("Momenta"));
  //  int NG(momenta->lens[1]);
  cc4s::Vector<> *cartesianMomenta(new cc4s::Vector<>[NG]);
  momenta->read_all(&cartesianMomenta[0][0]);

  // FIXME: give direct or reciprocal grid and calaculate all properties,
  // such as a,b,c and omega from it.

  // Detecting whether a grid is half or full by summming up all
  // G vectors. If it is full and symmetric, the length of the sum
  // should be equal to 0. But in case where the grid is not symmetric,
  // this method will fail.
  cc4s::Vector<> check_grid;
  for (int g(0); g < NG; ++g){
    check_grid+=cartesianMomenta[g];
  }
  LOG(1, "interpolation3D")<< "Length of the sum of all G vectors="
    << check_grid.length() << std::endl;
  if (check_grid.length() > 1e-10){
    LOG(1, "interpolation3D")<<
      "Working with half G grid, completing the other half..."<< std::endl;
    LOG(1, "interpolation3D")<<
      "!!! Always check if you are indeed working with half G grid !!!"
      << std::endl;
    cartesianGrid = new Momentum[(2*NG-1)];
    cartesianGrid[0] = Momentum(
      cartesianMomenta[0], 0.5*structureFactors[0], VofG[0]
    );
    for (int g(1); g<NG; ++g) {
      cartesianGrid[g] = Momentum(
        cartesianMomenta[g], 0.5*structureFactors[g], VofG[g]
      );
     cartesianGrid[(g+NG-1)] = Momentum(
        cartesianMomenta[g]*(-1.), 0.5*structureFactors[g], VofG[g]
      );
    }
    NG = (2*NG-1);
  }
  else {
    LOG(1, "interpolation3D")<< "Working with full G grid." << std::endl;
    LOG(1, "interpolation3D")<<
      "!!! Always check if you are indeed working with full G grid !!!"
      << std::endl;
    cartesianGrid = new Momentum[NG];
    for (int g(0); g<NG; ++g)
      cartesianGrid[g] = Momentum(
        cartesianMomenta[g], structureFactors[g], VofG[g]
      );
  }

  // sort according to vector length.
  std::sort(cartesianGrid, &cartesianGrid[NG], Momentum::sortByLength);

  // get the 3 unit vectors;
  cc4s::Vector<> a(cartesianGrid[1].v);

  // GC is the shortest vector.
  if (isArgumentGiven("shortestGvector")) {
    GC = getRealArgument("shortestGvector");
  }
  else {
    GC = a.length();
  }
  // determine a small length on the scale of the cell
  const double epsilon(1e-9*a.length());

  LOG(2, "GridSearch") << "b1=#1" << std::endl;
  int j=2;
  //a and b should not be parallel;
  while ((a.cross(cartesianGrid[j].v)).length() < epsilon) ++j;
  cc4s::Vector<> b(cartesianGrid[j].v);
  LOG(2, "GridSearch") << "b2=#" << j << std::endl;
  ++j;
  //a, b and c should not be on the same plane;
  while (abs((a.cross(b)).dot(cartesianGrid[j].v)) < epsilon) ++j;
  cc4s::Vector<> c(cartesianGrid[j].v);
  LOG(2, "GridSearch") << "b3=#" << j << std::endl;

  // print the basis vectors
  LOG(2, "GridSearch") << "b1=" << a << std::endl;
  LOG(2, "GridSearch") << "b2=" << b << std::endl;
  LOG(2, "GridSearch") << "b3=" << c << std::endl;

  //construct the transformation matrix
  std::vector<Vector<>> T(3);
  double Omega((a.cross(b)).dot(c));
  T[0] = b.cross(c)/Omega;
  T[1] = c.cross(a)/Omega;
  T[2] = a.cross(b)/Omega;

  std::vector<Vector<>> ReciprocalLattice(3);
  ReciprocalLattice[0] = a*M_PI*2.;
  ReciprocalLattice[1] = b*M_PI*2.;
  ReciprocalLattice[2] = c*M_PI*2.;

  auto ctfReciprocalLattice(
    new CTF::Tensor<>(2, std::vector<int>({3,3}).data())
  );
  auto ctfRealLattice(
    new CTF::Tensor<>(2, std::vector<int>({3,3}).data())
  );

  std::vector<int64_t> indices(ctfReciprocalLattice->wrld->rank == 0 ? 3*3 : 0);
  for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; }
  ctfReciprocalLattice->write(
    indices.size(), indices.data(), ReciprocalLattice.data()->coordinate
  );
  ctfRealLattice->write(
    indices.size(), indices.data(), T.data()->coordinate
  );

  if (isArgumentGiven("ReciprocalLattice")) {
    allocatedTensorArgument<>("ReciprocalLattice", ctfReciprocalLattice);
  }

  if (isArgumentGiven("RealLattice")) {
    allocatedTensorArgument<>("RealLattice", ctfRealLattice);
  }

  // determine bounding box in direct (reciprocal) coordinates
  Vector<> directMin, directMax;
  for (int g(0); g < NG; ++g) {
    for (int d(0); d < 3; ++d) {
      double directComponent(T[d].dot(cartesianGrid[g].v));
      directMin[d] = std::min(directMin[d], directComponent);
      directMax[d] = std::max(directMax[d], directComponent);
    }
  }
  LOG(2, "FiniteSizeInterpolation") << "directMin=" << directMin <<
    ", directMax=" << directMax << std::endl;

  // build grid for the entire bounding box
  Vector<int> boxDimensions, boxOrigin;
  int64_t boxSize(1);
  for (int d(0); d < 3; ++d) {
    boxSize *=
      boxDimensions[d] = std::floor(directMax[d] - directMin[d] + 1.5);
    boxOrigin[d] = std::floor(directMin[d] + 0.5);
  }
  LOG(2, "FiniteSizeInterpolation") << "boxOrigin=" << boxOrigin <<
    " boxDimensions=" << boxDimensions << std::endl;

  // allocate and initialize regular grid
  double *regularSG(new double[boxSize]);
  for (int64_t g(0); g < boxSize; ++g) regularSG[g] = 0;
  // enter known SG values
  for (int g(0); g < NG; ++g) {
    int64_t index(0);
    Vector<> directG;
    for (int d(2); d >= 0; --d) {
      directG[d] = T[d].dot(cartesianGrid[g].v);
      index *= boxDimensions[d];
      index += std::floor(directG[d] + 0.5) - boxOrigin[d];
    }
    if (regularSG[index] != 0.0) {
      LOG(2, "FiniteSizeInterpolation") <<
        "Overwriting previous grid value G_direct=" << directG <<
        ", index=" << index << std::endl;
    }
    regularSG[index] = cartesianGrid[g].s;
  }

  // check number of points in the interior and on the boundary
  int64_t interiorPointsCount(0), boundaryPointsCount(0);
  for (int z(1); z < boxDimensions[2]-1; ++z) {
    for (int y(1); y < boxDimensions[1]-1; ++y) {
      for (int x(1); x < boxDimensions[0]-1; ++x) {
        int64_t index(x + boxDimensions[0] * (y + boxDimensions[1]*z));
        bool inside(true);
        for (int dz(-1); dz <= 1; ++dz) {
          for (int dy(-1); dy <= 1; ++dy) {
            for (int dx(-1); dx <= 1; ++dx) {
              int64_t offset(dx + boxDimensions[0]*(dy + boxDimensions[1]*dz));
              inside &= regularSG[index+offset] != 0.0;
              if (!inside) break;
            }
            if (!inside) break;
          }
          if (!inside) break;
        }
        interiorPointsCount += inside ? 1 : 0;
        boundaryPointsCount += regularSG[index] != 0.0 && !inside ? 1 : 0;
      }
    }
  }
  LOG(2, "FiniteSizeInterpolation") << "Number of momentum points inside cutoff=" <<
    interiorPointsCount << ", Number of momentum points on boundary=" <<
    boundaryPointsCount << std::endl;

  // create trilinear or tricubic interpolator
  // TODO: use factory to select different interpolators, similar to mixers
  //int fTricubic(getIntegerArgument("fTricubic", 1));
  gte::IntpTricubic3<double> interpolatedSG(
  boxDimensions[0], boxDimensions[1], boxDimensions[2],
  boxOrigin[0], 1, boxOrigin[1], 1, boxOrigin[2], 1,
  regularSG,
  true
  );
  /**
    gte::IntpTrilinear3<double> interpolatedSG(
    boxDimensions[0], boxDimensions[1], boxDimensions[2],
    boxOrigin[0], 1, boxOrigin[1], 1, boxOrigin[2], 1,
    regularSG
    );
  **/
  // spherically sample
  int64_t N(
    getIntegerArgument("gridPointsFibonacci",DEFAULT_NUM_FIBONACCI)
  );
  double lastLength(-1);
  averageSGs.clear(); GLengths.clear();
  meanErrorSG.clear();
  double maxlength(cartesianGrid[0].l);
  for (int g(1); g < NG; ++g) {
    maxlength = std::max(maxlength,cartesianGrid[g].l);
  }
  int num=1000;   // 1000 gridpoints for G=0..G_{max}
  for (int g(0); g < num; ++g) {
    double length(maxlength/1000.*double(g));
    if (abs(length - lastLength) > 1e-3) {
      constructFibonacciGrid(length,N);
      double sumSG(0.);
      for (int f(0); f < N; ++f) {
        Vector<> directG;
        for (int d(0); d < 3; ++d) {
          directG[d] = T[d].dot(fibonacciGrid[f].v);
        }
        // lookup interpolated value in direct coordinates
        sumSG += interpolatedSG(directG[0], directG[1], directG[2]);
      }
      sumSG /= N;
      averageSGs.push_back(sumSG);
      GLengths.push_back(length);
      double meanError(0.);
      for (int f(0); f < N; ++f) {
        Vector<> directG;
        for (int d(0); d < 3; ++d) {
          directG[d] = T[d].dot(fibonacciGrid[f].v);
        }
	double interpol(interpolatedSG(directG[0], directG[1], directG[2]));
	meanError += std::abs(interpol - sumSG);
      }
      meanErrorSG.push_back(meanError/N);
      lastLength = length;
    }
  }

  //  for (int g(0); g<num; ++g){
  //  LOG(2,"sphericalAv") << GLengths[g] << " " << averageSGs[g] << " " << meanErrorSG[g] << std::endl;
  //  }

  //Define the 3D zone close to the Gamma point which needed to be
  //integrated over. Find the vectors which define it. Small BZ
  //  for (int t(0); t < 20; t++){
  //    LOG(0,"G vectors by length") << cartesianGrid[t].v << std::endl;
  //  }

  std::vector<Vector<>> smallBZ;
  smallBZ.push_back(cartesianGrid[1].v);
  for (int t(2); t<NG; t++){
    if (IsInSmallBZ(cartesianGrid[t].v, 1., smallBZ)){
      smallBZ.push_back(cartesianGrid[t].v);
    }
  }

  LOG(1,"interpolation3D") << "Number of basis vectors of small BZ="
    << smallBZ.size() << std::endl;
  for (std::vector<int>::size_type i = 0; i != smallBZ.size(); i++){
    LOG(2,"interpolation3D") << "smallBZ basis vector: " << smallBZ[i] << std::endl;
  }

  //integration in 3D
  double constantFactor(getRealArgument("constantFactor"));
  //cutOffRadius is set to a big default value 100 to ensure
  //the integration is over the whole G grid.
  double cutOffRadius(getRealArgument("cutOffRadius", 100));
  int N0(51), N1(51), N2(51);
  inter3D = 0.;
  sum3D   = 0.;
  int countNO(0);
  int countNOg(0);
  std::vector<Vector<>> gridWithinRadius;
  for (int i(0); i < NG; ++i) {
    if (cartesianGrid[i].l < cutOffRadius) {
      gridWithinRadius.push_back(cartesianGrid[i].v);
      if (cartesianGrid[i].l > epsilon) {
        sum3D += constantFactor/cartesianGrid[i].l/cartesianGrid[i].l
                 *cartesianGrid[i].s;
      }
    }
  }

  MpiCommunicator communicator(
    Cc4s::world->rank, Cc4s::world->np, Cc4s::world->comm
  );
  communicator.barrier();
  for (int t0(-N0); t0 <= N0; ++t0) {
    for (int t1(-N1); t1 <= N1; ++t1) {
      for (int t2(-N2); t2 <= N2; ++t2) {
        if ((std::abs(t0+t1+t2)) % Cc4s::world->np == Cc4s::world->rank) {
          Vector<double> directg;
          Vector<double> ga(((a/double(N0))*double(t0)));
          Vector<double> gb(((b/double(N1))*double(t1)));
          Vector<double> gc(((c/double(N2))*double(t2)));
          Vector<double> g(ga+gb+gc);
          //for each g that is within smallBZ, add its contribution
          //and that of all its
          //periodic images that differ only in a reciprocal lattice
          //to inter3D
          if (IsInSmallBZ(g, 2, smallBZ)){
            countNOg++;
            for (std::vector<int>::size_type i = 0; i != gridWithinRadius.size(); ++i) {
              //reset g to the vector that is within the first smallBZ.
              g=ga+gb+gc;
              //add an reciprocal lattice vector to g to get a periodic image of it.
              g += gridWithinRadius[i];
              for (int d(0); d <3; ++d){
                directg[d]=T[d].dot(g);
              }
              countNO++;
	      double interpol;
              if (g.length() > epsilon) {
		interpol = interpolatedSG(directg[0], directg[1], directg[2]);
                inter3D += interpol * constantFactor/g.length()/g.length();
	      }
	    }
          }
        }
      }
    }
  }


  double totalInter3D(0);
  int totalCountNOg(0);
  communicator.allReduce(inter3D, totalInter3D);
  communicator.allReduce(countNOg, totalCountNOg);

  LOG(2,   "integration3D") << "countNOg= " << totalCountNOg << std::endl;
  LOG(2, "interpolation3D") <<    "sum3D= " <<         sum3D << std::endl;
  inter3D=totalInter3D/totalCountNOg;
  LOG(2, "interpolation3D") << "Number of points in summation=" << countNO << std::endl;
}

void FiniteSizeCorrection::dryInterpolation3D() {
}

//scale=1 is used to search for the vectors which
//define smallBZ; scale = 2 is used to tell if a vector is
//within the smallBZ or not. For a vector that is within the smallBZ,
//its projection on any vectors which define smallBZ must be less
//than 1/2.
bool FiniteSizeCorrection::IsInSmallBZ(
  Vector<> point, double scale, std::vector<Vector<>> smallBZ
)
{
  std::vector<int>::size_type countVector(0);
  for (std::vector<int>::size_type i = 0; i != smallBZ.size(); i++){
    // FIXME: use an epsilon instead of 1e-9
    double epsilon(smallBZ[0].length() * 1e-10);
    if (
      abs(abs(smallBZ[i].dot(point))/smallBZ[i].length()/smallBZ[i].length()
       *scale -1.0) < epsilon
       || abs(smallBZ[i].dot(point))/smallBZ[i].length()/smallBZ[i].length()*scale
       -1.0 < -epsilon
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

double FiniteSizeCorrection::integrate(
  cc4s::Inter1D<double> Int1d,
  double start, double end, int steps
){
  double s = 0;
  double h = (end-start)/steps;
  for (int i = 0; i < steps; ++i)
  s += simpson(Int1d, start + h*i, h);
  return h*s;
}

double FiniteSizeCorrection::simpson(
  cc4s::Inter1D<double> Int1d,
  double x, double h
){
  return (SGxVG(Int1d, x) + 4*SGxVG(Int1d, x+h/2.) + SGxVG(Int1d, x+h))/6.;
}

double FiniteSizeCorrection::SGxVG(
  cc4s::Inter1D<double> Int1d, double x
){
  return (x > 0. && x<GC) ? (cos(x/GC*M_PI)+1)*1./2/x/x*Int1d.getValue(x)*x*x : 0.;
}

void FiniteSizeCorrection::calculateFiniteSizeCorrection() {
  cc4s::Inter1D<double> Int1d(
    GLengths.size(), GLengths.data(), averageSGs.data()
  );
  Int1d.cubicSpline(0., 0., "M");
  double x=0.;
  for (int i(1); i<1000; i++){
    LOG(2, "Interpolation") << x << " " << Int1d.getValue(x) << std::endl;
    x = i*0.001;
  }
  for (unsigned int i(0); i < GLengths.size(); ++i){
    LOG(2, "StructureFactor") << GLengths[i] << " " << averageSGs[i] << std::endl;
  }
  
  int kpoints(getIntegerArgument("kpoints",1));
  double volume(getRealArgument("volume"));
  double constantFactor(getRealArgument("constantFactor"));

  double r1 = integrate(Int1d, 0.0, GC, 1000)*constantFactor*volume*kpoints*4*M_PI;
  double  sumSGVG(0.);

  // we assume the first entry of cartesianGrid corresponds to G=0
  for (int d(1); d < NG; ++d) {
    sumSGVG += cartesianGrid[d].vg * cartesianGrid[d].s;
  }

  LOG(0, "FiniteSize") << "Uncorrected e=" << sumSGVG << std::endl;
  LOG(0, "FiniteSize") << "Corrected e=" << sumSGVG+inter3D-sum3D << std::endl;

  {
    std::stringstream stream;
    stream << std::setprecision(10) << "Summation within cutoff radius=" << sum3D << std::endl;
    LOG(1,"FiniteSize") << stream.str();
  }
  {
    std::stringstream stream;
    stream << std::setprecision(10) << "Integral within cutoff radius="<< inter3D << std::endl;
    LOG(1,"FiniteSize") << stream.str();
  }

  {
    std::stringstream stream;
    stream << std::setprecision(10) << "Spherical Averaging Correction=" << r1 << std::endl;
    LOG(1,"FiniteSize") << stream.str();
  }
  {
    std::stringstream stream;
    stream << std::setprecision(10) << "Spherical Averaging Corrected="<< sumSGVG + r1 << std::endl;
    LOG(1,"FiniteSize") << stream.str();
  }

  setRealArgument("CorrectedEnergy"  , sumSGVG+inter3D-sum3D);
}

void FiniteSizeCorrection::dryCalculateFiniteSizeCorrection() {
  getIntegerArgument("kpoints", 1);
  getRealArgument("volume");
  getRealArgument("constantFactor");

  setRealArgument("CorrectedEnergy", 0.0);
}


