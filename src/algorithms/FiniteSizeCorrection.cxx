#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#define _USE_MATH_DEFINES
#include <math.h>

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
  calculateStructureFactor();
  constructFibonacciGrid();
  calculateFiniteSizeCorrection();
}

void FiniteSizeCorrection::calculateStructureFactor() {
//Definition of the variables
  Tensor<complex> *GammaGai(
	getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );

// local allocation
//  int a(7);
// heap allocation (survive after function return)
//  int *b(new int(7));

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

  //Get Tabij
  Tensor<> *realTabij(getTensorArgument("DoublesAmplitudes"));
  Tensor<complex> Tabij(
    4, realTabij->lens, realTabij->sym, *realTabij->wrld, "Tabij"
  );
  toComplexTensor(*realTabij, Tabij);

  //construct SG
  NG = CGai.lens[0];
  Vector<complex> *SG(new Vector<complex>(NG, *CGai.wrld, "SG"));
  (*SG)["G"] =   2.0 * conjCGai["Gai"] * CGai["Gbj"] * Tabij["abij"];
// BUG: the following line yields wrong sign:
//  (*SG)["G"] -= 1.0 * conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
  (*SG)["G"] += -1.0 * conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
  Vector<> *realSG(new Vector<>(NG, *CGai.wrld, "realSG"));
  fromComplexTensor(*SG, *realSG);
  allocatedTensorArgument<>("StructureFactor", realSG);
  //Get EMp2
  //Scalar<complex> EMp2(*CGai.wrld);
  //EMp2[""] = (*SG)["G"] * VG["G"];
  Scalar<> EMp2(*CGai.wrld);
  EMp2[""] = (*realSG)["G"] * (*realVG)["G"];
  double DEMp2(std::real(EMp2.get_val()));
  setRealArgument("EMp2", DEMp2);  

  allocatedTensorArgument<>("VG", realVG);

  structureFactors = new double[NG];
  realSG->read_all(structureFactors);
}

void FiniteSizeCorrection::calculateFiniteSizeCorrection() {
  // ...
  for (int g(0); g < NG; ++g) {
    LOG(1, "FiniteSizeCorrection") << structureFactors[g] << std::endl;
  }
}

void FiniteSizeCorrection::constructFibonacciGrid(int N, double R) {
  //This function construct a Fibonacci grid on a sphere with a certain radius.
  //Returns a vector of vectors: {x,y,z}
  //The N should be fixed and R should be a vector which is selected by another 
  //function which determines the R's
  double inc = M_PI * (3 - std::sqrt(5));
  double off = 2. / N;
  double r2d = 180./ M_PI;
  double x, y, z, r, phi, theta;
  new std::vector < vector<double>> fibonacciGrid;
  std::vector <double> coordinate;
  for (int k(0); k < N; ++k) {
    y = R* (k*off - 1. + 0.5*off);
    r = std::sqrt(R**2 - y*y);
    phi = k * inc;
    x = std::cos(phi)*r;
    z = std::sin(phi)*r;
    coordinate.push_back(x);
    coordinate.push_back(y);
    coordinate.push_back(z);
    fibonacciGrid.push_back(coordinate);
  }
}

