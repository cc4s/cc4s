#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Vector.hpp>
#include <math/Interpolation.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

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
  interpolation3D();
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
  CTF::Vector<complex> *SG(new CTF::Vector<complex>(NG, *CGai.wrld, "SG"));
  (*SG)["G"] =   2.0 * conjCGai["Gai"] * CGai["Gbj"] * Tabij["abij"];
// BUG: the following line yields wrong sign:
//  (*SG)["G"] -= 1.0 * conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
  (*SG)["G"] += -1.0 * conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
  CTF::Vector<> *realSG(new CTF::Vector<>(NG, *CGai.wrld, "realSG"));
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
  //for (int g(0); g < NG; ++g) {
  //  LOG(1, "FiniteSizeCorrection") << structureFactors[g] << std::endl;
  //}
}

void FiniteSizeCorrection::constructFibonacciGrid() {
  //This function construct a Fibonacci grid on a sphere with a certain radius.
  //Returns a vector of vectors: {x,y,z}
  //The N should be fixed and R should be a vector which is selected by another 
  //function which determines the R's
  N = 128;
  double R = 1.0;
  double inc = M_PI * (3 - std::sqrt(5));
  //cc4s::Vector<> *fibonacciGrid(new cc4s::Vector<>[N]);
  std::vector<std::vector<double>> fibonacciGrid(N, std::vector<double>(3));
  for (int k(0); k < N; ++k) {
    double z((2.0*k+1)/N - 1.0);
    double r(R * std::sqrt(1.0 - z*z));
    double phi(k * inc);
    fibonacciGrid[k][0] = r * std::cos(phi);
    fibonacciGrid[k][1] = r * std::sin(phi);
    fibonacciGrid[k][2] = R * z;
    //LOG(1, "FibonacciGrid") << z << "; " << fibonacciGrid[k] << std::endl;
  }
//  LOG(1, "FibonacciGrid") << fibonacciGrid[0].approximately(fibonacciGrid[1]) << std::endl;
//  LOG(1, "FibonacciGrid") << fibonacciGrid[1].approximately(fibonacciGrid[1]) << std::endl;
}

void FiniteSizeCorrection::interpolation3D() {
  Tensor<> *momenta(getTensorArgument<>("Momenta"));
  cc4s::Vector<> *regularGrid(new cc4s::Vector<>[NG]);
// or alternatively:
//  std::vector<cc4s::Vector<>> regularGrid(NG);
  momenta->read_all(&(regularGrid[0][0]));
  class Momentum {
  public:
    cc4s::Vector<> v;
    double s;
    double l;
    Momentum(): s(0.0), l(0.0) {
    }
    Momentum(cc4s::Vector<> v_, double s_) {
      v = v_; 
      s = s_;
      l = v_.length();
    }
    bool operator < (Momentum const &m) {
      return v < m.v;
    }
  };
  Momentum *momentumGrid(new Momentum[NG]);
  for (int g(0); g<NG; ++g) {
    momentumGrid[g] = Momentum(regularGrid[g], structureFactors[g]);
  }
  
  std::sort(momentumGrid, &momentumGrid[NG]);
  //look for the lengths of unit vectors in each direction
  cc4s::Vector<double,NG> x, y, z;
  for (int t(0); t<NG; ++t){
    x[t] = std::abs(regularGrid[t][0]);
    y[t] = std::abs(regularGrid[t][1]);
    z[t] = std::abs(regularGrid[t][2]);
  }
  std::sort(x[0], x[NG]);
  std::sort(y[0], y[NG]);
  std::sort(z[0], z[NG]);
  //for (int g(0); g<NG; ++g) {
  //  LOG(1, "Sorted")  << "momentumGrid[" << g << "]=" << momentumGrid[g].v 
  //  << ",l " << momentumGrid[g].l << ", " << momentumGrid[g].s << std::endl;
  //}

    //LOG(1, "test") << "length[" << g << "]=" << regularGrid[g].length() << std::endl;
    //LOG(1, "Unsorted")  << "momentumGrid[" << g << "]=" << momentumGrid[g].v
    //<< ",l " << momentumGrid[g].l << ", " << momentumGrid[g].s << std::endl;
 // for (int n(0); n < NG; ++n){
 //   onRegularGrid[n].initial(regularGrid[n], structureFactors[n]);
 // }
 // cout << onRegularGrid.coordinates
 // 
 // for (int t(0); t < N; ++t){
 //   for (int d(0); d < NG; ++d){
 //      
 //   }
 // }
  //double x[3]={0.3,0.4, 0.7};
  //double v[8]={0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  //cc4s::Interpolation<double> Intp;
  //double tmp = Intp.Trilinear(x,v);
  //LOG(1, "linear") << tmp << std::endl;
  //std::cout << "Linear" << cc4s::Interpolation<>::Linear(x, v) <<std::endl;
  //remove duplicated grid points
  //for (int d(0); d < NG; ++d) {
  //  LOG(1, "regularGrid") << regularGrid[d] << std::endl;
  //}
}

