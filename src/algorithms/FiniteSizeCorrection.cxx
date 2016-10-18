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


class FiniteSizeCorrection::Momentum {
  public:
    cc4s::Vector<> v;
    double s;
    double l;
    Momentum(): s(0.0), l(0.0) {
    }
    Momentum(cc4s::Vector<> v_, double s_=0.) {
      v = v_; 
      s = s_;
      l = v_.length();
    }
    double locate(Momentum *m, int const n) {
      cc4s::Vector<> u(v);
      if (v[3] < 0.) u= v*(-1.);
      for (int d(0); d < n; ++d) {
       // LOG(1,"locating") << "d=#" << d <<std::endl;
        if (u.approximately(m[d].v)) return m[d].s;
      }  
      return 0;
    }
    static bool sortbyl (Momentum const &n, Momentum const &m) {
      return n.l < m.l;
    }
    static bool sortbyv (Momentum const &n, Momentum const &m) {
      return n.v < m.v;
    }
};



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
  double R = 0.2;
  double inc = M_PI * (3 - std::sqrt(5));
  fibonacciGrid = new Momentum[N];
  //std::vector<std::vector<double>> fibonacciGrid(N, std::vector<double>(3));
  for (int k(0); k < N; ++k) {
    double z((2.0*k+1)/N - 1.0);
    double r(R * std::sqrt(1.0 - z*z));
    double phi(k * inc);
    fibonacciGrid[k].v[0] = r * std::cos(phi);
    fibonacciGrid[k].v[1] = r * std::sin(phi);
    fibonacciGrid[k].v[2] = R * z;
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
  Momentum *momentumGrid(new Momentum[NG]);
  for (int g(0); g<NG; ++g) {
    momentumGrid[g] = Momentum(regularGrid[g], structureFactors[g]);
  }

  //test sort by v
  //std::sort(momentumGrid, &momentumGrid[NG], Momentum::sortbyv);
  ////
  //for (int g(0); g<NG; ++g) {
  //  LOG(1, "Sortedbyv")  << "momentumGrid[" << g << "]=" << momentumGrid[g].v 
  //  << ",l " << momentumGrid[g].l << ", " << momentumGrid[g].s << std::endl;
  //}
  //sort according to vector length. 
  
  std::sort(momentumGrid, &momentumGrid[NG], Momentum::sortbyl);
  
  for (int g(0); g<NG; ++g) {
    LOG(1, "Sortedbyl")  << "momentumGrid[" << g << "]=" << momentumGrid[g].v 
    << ",l " << momentumGrid[g].l << ", " << momentumGrid[g].s << std::endl;
  }

  //get the 3 unit vectors;
  cc4s::Vector<> a(momentumGrid[1].v);
  LOG(1, "GridSearch") << "b1=#1" << std::endl;
  //the 0th element is 0, avoid it.
  int j=2;
  //a and b should not be parallel;
  while ((a.cross(momentumGrid[j].v)).length() < 1e-8) ++j;
  cc4s::Vector<> b(momentumGrid[j].v);
  LOG(1, "GridSearch") << "b2=#" << j << std::endl;
  ++j;
  //a, b and c should not be on the same plane;
  while (abs((a.cross(b)).dot(momentumGrid[j].v)) < 1e-8) ++j;
  cc4s::Vector<> c(momentumGrid[j].v);
  LOG(1, "GridSearch") << "b3=#" << j << std::endl;
  LOG(1, "GridSearch") << "b1=" << a << std::endl;
  LOG(1, "GridSearch") << "b2=" << b << std::endl;
  LOG(1, "GridSearch") << "b3=" << c << std::endl;
  
  //construct the transformation matrix  
  cc4s::Vector<> *T(new cc4s::Vector<>[3]);
  double Omega((a.cross(b)).dot(c));
  T[0] = b.cross(c)/Omega;
  T[1] = c.cross(a)/Omega;
  T[2] = a.cross(b)/Omega;
  double x, y, z;
  for (int d(0); d<NG; ++d){
    x = T[0].dot(momentumGrid[d].v);
    y = T[1].dot(momentumGrid[d].v);
    z = T[2].dot(momentumGrid[d].v);
    momentumGrid[d].v[0] = (abs(x) < 1e-8) ? 0 : x;
    momentumGrid[d].v[1] = (abs(y) < 1e-8) ? 0 : y;
    momentumGrid[d].v[2] = (abs(z) < 1e-8) ? 0 : z;
   // LOG(1, "Transformed")  << "momentumGrid[" << d << "]=" << momentumGrid[d].v 
   // << ",l " << momentumGrid[d].l << ", " << momentumGrid[d].s << std::endl;
  }
  //for (int d(0); d<NG; ++d){
  //  LOG(1, "Locate")  << momentumGrid[d].locate(momentumGrid, NG) 
  //    << " momentumGrid["<< d << "]" << momentumGrid[d].s << std::endl;
  //}
 
  for (int d(0); d<N; ++d){
    x = T[0].dot(fibonacciGrid[d].v);
    y = T[1].dot(fibonacciGrid[d].v);
    z = T[2].dot(fibonacciGrid[d].v);
    fibonacciGrid[d].v[0] = (abs(x) < 1e-8) ? 0 : x;
    fibonacciGrid[d].v[1] = (abs(y) < 1e-8) ? 0 : y;
    fibonacciGrid[d].v[2] = (abs(z) < 1e-8) ? 0 : z;
    //LOG(1, "Transformed")  << "fibonacciGrid[" << d << "]=" << fibonacciGrid[d].v 
    //<< ",l " << fibonacciGrid[d].l << ", " << fibonacciGrid[d].s << std::endl;
  }
  //sort momentumGrid w.r.t. x, y, z in order to find the SG
  std::sort(momentumGrid, &momentumGrid[NG], Momentum::sortbyv);
  for (int d(0); d<NG; ++d){
    LOG(1, "Sortbyv")  << "momentumGrid[" << d << "]=" << momentumGrid[d].v 
    << ",l " << momentumGrid[d].l << ", " << momentumGrid[d].s << std::endl;
  }
  //Trilinear interpolation on each point
  Momentum vertex[8];
  for (int d(0); d <N; ++d) {
    int xmin=std::floor(fibonacciGrid[d].v[0]);
    int xmax=std::ceil(fibonacciGrid[d].v[0]);
    int ymin=std::floor(fibonacciGrid[d].v[1]);
    int ymax=std::ceil(fibonacciGrid[d].v[1]);
    int zmin=std::floor(fibonacciGrid[d].v[2]);
    int zmax=std::ceil(fibonacciGrid[d].v[2]);
    vertex[0].v[0] = xmax;
    vertex[0].v[1] = ymin;
    vertex[0].v[2] = zmin;
    vertex[1].v[0] = xmax;
    vertex[1].v[1] = ymax;
    vertex[1].v[2] = zmin;
    vertex[2].v[0] = xmax;
    vertex[2].v[1] = ymax;
    vertex[2].v[2] = zmax;
    vertex[3].v[0] = xmax;
    vertex[3].v[1] = ymin;
    vertex[3].v[2] = zmax;
    vertex[4].v[0] = xmin;
    vertex[4].v[1] = ymin;
    vertex[4].v[2] = zmin;
    vertex[5].v[0] = xmin;
    vertex[5].v[1] = ymax;
    vertex[5].v[2] = zmin;
    vertex[6].v[0] = xmin;
    vertex[6].v[1] = ymax;
    vertex[6].v[2] = zmax;
    vertex[7].v[0] = xmin;
    vertex[7].v[1] = ymin;
    vertex[7].v[2] = zmax;
    double x[3] = {
      fibonacciGrid[d].v[0]-xmin, fibonacciGrid[d].v[1]-ymin,
      fibonacciGrid[d].v[2]-zmin
                  };
    double v[8] = {
      vertex[0].locate(momentumGrid,NG),vertex[1].locate(momentumGrid,NG),
      vertex[2].locate(momentumGrid,NG),vertex[3].locate(momentumGrid,NG),
      vertex[4].locate(momentumGrid,NG),vertex[5].locate(momentumGrid,NG),
      vertex[6].locate(momentumGrid,NG),vertex[7].locate(momentumGrid,NG)
                  };
    for (int t(0); t< 8; ++t) {
      LOG(1, "vertex") << v[t] << std::endl; 
    }
    
    cc4s::Interpolation<double> intp;
    fibonacciGrid[d].s = intp.Trilinear(x,v);
    LOG(1, "linear") << fibonacciGrid[d].s << " v="
      << fibonacciGrid[d].v << std::endl;
  } 

  

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

