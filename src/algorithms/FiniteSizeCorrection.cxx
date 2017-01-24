#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Vector.hpp>
#include <math/Interpolation.hpp>
#include <tcc/DryTensor.hpp>
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
  //constructFibonacciGrid();
  interpolation3D();
  calculateFiniteSizeCorrection();
}


class FiniteSizeCorrection::Momentum {
  public:
    cc4s::Vector<> v;
    double s;
    double l;
    double vg;
    Momentum(): s(0.0), l(0.0) ,vg(0.) {
    }
    Momentum(cc4s::Vector<> v_, double s_=0., double vg_=0.) {
      v = v_; 
      s = s_;
      l = v_.length();
      vg = vg_;
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
  Scalar<> EMp2(*CGai.wrld);
  EMp2[""] = (*realSG)["G"] * (*realVG)["G"];
  double DEMp2(std::real(EMp2.get_val()));
  setRealArgument("EMp2", DEMp2);  

  allocatedTensorArgument<>("VG", realVG);
  VofG = new double[NG];
  realVG->read_all(VofG);
  structureFactors = new double[NG];
  realSG->read_all(structureFactors);
}


void FiniteSizeCorrection::constructFibonacciGrid(double R) {
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
  cc4s::Vector<> *regularGrid(new cc4s::Vector<>[NG]);
  momenta->read_all(&(regularGrid[0][0]));
  momentumGrid = new Momentum[2*NG];
  for (int g(0); g<NG; ++g) {
    momentumGrid[g] = Momentum(regularGrid[g], structureFactors[g], VofG[g]);
    momentumGrid[g+NG] = Momentum(regularGrid[g]*(-1.), structureFactors[g], VofG[g]);
  }

  //sort according to vector length. 
  
  std::sort(momentumGrid, &momentumGrid[2*NG], Momentum::sortbyl);
  

  //get the 3 unit vectors;
  cc4s::Vector<> a(momentumGrid[2].v);

  // GC is the shortest vector.
  if (isArgumentGiven("shortestGvector")) {
    GC = getRealArgument("shortestGvector");
  }
  else {
    GC = a.length();
  }

  LOG(1, "GridSearch") << "b1=#2" << std::endl;
  //the 0th and 1st elements are 0, avoid it.
  int j=3;
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
  for (int d(0); d<2*NG; ++d){
    x = T[0].dot(momentumGrid[d].v);
    y = T[1].dot(momentumGrid[d].v);
    z = T[2].dot(momentumGrid[d].v);
    momentumGrid[d].v[0] = (abs(x) < 1e-8) ? 0 : x;
    momentumGrid[d].v[1] = (abs(y) < 1e-8) ? 0 : y;
    momentumGrid[d].v[2] = (abs(z) < 1e-8) ? 0 : z;
  }
 
  //Determine the radii at which to construct the fibonacciGrids.
  std::sort(regularGrid, &regularGrid[NG], Vector<double,3>::sortByLength);
  numBins=1;
  for (int d(1); d < NG; ++d) {
    if (abs(regularGrid[d].length()-regularGrid[d-1].length()) < 1e-3) continue;
    else ++numBins;
  }
  aveSG = new double[numBins];
  lengthG = new double[numBins];
  aveSG[0]=0.;
  lengthG[0]=0.;
  numBins = 1;
  for (int d(1); d < NG; ++d) {
    if (abs(regularGrid[d].length()-regularGrid[d-1].length()) < 1e-3) 
      continue;
    constructFibonacciGrid(regularGrid[d].length());
    for (int g(0); g<N; ++g){
      x = T[0].dot(fibonacciGrid[g].v);
      y = T[1].dot(fibonacciGrid[g].v);
      z = T[2].dot(fibonacciGrid[g].v);
      fibonacciGrid[g].v[0] = (abs(x) < 1e-8) ? 0 : x;
      fibonacciGrid[g].v[1] = (abs(y) < 1e-8) ? 0 : y;
      fibonacciGrid[g].v[2] = (abs(z) < 1e-8) ? 0 : z;
    }

    //Trilinear interpolation on each point
    Momentum vertex[8];
    double average=0.;
    for (int t(0); t < N; ++t) {
      int xmin=std::floor(fibonacciGrid[t].v[0]);
      int xmax=std::ceil(fibonacciGrid[t].v[0]);
      int ymin=std::floor(fibonacciGrid[t].v[1]);
      int ymax=std::ceil(fibonacciGrid[t].v[1]);
      int zmin=std::floor(fibonacciGrid[t].v[2]);
      int zmax=std::ceil(fibonacciGrid[t].v[2]);
      vertex[0].v[0] = xmin;
      vertex[0].v[1] = ymin;
      vertex[0].v[2] = zmin;
      vertex[1].v[0] = xmin;
      vertex[1].v[1] = ymax;
      vertex[1].v[2] = zmin;
      vertex[2].v[0] = xmin;
      vertex[2].v[1] = ymin;
      vertex[2].v[2] = zmax;
      vertex[3].v[0] = xmin;
      vertex[3].v[1] = ymax;
      vertex[3].v[2] = zmax;
      vertex[4].v[0] = xmax;
      vertex[4].v[1] = ymin;
      vertex[4].v[2] = zmin;
      vertex[5].v[0] = xmax;
      vertex[5].v[1] = ymax;
      vertex[5].v[2] = zmin;
      vertex[6].v[0] = xmax;
      vertex[6].v[1] = ymin;
      vertex[6].v[2] = zmax;
      vertex[7].v[0] = xmax;
      vertex[7].v[1] = ymax;
      vertex[7].v[2] = zmax;

      double x[3] = {
        fibonacciGrid[t].v[0]-xmin, fibonacciGrid[t].v[1]-ymin,
        fibonacciGrid[t].v[2]-zmin
                    };

      double v[8] = {
        vertex[0].locate(momentumGrid,2*NG),vertex[1].locate(momentumGrid,2*NG),
        vertex[2].locate(momentumGrid,2*NG),vertex[3].locate(momentumGrid,2*NG),
        vertex[4].locate(momentumGrid,2*NG),vertex[5].locate(momentumGrid,2*NG),
        vertex[6].locate(momentumGrid,2*NG),vertex[7].locate(momentumGrid,2*NG)
                    };
      
      cc4s::Inter3D<double> intp;
      fibonacciGrid[t].s = intp.Trilinear(x,v);
      average += fibonacciGrid[t].s;
    }
    average = average / N; 
    aveSG[numBins] = average;
    lengthG[numBins] =  regularGrid[d].length();
    numBins++;
    
  }  
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
  cc4s::Inter1D<double> Int1d(numBins, lengthG, aveSG);
  //double xx[200], yy[200];
  //for (int j(0); j < 200; j++) {
  //  xx[j] = j*0.01;
  //  yy[j] = sin(xx[j]);
  //}
  //cc4s::Inter1D<double> Int1d(200, xx, yy);
  Int1d.cubicSpline(0., 0., "M");
  double x=0.;
  for (int i(1); i<1000; i++){
    LOG(1, "IntTest") << x << " " << Int1d.getValue(x) << std::endl;
    x = i*0.001;
  }
  for (int i(0); i<numBins; i++){
    LOG(1, "SGTest") << lengthG[i] << " " << aveSG[i] << std::endl;
  }
  int kpoints(getIntegerArgument("kpoints"));
  double volume(getRealArgument("volume"));
  double constantFactor(getRealArgument("constantFactor"));
  double r1 = integrate(Int1d, 0.0, GC, 1000)*constantFactor*volume*kpoints*4*M_PI;
  double sumSGVG(0.);
  double sumSGVG1(0.);
  for (int d(0); d < NG; ++d){
    sumSGVG1 += VofG[d] * structureFactors[d];
  }
  LOG(1,"integrate") << r1 << " sum= " << sumSGVG << " sum1= "
   << sumSGVG1<< " GC=" << GC << std::endl;
}
