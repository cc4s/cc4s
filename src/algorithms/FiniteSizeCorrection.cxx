#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Vector.hpp>
#include <math/Interpolation.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <iostream>
// FIXME: use common way for math constants
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
    Momentum(): s(0.0), l(0.0), vg(0.) {
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

    static bool sortByLength (Momentum const &n, Momentum const &m) {
      return n.l < m.l;
    }
    static bool sortByVector (Momentum const &n, Momentum const &m) {
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
  //Get energy from amplitudes
  Scalar<> EMp2(*CGai.wrld);
  EMp2[""] = (*realSG)["G"] * (*realVG)["G"];
  double DEMp2(std::real(EMp2.get_val()));
  setRealArgument("EnergyFromAmplitudes", DEMp2);

  //  allocatedTensorArgument<>("VG", realVG);
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
  cc4s::Vector<> *cartesianMomenta(new cc4s::Vector<>[NG]);
  momenta->read_all(&cartesianMomenta[0][0]);

  // FIXME: give direct or reciprocal grid and calaculate all properties
  // from it.
  cartesianGrid = new Momentum[2*NG];
  // complete momentum grid in a Gamma only calculation
  for (int g(0); g<NG; ++g) {
    cartesianGrid[g] = Momentum(
      cartesianMomenta[g], structureFactors[g], VofG[g]
    );
    cartesianGrid[g+NG] = Momentum(
      cartesianMomenta[g]*(-1.), structureFactors[g], VofG[g]
    );
  }

  // sort according to vector length.
  std::sort(cartesianGrid, &cartesianGrid[2*NG], Momentum::sortByLength);

  // get the 3 unit vectors;
  cc4s::Vector<> a(cartesianGrid[2].v);

  // GC is the shortest vector.
  if (isArgumentGiven("shortestGvector")) {
    GC = getRealArgument("shortestGvector");
  }
  else {
    GC = a.length();
  }

  LOG(2, "GridSearch") << "b1=#2" << std::endl;
  // the 0th and 1st elements are 0, avoid it.
  int j=3;
  //a and b should not be parallel;
  while ((a.cross(cartesianGrid[j].v)).length() < 1e-8) ++j;
  cc4s::Vector<> b(cartesianGrid[j].v);
  LOG(2, "GridSearch") << "b2=#" << j << std::endl;
  ++j;
  //a, b and c should not be on the same plane;
  while (abs((a.cross(b)).dot(cartesianGrid[j].v)) < 1e-8) ++j;
  cc4s::Vector<> c(cartesianGrid[j].v);
  LOG(2, "GridSearch") << "b3=#" << j << std::endl;

  // print the basis vectors
  LOG(2, "GridSearch") << "b1=" << a << std::endl;
  LOG(2, "GridSearch") << "b2=" << b << std::endl;
  LOG(2, "GridSearch") << "b3=" << c << std::endl;

  //construct the transformation matrix
  cc4s::Vector<> *T(new cc4s::Vector<>[3]);
  double Omega((a.cross(b)).dot(c));
  T[0] = b.cross(c)/Omega;
  T[1] = c.cross(a)/Omega;
  T[2] = a.cross(b)/Omega;

  // determine bounding box in direct (reciprocal) coordinates
  Vector<> directMin, directMax;
  for (int g(0); g < 2*NG; ++g) {
    for (int d(0); d < 3; ++d) {
      double directComponent(T[d].dot(cartesianGrid[g].v));
      directMin[d] = std::min(directMin[d], directComponent);
      directMax[d] = std::max(directMax[d], directComponent);
    }
  }

  // build grid for the entire bounding box
  int boxDimensions[3], boxOrigin[3];
  int64_t boxSize(1);
  for (int d(0); d < 3; ++d) {
    boxSize *=
      boxDimensions[d] = static_cast<int>(directMax[d] - directMin[d] + 0.5);
    boxOrigin[d] = static_cast<int>(directMin[d] + 0.5);
  }
  // allocate and initialize regular grid
  double *regularSG(new double[boxSize]);
  for (int64_t g(0); g < boxSize; ++g) regularSG[g] = 0;
  // enter known Sq values
  for (int g(0); g < 2*NG; ++g) {
    int64_t index(0);
    for (int d(2); d >= 0; --d) {
      int directComponent(static_cast<int>(T[d].dot(cartesianGrid[g].v) + 0.5));
      index *= boxDimensions[d];
      index += directComponent - boxOrigin[d];
    }
    regularSG[index] = cartesianGrid[g].s;
  }

  // TODO: create trilinear or tricubic interpolator

  // spherically sample
  double lastLength(0);
  averageSGs.clear(); GLengths.clear();
  for (int g(1); g < 2*NG; ++g) {
    double length(cartesianGrid[g].l);
    if (abs(length - lastLength) > 1e-3) {
      constructFibonacciGrid(length);
      double averageSG(0);
      for (int f(0); f < N; ++f) {
        Vector<> sampledDirect;
        for (int d(0); d < 3; ++d) {
          sampledDirect[d] = T[d].dot(fibonacciGrid[f].v);
        }
        // TODO: lookup interpolated value in direct coordinates
        averageSG += 0;
      }
      averageSGs.push_back(averageSG / N);
      GLengths.push_back(length);
    }
    lastLength = length;
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
  cc4s::Inter1D<double> Int1d(
    GLengths.size(), GLengths.data(), averageSGs.data()
  );
  //double xx[200], yy[200];
  //for (int j(0); j < 200; j++) {
  //  xx[j] = j*0.01;
  //  yy[j] = sin(xx[j]);
  //}
  //cc4s::Inter1D<double> Int1d(200, xx, yy);
  Int1d.cubicSpline(0., 0., "M");
  double x=0.;
  for (int i(1); i<1000; i++){
    LOG(2, "Interpolation") << x << " " << Int1d.getValue(x) << std::endl;
    x = i*0.001;
  }
  for (unsigned int i(0); i < GLengths.size(); ++i){
    LOG(2, "StructureFactor") << GLengths[i] << " " << averageSGs[i] << std::endl;
  }
  int kpoints(getIntegerArgument("kpoints", 1));
  double volume(getRealArgument("volume"));
  double constantFactor(getRealArgument("constantFactor"));
  // the factor 2 is only needed when half of the G grid is provided (this is the case
  // for a Gamma point calculation in vasp).
  double r1 = 2*integrate(Int1d, 0.0, GC, 1000)*constantFactor*volume*kpoints*4*M_PI;
  double  sumSGVG(0.);
  for (int d(0); d < NG; ++d){
    sumSGVG += VofG[d] * structureFactors[d];
  }
  LOG(1,"FiniteSize") << "Uncorrected e="  << sumSGVG       << std::endl;
  LOG(1,"FiniteSize") << "Correction  e="  << r1            << std::endl;
  LOG(1,"FiniteSize") << "Corrected   e="  << sumSGVG + r1  << std::endl;

  setRealArgument("CorrectedEnergy"  , sumSGVG + r1);
}
