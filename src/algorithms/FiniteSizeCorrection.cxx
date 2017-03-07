#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Vector.hpp>
#include <math/Interpolation.hpp>
#include <gte/TricubicInterpolation.hpp>
#include <gte/TrilinearInterpolation.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <iostream>
// FIXME: use common way for math constants
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "mpi.h"

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
  }
  else {
    LOG(0,"FiniteSize") << "Calculating structure factor" << std::endl;
    calculateStructureFactor();
  }
  //constructFibonacciGrid();
  LOG(0,"FiniteSize") << "Interpolating and integrating" << std::endl;
  interpolation3D();
  LOG(0,"FiniteSize") << "Caclulating finite size correction" << std::endl;
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

void FiniteSizeCorrection::readFromFile(){
  LOG(0,"readFromFile") << "reading " << std::endl;
  Tensor<> *realVG(getTensorArgument<>("CoulombKernel"));
  Tensor<> *realSG(getTensorArgument<>("StructureFactor"));
  LOG(1,"readFromFile") << "success\n Loading into Vectors " << std::endl;
  NG=realVG->lens[0];
  VofG = new double[NG];
  realVG->read_all(VofG);
  LOG(1,"readFromFile") << "VofG Finished" << std::endl;
  structureFactors = new double[NG];
  realSG->read_all(structureFactors);
  LOG(1,"readFromFile") << "Finished" << std::endl;
}

void FiniteSizeCorrection::calculateStructureFactor() {

//Definition of the variables
  Tensor<complex> *GammaFai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );

  Tensor<complex> *GammaGai;

  if (isArgumentGiven("CoulombVertexSingularVectors")) {
    Tensor<complex> *UGF(
      getTensorArgument<complex>("CoulombVertexSingularVectors")
    );
    int lens[]= {UGF->lens[0], GammaFai->lens[1], GammaFai->lens[2]};
    GammaGai = new Tensor<complex>(
      3, lens, GammaFai->sym, *GammaFai->wrld, "GammaGai"
    );
    (*GammaGai)["Gai"] = (*GammaFai)["Fai"] * (*UGF)["GF"];
  } else {
    GammaGai = GammaFai;
  }

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

  if (isArgumentGiven("SinglesAmplitudes") ) {
    //Get Tai
    Tensor<> *realTai(getTensorArgument("SinglesAmplitudes"));
    Tensor<complex> Tai(
      2, realTai->lens, realTai->sym, *realTai->wrld, "Tai"
    );
    toComplexTensor(*realTai, Tai);
    Tabij["abij"] += Tai["ai"] * Tai["bj"];
  }

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
  //  Scalar<> EMp2(*CGai.wrld);
  //  EMp2[""] = (*realSG)["G"] * (*realVG)["G"];
  //  double DEMp2(std::real(EMp2.get_val()));
  //  setRealArgument("EnergyFromAmplitudes", DEMp2);

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

  // FIXME: give direct or reciprocal grid and calaculate all properties,
  // such as a,b,c and omega from it.
  cartesianGrid = new Momentum[2*NG];
  // complete momentum grid in a Gamma only calculation
  // TODO: detect if grid needs to be doubled
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
  // determine a small length on the scale of the cell
  const double epsilon(1e-9*a.length());

  LOG(2, "GridSearch") << "b1=#2" << std::endl;
  // the 0th and 1st elements are 0, avoid it.
  int j=3;
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
  for (int g(0); g < 2*NG; ++g) {
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
  double lastLength(-1);
  averageSGs.clear(); GLengths.clear();
  for (int g(0); g < 2*NG; ++g) {
    double length(cartesianGrid[g].l);
    if (abs(length - lastLength) > 1e-3) {
      constructFibonacciGrid(length);
      double sumSG(0);
      // TODO: use parameter instead of fixed Fibonacci grid size
      for (int f(0); f < N; ++f) {
        Vector<> directG;
        for (int d(0); d < 3; ++d) {
          directG[d] = T[d].dot(fibonacciGrid[f].v);
        }
        // lookup interpolated value in direct coordinates
        sumSG += interpolatedSG(directG[0], directG[1], directG[2]);
      }
      averageSGs.push_back(sumSG / N);
      GLengths.push_back(length);
      lastLength = length;
    }
  }

  //Define the 3D zone close to the Gamma point which needed to be integrated over. Find the vectors which define it. Small BZ
  //  for (int t(0); t < 20; t++){
  //    LOG(0,"G vectors by length") << cartesianGrid[t].v << std::endl;
  //  }

  std::vector<Vector<>> smallBZ;
  smallBZ.push_back(cartesianGrid[2].v);
  for (int t(3); t<2*NG; t++){
    if (IsInSmallBZ(cartesianGrid[t].v, 1., smallBZ)){
      smallBZ.push_back(cartesianGrid[t].v);
    }
  }

  LOG(0,"FiniteSize") << "Size of vectors="  << smallBZ.size() << std::endl;
  for (std::vector<int>::size_type i = 0; i != smallBZ.size(); i++){
    LOG(2,"interpolation3D") << "smallBZ basis vector: " << smallBZ[i] << std::endl;
  }

  //integration in 3D
  double constantFactor(getRealArgument("constantFactor"));
  double cutOffRadius(getRealArgument("cutOffRadius", 1e-5));
  int N0(51), N1(51), N2(51);
  inter3D = 0.;
  sum3D   = 0.;
  int countNO(0);
  int countNOg(0);
  std::vector<Vector<>> gridWithinRadius;
  for (int i(0); i < 2*NG; ++i) {
    if (cartesianGrid[i].l < cutOffRadius) {
      gridWithinRadius.push_back(cartesianGrid[i].v);
      if (cartesianGrid[i].l > epsilon) {
        sum3D += constantFactor/cartesianGrid[i].l/cartesianGrid[i].l
                 *cartesianGrid[i].s;
      }
    }
  }

  MPI_Barrier(Cc4s::world->comm);
  for (int t0(-N0); t0 <= N0; ++t0) {
    for (int t1(-N1); t1 <= N1; ++t1) {
      for (int t2(-N2); t2 <= N2; ++t2) {
        if ((std::abs(t0+t1+t2)) % Cc4s::world->np == Cc4s::world->rank) {
          Vector<double> directg;
          Vector<double> ga(((a/double(N0))*double(t0)));
          Vector<double> gb(((b/double(N1))*double(t1)));
          Vector<double> gc(((c/double(N2))*double(t2)));
          Vector<double> g(ga+gb+gc);
          if (IsInSmallBZ(g, 2, smallBZ)){
            countNOg++;
            for (std::vector<int>::size_type i = 0; i != gridWithinRadius.size(); ++i) {
              g=ga+gb+gc;
              g += gridWithinRadius[i];
              for (int d(0); d <3; ++d){
                directg[d]=T[d].dot(g);
              }
              countNO++;
              if (g.length() > epsilon) {
                inter3D += interpolatedSG(
                  directg[0], directg[1],
                  directg[2]
                ) * constantFactor/g.length()/g.length();
              }
            }
          }
        }
      }
    }
  }
  double totalInter3D(0);
  int totalCountNOg(0);
  MPI_Allreduce(
    &inter3D, &totalInter3D, 1, MPI_DOUBLE, MPI_SUM, Cc4s::world->comm
  );
  MPI_Allreduce(
    &countNOg, &totalCountNOg, 1, MPI_INT, MPI_SUM, Cc4s::world->comm
  );
  MPI_Barrier(Cc4s::world->comm);

  LOG(2,"integration3D") << "countNOg= " << totalCountNOg <<std::endl;
  sum3D = sum3D/2.; //Both readFromFile=0 and 1 needs to be divided by 2
  LOG(2, "interpolation3D") << "sum3D= " << sum3D << std::endl;
  inter3D=totalInter3D/totalCountNOg;
  LOG(2,"interpolation3D") << "Number of points in summation="<< countNO << std::endl;
}

bool FiniteSizeCorrection::IsInSmallBZ(
  Vector<> point, double scale, std::vector<Vector<>> smallBZ
)
{
  std::vector<int>::size_type countVector(0);
  for (std::vector<int>::size_type i = 0; i != smallBZ.size(); i++){
    // FIXME: use an epsilon instead of 1e-9
    if (
      abs(abs(smallBZ[i].dot(point))/smallBZ[i].length()/smallBZ[i].length()
       *scale -1.0) < 1e-9 ||
      abs(smallBZ[i].dot(point))/smallBZ[i].length()/smallBZ[i].length()*scale
       -1.0 < -1e-9
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
  int kpoints(getIntegerArgument("kpoints", 1));
  double volume(getRealArgument("volume"));
  double constantFactor(getRealArgument("constantFactor"));
  // the factor 1/2 is only needed when half of the G grid is provided (this is the case
  // for a Gamma point calculation in vasp).
  double r1 = 1./2.*integrate(Int1d, 0.0, GC, 1000)*constantFactor*volume*kpoints*4*M_PI;
  double  sumSGVG(0.);
  int fReadFromFile(getIntegerArgument("fReadFromFile", 0));
  if (fReadFromFile == 1) r1=r1*2.;
  //TODO:Exclude the integrated region from the summation.
  for (int d(0); d < NG; ++d){
    sumSGVG += VofG[d] * structureFactors[d];
    }
  inter3D = inter3D/2.;//both read or not from file should divide by 2

  LOG(0,"FiniteSize") << "Uncorrected e="  << std::setprecision(10) << sumSGVG
    << std::endl;
  LOG(0,"FiniteSize") << "Corrected   e= " << std::setprecision(10) 
    << sumSGVG+inter3D-sum3D << std::endl;

  LOG(1,"FiniteSize") << "Integral within cutoff radius= " 
    << std::setprecision(10) << inter3D    << std::endl;
  LOG(1,"FiniteSize") << "Sumation within cutoff radius= " 
    << std::setprecision(10) << sum3D      << std::endl;

  LOG(1,"FiniteSize") << "Spherical Averaging Correction  e="
    << std::setprecision(10)  << r1        << std::endl;
  LOG(1,"FiniteSize") << "Spherical Averaging Corrected   e="
    << std::setprecision(10)  << sumSGVG + r1 << std::endl;

  setRealArgument("CorrectedEnergy"  , sumSGVG+inter3D-sum3D);
}
