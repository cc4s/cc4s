#include <algorithms/RpaApxEnergy.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

#include <ctf.hpp>
#include <algorithm>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(RpaApxEnergy);

RpaApxEnergy::RpaApxEnergy(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

RpaApxEnergy::~RpaApxEnergy() {
}

namespace cc4s {
  // provides vector space structure for imaginary frequency grids
  class FrequencyGrid {
  public:
    FrequencyGrid(const int N): nus(N), ws(N) { }
    FrequencyGrid(const FrequencyGrid &g): nus(g.nus), ws(g.ws) { }
    FrequencyGrid &operator =(const FrequencyGrid &g) {
      nus = g.nus; ws = g.ws;
      return *this;
    }
    FrequencyGrid &operator +=(const FrequencyGrid &g) {
      for (size_t n(0); n < nus.size(); ++n) {
        nus[n] += g.nus[n];
        ws[n] += g.ws[n];
      }
      return *this;
    }
    FrequencyGrid &operator -=(const FrequencyGrid &g) {
      for (size_t n(0); n < nus.size(); ++n) {
        nus[n] -= g.nus[n];
        ws[n] -= g.ws[n];
      }
      return *this;
    }
    FrequencyGrid &operator *=(double s) {
      for (size_t n(0); n < nus.size(); ++n) {
        nus[n] *= s;
        ws[n] *= s;
      }
      return *this;
    }
    FrequencyGrid &operator /=(double s) {
      for (size_t n(0); n < nus.size(); ++n) {
        nus[n] /= s;
        ws[n] /= s;
      }
      return *this;
    }
    FrequencyGrid &operator -() {
      for (size_t n(0); n < nus.size(); ++n) {
        nus[n] = -nus[n];
        ws[n] = -ws[n];
      }
      return *this;
    }
    double dot(const FrequencyGrid &g) {
      double result(0.0);
      for (size_t n(0); n < nus.size(); ++n) {
        result += nus[n] * g.nus[n];
        result += ws[n] * g.ws[n];
      }
      return result;
    }
    std::vector<double> nus, ws;
  };

  inline FrequencyGrid operator +(
    const FrequencyGrid &a, const FrequencyGrid &b
  ) {
    FrequencyGrid result(a);
    return result += b;
  }

  inline FrequencyGrid operator -(
    const FrequencyGrid &a, const FrequencyGrid &b
  ) {
    FrequencyGrid result(a);
    return result -= b;
  }

  inline FrequencyGrid operator *(
    const double s, const FrequencyGrid &g
  ) {
    FrequencyGrid result(g);
    return result *= s;
  }

  inline FrequencyGrid operator /(
    const FrequencyGrid &g, const double s
  ) {
    FrequencyGrid result(g);
    return result /= s;
  }

  class FrequencyGridOptimizer {
  public:
    FrequencyGridOptimizer(
      const int N, CTF::Tensor<double> epsi, CTF::Tensor<double> epsa
    ):
      grid(N)
    {
      int No(epsi.lens[0]);
      int Nv(epsa.lens[0]);
      Dai = new CTF::Tensor<double>(2, std::vector<int>({Nv,No}).data());
      (*Dai)["ai"]  = epsa["a"];
      (*Dai)["ai"] -= epsi["i"];

      // take N quantiles as initial estimates for grid points
      std::vector<double> deltas(No*Nv);
      Dai->read_all(deltas.data(), true);
      std::sort(deltas.begin(), deltas.end());
      double lastNu(0);
      for (size_t n(0); n < grid.ws.size(); ++n) {
        grid.nus[n] = deltas[No*Nv*n/N + 0*No*Nv/(2*N)];
        grid.ws[n] = grid.nus[n] - lastNu;
        lastNu = grid.nus[n];
      }
      Eai = new CTF::Tensor<double>(false, *Dai);
    }

    void optimize(const int stepCount) {
      LOG(1, "RPA") << "optimizing grid" << std::endl;
      double E(getError(grid));
      FrequencyGrid lastDelta(-getGradient(grid));
      FrequencyGrid lastDirection(lastDelta);
      E = lineSearch(grid, lastDirection);
      NEW_FILE("E.dat");
      NEW_FILE("nu.dat");
      NEW_FILE("w.dat");
      for (int m(0); m < stepCount; ++m) {
        FrequencyGrid Delta(-getGradient(grid));
        double beta(
          std::max(0.0, Delta.dot(Delta-lastDelta) / lastDelta.dot(lastDelta))
        );
        FrequencyGrid direction(Delta + beta*lastDirection);
        E = lineSearch(grid, direction);
        FILE("E.dat") << m << " " << E << std::endl;
        FILE("nu.dat") << m;
        FILE("w.dat") << m;
        for (size_t n(0); n < grid.nus.size(); ++n) {
          FILE("nu.dat") << " " << grid.nus[n];
          FILE("w.dat") << " " << grid.ws[n];
        }
        FILE("nu.dat") << std::endl;
        FILE("w.dat") << std::endl;
        LOG(1, "RPA") << "error=" << E << ", beta=" << beta << std::endl;
        lastDirection = direction;
        lastDelta = Delta;
      }
    }

    double lineSearch(FrequencyGrid &grid, const FrequencyGrid &direction) {
      // starting vectors
      double E0(getError(grid));
      double E1(getError(grid + direction));
      double Em;
      double alpham(0.5);
      double alpha1(1.0);
      // increase alpha until error also increases
      while (E1 <= E0) {
        alpham = alpha1;
        alpha1 *= 2.0;
        Em = E1;
        E1 = getError(grid + alpha1*direction);
      } while (E1 <= E0);
      double alpha0(0.0);
      // interval bisection to machine precision
      while (alpham-alpha0 > 1e-16 && alpha1-alpham > 1e-16) {
        if (E0 < E1) {
          alpha1 = alpham;
          E1 = Em;
        } else {
          alpha0 = alpham;
          E0 = Em;
        }
        alpham = 0.5*(alpha0+alpha1);
        Em = getError(grid + alpham*direction);
      }
      grid += alpham*direction;
      return Em;
    }

    void testGradient(const double stepSize) {
      double E0(getError(grid));
      FrequencyGrid grad(getGradient(grid));
      FrequencyGrid numGrad(grid);
      for (size_t n(0); n < grid.ws.size(); ++n) {
        grid.ws[n] += stepSize;
        double E1(getError(grid));
        grid.ws[n] -= stepSize;
        numGrad.ws[n] = (E1-E0) / stepSize;
        LOG(1, "RPA") << "dE/dw[" << n << "] = " <<
          grad.ws[n] << " ~ " << numGrad.ws[n] << std::endl;
      }
      for (size_t n(0); n < grid.ws.size(); ++n) {
        grid.nus[n] += stepSize;
        double E1(getError(grid));
        grid.nus[n] -= stepSize;
        numGrad.nus[n] = (E1-E0) / stepSize;
        LOG(1, "RPA") << "dE/dnu[" << n << "] = " <<
          grad.nus[n] << " ~ " << numGrad.nus[n] << std::endl;
      }
    }

    static double propagator(
      const double delta, const double nu, const int e = 2
    ) {
      return 4*delta*delta / std::pow(delta*delta + nu*nu, e);
    }

    double getError(const FrequencyGrid &grid) {
      // error of numerical quadrature with current grid wn[n] & nus[n] from
      // analytic value 1/(eps_a-eps_i) for each ai
      CTF::Transform<double, double>(
        std::function<void(double, double &)>(
          [&grid](double delta, double &error) {
            error = 0.0;
            // use quadrature to integrate Mp2 propagator
            for (size_t n(0); n < grid.ws.size(); ++n) {
              error += grid.ws[n] * propagator(delta, grid.nus[n]);
            }
            error /= Pi();
            // subtract analytic result = 1/delta
            error -= 1 / delta;
          }
        )
      ) (
        (*Dai)["ai"], (*Eai)["ai"]
      );

      // get frobenius norm of errors, which is to be optimized
      double E(Eai->norm2());
      return 0.5*E*E;
    }

    // expects a call of getError with the same grid first
    FrequencyGrid getGradient(const FrequencyGrid &grid) {
      FrequencyGrid gradGrid(grid);
      for (size_t n(0); n < grid.ws.size(); ++n) {
        // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
        // integration weight ws[n] for each ai
        auto DWai(*Dai);
        CTF::Transform<double, double>(
          std::function<void(double, double &)>(
            [&grid, n](double error, double &delta) {
              double diff(error * propagator(delta, grid.nus[n]) / Pi());
              delta = diff;
            }
          )
        ) (
          (*Eai)["ai"], DWai["ai"]
        );
        // sum contributions from all energy gaps eps_a-eps_i
        CTF::Scalar<double> DW;
        DW[""] = DWai["ai"];
        gradGrid.ws[n] = DW.get_val();

        // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
        // integration point nus[n] for each ai
        auto DNuai(*Dai);
        CTF::Transform<double, double>(
          std::function<void(double, double &)>(
            [&grid, n](double error, double &delta) {
              double diff(
                error *
                (-2) * grid.ws[n] * propagator(delta, grid.nus[n], 3) *
                2*grid.nus[n]
              );
              delta = diff / Pi();
            }
          )
        ) (
          (*Eai)["ai"], DNuai["ai"]
        );
        // sum contributions from all energy gaps eps_a-eps_i
        CTF::Scalar<double> DNu;
        DNu[""] = DNuai["ai"];
        gradGrid.nus[n] = DNu.get_val();
      }
      return gradGrid;
    }

    FrequencyGrid grid;

    // Dai = eps_a-eps_i for each ai
    CTF::Tensor<double> *Dai;

    // error of numerical quadrature of current grid wn[n] & nus[n] from
    // analytic value 1/(eps_a-eps_i) for each ai.
    // The current grid is the argument of last getError(grid).
    CTF::Tensor<double> *Eai;
  };
}

void RpaApxEnergy::computeFrequencyGrid() {
  FrequencyGridOptimizer optimizer(
    getIntegerArgument("imaginaryFrequencies", 6),
    *getTensorArgument<>("HoleEigenEnergies"),
    *getTensorArgument<>("ParticleEigenEnergies")
  );
  optimizer.optimize(getIntegerArgument("stepCount", 1024));
}

void RpaApxEnergy::run() {
  computeFrequencyGrid();

  // create tcc infrastructure
  typedef CtfMachineTensor<complex> MachineTensor;
  auto machineTensorFactory(MachineTensor::Factory::create());
  auto tcc(tcc::Tcc<complex>::create(machineTensorFactory));

  // read the Coulomb vertex GammaGqr
  auto GammaFqr(getTensorArgument<complex>("CoulombVertex"));

  // read the Particle/Hole Eigenenergies
  auto epsi(getTensorArgument<>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // get index ranges for particles and holes
  int NF(GammaFqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(GammaFqr->lens[1]);

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int FiaStart[] = {0, iStart,aStart};
  int FiaEnd[]   = {NF,iEnd,  aEnd};
  int FaiStart[] = {0, aStart,iStart};
  int FaiEnd[]   = {NF,aEnd,  iEnd};
  auto ctfGammaFia(GammaFqr->slice(FiaStart, FiaEnd));
  auto ctfGammaFai(GammaFqr->slice(FaiStart, FaiEnd));
  auto ctfConjGammaFia(ctfGammaFia);
  conjugate(ctfConjGammaFia); 
  auto ctfConjGammaFai(ctfGammaFai);
  conjugate(ctfConjGammaFai);
  auto GammaFia(tcc->createTensor(MachineTensor::create(ctfGammaFia)));
  auto GammaFai(tcc->createTensor(MachineTensor::create(ctfGammaFai)));
  auto conjGammaFia(tcc->createTensor(MachineTensor::create(ctfConjGammaFia)));
  auto conjGammaFai(tcc->createTensor(MachineTensor::create(ctfConjGammaFai)));

  auto realNun(getTensorArgument("ImaginaryFrequencyPoints"));
  int Nn(realNun->lens[0]);

  auto realWn(getTensorArgument("ImaginaryFrequencyWeights"));
  CTF::Tensor<complex> complexWn(1, &Nn, *realWn->wrld);
  toComplexTensor(*realWn, complexWn);

  CTF::Tensor<complex> Dain(3, std::vector<int>({Nv,No,Nn}).data());
  Dain["ain"]  = (*epsa)["a"];
  Dain["ain"] -= (*epsi)["i"];
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double nu, complex &d){ d = 1.0 / (d - complex(0,1)*nu); }
    )
  ) (
    (*realNun)["n"], Dain["ain"]
  );
  CTF::Tensor<complex> conjDain(Dain);
  conjugate(conjDain);

  auto Pain(tcc->createTensor(MachineTensor::create(Dain)));
  auto conjPain(tcc->createTensor(MachineTensor::create(conjDain)));

  auto energy(tcc->createTensor(std::vector<int>(), "energy"));
  auto chiVFGn(tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chiV"));
  tcc->compile(
    (
      (*chiVFGn)["FGn"] <<=
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      (*chiVFGn)["FGn"] +=
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"]
    )
  )->execute();

  CTF::Scalar<complex> ctfEnergy;
  ctfEnergy[""] = energy->getMachineTensor<MachineTensor>()->tensor[""];
  LOG(0, "RPA") << std::real(ctfEnergy.get_val());
}

