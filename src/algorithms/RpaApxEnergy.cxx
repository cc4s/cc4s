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
  class FrequencyGridOptimizer {
  public:
    FrequencyGridOptimizer(
      const int N, CTF::Tensor<double> epsi, CTF::Tensor<double> epsa
    ):
      ws(N), nus(N), dws(N), dnus(N)
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
      for (size_t n(0); n < ws.size(); ++n) {
        nus[n] = deltas[No*Nv*n/N + No*Nv/(2*N)];
        ws[n] = nus[n] - lastNu;
        lastNu = nus[n];
      }
    }

    static double propagator(
      const double delta, const double nu, const int e = 2
    ) {
      return 4*delta*delta / std::pow(delta*delta + nu*nu, e);
    }

    void optimize(const int stepCount, const double stepSize) {
      NEW_FILE("L.dat");
      for (int m(0); m < stepCount; ++m) {
        // error of numerical quadrature with current grid wn[n] & nus[n] from
        // analytic value 1/(eps_a-eps_i) for each ai
        auto Eai(*Dai);
        CTF::Transform<double>(
          std::function<void(double &)>(
            [this](double &delta) {
              double error(0);
              // use quadrature to integrate Mp2 propagator
              for (size_t n(0); n < ws.size(); ++n) {
                error += ws[n] * propagator(delta, nus[n]);
              }
              error /= Pi();
              error -= 1 / delta;
              delta = error;
            }
          )
        ) (
          Eai["ai"]
        );

        // get frobenius norm of errors, which is to be optimized
        double L(Eai.norm2()*0.5);
        FILE("L.dat") << m << " " << L;

        for (size_t n(0); n < ws.size(); ++n) {
          // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
          // integration weight ws[n] for each ai
          auto DWai(*Dai);
          CTF::Transform<double, double>(
            std::function<void(double, double &)>(
              [this, n](double error, double &delta) {
                double diff(error * propagator(delta, nus[n]));
                delta = diff / Pi();
              }
            )
          ) (
            Eai["ai"], DWai["ai"]
          );
          // sum contributions from all energy gaps eps_a-eps_i
          CTF::Scalar<double> DW;
          DW[""] = DWai["ai"];
          dws[n] = DW.get_val();

          // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
          // integration point nus[n] for each ai
          auto DNuai(*Dai);
          CTF::Transform<double, double>(
            std::function<void(double, double &)>(
              [this, n](double error, double &delta) {
                double diff(
                  error *
                  (-2) * ws[n] * propagator(delta, nus[n], 3) *
                  2*nus[n]
                );
                delta = diff / Pi();
              }
            )
          ) (
            Eai["ai"], DNuai["ai"]
          );
          // sum contributions from all energy gaps eps_a-eps_i
          CTF::Scalar<double> DNu;
          DNu[""] = DNuai["ai"];
          dnus[n] = DNu.get_val();

          // update grid
          ws[n] -= stepSize * dws[n];
          nus[n] -= stepSize * dnus[n];
        }

        for (size_t n(0); n < ws.size(); ++n) {
          FILE("L.dat") << " " << ws[n];
        }
        for (size_t n(0); n < ws.size(); ++n) {
          FILE("L.dat") << " " << nus[n];
        }
        FILE("L.dat") << std::endl;
      }
    }

    // current frequency integration grid
    std::vector<double> ws;
    std::vector<double> nus;
    // gradient of sum_ai(0.5*error(eps_a-eps_i)^2) w.r.t. ws and nus
    std::vector<double> dws;
    std::vector<double> dnus;

    CTF::Tensor<double> *Dai;
  };
}

void RpaApxEnergy::computeFrequencyGrid() {
  FrequencyGridOptimizer optimizer(
    getIntegerArgument("imaginaryFrequencies", 6),
    *getTensorArgument<>("HoleEigenEnergies"),
    *getTensorArgument<>("ParticleEigenEnergies")
  );
  optimizer.optimize(
    getIntegerArgument("stepCount", 1024),
    getRealArgument("stepSize", 1.0)
  );
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

