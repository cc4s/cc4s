#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/MultiCombinations.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalMp2EnergyFromCoulombIntegrals);

ThermalMp2EnergyFromCoulombIntegrals::ThermalMp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalMp2EnergyFromCoulombIntegrals::~ThermalMp2EnergyFromCoulombIntegrals() {
}

void ThermalMp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
//  Tensor<> *Ni(getTensorArgument("ThermalHoleOccupancies"));
//  Tensor<> *Na(getTensorArgument("ThermalParticleOccupancies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  beta = 1/getRealArgument("Temperature");

  // compute \Delta^{ab}_{ij} = eps_a + eps_b - eps_i - eps_j
  Dabij = NEW(Tensor<>, false, *Vabij);
  (*Dabij)["abij"] =  (*epsa)["a"];
  (*Dabij)["abij"] += (*epsa)["b"];
  (*Dabij)["abij"] -= (*epsi)["i"];
  (*Dabij)["abij"] -= (*epsi)["j"];

/* tests:
  LOG(1, "FT-MP2") <<
    MultiCombinations(5,3).multinomial(std::vector<unsigned int>({2,0,1,0,0}))
    << std::endl;
  testDLogZMp2(0);
  testDLogZMp2(1);
  testDLogZMp2(2);
  testDLogZMp2(3);
  testDLogZHf(0);
  testDLogZHf(1);
  testDLogZHf(2);
  testDLogZHf(3);
  testDLogZH0(0);
  testDLogZH0(1);
  testDLogZH0(2);
  testDLogZHf(0, D_MU);
  testDLogZHf(1, D_MU);
  testDLogZMp2(0, D_MU);
  testDLogZMp2(1, D_MU);
*/

  // TODO: compute Helmholtz free energy from
  // F0 + Omega_1 + Omega_2 - 1/2 (dOmega_1/dmu)^2 / (d^2Omega_0/dmu^2)
  // [Fetta, Walecka, Eq.(30.69)]

  real Omega0(-getDLogZH0(0, D_MU)/beta);
  real Omega1(-getDLogZHf(0, D_MU)/beta);
  real Omega2(-getDLogZMp2(0, D_MU)/beta);

  real N0_0(-getDLogZH0(1, D_MU)/beta);
  real N1_0(-getDLogZHf(1, D_MU)/beta);
  real N2_0(-getDLogZMp2(1, D_MU)/beta);

  real N0_1(-getDLogZH0(2, D_MU)/beta);
  real N1_1(-getDLogZHf(2, D_MU)/beta);

  real N0_2(-getDLogZH0(3, D_MU)/beta);

  real mu1(-N1_0/N0_1);
  real mu2(-(N2_0+N1_1*mu1+0.5*N0_2*mu1*mu1)/N0_1);

  real Fc(Omega1 + Omega2 - N1_0*mu1 - 0.5*N0_1*mu1*mu1);
  EMIT() << YAML::Key << "Omega0" << YAML::Value << Omega0;
  EMIT() << YAML::Key << "Omega1" << YAML::Value << Omega1;
  EMIT() << YAML::Key << "Omega2" << YAML::Value << Omega2;
  EMIT() << YAML::Key << "N0_0" << YAML::Value << N0_0;
  EMIT() << YAML::Key << "N1_0" << YAML::Value << N1_0;
  EMIT() << YAML::Key << "N2_0" << YAML::Value << N2_0;
  EMIT() << YAML::Key << "N0_1" << YAML::Value << N0_1;
  EMIT() << YAML::Key << "N1_1" << YAML::Value << N1_1;
  EMIT() << YAML::Key << "N0_2" << YAML::Value << N0_2;
  EMIT() << YAML::Key << "mu1" << YAML::Value << mu1;
  EMIT() << YAML::Key << "mu2" << YAML::Value << mu2;

  setRealArgument("ThermalFreeEnergy", Fc);
}

void ThermalMp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("ThermalPPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}


/**
 * \brief Computes the nth derivative of the MP2 contribution to log(Z)
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * w.r.t. (-beta)
 **/
cc4s::real ThermalMp2EnergyFromCoulombIntegrals::getDLogZMp2(
  const unsigned int n, const bool dbeta
) {
  Tensor<> Tabij(false, *getTensorArgument("ThermalPPHHCoulombIntegrals"));
  if (n == 0) {
    // all derivative degrees are zero
    addLogZMp2Amplitudes(Tabij, std::vector<unsigned int>(5));
  } else {
    if (dbeta) {
      LOG(2, "FT-MP2") << "computing MP2 contribution of d^"
        << n << " log(Z(beta)) / (-dbeta)^" << n << " ..." << std::endl;
      // go through all possiblities to distribute n derivatives among 5 terms
      for (auto &derivedTerms: MultiCombinations(5, n)) {
        // convert to degree of derivative for each term
        std::vector<unsigned int> degrees(5);
        for (auto term: derivedTerms) {
          ++degrees[term];
        };
        real multiplicity( MultiCombinations(5, n).multinomial(degrees) );
        addLogZMp2Amplitudes(Tabij, degrees, dbeta, multiplicity);
      }
    } else {
      LOG(2, "FT-MP2") << "computing MP2 contribution of d^"
        << n << " log(Z(beta)) / (beta*dmu)^" << n << " ..." << std::endl;
      // go through all possiblities to distribute n derivatives among 4 terms
      for (auto &derivedTerms: MultiCombinations(4, n)) {
        // convert to degree of derivative for each term
        std::vector<unsigned int> derivedDegrees(4);
        for (auto term: derivedTerms) {
          ++derivedDegrees[term];
        };
        // the first term will not be derived
        std::vector<unsigned int> degrees(1);
        degrees.insert(
          degrees.end(), derivedDegrees.begin(), derivedDegrees.end()
        );
        real multiplicity( MultiCombinations(4,n).multinomial(derivedDegrees) );
        addLogZMp2Amplitudes(Tabij, degrees, dbeta, multiplicity);
      }
    }
  }

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  energy[""] = (+0.5) * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real direct( energy.get_val() );
  energy[""] = (-0.5) * spins * Tabij["abij"] * (*Vabij)["abji"];
  real exchange( energy.get_val() );
  writeContribution("MP2 direct", n, dbeta, direct);
  writeContribution("MP2 exchange", n, dbeta, exchange);
  return direct+exchange;
}

/**
 * \brief Adds the contribution of the thermal MP2 amplitudes
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * to the tensor Tabij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 5 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addLogZMp2Amplitudes(
  Tensor<> &Tabij,
  const std::vector<unsigned int> &degrees, const bool dbeta,
  const real multiplicity
) {
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  // start with Vabij
  Tensor<> tabij(*Vabij);
  // apply the derivative of the 5 terms of the respective degrees:
  // Tabij *= f^a = 1/(1+exp(-eps_a*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[1], dbeta)
    )
  ) (
    (*epsa)["a"], tabij["abij"]
  );
  // Tabij *= f^b
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[2], dbeta)
    )
  ) (
    (*epsa)["b"], tabij["abij"]
  );
  // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[3], dbeta)
    )
  ) (
    (*epsi)["i"], tabij["abij"]
  );
  // Tabij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[4], dbeta)
    )
  ) (
    (*epsi)["j"], tabij["abij"]
  );
  // Tabij *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalMp2Propagation<>(beta, degrees[0])
    )
  ) (
    (*Dabij)["abij"], tabij["abij"]
  );
  // add this contribution
  Tabij["abij"] += multiplicity * tabij["abij"];
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ (-1) \int_0^\beta{\rm d}\tau_1\,f_i\, f_j \f$
 * and returns the result in the tensor Tij.
 **/
cc4s::real ThermalMp2EnergyFromCoulombIntegrals::getDLogZHf(
  const unsigned int n, const bool dbeta
) {
  Tensor<> Tij(2, getTensorArgument("ThermalHHHHCoulombIntegrals")->lens);
  if (n == 0) {
    // all derivative degrees are zero
    addLogZHfAmplitudes(Tij, std::vector<unsigned int>(3));
  } else {
    if (dbeta) {
      LOG(2, "FT-MP2") << "computing HF contribution of d^"
        << n << " log(Z(beta)) / d(-beta)^" << n << " ..." << std::endl;
      // go through all possiblities to distribute n derivatives among 3 terms
      for (auto &derivedTerms: MultiCombinations(3, n)) {
        // convert to degree of derivative for each term
        std::vector<unsigned int> degrees(3);
        for (auto term: derivedTerms) {
          ++degrees[term];
        };
        real multiplicity( MultiCombinations(3, n).multinomial(degrees) );
        addLogZHfAmplitudes(Tij, degrees, dbeta, multiplicity);
      }
    } else {
      LOG(2, "FT-MP2") << "computing HF contribution of d^"
        << n << " log(Z(beta)) / (beta*dmu)^" << n << " ..." << std::endl;
      // go through all possiblities to distribute n derivatives among 2 terms
      for (auto &derivedTerms: MultiCombinations(2, n)) {
        // convert to degree of derivative for each term
        std::vector<unsigned int> derivedDegrees(2);
        for (auto term: derivedTerms) {
          ++derivedDegrees[term];
        };
        // the first term will not be derived
        std::vector<unsigned int> degrees(1);
        degrees.insert(
          degrees.end(), derivedDegrees.begin(), derivedDegrees.end()
        );
        real multiplicity( MultiCombinations(2,n).multinomial(derivedDegrees) );
        addLogZHfAmplitudes(Tij, degrees, dbeta, multiplicity);
      }
    }
  }

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  // we compute -Veff + HF = -2*HF + HF = -HF, so the sign is negative
// TODO: what about Hartree term?
  energy[""] = (-1.0) * (+0.5) * spins * spins * Tij["ij"] * (*Vijkl)["ijij"];
  real direct( energy.get_val() );
  energy[""] = (-1.0) * (-0.5) * spins * Tij["ij"] * (*Vijkl)["ijji"];
  real exchange( energy.get_val() );

  writeContribution("Hartree", n, dbeta, direct);
  writeContribution("Exchange", n, dbeta, exchange);
  return direct+exchange;
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ (-1) \int_0^\beta{\rm d}\tau_1\,f_i\, f_j \f$
 * and returns the result in the tensor Tij.
 **/
cc4s::real ThermalMp2EnergyFromCoulombIntegrals::getDLogZHfDEps(
) {
  Tensor<> Tij(2, getTensorArgument("ThermalHHHHCoulombIntegrals")->lens);
  Tensor<> *ni(getTensorArgument("ThermalHoleOccupancies"));
  Tensor<> *na(getTensorArgument("ThermalParticleOccupancies"));
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Tensor<> *depsdmu(new Tensor<>(*ni));
  // FIXME: assumes No=Nv=Np
  (*depsdmu)["k"] =  beta * spins * (*Vijkl)["kjkj"] * (*ni)["j"] * (*na)["j"];
  (*depsdmu)["k"] -= beta * (*Vijkl)["kjjk"] * (*ni)["j"] * (*na)["j"];
  // write deps/dmu
  allocatedTensorArgument("ThermalDEpsDMu", depsdmu);

  Scalar<> energy;
  // we compute -Veff + HF = -2*HF + HF = -HF, so the sign is negative
  energy[""] = (+beta) * spins * spins * (*depsdmu)["j"] * (*Vijkl)["ijij"]
    * (*na)["j"] * (*ni)["i"] * (*ni)["j"];
  real direct( energy.get_val() );
  energy[""] = (-beta) * spins * (*depsdmu)["j"] * (*Vijkl)["ijji"]
    * (*na)["j"] * (*ni)["i"] * (*ni)["j"];
  real exchange( energy.get_val() );

  writeContribution("Hartree deps/dmu", 1, D_MU, direct);
  writeContribution("Exchange deps/dmu", 1, D_MU, exchange);
  return direct+exchange;
}


/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ \int_0^\beta{\rm d}\tau\, f_i\, f_j \f$
 * to the tensor Tij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 3 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addLogZHfAmplitudes(
  Tensor<> &Tij,
  const std::vector<unsigned int> &degrees,
  const bool dbeta,
  const real multiplicity
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> tij(false, Tij);
  // start with 1

  // the first term is (-1) int_0^beta dtau = -beta
  switch (degrees[0]) {
  case 0:
    // degrees[0] will always be zero for d/dmu
    tij["ij"] = -beta;
    break;
  case 1:
    tij["ij"] = dbeta ? 1 : 0;
    break;
  default:
    // all higher derivatives of -beta w.r.t. (-beta) are zero
    // add nothing to Tabij and return
    return;
  }
  // apply the derivative of the other 2 terms of the respective degrees:
  // Tij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[1], dbeta)
    )
  ) (
    (*epsi)["i"], tij["ij"]
  );
  // Tij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[2], dbeta)
    )
  ) (
    (*epsi)["j"], tij["ij"]
  );
  // add this contribution
  Tij["ij"] += multiplicity * tij["ij"];
}

cc4s::real ThermalMp2EnergyFromCoulombIntegrals::getDLogZH0(
  const unsigned int n, const bool dbeta
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> Ti(false, *epsi);
  Ti["i"] = 1;
  switch (n) {
  case 0:
    Transform<real, real>(
      std::function<void(real, real &)>(
        [this](real eps, real &T) {
          T = std::log(1 + std::exp(-beta*eps));
        }
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    break;
  default:
    // the first derivative is eps_i*f_i (dbeta) or f_i (dmu),
    // use the ThermalContraction clas providing arbitrary derivatives of f_i
    Transform<real, real>(
      std::function<void(real, real &)>(
        ThermalContraction<>(beta, false, n-1, dbeta)
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    if (dbeta) {
      Ti["i"] *= (*epsi)["i"];
    }
    break;
  }
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> dLogZH0;
  // sum over all relevant (thermal hole) states
  dLogZH0[""] = spins * Ti["i"];
  writeContribution("single body reference", n, dbeta, dLogZH0.get_val());
  return dLogZH0.get_val();
}

void ThermalMp2EnergyFromCoulombIntegrals::writeContribution(
  const std::string &contribution,
  const unsigned int n, const bool dbeta, const real dLogZ
) {
  std::stringstream term;
  real m(dLogZ);
  switch (n) {
  case 0:
    term << "grand potential from " << contribution;
    // convert log(Z) to free energy for writing into log
    m *= -getRealArgument("Temperature");
    break;
  default:
    term << n << ".central moment of " << (dbeta ? "H" : "N") << " from "
      << contribution;
    break;
  }
  LOG(0, "FT-MP2") << term.str() << "=" << m << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDLogZMp2(
  const unsigned int n, const bool dbeta
) {
  real dLogZ, exactDLogZ;
  std::string d;
  if (dbeta) {
    real b(beta), db(0.00001);
    beta = b+db;
    dLogZ = getDLogZMp2(n);
    beta = b-db;
    dLogZ -= getDLogZMp2(n);
    beta = b;
    dLogZ /= -2*db;
    exactDLogZ = getDLogZMp2(n+1);
    d = "-dbeta";
  } else {
    Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
    real dmu( 0.00001 );
    (*epsi)["i"] -= dmu; (*epsa)["a"] -= dmu;
    dLogZ = getDLogZMp2(n, dbeta);
    (*epsi)["i"] -= -2.0*dmu; (*epsa)["a"] -= -2.0*dmu;
    dLogZ -= getDLogZMp2(n, dbeta);
    (*epsi)["i"] -= dmu; (*epsa)["a"] -= dmu;
    dLogZ /= 2*beta*dmu;
    exactDLogZ = getDLogZMp2(n+1, dbeta);
    d = "beta*dmu";
  }
  LOG(1, "FT-MP2") << "MP2 numerical(d^" << (n+1) << " log(Z(beta)) / (" << d
    << ")^"  << (n+1) << " = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "MP2     exact(d^" << (n+1) << " log(Z(beta)) / (" << d
    << ")^"  << (n+1) << " = " << exactDLogZ << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDLogZHf(
  const unsigned int n, const bool dbeta
) {
  real dLogZ, exactDLogZ;
  std::string d;
  if (dbeta) {
    real b(beta), db(0.00001);
    beta = b+db;
    dLogZ = getDLogZHf(n, dbeta);
    beta = b-db;
    dLogZ -= getDLogZHf(n, dbeta);
    beta = b;
    dLogZ /= -2*db;
    exactDLogZ = getDLogZHf(n+1, dbeta);
    d = "-dbeta";
  } else {
    Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
    real dmu( 0.00001 );
    (*epsi)["i"] -= dmu; (*epsa)["a"] -= dmu;
    dLogZ = getDLogZHf(n, dbeta);
    (*epsi)["i"] -= -2.0*dmu; (*epsa)["a"] -= -2.0*dmu;
    dLogZ -= getDLogZHf(n, dbeta);
    (*epsi)["i"] -= dmu; (*epsa)["a"] -= dmu;
    dLogZ /= 2*beta*dmu;
    exactDLogZ = getDLogZHf(n+1, dbeta);
    d = "beta*dmu";
  }
  LOG(1, "FT-MP2") << "HF numerical(d^" << (n+1) << " log(Z(beta)) / (" << d
    << ")^"  << (n+1) << " = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "HF     exact(d^" << (n+1) << " log(Z(beta)) / (" << d
    << ")^"  << (n+1) << " = " << exactDLogZ << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDLogZH0(
  const unsigned int n, const bool dbeta
) {
  real b(beta), db(0.00001);
  beta = b+db;
  real dLogZ( getDLogZH0(n) );
  beta = b-db;
  dLogZ -= getDLogZH0(n);
  beta = b;
  dLogZ /= -2*db;
  real exactDLogZ( getDLogZH0(n+1) );
  LOG(1, "FT-MP2") << "H0 numerical(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "H0     exact(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << exactDLogZ << std::endl;
}

