/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/

// compute Helmholtz free energy up to second order in the Coulomb perturbation
// according to
// F0 + Omega_1 + Omega_2 - 1/2 (dOmega_1/dmu)^2 / (d^2Omega_0/dmu^2)
// at mu=mu0 where mu0 is the Hartree--Fock chemical potential,
// following
// [Kohn, Luttinger, PR (1960), Eq.(18)] and [Fetta, Walecka, Eq.(30.69)]

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
  Dai = NEW(Tensor<>, 2, &Vabij->lens[1]);
  (*Dai)["ai"] =  (*epsa)["a"];
  (*Dai)["ai"] -= (*epsi)["i"];

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

  if (isArgumentGiven("chemicalPotentialShift")) {
    deltaMu = getRealArgument("chemicalPotentialShift");
    shiftedChemicalPotential();
  } else {
    deltaMu = 0.0;
    expandedChemicalPotential();
  }
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

void ThermalMp2EnergyFromCoulombIntegrals::shiftedChemicalPotential() {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;

  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> Nc(1, epsa->lens);
  Tensor<> Nk(1, epsi->lens);
  Tensor<> nk(1, epsi->lens);
  // particles in perturbation: contraction with shifted chemical potential
  Nc["c"] = 1.0;
  // Nc *= f^c = 1/(1+exp(-(eps_c-deltaMu)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, true)
    )
  ) (
    (*epsa)["c"], Nc["c"]
  );
  // holes in perturbation: contraction with shifted chemical potential
  Nk["k"] = 1.0;
  // Nk *= f_k = 1/(1+exp(+(eps_k-deltaMu)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, false)
    )
  ) (
    (*epsi)["k"], Nk["k"]
  );
  // terms in effective potential: contraction with Hartree--Fock chemical pot.
  nk["k"] = 1.0;
  // nk *= f_k = 1/(1+exp(+(eps_k-0)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, 0.0, false)
    )
  ) (
    (*epsi)["k"], nk["k"]
  );

  // zeroth order
  real Omega0(-getDLogZH0(0, D_BETA)/beta);

  // first order
  // Hartree and exchange term, use shifted occupancies Nk
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  energy[""] = (+0.5) * spins * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijij"];
  real ED1( energy.get_val() );
  energy[""] = (-0.5) * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijji"];
  real EX1( energy.get_val() );
  // minus effective potential, use Hartree--Fock occupancies nk
  energy[""] = (-1.0) * spins * spins * nk["i"] * nk["j"] * (*Vijkl)["ijij"];
  real EE1( energy.get_val() );
  energy[""] = (+1.0) * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijji"];

  // second order:
  Tensor<> *Vaijk(getTensorArgument("ThermalPHHHCoulombIntegrals"));
  // operator for singles
  Tensor<> Fai(2, Vaijk->lens);
  // contraction weight = weight of perturation - weight of effective pot.
  Nk["k"] -= nk["k"];
  // setup singles operator with difference of perturbation - effective pot.
  Fai["ai"] =  (+1.0) * spins * (*Vaijk)["akik"] * Nk["k"];
  Fai["ai"] += (-1.0) * (*Vaijk)["akki"] * Nk["k"];
  // restore normal weight of pertrubation holes for further calculations
  Nk["k"] += nk["k"];
  // singles:
  // start with Fai
  Tensor<> Tai(Fai);
  // Tai *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalMp2Propagation<>(beta, 0))
  ) (
    (*Dai)["ai"], Tai["ai"]
  );
  // no symmetry, one loop, one hole contracted: +1.0 * spins
  energy[""] = (+1.0) * spins * Tai["ai"] * Fai["ai"] * Nk["i"] * Nc["c"];
  real ES2( -energy.get_val()/beta );

  // doubles:  
  // start with Vabij
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> Tabij(false, *Vabij);
  Tabij["abij"] = (*Vabij)["abij"] * Nc["a"] * Nc["b"] * Nk["i"] * Nk["j"];
  // Tabij *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalMp2Propagation<>(beta, 0))
  ) (
    (*Dabij)["abij"], Tabij["abij"]
  );
  energy[""] = (+0.5) * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real ED2( -energy.get_val()/beta );
  energy[""] = (-0.5) * spins * Tabij["abij"] * (*Vabij)["abji"];
  real EX2( -energy.get_val()/beta );

  real FHf(Omega0+ED1+EX1+EE1);
  real Fc(ES2+ED2+EX2);
  EMIT() << YAML::Key << "Omega0" << YAML::Value << Omega0;
  EMIT() << YAML::Key << "D1" << YAML::Value << ED1;
  EMIT() << YAML::Key << "X1" << YAML::Value << EX1;
  EMIT() << YAML::Key << "eff1" << YAML::Value << EE1;
  EMIT() << YAML::Key << "S2" << YAML::Value << ES2;
  EMIT() << YAML::Key << "D2" << YAML::Value << ED2;
  EMIT() << YAML::Key << "X2" << YAML::Value << EX2;
  EMIT() << YAML::Key << "Hartree-Fock-free-energy" << YAML::Value << FHf;
  EMIT() << YAML::Key << "correlation-free-energy" << YAML::Value << Fc;
}

void ThermalMp2EnergyFromCoulombIntegrals::expandedChemicalPotential() {
  // expand Omega = Omega0 + Omega1 + Omega2 + ... in orders of perturbation
  real Omega0(-getDLogZH0(0, D_MU)/beta);
// FIXME: Omega1(mu) = -eff + HF + XG = -1*(HF+XG)(mu0) + 0.5*(HF+XG)(mu)
// i.e. Omega1(mu0) = -0.5*(HF+XG)(mu0) but dOmega1/dmu(mu0) has opposite sign
  real Omega1(-getDLogZHf(0, D_MU)/beta);
  real Omega2(-getDLogZMp2(0, D_MU)/beta);

  // expand N = N0 + N1 + N2 + ... in orders of perturbation
  real N0_0(getDLogZH0(1, D_MU));
// FIXME: Omega1(mu) = -eff + HF + XG = -1*(HF+XG)(mu0) + 0.5*(HF+XG)(mu)
// i.e. Omega1(mu0) = -0.5*(HF+XG)(mu0) but dOmega1/dmu(mu0) has opposite sign
  real N1_0(-getDLogZHf(1, D_MU));
  real N2_0(getDLogZMp2(1, D_MU));

  // get derivatives of N0, N1
  real N0_1(getDLogZH0(2, D_MU)*beta);
// FIXME: Omega1(mu) = -eff + HF + XG = -1*(HF+XG)(mu0) + 0.5*(HF+XG)(mu)
// i.e. Omega1(mu0) = -0.5*(HF+XG)(mu0) but dOmega1/dmu(mu0) has opposite sign
  real N1_1(-getDLogZHf(2, D_MU)*beta);

  real N0_2(getDLogZH0(3, D_MU)*beta*beta);

  // expand mu = mu0 + mu1 + mu2 + ... in orders of perturbation
  real mu0(getRealArgument("ChemicalPotential"));
  real mu1(0), mu2(0);
  if (std::abs(N0_1) > 1e-7) {
    // determine mu1 from 1.order serires expansion of N around mu0 including
    // only terms of 0th or 1st order in perturbation:
    // N_fixed = N0 + N1 + dN0/dmu * (mu-mu0) with mu=mu0+m1=mu1, thus
    mu1 = -N1_0/N0_1;
    // similarly, determine mu2 from 2.order series expansion of N around mu0
    // including only terms through 2nd order in perturbation:
    // N_fixed = N0+N1+N2 + d(N0+N1)/dmu *(mu-mu0) + 1/2 d^2N0/dmu^2 (mu-mu0)^2
    // leading to
    mu2 = -(N2_0+N1_1*mu1+0.5*N0_2*mu1*mu1) / N0_1;
  }

  // use series expansion of Omega around Omega0 with above terms
  // and dOmegai/dmu = -Ni then
  // write Helmholtz free energy F = Omega(mu) + (mu0+mu1+mu2+...)*N_fixed
  // and split into Hartree--Fock and correlation contribution
  real FHf(Omega0 + mu0*N0_0 + Omega1);
  real Fc(Omega2 - N1_0*mu1 - 0.5*N0_1*mu1*mu1);

  // TODO: implement mixed derivatives e.g. d^2logZ/(dmu dbeta)
  // to evaluate E around mu0:
  // E(mu) = E0+E1+E2 + d(E0+E1)/dmu (mu-mu0) + 1/2 d^2E0/dmu^2 *(mu-mu0)^2

  EMIT() << YAML::Key << "Omega0" << YAML::Value << Omega0;
  EMIT() << YAML::Key << "Omega1" << YAML::Value << Omega1;
  EMIT() << YAML::Key << "Omega2" << YAML::Value << Omega2;
  EMIT() << YAML::Key << "N0_0" << YAML::Value << N0_0;
  EMIT() << YAML::Key << "N1_0" << YAML::Value << N1_0;
  EMIT() << YAML::Key << "N2_0" << YAML::Value << N2_0;
  EMIT() << YAML::Key << "N0_1" << YAML::Value << N0_1;
  EMIT() << YAML::Key << "N1_1" << YAML::Value << N1_1;
  EMIT() << YAML::Key << "N0_2" << YAML::Value << N0_2;
  EMIT() << YAML::Key << "mu0" << YAML::Value << mu0;
  EMIT() << YAML::Key << "mu1" << YAML::Value << mu1;
  EMIT() << YAML::Key << "mu2" << YAML::Value << mu2;
  EMIT() << YAML::Key << "Hartree-Fock-free-energy" << YAML::Value << FHf;
  EMIT() << YAML::Key << "correlation-free-energy" << YAML::Value << Fc;

  setRealArgument("ThermalFreeEnergy", FHf+Fc);
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
  const std::vector<unsigned int> &degrees,
  const bool dbeta,
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
      ThermalContraction<>(beta, deltaMu, true, degrees[1], dbeta)
    )
  ) (
    (*epsa)["a"], tabij["abij"]
  );
  // Tabij *= f^b
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, true, degrees[2], dbeta)
    )
  ) (
    (*epsa)["b"], tabij["abij"]
  );
  // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, false, degrees[3], dbeta)
    )
  ) (
    (*epsi)["i"], tabij["abij"]
  );
  // Tabij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, false, degrees[4], dbeta)
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
  Tensor<> *VHijkl(Vijkl);
  if (isArgumentGiven("ThermalHartreeCoulombIntegrals")) {
    VHijkl = getTensorArgument("ThermalHartreeCoulombIntegrals");
  }
  // we compute -Veff + HF = -2*HF + HF = -HF, so the sign is negative
  energy[""] = (-1.0) * (+0.5) * spins * spins * Tij["ij"] * (*VHijkl)["ijij"];
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
      ThermalContraction<>(beta, deltaMu, false, degrees[1], dbeta)
    )
  ) (
    (*epsi)["i"], tij["ij"]
  );
  // Tij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, false, degrees[2], dbeta)
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
          T = std::log(1 + std::exp(-beta*(eps-deltaMu)));
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
        ThermalContraction<>(beta, deltaMu, false, n-1, dbeta)
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

