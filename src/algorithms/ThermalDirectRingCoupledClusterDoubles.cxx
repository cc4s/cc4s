#include <algorithms/ThermalDirectRingCoupledClusterDoubles.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalDirectRingCoupledClusterDoubles);

ThermalDirectRingCoupledClusterDoubles::ThermalDirectRingCoupledClusterDoubles(
  std::vector<Argument> const &argumentList
): ThermalClusterDoublesAlgorithm(argumentList) {
}

ThermalDirectRingCoupledClusterDoubles::
  ~ThermalDirectRingCoupledClusterDoubles(
) {
}

void ThermalDirectRingCoupledClusterDoubles::update(int n) {
  // Read Vabij
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // use double buffering for amplitudes
  Tensor<> *thisTabij(Tabij[n&1]);
  Tensor<> *nextTabij(Tabij[(n^1)&1]);

  // get the occupancies for contractions
  Tensor<> *Ni(getTensorArgument("ThermalHoleOccupancies"));
  Tensor<> *Na(getTensorArgument("ThermalParticleOccupancies"));

  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(),
    abbreviation.begin(), ::toupper
  );


  // I. Start with contribution where the two particle hole pairs are either
  // both incomming (the closing and the quadratic case)
  // or both outgoing (the Coulomb case).
  Tensor<> PVabij(*Vabij);
  PVabij.set_name("PVabij");
  // propgate the connected states (i.e. integrated over the interval)
  PPHHImaginaryTimePropagation pphhPropagation(beta / samples);
  Transform<>(std::function<void(double, double, double &)>(pphhPropagation)) (
    (*Dai)["ai"], (*Dai)["bj"], PVabij["abij"]
  );

  // I.A. close the current amplitudes and compute the contributions
  // to the updated energy
  // Note that contractions need to involve the thermal occupancies
  (*directEnergy)[""] += 2.0 * (*thisTabij)["abij"] *
    (*Na)["a"] * (*Na)["b"] * (*Ni)["i"] * (*Ni)["j"] * PVabij["abij"];
  // TODO: think of how to factor that...
  (*exchangeEnergy)[""] += -1.0 * (*thisTabij)["abji"] *
    (*Na)["a"] * (*Na)["b"] * (*Ni)["i"] * (*Ni)["j"] * PVabij["abij"];


  // II. calculate all contributions to the next amplitudes having
  // exactly one occurrance of interaction H_1

  // II.A Start with the contributions that can use the same propagated
  // particle hole pairs PVabij

  // II.A.1 The Coulomb contribution
  (*nextTabij)["abij"] = PVabij["cdkl"];

  // II.B The other terms require the current amplitudes where one particle
  // hole pair propgates freely from the start to the end of the interval.
  // Propagate the left pair accordingly.
  FreePHImaginaryTimePropagation phPropagation(beta / samples);
  Transform<>(std::function<void(double, double &)>(phPropagation)) (
    (*Dai)["ai"], (*thisTabij)["abij"]
  );

  // II.B.1 The quadratic contribution. Note that the left pair of the
  // amplitudes is prepated for free propagation.
  (*nextTabij)["abij"] +=
    (*thisTabij)["acik"] *
      (*Na)["c"] * (*Ni)["k"] * PVabij["cdkl"] * (*Na)["d"] * (*Ni)["l"] *
      (*thisTabij)["bdjl"];

  // II.C Now, add contributions where the one particle hole propgates to,
  // the other from the Coulomb interaction. Wlog, the left one propgates to it.
  PVabij["abij"] = (*Vabij)["abij"];
  HPPHImaginaryTimePropagation hpphPropagation(beta / samples);
  Transform<>(std::function<void(double, double, double &)>(hpphPropagation)) (
    (*Dai)["ai"], (*Dai)["bj"], PVabij["abij"]
  );

  // II.C.1 The linear term with V contracted on the right side
  (*nextTabij)["abij"] +=
    (*thisTabij)["acik"] * (*Na)["c"] * (*Ni)["k"] * PVabij["cdkl"];

  // II.C.2 the linear term with V contracted on the left side.
  // Note that we need to reverse the order in PVabij and thisTabij
  // since their left side is prepared for contraction or free propagation,
  // respectively.
  (*nextTabij)["abij"] +=
      PVabij["dali"] * (*Na)["d"] * (*Ni)["l"] * (*thisTabij)["bdjl"];


  // III. finally, also propagate the right particle hole pair of the current
  // amplitudes freely
  Transform<>(std::function<void(double, double &)>(phPropagation)) (
    (*Dai)["bj"], (*thisTabij)["abij"]
  );
  // and add it to the nextTabij. These are the contributions without an
  // occurrance of the interaction H_1 within the current interval
  (*nextTabij)["abij"] += (*thisTabij)["abij"];
}


void ThermalDirectRingCoupledClusterDoubles::dryUpdate() {
  // Read the DRCCD amplitudes Tabij
  //DryTensor<> *Tabij(
  getTensorArgument<double, DryTensor<double>>("DrccdDoublesAmplitudes");
  //);

  // Read the Coulomb Integrals Vabij
  DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));

  // Allocate Tensors for T2 amplitudes
  DryTensor<> Rabij(*Vabij);

  // Define intermediates
  DryTensor<> Cabij(*Vabij);

  DryTensor<> Dabij(*Vabij);
}

