#include <algorithms/CcsdEnergyFromCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEnergyFromCoulombIntegrals);

CcsdEnergyFromCoulombIntegrals::CcsdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}

CcsdEnergyFromCoulombIntegrals::~CcsdEnergyFromCoulombIntegrals() {
}

void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
) {
  iterate<double>(i, TaiMixer, TabijMixer);
}

void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
) {
  iterate<complex>(i, TaiMixer, TabijMixer);
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////
template <typename F>
void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<F> *TaiMixer, Mixer<F> *TabijMixer
) {
}

