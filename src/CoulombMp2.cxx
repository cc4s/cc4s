#include <CoulombMp2.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

CoulombMp2::CoulombMp2(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
  
}

CoulombMp2::~CoulombMp2() {
}

/**
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
 */
void CoulombMp2::run() {
  Tensor<> *iEps(getTensorArgument("iEps"));
  Tensor<> *aEps(getTensorArgument("aEps"));
  Tensor<> *vabij(getTensorArgument("vabij"));
 
  Tensor<> Dabij(vabij);
  Tensor<> Tabij(vabij);

  Dabij["abij"] += (*iEps)["i"];
  Dabij["abij"] += (*iEps)["j"];
  Dabij["abij"] -= (*aEps)["a"];
  Dabij["abij"] -= (*aEps)["b"];
  // NOTE: ctf double counts if lhs tensor is SH,SH
  Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*vabij),"abij", Dabij,"abij", 0.0,"abij", fDivide);

  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * (*vabij)["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * (*vabij)["abij"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;
  LOG(0) << "e=" << e << std::endl;

  setRealArgument("mp2Energy", e);
}

