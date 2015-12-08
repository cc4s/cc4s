#include <CoulombRpa.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

CoulombRpa::CoulombRpa(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  
}

CoulombRpa::~CoulombRpa() {
}

/**
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
 */
void CoulombRpa::run() {
  Tensor<> *vabij(getTensorArgument("vabij"));

  int nv(vabij->lens[0]);
  int no(vabij->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *tabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("tabij", tabij);
  
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) << "Solving RPA Amplitude Equations:" << std::endl;
 
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration : " << i << std::endl;
    iterateCoulombRpa();
    energy[""] = 2.0 * (*tabij)["abij"] * (*vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*tabij)["abji"] * (*vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
  }

  setRealArgument("rpaEnergy", e);
}

void CoulombRpa::iterateCoulombRpa() {
  // get tensors
  Tensor<> *iEps(getTensorArgument("iEps"));
  Tensor<> *aEps(getTensorArgument("aEps"));
  Tensor<> *vabij(getTensorArgument("vabij"));
  Tensor<> *tabij(getTensorArgument("tabij"));
  Tensor<> Rabij(*vabij);
  Tensor<> Cabij(*vabij);
  Tensor<> Dabij(*vabij);

  Rabij["abij"] = (*vabij)["abij"];
  Rabij["abij"] += 2.0 * (*vabij)["acik"] * (*tabij)["cbkj"];
  Cabij["abij"] =  2.0 * (*vabij)["cbkj"] * (*tabij)["acik"];
  Rabij["abij"] += Cabij["abij"];
  Rabij["abij"] += 2.0 * Cabij["acik"] * (*tabij)["cbkj"];

  Dabij["abij"] += (*iEps)["i"];
  Dabij["abij"] += (*iEps)["j"];
  Dabij["abij"] -= (*aEps)["a"];
  Dabij["abij"] -= (*aEps)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
  Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  tabij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
}

