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
  TensorData<> *iEpsData(getTensorDataArgument("iEps"));
  TensorData<> *aEpsData(getTensorDataArgument("aEps"));
  TensorData<> *vabijData(getTensorDataArgument("vabij"));
  TensorData<> *tabijData(getTensorDataArgument("tabij"));
 
  int nv(vabijData->value->lens[0]);
  int no(vabijData->value->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  tabijData->value = new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij");
 
  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  LOG(0) << "Solving RPA Amplitude Equations:" << std::endl;
 
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration : " << i << std::endl;
    iterateCoulombRpa();
    energy[""] = 2.0 * (*tabijData->value)["abij"] * (*vabijData->value)["abij"];
    dire = energy.get_val();
    energy[""] = (*tabijData->value)["abji"] * (*vabijData->value)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
  }

}

void CoulombRpa::iterateCoulombRpa() {
  // Define Tensors
  TensorData<> *iEpsData(getTensorDataArgument("iEps"));
  TensorData<> *aEpsData(getTensorDataArgument("aEps"));
  TensorData<> *vabijData(getTensorDataArgument("vabij"));
  TensorData<> *tabijData(getTensorDataArgument("tabij"));
  Tensor<> Rabij(*vabijData->value);
  Tensor<> Cabij(*vabijData->value);
  Tensor<> Dabij(*vabijData->value);

  LOG(0) << "Setting Rabij." << std::endl;
  Rabij["abij"] = (*vabijData->value)["abij"];
//  Cabij["abij"] = 2.0 * (*vabijData->value)["acik"] * (*tabijData->value)["cbkj"];
  LOG(0) << "Contracting Vabij with Tabij." << std::endl;
  Cabij["abij"] =  2.0 * (*vabijData->value)["cbkj"] * (*tabijData->value)["acik"];
  LOG(0) << "done." << std::endl;
  Rabij["abij"] += Cabij["abij"];
  Rabij["abij"] += Cabij["baji"];
  LOG(0) << "Contracting Cabij with Tabij." << std::endl;
  Rabij["abij"] += 2.0 * Cabij["acik"] * (*tabijData->value)["cbkj"];
  LOG(0) << "done." << std::endl;

  LOG(0) << "Calculating Dabij." << std::endl;
  Dabij["abij"] += (*iEpsData->value)["i"];
  Dabij["abij"] += (*iEpsData->value)["j"];
  Dabij["abij"] -= (*aEpsData->value)["a"];
  Dabij["abij"] -= (*aEpsData->value)["b"];
  LOG(0) << "done." << std::endl;
  // NOTE: ctf double counts if lhs tensor is SH,SH
  //Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  tabijData->value->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);

}


