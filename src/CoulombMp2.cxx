#include <CoulombMp2.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

CoulombMp2::CoulombMp2(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  
}

CoulombMp2::~CoulombMp2() {
}

/**
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
 */
void CoulombMp2::run() {
  TensorData<> *iEpsData(getTensorDataArgument("iEps"));
  TensorData<> *aEpsData(getTensorDataArgument("aEps"));
  TensorData<> *vabijData(getTensorDataArgument("vabij"));
  // allocate
  //int nv(vabijData->value->lens[1]);
  //int no(vabijData->value->lens[3]);
  //int nv(vabijData->value->lens[1]);
  //int no(vabijData->value->lens[3]);
  //int lens[] = { nv, nv, no, no };
  //int syms[] = { NS, NS, NS, NS };
 
  Tensor<> Dabij(vabijData->value);
  Tensor<> Tabij(vabijData->value);
  //Tensor<> Dabij(*vabijData->value);
  //Tensor<> Tabij(*vabijData->value);

  Dabij["abij"] += (*iEpsData->value)["i"];
  Dabij["abij"] += (*iEpsData->value)["j"];
  Dabij["abij"] -= (*aEpsData->value)["a"];
  Dabij["abij"] -= (*aEpsData->value)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
  Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*vabijData->value),"abij", Dabij,"abij", 0.0,"abij", fDivide);

  
  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * (*vabijData->value)["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * (*vabijData->value)["abij"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;
  LOG(0) << "e=" << e << std::endl;

//  int lens[] = { nv, nv, no, no };
//  int syms[] = { NS, NS, NS, NS };
//  vabijData->value = new Tensor<>(4, lens, syms, *Cc4s::world, "Cabij");
 
//  (*vabijData->value)["abij"] =  (*aiCoulombVertexRealData->value)["gai"]*
//                                 (*aiCoulombVertexRealData->value)["gbj"];
// read from tensors: aiCoulombVertexImagData->value
// allocate and write to tensor vabijData->value
//  (*vabijData->value)["i.."] = (*
}

