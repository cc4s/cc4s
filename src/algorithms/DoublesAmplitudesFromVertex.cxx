#include <algorithms/DoublesAmplitudesFromVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(DoublesAmplitudesFromVertex);

DoublesAmplitudesFromVertex::DoublesAmplitudesFromVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DoublesAmplitudesFromVertex::~DoublesAmplitudesFromVertex() {
}

void DoublesAmplitudesFromVertex::run() {
  // read the amplitudes vertex YLai
  Tensor<complex> *YLai(getTensorArgument<complex>("DoublesAmplitudesVertex"));

  // get Nv,No
  int Nv(YLai->lens[1]);
  int No(YLai->lens[2]);

  // allocate amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<> *Tabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Tabij"));

  // split YLai into real and imaginary parts
  Tensor<> realYLai(3, YLai->lens, YLai->sym, *YLai->wrld, "realYLai");
  Tensor<> imagYLai(3, YLai->lens, YLai->sym, *YLai->wrld, "imagYLai");
  fromComplexTensor(*YLai, realYLai, imagYLai);

  (*Tabij)["abij"]  = realYLai["Lai"] * realYLai["Lbj"];
  (*Tabij)["abij"] -= imagYLai["Lai"] * imagYLai["Lbj"];

  allocatedTensorArgument("DoublesAmplitudes", Tabij);
}

void DoublesAmplitudesFromVertex::dryRun() {
  // read the Coulomb vertex GammaGqr
  DryTensor<complex> *YLai(
    getTensorArgument<complex, DryTensor<complex>>("DoublesAmplitudesVertex")
  );

  // get Nv,No
  int Nv(YLai->lens[1]);
  int No(YLai->lens[2]);

  // allocate amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> *Tabij(new DryTensor<>(4, vvoo, syms, SOURCE_LOCATION));
  Tabij->use();

  // split YLai into real and imaginary parts
  DryTensor<> realYLai(3, YLai->lens.data(), YLai->syms.data(),SOURCE_LOCATION);
  DryTensor<> imagYLai(3, YLai->lens.data(), YLai->syms.data(),SOURCE_LOCATION);
}

