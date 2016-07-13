#include <algorithms/ReduceCoulombVertex.hpp>
#include <util/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ReduceCoulombVertex);

ReduceCoulombVertex::ReduceCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ReduceCoulombVertex::~ReduceCoulombVertex() {
}

void ReduceCoulombVertex::run() {
  Tensor<complex> *GammaGqr(
    getTensorArgument<complex>("CoulombVertex")
  );
  Tensor<complex> *UGg(
    getTensorArgument<complex>("EnergyMatrixTransform")
  );
  int lens[] = { UGg->lens[1], GammaGqr->lens[1], GammaGqr->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *Gammagqr = new Tensor<complex>(
    3, lens, syms, *GammaGqr->wrld, "Gammagqr"
  );
  allocatedTensorArgument<complex>(
    "ReducedCoulombVertex", Gammagqr
  );
  (*Gammagqr)["gqr"] = (*GammaGqr)["Gqr"] * (*UGg)["Gg"];
}

void ReduceCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex, DryTensor<complex>>("CoulombVertex")
  );
  DryTensor<complex> *UGg(
    getTensorArgument<complex, DryTensor<complex>>("EnergyMatrixTransform")
  );
  int lens[] = { UGg->lens[1], GammaGqr->lens[1], GammaGqr->lens[2] };
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *Gammagqr = new DryTensor<complex>(3, lens, syms);
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ReducedCoulombVertex", Gammagqr
  );
}

