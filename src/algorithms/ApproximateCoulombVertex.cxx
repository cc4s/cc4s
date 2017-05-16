#include <algorithms/ApproximateCoulombVertex.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ApproximateCoulombVertex);

ApproximateCoulombVertex::ApproximateCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ApproximateCoulombVertex::~ApproximateCoulombVertex() {
}

void ApproximateCoulombVertex::run() {
  Tensor<complex> *GammaGqr(
    getTensorArgument<complex>("FullCoulombVertex")
  );
  Tensor<complex> *UGF(
    getTensorArgument<complex>("CoulombVertexSingularVectors")
  );
  Tensor<complex> UTGF(*UGF);
  conjugate(UTGF);
  int lens[] = { UGF->lens[1], GammaGqr->lens[1], GammaGqr->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *GammaFqr = new Tensor<complex>(
    3, lens, syms, *GammaGqr->wrld, "GammaFqr"
  );
  allocatedTensorArgument<complex>(
    "CoulombVertex", GammaFqr
  );
  (*GammaFqr)["Fqr"] = (*GammaGqr)["Gqr"] * UTGF["GF"];
}

void ApproximateCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex, DryTensor<complex>>("FullCoulombVertex")
  );
  DryTensor<complex> *UGF(
    getTensorArgument<complex, DryTensor<complex>>(
      "CoulombVertexSingularVectors"
    )
  );
  DryTensor<complex> UTGF(*UGF);
  int lens[] = { UGF->lens[1], GammaGqr->lens[1], GammaGqr->lens[2] };
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *GammaFqr = new DryTensor<complex>(
    3, lens, syms, SOURCE_LOCATION
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "CoulombVertex", GammaFqr
  );
}

