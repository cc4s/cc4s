#include <CtfTest.hpp>
#include <util/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(CtfTest);

CtfTest::CtfTest(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CtfTest::~CtfTest() {
}

/**
 * \brief Calculates Coulomb integrals from aiCoulombVertexReal/Imag
 */
void CtfTest::run() {
  int64_t no(getIntegerArgument("no"));
  int64_t nv(getIntegerArgument("nv"));
  Matrix<> a(no, nv, NS, *Cc4s::world, "Aia");
  Matrix<> b(nv, no, NS, *Cc4s::world, "Bai");
  Matrix<> c(no, no, NS, *Cc4s::world, "Bai");
  setRandomTensor(a);
  setRandomTensor(b);
  LOG(2) << "Starting..." << std::endl;
  c["ij"] = a["ia"] * b["aj"];
  LOG(2) << "DONE." << std::endl;
}

