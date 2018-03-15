#include <algorithms/RealStructure.hpp>

#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(RealStructure);

RealStructure::RealStructure(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

RealStructure::~RealStructure() {
}

void RealStructure::run() {

  Tensor<> *momenta(getTensorArgument<>("Momenta"));
  Tensor<> *StructureFactor(getTensorArgument<>("StructureFactor"));

  cartesianMomenta = NEW(cc4s::Vector<>)
  int NG(momenta->lens[0]);
  momenta->read_all(&cartesianMomenta[0][0]);

  cc4s::Vector<> check_grid;
  for (int g(0); g<NG; ++g){
    check_grid+=cartisianMomenta[g];
  }
  if (check_grid.length() > 1e-10){
    LOG(1,"RealStructure") << "Alles gute.";
  }

}


