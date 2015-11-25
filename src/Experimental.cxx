#include "Experimental.hpp"
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

void Experimental::readRandom(Tensor<> *tensor, int seed) {
  int64_t indicesCount;
  int64_t *indices;
  double *values;
  if (world->rank == 0) {
    std::cout << "Fetching randomly " << tensor->get_name() << "...";
  }
  tensor->read_local(&indicesCount, &indices, &values);
  for (int64_t j(0); j < indicesCount; ++j) {
    values[j] = ((indices[j]*13 + seed)%13077)/13077. -.5;
    values[j] = ((indices[j]*13 + seed+1)%13077)/13077. -.5;
  }
  tensor->write(indicesCount, indices, values);
  free(indices); free(values);
  if (world->rank == 0) std::cout << " OK" << std::endl;
}

void Experimental::testSymmetries() {
  int symsAS[] = {AS, NS};
  int symsNS[] = {NS, NS};
  int lens[] = {3, 3};
  Tensor<> a(2, lens, symsAS, *world, "a");
  Tensor<> n(2, lens, symsNS, *world, "n");
  Tensor<> ns(2, lens, symsNS, *world, "ns");
  double givenValues[] = {
     0.0,  1.0, 2.0,
    -1.0,  0.0, 5.0,
    -2.0, -1.0, 0.0
  };
  double nonsymmetricValues[] = {
     1.0,  2.0, 3.0,
     4.0,  5.0, 6.0,
     7.0,  8.0, 9.0
  };
//  int64_t givenIndices[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  int64_t indicesCount, *indices;
  double *values;
  ns.read_local(&indicesCount, &indices, &values);
  for (int i(0); i < indicesCount; ++i) {
    values[i] = nonsymmetricValues[indices[i]];
  }
  ns.write(indicesCount, indices, values);
  free(indices); free(values);
  n.read_local(&indicesCount, &indices, &values);
  for (int i(0); i < indicesCount; ++i) {
    values[i] = givenValues[indices[i]];
  }
  n.write(indicesCount, indices, values);
  free(indices); free(values);

// test reading in different indices at different ranks:
// works as expected here only for 3 or more processors
//  n.write(3, &givenIndices[world->rank*3], &givenValues[world->rank*3]);
  // BUG: the tensor is internally (anti-)symmetrized by
  // a["ij"] = n["ij"] +(-) n["ji"]
/*
  a["ij"] = 0.5*n["ij"];
  a.read_all(&indicesCount, &values, true);
  for (int i(0); i < indicesCount; ++i) {
    if (world->rank == 0) {
      std::cout << " " << values[i];
      if (i % 3 == 2) std::cout << std::endl;
    }
  }
  free(values);
*/
  ns["ij"] -= ns["ji"];
  a["ij"] = 0.5 * ns["ij"];
  // test AS matrix * NS matrix
  n["ij"] = a["ik"] * n["kj"];
  // works as expected
  n.read_all(&indicesCount, &values, true);
  for (int i(0); i < indicesCount; ++i) {
    if (world->rank == 0) {
      std::cout << " " << values[i];
      if (i % 3 == 2) std::cout << std::endl;
    }
  }
  free(values);  

  // test slicing in (anti-)symmetrical tensor
  int sliceLens[] = {2, 2};
  Tensor<> s(2, sliceLens, symsNS, *world, "s");
  double sliceValues[] = {
    7.0, 11.0,
    0.0, 13.0
  };
  s.read_local(&indicesCount, &indices, &values);
  for (int i(0); i < indicesCount; ++i) {
    values[i] = sliceValues[indices[i]];
  }
  s.write(indicesCount, indices, values);
  free(indices); free(values);
  int aMin[] = {0, 1}, aMax[] = {2, 3};
  int sMin[] = {0, 0}, sMax[] = {2, 2};
  a.slice(aMin, aMax, 1.0, s, sMin, sMax, 1.0);
  a.read_all(&indicesCount, &values, true);
  for (int i(0); i < indicesCount; ++i) {
    if (world->rank == 0) {
      std::cout << " " << values[i];
      if (i % 3 == 2) std::cout << std::endl;
    }
  }
  free(values);  
  // slicing on (anti-)symmetrical tensors works as expected
}
