/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Cc4s.hpp"
#include "Exception.hpp"
#include <ctf.hpp>
#include <iostream>

using namespace CTF;

Cc4s::Cc4s(
  CTF::World *world_, Options const &options
):
  world(world_),
  no(options.no),
  nv(options.nv),
  nG(options.nG),
  niter(options.niter)
{
  chiReal = new Chi(world, options, 0);
  chiImag = new Chi(world, options, 1);
  V = new CoulombIntegrals(chiReal, chiImag);
  T = new Amplitudes(V);
}

Cc4s::~Cc4s() {
  delete T;
  delete V;
  delete chiReal;
  delete chiImag;
}

void Cc4s::run() {
  // NOTE: should be (*V)["ijab"]
  Scalar<> energy(*world);
  double norm;
  std::cout << "doing ..." << std::endl;
  energy[""] =
    2.0*(*T)["abij"]*(*V)["abij"] -
    (*T)["abij"]*(*V)["abji"];
  if (world->rank == 0) {
    printf("e=%lf\n", energy.get_val());
  }
  for (int i(0); i < niter; ++i) {
    double d = MPI_Wtime();
    iterateAmplitudes();
    // NOTE: should be (*V)["ijab"]
    energy[""] =
      2.0*(*T)["abij"]*(*V)["abij"] -
      (*T)["abij"]*(*V)["abji"];
    norm = T->abij->norm2();
    if (world->rank == 0) {
      printf("%d: (%d nodes) in time = %lf, |T| = %lf\n",
          i+1, world->np, MPI_Wtime()-d, norm);
      printf("e=%lf\n", energy.get_val());
    }
  }
}


double divide(double a, double b) {
  return a / b;
}

void Cc4s::iterateAmplitudes() {
  int syms[] = {NS, NS, NS, NS};
  Tensor<> T21 = Tensor<>(4, T->abij->lens, syms, *world, "T21");
  T21["abij"] = (*T)["abij"] + .5*(*T)["ai"]*(*T)["bj"];

  Tensor<> tZabij = Tensor<>(4, T->abij->lens, syms, *world, "tZabij");
  tZabij["abij"] = (*V)["abij"];

/*
  for (int b(0); b < nv; b += no) {
//    for (int a(b); a < nv; a += no) {
    for (int a(0); a < nv; a += no) {
      Tensor<> Vxycd(V->getSlice(a, b));
      Vxycd.set_name("Vxycd");
      int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);

      int Tbegin[] = {0, 0, 0, 0};
      int lens[] = {na, nb, no, no};
//      int syms[] = {Vxycd.sym[0], NS, SY, NS};
      int syms[] = {NS, NS, NS, NS};
      Tensor<> Txyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
      Txyij["xyij"] = Vxycd["xyef"] * T21["efij"];

      int tzBegin[] = {a, b, 0, 0};
      int tzEnd[] = {a+na, b+nb, no, no};
      tZabij.slice(
        tzBegin,tzEnd,1.0, Txyij,Tbegin,lens,0.5
      );
    }
  }
*/
  tZabij["abij"] += 0.5 * (*V)["abef"] * T21["efij"];

//  int syms[] = {SY, NS, SY, NS};
  Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
  Dabij["abij"] += (*V)["i"];
  Dabij["abij"] += (*V)["j"];
  Dabij["abij"] -= (*V)["a"];
  Dabij["abij"] -= (*V)["b"];

  Bivar_Function<> fctr(&divide);
  T->abij->contract(1.0, tZabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
} 

void Cc4s::testSymmetries() {
  int symsAS[] = {AS, NS};
  int symsNS[] = {NS, NS};
  int lens[] = {3, 3};
  Tensor<> a(2, lens, symsAS, *world, "a");
  Tensor<> n(2, lens, symsNS, *world, "n");
  double givenValues[] = {
     0.0,  1.0, 2.0,
    -1.0,  0.0, 5.0,
    -2.0, -1.0, 0.0
  };
  int64_t givenIndices[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  int64_t indicesCount, *indices;
  double *values;
/*
  n.read_local(&indicesCount, &indices, &values);
  for (int i(0); i < indicesCount; ++i) {
    values[i] = givenValues[indices[i]];
  }
  n.write(indicesCount, indices, values);
  free(indices); free(values);
*/
// test reading in different indices at different ranks:
// works as expected here only for 3 or more processors
  n.write(3, &givenIndices[world->rank*3], &givenValues[world->rank*3]);
  // BUG: the tensor is internally (anti-)symmetrized by
  // a["ij"] = n["ij"] +(-) n["ji"]
  a["ij"] = 0.5*n["ij"];
  a.read_all(&indicesCount, &values, true);
  for (int i(0); i < indicesCount; ++i) {
    if (world->rank == 0) {
      std::cout << " " << values[i];
      if (i % 3 == 2) std::cout << std::endl;
    }
  }
  free(values);
  n["ij"] = a["ij"];
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


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  try {
    World *world = new World(argumentCount, arguments);
    Cc4s cc4s(world, Options(argumentCount, arguments));
//    cc4s.testSymmetries();
    cc4s.run();
  } catch (Exception *cause) {
    std::cout << cause->getMessage() << std::endl;
  }

  MPI_Finalize();
  return 0;
}

