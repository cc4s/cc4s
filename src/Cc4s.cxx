/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "cc4s.hpp"
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
  chi = new Chi(world, options);
  V = new CoulombIntegrals(chi);
  T = new Amplitudes(V);
}

Cc4s::~Cc4s() {
  delete T;
  delete V;
  delete chi;
}

void Cc4s::run() {
  // NOTE: should be V->get(IJAB)
  Scalar<> energy(*world);
  double norm;
  energy[""] = T->get(ABIJ)["abij"]*V->get(ABIJ)["abij"];
  if (world->rank == 0) {
    printf("e=%lf\n", energy.get_val());
  }
  for (int i(0); i < niter; ++i) {
    double d = MPI_Wtime();
    iterateAmplitudes();
    // NOTE: should be V->get(IJAB)
    energy[""] = T->get(ABIJ)["abij"]*V->get(ABIJ)["abij"];
    norm = T->get(ABIJ).norm2();
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
  Tensor<> T21 = Tensor<>(T->get(ABIJ));
  T21.set_name("T21");
  T21["abij"] += .5*T->get(AI)["ai"]*T->get(AI)["bj"];

  Tensor<> tZabij(V->get(ABIJ));
  tZabij.set_name("tZabij");

  for (int b(0); b < nv; b += no) {
//    for (int a(b); a < nv; a += no) {
    for (int a(0); a < nv; a += no) {
      Tensor<> Vxycd(V->getSlice(ABCD, a, b));
      Vxycd.set_name("Vxycd");
      int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);

      int Tbegin[] = {0, 0, 0, 0};
      int lens[] = {na, nb, no, no};
//      int syms[] = {Vxycd.sym[0], NS, SY, NS};
      int syms[] = {Vxycd.sym[0], NS, NS, NS};
      Tensor<> Txyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
      Txyij["xyij"] = Vxycd["xycd"] * T21["cdij"];

      int tzBegin[] = {a, b, 0, 0};
      int tzEnd[] = {a+na, b+nb, no, no};
      tZabij.slice(
        tzBegin,tzEnd,1.0, Txyij,Tbegin,lens,0.5
      );
    }
  }

//  tZabij["abij"] += 0.5*V->get(ABCD)["abef"]*T21["efij"];
//  add_Vxyef_T21efij(tZabij, T21);
  Tensor<> Dabij(4, V->get(ABIJ).lens, V->get(ABIJ).sym, *world);
  Dabij["abij"] += V->get(I)["i"];
  Dabij["abij"] += V->get(I)["j"];
  Dabij["abij"] -= V->get(A)["a"];
  Dabij["abij"] -= V->get(A)["b"];

  Bivar_Function<> fctr(&divide);
  T->get(ABIJ).contract(1.0, tZabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
} 

/**
 * \deprecated
 */
void Cc4s::add_Vxyef_T21efij(Tensor<> &Zabij, Tensor<> &T21) {
// for comparison:
  Zabij["abij"] += .5*V->get(ABCD)["abef"]*T21["efij"];
  return;
  for (int a = 0; a < nv; a += no) {
    // in case nv is not a multiple of no
    int na(std::min(nv-a, no));
    for (int b = 0; b < nv; b += no) {
      // in case nv is not a multiple of no
      int nb(std::min(nv-b, no));
      int abvv[] = {na,nb,nv,nv};
      int aboo[] = {na,nb,no,no};
      // NOTE: the sliced tensors are not symmetrical in the first two indices
      // except for a=b
      // TODO: respect symmetry in a,b
      int sym[] = {NS,NS,NS,NS};
      // slice begin and end of the 
      int Z_begin[] = {a, b, 0, 0};
      int Z_end[] = {a+na, b+nb, no, no};
      int slicedZ_begin[] = {0, 0, 0, 0};
      int slicedZ_end[] = {na, nb, no, no};
      // allocate tensor holding a slice of Vabcd
      Tensor<> slicedV(4, abvv, sym, *world, "slicedVabcd", 1);
      // fetch or recalculate the entires of this slice
      V->calculate_xycd(slicedV, a, b);
      // allocate temporary tensor holding the slice of Zabij for each a and b
      Tensor<> slicedZabij(4, aboo, sym, *world, "slicedZabij", 1);
      slicedZabij["abij"] = slicedV["abef"]*T21["efij"];
      // add half of the sliced Zabij to Zabij at the respective position
      Zabij.slice(
        Z_begin, Z_end, 1.0,
        slicedZabij, slicedZ_begin, slicedZ_end, 0.5
      );
    }
  }
}


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  try {
    World *world = new World(argumentCount, arguments);
    Cc4s cc4s(world, Options(argumentCount, arguments));
    cc4s.run();
  } catch (Exception *cause) {
    std::cout << cause->getMessage() << std::endl;
  }

  MPI_Finalize();
  return 0;
}

