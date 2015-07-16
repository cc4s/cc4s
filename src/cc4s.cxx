/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <ctf.hpp>
using namespace CTF;

#include "cc4s.hpp"

World *cc4s::world;
int cc4s::rank, cc4s::np, cc4s::no, cc4s::nv, cc4s::nG, cc4s::niter;
CoulombIntegrals *cc4s::V;
Tensor<> *cc4s::Tabij, *cc4s::Tai;
Chi *cc4s::chi;
bool cc4s::profile;


void cc4s::startup() {
  V = new CoulombIntegrals();
  V->calculate();
  // copy Vabij and Vai into the initial amplitudes
  Tabij = new Tensor<>(V->abij);
  Tai = new Tensor<>(V->ai);
}

void cc4s::cleanup() {
}

void cc4s::run() {
  Scalar<> energy(*world);
  double norm;
  // NOTE: should be V->ijab
  energy[""] = Tabij[0]["abij"]*V->abij[0]["abij"];
  if (rank == 0) {
    printf("e=%lf\n", energy.get_val());
  }
  for (int i(0); i < cc4s::niter; ++i) {
    double d = MPI_Wtime();
    iterateAmplitudes();
    // NOTE: should be V->ijab
    energy[""] = Tabij[0]["abij"]*V->abij[0]["abij"];
    norm = Tabij->norm2();
    if (rank == 0) {
      printf("%d: (%d nodes) in time = %lf, |T| = %lf\n",
          i+1, np, MPI_Wtime()-d, norm);
      printf("e=%lf\n", energy.get_val());
    }
  }
}


double divide(double a, double b) {
  return a / b;
}

void cc4s::iterateAmplitudes() {
  Tensor<> T21 = Tensor<>(*Tabij);
  T21["abij"] += .5*Tai[0]["ai"]*Tai[0]["bj"];

  Tensor<> tZabij(*V->abij);
  add_Vxyef_T21efij(tZabij, T21);

  Tensor<> Dabij(4, V->abij->lens, V->abij->sym, *cc4s::world);
  Dabij["abij"] += V->i[0]["i"];
  Dabij["abij"] += V->i[0]["j"];
  Dabij["abij"] -= V->a[0]["a"];
  Dabij["abij"] -= V->a[0]["b"];

  Bivar_Function<> fctr(&divide);
  Tabij->contract(1.0, tZabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
} 

void cc4s::add_Vxyef_T21efij(Tensor<> &Zabij, Tensor<> &T21) {
// for comparison:
//  Zabij["abij"] += .5*V->abcd[0]["abef"]*T21["efij"];
//  return;
  for (int a = 0; a < cc4s::nv; a += cc4s::no) {
    // in case nv is not a multiple of no
    int na(std::min(cc4s::nv-a, cc4s::no));
    for (int b = 0; b < cc4s::nv; b += cc4s::no) {
      // in case nv is not a multiple of no
      int nb(std::min(cc4s::nv-b, cc4s::no));
      int abvv[] = {na,nb,cc4s::nv,cc4s::nv};
      int aboo[] = {na,nb,cc4s::no,cc4s::no};
      // NOTE: the sliced tensors are not symmetrical in the first two indices
      // except for a=b
      // TODO: respect symmetry in a,b
      int sym[] = {NS,NS,NS,NS};
      // slice begin and end of the 
      int Z_begin[] = {a, b, 0, 0};
      int Z_end[] = {a+na, b+nb, cc4s::no, cc4s::no};
      int slicedZ_begin[] = {0, 0, 0, 0};
      int slicedZ_end[] = {na, nb, cc4s::no, cc4s::no};
      // allocate tensor holding a slice of Vabcd
      Tensor<> slicedV(4, abvv, sym, *cc4s::world, "slicedVabcd", 1);
      // fetch or recalculate the entires of this slice
      V->calculate_xycd(slicedV, a, b);
      // allocate temporary tensor holding the slice of Zabij for each a and b
      Tensor<> slicedZabij(4, aboo, sym, *cc4s::world, "slicedZabij", 1);
      slicedZabij["abij"] = slicedV["abef"]*T21["efij"];
      // add half of the sliced Zabij to Zabij at the respective position
      Zabij.slice(
        Z_begin, Z_end, 1.0,
        slicedZabij, slicedZ_begin, slicedZ_end, 0.5
      );
    }
  }
}


char* getCmdOption(char ** begin,
                   char ** end,
                   const   std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end){
    return *itr;
  }
  return 0;
}


int main(int argc, char **argv){
  int const in_num = argc;
  char **input_str = argv;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &cc4s::rank);
  MPI_Comm_size(MPI_COMM_WORLD, &cc4s::np);

  if (getCmdOption(input_str, input_str+in_num, "-nG")){
    cc4s::no = atoi(getCmdOption(input_str, input_str+in_num, "-nG"));
    if (cc4s::nG < 0) cc4s::nG = 10;
  } else cc4s::nG = 10;
  if (getCmdOption(input_str, input_str+in_num, "-no")){
    cc4s::no = atoi(getCmdOption(input_str, input_str+in_num, "-no"));
    if (cc4s::no < 0) cc4s::no = 4;
  } else cc4s::no = 4;
  if (getCmdOption(input_str, input_str+in_num, "-nv")){
    cc4s::nv = atoi(getCmdOption(input_str, input_str+in_num, "-nv"));
    if (cc4s::nv < 0) cc4s::nv = 6;
  } else cc4s::nv = 6;
  if (getCmdOption(input_str, input_str+in_num, "-niter")){
    cc4s::niter = atoi(getCmdOption(input_str, input_str+in_num, "-niter"));
    if (cc4s::niter < 0) cc4s::niter = 1;
  } else cc4s::niter = 1;

  cc4s::profile = false;
  cc4s::world = new World(argc, argv);

  cc4s::startup();
  cc4s::run();
  cc4s::cleanup();

  MPI_Finalize();
  return 0;
}

