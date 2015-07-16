/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <ctf.hpp>
using namespace CTF;

// #include "CoulombIntegrals.hpp"


double divide(double a, double b) {
  return a / b;
}

class cc4s {
  public:
    static World *world;
    static int rank, np, no, nv, nG;
    static CoulombIntegrals *V;
    static Chi *chi;
    static profile(false);

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts a total index of a tensor entry into the positions
     * along each dimension.
     */
    // FIXME: cannot handle out of bound indicies
    static void from_index(Tensor<> const &t, int64_t index, int *pos) {
      for (int d(0); d < t.dim; ++d) {
        pos[d] = index % (int64_t)t.lens[d];
        index /= d.lens[d];
      }
    }

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts given positions along each dimension into a total index.
     */
    static int64_t to_index(Tensor<> const &t, int const *pos) {
      int64_t index(0);
      for (int d(d.dim-1); d >= 0; --d) {
        index *= d.lens[d];
        index += pos[d];
      }
      return index;
    }

};

class Chi {
  public:
    Tensor<> *ab, *ai, *ij:

    Chi() {
      // keep the chi tensors in the memory for now
      {
        int lens[] = {cc4s:nG,cc4s:nv, cc4s:nv};
        int syms[] = {NS, SY, NS};
        ab = new Tensor<>(3, lens, syms, *cc4s::world, "Xab", cc4s::profile);
      }
      {
        int lens[] = {cc4s:nG,cc4s:nv, cc4s:no};
        int syms[] = {NS, NS, NS};
        ai = new Tensor<>(3, lens, syms, *cc4s::world, "Xai", cc4s::profile);
      }
      {
        int lens[] = {cc4s:nG,cc4s:no, cc4s:no};
        int smys[] = {NS, NS, NS};
        ij = new Tensor<>(3, lens, smys, *cc4s::world, "Xij", cc4s::profile);
      }
    }
    ~Chi() {
      delete ab, ai, ij;
    }
};

class CoulombIntegrals {
  public:
    // x,y,z,w denotes indices in slices of virtual indices
    Tensor<> *a, *i, *ai, *abij, *xycd;
    // for testing
    Tensor<> *abcd;

    CoulombIntegrals() {
      {
        int lens[] = {cc4s::nv, cc4s::no};
        int smys[] = {NS, NS}
        ai = new Tensor<>(2, lens, syms, *cc4s::world, "Vai", cc4s::profile)
      {
        // sliced indices are of length no:
        int lens[] = {cc4s::no, cc4s::no, cc4s::nv, cc4s::nv};
        int smys[] = {NS, NS, SY, NS};
        xycd = new Tensor<>(4, lens, syms, *cc4s::world, "Vxycd",cc4s::profile);
      }
    }

    void calculate() {
      int64_t j, indicesCount, *indices;
      double *values;
      Tensor<> *tensors[] = {ai};

      for (int i(0); i < sizeof(tensors)/sizeof(*tensors); ++i) {
        tensors[i]->read_local(&indicesCount, &indices, &values);
        for (j(0); j < indicesCount; ++j) {
          values[j] = ((indices[j]*16+i)%13077)/13077. -.5;
        }
        tensors[i]->write(indicesCount, indices, values);
        free(indices); free(values);
      }
    }
};


void fetch_slicedV(Tensor<> &slicedV, int a, int b) {
  int64_t indicesCount, *indices;
  double *values;
  // TODO: make accessible
  int nv = slicedV.lens[2];
  // lengths in the full Vabcd tensor
  int absoluteLengths[4] = {nv, nv, nv, nv};
  int positions[4];
  int64_t absoluteIndex;
  slicedV.read_local(&indicesCount, &indices, &values);
  for (int j(0); j < indicesCount; ++j) {
    // get position within slice 
    from_index(indices[j], 4, slicedV.lens, positions);
    // compute absolute position in Vabcd
    positions[0] += a;
    positions[1] += b;
    // get index of absolute position
    absoluteIndex = to_index(4, absoluteLengths, positions);
    // use it for calculating the respective element
    // NOTE: Vabcd is the array #6:
    values[j] = ((absoluteIndex*16+6)%13077)/13077. -.5;
  }
  slicedV.write(indicesCount, indices, values);
  free(indices), free(values);
}

void add_Vabef_T21efij(Tensor<> &Zabij, Integrals &V, Tensor<> &T21) {
// for comparison:
//  Zabij["abij"] += .5*V["abef"]*T21["efij"];
//  return;
  for (int a = 0; a < V.nv; a += V.no) {
    // in case nv is not a multiple of no
    int na(std::min(V.nv-a, V.no));
    for (int b = 0; b < V.nv; b += V.no) {
      // in case nv is not a multiple of no
      int nb(std::min(V.nv-b, V.no));
      int abvv[] = {na,nb,V.nv,V.nv};
      int aboo[] = {na,nb,V.no,V.no};
      // NOTE: the sliced tensors are not symmetrical in the first two indices
      // except for a=b
      // TODO: respect symmetry in a,b
      int sym[] = {NS,NS,NS,NS};
      // slice begin and end of the 
      int Z_begin[] = {a, b, 0, 0};
      int Z_end[] = {a+na, b+nb, V.no, V.no};
      int slicedZ_begin[] = {0, 0, 0, 0};
      int slicedZ_end[] = {na, nb, V.no, V.no};
      // allocate tensor holding a slice of Vabcd
      Tensor<> slicedV(4, abvv, sym, *V.dw, "slicedVabcd", 1);
      // fetch or recalculate the entires of this slice
      fetch_slicedV(slicedV, a, b);
      // allocate temporary tensor holding the slice of Zabij for each a and b
      Tensor<> slicedZabij(4, aboo, sym, *V.dw, "slicedZabij", 1);
      slicedZabij["abij"] = slicedV["abef"]*T21["efij"];
      // add half of the sliced Zabij to Zabij at the respective position
      Zabij.slice(
        Z_begin, Z_end, 1.0,
        slicedZabij, slicedZ_begin, slicedZ_end, 0.5
      );
    }
  }
}

void ccsd::iterate() {
  Tensor<> T21 = Tensor<>(T.abij);
  T21["abij"] += .5*T.ai["ai"]*T.ai["bj"];

  Tensor<> tZabij(V.abij);
  add_Vabef_T21efij(tZabij, V, T21);

  Tensor<> Dabij(4, V.abij->lens, V.abij->sym, *cc4s::world);
  Dabij["abij"] += V.i["i"];
  Dabij["abij"] += V.i["j"];
  Dabij["abij"] -= V.a["a"];
  Dabij["abij"] -= V.a["b"];

  Bivar_Function<> fctr(&divide);
  T.abij->contract(1.0, tZabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
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
  int rank, np, niter, no, nv, sched_nparts, i;
  int const in_num = argc;
  char **input_str = argv;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &cc4s::rank);
  MPI_Comm_size(MPI_COMM_WORLD, &cc4s::np);

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

  cc4s::world = new World(argc, argv);
    Scalar<> energy(dw);
    {
      Integrals V(no, nv, dw);
      V.fill_rand();
      Amplitudes T(no, nv, dw);
      T.fill_rand();
      Timer_epoch tccsd("CCSD");
      tccsd.begin();
      double norm;
      energy[""] = T["abij"]*V["ijab"];
      if (rank == 0) {
        printf("e=%lf\n", energy.get_val());
      }
      for (i=0; i<niter; i++){
        double d = MPI_Wtime();
        ccsd(V,T,sched_nparts);
        energy[""] = T["abij"]*V["ijab"];
        norm = T.abij->norm2();
        if (rank == 0) {
          printf("%d: (%d nodes) in time = %lf, |T| = %lf\n",
              i+1, np, MPI_Wtime()-d, norm);
          printf("e=%lf\n", energy.get_val());
        }
      }
      tccsd.end();
    }

  MPI_Finalize();
  return 0;
}

