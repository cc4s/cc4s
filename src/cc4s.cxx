/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <ctf.hpp>
using namespace CTF;

// #include "CoulombIntegrals.hpp"


double divide(double a, double b) {
  return a / b;
}

class CoulombIntegrals;
class Chi;

class cc4s {
  public:
    static World *world;
    static int rank, np, no, nv, nG, niter;
    static CoulombIntegrals *V;
    static Tensor<> *Tabij, *Tai;
    static Chi *chi;
    static bool profile;

    static void startup();
    static void cleanup();
    static void run();

  protected:
    static void iterateAmplitudes();
    static void add_Vxyef_T21efij(Tensor<> &Zabij, Tensor<> &T21);

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts a total index of a tensor entry into the positions
     * along each dimension.
     */
    // FIXME: cannot handle out of bound indicies
    static void from_index(Tensor<> const &t, int64_t index, int *pos) {
      for (int d(0); d < t.order; ++d) {
        pos[d] = index % (int64_t)t.lens[d];
        index /= t.lens[d];
      }
    }

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts given positions along each dimension into a total index.
     */
    static int64_t to_index(int dim, int const *lens, int const *pos) {
      int64_t index(0);
      for (int d(dim-1); d >= 0; --d) {
        index *= lens[d];
        index += pos[d];
      }
      return index;
    }
  friend class Chi;
  friend class CoulombIntegrals;
};

class Chi {
  public:
    Tensor<> *ab, *ai, *ij;

    Chi() {
      // keep the chi tensors in the memory for now
      {
        int lens[] = {cc4s::nG, cc4s::nv, cc4s::nv};
        int syms[] = {NS, SY, NS};
        ab = new Tensor<>(3, lens, syms, *cc4s::world, "Xab", cc4s::profile);
      }
      {
        int lens[] = {cc4s::nG, cc4s::nv, cc4s::no};
        int syms[] = {NS, NS, NS};
        ai = new Tensor<>(3, lens, syms, *cc4s::world, "Xai", cc4s::profile);
      }
      {
        int lens[] = {cc4s::nG, cc4s::no, cc4s::no};
        int smys[] = {NS, NS, NS};
        ij = new Tensor<>(3, lens, smys, *cc4s::world, "Xij", cc4s::profile);
      }
    }
    ~Chi() {
      delete ab; delete ai; delete ij;
    }
};

class CoulombIntegrals {
  public:
    Tensor<> *a, *i, *ai, *abij;
// NOTE: only for testing
    Tensor<> *abcd;

    CoulombIntegrals() {
      a = new Vector<>(cc4s::nv, *cc4s::world, "Va", cc4s::profile);
      i = new Vector<>(cc4s::no, *cc4s::world, "Vi", cc4s::profile);
      {
        int lens[] = {cc4s::nv, cc4s::no};
        int syms[] = {NS, NS};
        ai = new Tensor<>(2, lens, syms, *cc4s::world, "Vai", cc4s::profile);
      }
      {
        int lens[] = {cc4s::nv, cc4s::nv, cc4s::no, cc4s::no};
        int syms[] = {SY, NS, SY, NS};
        abij = new Tensor<>(4, lens, syms, *cc4s::world, "Vabij",cc4s::profile);
      }
// NOTE: only for testing
      {
        int lens[] = {cc4s::nv, cc4s::nv, cc4s::nv, cc4s::nv};
        int syms[] = {SY, NS, SY, NS};
        abcd = new Tensor<>(4, lens, syms, *cc4s::world, "Vabcd",cc4s::profile);
      }
    }

    void calculate() {
      int64_t indicesCount, *indices;
      double *values;
      Tensor<> *tensors[] = {
        a, i, ai, abij
// NOTE: only for testing
        , abcd
      };

      for (int i(0); i < (int)sizeof(tensors)/(int)sizeof(*tensors); ++i) {
        tensors[i]->read_local(&indicesCount, &indices, &values);
        for (int64_t j(0); j < indicesCount; ++j) {
          values[j] = ((indices[j]*16+i)%13077)/13077. -.5;
        }
        tensors[i]->write(indicesCount, indices, values);
        free(indices); free(values);
      }
    }

    void calculate_xycd(Tensor<> &xycd, int a, int b) {
      int64_t indicesCount, *indices;
      double *values;
      // lengths in the full Vabcd tensor
      int absoluteLengths[] = {cc4s::nv, cc4s::nv, cc4s::nv, cc4s::nv};
      int positions[4];
      int64_t absoluteIndex;
      xycd.read_local(&indicesCount, &indices, &values);
      for (int j(0); j < indicesCount; ++j) {
        // get position within slice 
        cc4s::from_index(xycd, indices[j], positions);
        // compute absolute position in Vabcd
        positions[0] += a;
        positions[1] += b;
        // get index of absolute position
        absoluteIndex = cc4s::to_index(4, absoluteLengths, positions);
        // use it for calculating the respective element
        // NOTE: Vabcd is the array #4:
        values[j] = ((absoluteIndex*16+4)%13077)/13077. -.5;
      }
      xycd.write(indicesCount, indices, values);
      free(indices), free(values);
    }
};

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

void cc4s::iterateAmplitudes() {
  Tensor<> T21 = Tensor<>(Tabij);
  T21["abij"] += .5*Tai[0]["ai"]*Tai[0]["bj"];

  Tensor<> tZabij(V->abij);
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
  Zabij["abij"] += .5*V.abcd["abef"]*T21["efij"];
  return;
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

