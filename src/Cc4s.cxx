/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Cc4s.hpp"
#include "Exception.hpp"
#include <ctf.hpp>
#include <iostream>
#include <fstream>

using namespace CTF;

Cc4s::Cc4s(
  CTF::World *world_, Options *options_
):
  world(world_),
  options(options_)
{
}

Cc4s::~Cc4s() {
  delete T;
  delete V;
  delete chiReal;
  delete chiImag;
}


void Cc4s::run() {
  // Read from disk
  readFTOD();

  Scalar<> energy(*world);
  double e, dire, exce, norm;
  // NOTE: should be (*V)["ijab"]
  energy[""] = 0.25 * (*T)["abij"]*(*V)["abij"];
  e = energy.get_val();
  if (world->rank == 0) {
    std::cout << "e=" << e << std::endl;
  }
  for (int i(0); i < options->niter; ++i) {
    double d = MPI_Wtime();
    iterateMp2();
    // NOTE: should be (*V)["ijab"]
    energy[""] = (*T)["abij"]*(*V)["abij"];
    dire = energy.get_val();
    energy[""] = (*T)["abji"]*(*V)["abij"];
    exce = -0.5 * energy.get_val();
    e = dire + exce;
    norm = T->abij->norm2();
    if (world->rank == 0) {
      std::cout << i+1 << ": on " << world->np << " node(s) in time " <<
        (MPI_Wtime()-d) << ", |T| = " << norm << std::endl;
      std::cout << "e=" << e << std::endl;
    }
  }
}


double divide(double a, double b) {
  return a / b;
}

void Cc4s::iterateMp2() {
  {
    int syms[] = {NS, NS, NS, NS};
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = 0.5 * Dabij["abij"];

    Bivar_Function<> fctr(&divide);
    T->abij->contract(1.0, *V->abij, "abij", Dabij, "abij", 0.0, "abij", fctr);
  }
}

void Cc4s::iterateCcsd() {
  int no = options->no;
  int nv = options->nv;
  Tensor<> T21 = Tensor<>(T->abij);
  // NOTE: ctf double counts if lhs tensor is AS
  T21["abij"] += 0.5 * (*T)["ai"] * (*T)["bj"];
  Tensor<> tZabij = Tensor<>(V->abij);

  if (!V->abcd) {
    for (int b(0); b < options->nv; b += options->no) {
  //    for (int a(b); a < nv; a += no) {
      for (int a(0); a < nv; a += no) {
        if (world->rank == 0) {
          std::cout << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
        }
        Tensor<> Vxycd(V->getSlice(a, b));
        Vxycd.set_name("Vxycd");
        int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);
        int Tbegin[] = {0, 0, 0, 0};
        int lens[] = {na, nb, no, no};
//        int syms[] = {Vxycd.sym[0], NS, AS, NS};
        int syms[] = {NS, NS, NS, NS};
        Tensor<> Txyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
        Txyij["xyij"] = Vxycd["xyef"] * T21["efij"];

        int tzBegin[] = {a, b, 0, 0};
        int tzEnd[] = {a+na, b+nb, options->no, options->no};
        tZabij.slice(
          tzBegin,tzEnd,1.0, Txyij,Tbegin,lens,0.5
        );
	//TODO: if a!=b remove double counting of ctf, use b>=a loop
      }
    }
  } else {
    tZabij["abij"] += 0.5 * (*V)["abef"] * T21["efij"];
  }

  {
    int syms[] = {SH, NS, SH, NS};
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = 0.5 * Dabij["abij"];

    Bivar_Function<> fctr(&divide);
    T->abij->contract(1.0, tZabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
  }
} 


/**
 * \brief Reads the Fourier transformed overlap densities from disk.
 */
void Cc4s::readFTOD() {
  if (world->rank == 0) {
    std::cout <<
      "Reading Fourier transformed overlap densities from FTODDUMP...";
  }
  std::ifstream file("FTODDUMP");
  std::string line;
  // read a comment line
  std::getline(file, line);
  // NOTE: currently unused
  int nSpins, nk;
  // read the size data
  std::getline(file, line);
  std::stringstream lineStream(line);
  lineStream >> options->no >> options->nv >> options->nG >> nSpins >> nk;

  // allocate chi and Coulomb integral tensors
  chiReal = new Chi(world, options);
  chiImag = new Chi(world, options);
  V = new CoulombIntegrals(chiReal, chiImag, world, options);

  // allocate local indices and values of the chi tensors
  int64_t nG(options->nG); int np(options->no+options->nv);
  // distributed along g: round up g/nprocs 
  int64_t maxValuesCount((nG+world->np-1)/world->np * np*np);
  // each process can have at most maxValuesCount entires
  double *reals(new double[maxValuesCount]);
  double *imags(new double[maxValuesCount]);
  int64_t *indices(new int64_t[maxValuesCount]);
  int64_t valuesCount(0);
  // allocate local indices and values of eigenenergies
  double *iValues(new double[options->no]);
  double *aValues(new double[options->nv]);
  int64_t *iIndices(new int64_t[options->no]);
  int64_t *aIndices(new int64_t[options->nv]);

  // read another comment line
  std::getline(file, line);
  // start reading the file line by line
  while (std::getline(file, line)) {
    std::stringstream lineStream(line);
    double real, imag;
    int g, p, q, spin;
    // parse the line
    lineStream >> real >> imag >> g >> p >> q >> spin;
    if (g > 0) {
      // chi_q^p(g), spin is ignored
      // g, p and q are zero based
      --g; --p; --q;
      // distribed along g: current rank is responsible, only
      if (g % world->np == world->rank) {
        reals[valuesCount] = real;
        imags[valuesCount] = imag;
        indices[valuesCount] = g + nG*(p + np*q);
        ++valuesCount;
      }
    } else {
      // eigenenergy with eps_p, all other indices are to be ignored
      // they are written only on the root
      if (world->rank == 0) {
        --p;
        if (p < options->no) {
          iValues[p] = real;
          iIndices[p] = p;
        } else {
          aValues[p-options->no] = real;
          aIndices[p-options->no] = p-options->no;
        }
      }
    }
  }
  chiReal->gpq->write(valuesCount, indices, reals);
  chiImag->gpq->write(valuesCount, indices, imags);
  int64_t iValuesCount(world->rank == 0 ? options->no : 0);
  int64_t aValuesCount(world->rank == 0 ? options->nv : 0);
  V->i->write(iValuesCount, iIndices, iValues);
  V->a->write(aValuesCount, aIndices, aValues);
  delete[] indices; delete[] reals; delete[] imags;
  delete[] iIndices; delete[] aIndices; delete[] iValues; delete[] aValues;
  file.close();
  if (world->rank == 0) {
    std::cout << " OK" << std::endl;
  }

  double realNorm = chiReal->gpq->norm2();
  double imagNorm = chiImag->gpq->norm2();
  double iNorm = V->i->norm2();
  double aNorm = V->a->norm2();
  if (world->rank == 0) {
    std::cout <<
      "2-Norm of FTOD = (" << realNorm << "," << imagNorm << ")" << std::endl;
    std::cout <<
      "2-Norm of (eps_i,eps_a) = (" << iNorm << "," << aNorm << ")" << std::endl;
  }
  // calculate Coulomb integrals from Fourier transformed overlap densities
  V->fetch();
  // write V(1,1,1,1) for testing
  int64_t readIndices[] = { 0 };
  double readValues[] = { 0.0 };
  V->ijkl->read(1l, readIndices, readValues);
  if (world->rank == 0) {
    std::cout << "V(1,1,1,1) = " << readValues[0] << std::endl;
  }

  // allocate and calculate the intial amplitudes
  T = new Amplitudes(V, world, options);
}


void Cc4s::testSymmetries() {
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


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  try {
    World *world = new World(argumentCount, arguments);
    Options *options = new Options(argumentCount, arguments); 
    Cc4s cc4s(world, options);
//    cc4s.testSymmetries();
    cc4s.run();
  } catch (Exception *cause) {
    std::cout << cause->getMessage() << std::endl;
  }

  MPI_Finalize();
  return 0;
}

