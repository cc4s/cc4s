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
  (*T)["abij"] = 0.0;
  //energy[""] = 0.25 * (*T)["abij"]*(*V)["abij"];
  e = energy.get_val();
  if (world->rank == 0) {
    std::cout << "e=" << e << std::endl;
  }
  for (int i(0); i < options->niter; ++i) {
    double d = MPI_Wtime();
//    iterateMp2();
    iterateRccd();
    // NOTE: should be (*V)["ijab"]
    energy[""] = 2.0 * (*T)["abij"] * (*V)["abij"];
    dire = energy.get_val();
    energy[""] = (*T)["abji"] * (*V)["abij"];
    exce = -1.0 * energy.get_val();
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

void Cc4s::iterateRpa() {
  {
    int syms[] = {NS, NS, NS, NS};
    // Define Tensors
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    //Allocate Tensors for RPA
    Tensor<> Rabij = Tensor<>(V->abij);
    Tensor<> Cabij = Tensor<>(V->abij);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);

    if (world->rank == 0) {
      std::cout << "Solving RPA Amplitude Equations:" << std::endl;
    }


    Rabij["abij"] = (*V)["abij"];
    Rabij["abij"] += 2.0 * (*V)["acik"] * (*T)["cbkj"];
    Cabij["abij"] =  2.0 * (*V)["cbkj"] * (*T)["acik"];
    Rabij["abij"] += Cabij["abij"];
    //(*V)["cbkj"]*(*T)["acjk"];
    Rabij["abij"] += 2.0 * Cabij["acik"] * (*T)["cbkj"];


    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fctr(&divide);
    T->abij->contract(1.0, Rabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
  }
}

void Cc4s::iterateRccd() {
  {
    int syms[] = {NS, NS, NS, NS};
    // Define Tensors
    //Allocate Tensors for T1 amplitude equations
    Tensor<> Dai = Tensor<>(V->ai);
    Tensor<> Rai = Tensor<>(T->ai);
    //Allocate Tensors for T2 amplitudes
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Tensor<> Rabij = Tensor<>(V->abij);
    Tensor<> Fba = Tensor<>(V->ab);
    Tensor<> Fji = Tensor<>(V->ij);
    Tensor<> Fai = Tensor<>(V->ai);
    //intermediates
    Tensor<> Cai = Tensor<>(T->ai);
    Tensor<> Lac = Tensor<>(V->ab);
    Tensor<> Kac = Tensor<>(V->ab);
    Tensor<> Lki = Tensor<>(V->ij);
    Tensor<> Kki = Tensor<>(V->ij);
    Tensor<> Cklij = Tensor<>(V->ijkl);
    Tensor<> Cakic = Tensor<>(V->aijb);
    Tensor<> Cakci = Tensor<>(V->aibj);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);


    if (world->rank == 0) {
      std::cout << "Solving restricted T2 CCD Amplitude Equations:" << std::endl;
    }



//Build Kac
//    Kac["ac"] = Fba["ac"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["adkl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["adkl"];

//Build Lac
    Lac["ac"] = Kac["ac"];

//Build Kki
//    Kki["ki"] = Fji["ki"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["cdil"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["cdil"];

//Build Lki
    Lki["ki"] = Kki["ki"];

//  Contract L_ac with T2 Amplitudes
    Rabij["abij"] = Lac["ac"] * (*T)["cbij"];

//  Contract L_ki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*T)["abkj"];

//Build C_akic
    Cakic["akic"] = (*V)["akic"];
    Cakic["akic"] -= 0.5 * (*V)["lkdc"] * (*T)["dail"];
    Cakic["akic"] += (*V)["lkdc"] * (*T)["adil"];
    Cakic["akic"] -= 0.5 * (*V)["lkcd"] * (*T)["adil"];

//Build C_akci
    Cakci["akci"] = (*V)["akci"];
    Cakci["akci"] -= 0.5 * (*V)["lkcd"] * (*T)["dail"];

//  Contract C_akic and C_akci intermediates with T2 amplitudes

    Rabij["abij"] += 2.0 * Cakic["akic"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakic["akic"] * (*T)["bckj"];

    Rabij["abij"] -= Cakci["akci"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakci["bkci"] * (*T)["ackj"];

// Symmetrize Rabij by applying permutation operator
    // to save memory we use Cakci as intermediate for the permutation operator 
    Cakci["aibj"] = Rabij["abij"];
    Rabij["abij"] += Cakci["bjai"]; 

//////////////////////////////////////////////////////////////////////
// Now add all terms to Rabij that do not need to be symmetrized with
// the permutation operator
//////////////////////////////////////////////////////////////////////

//  Rabij are the Tabij amplitudes for the next iteration and need to be build
    Rabij["abij"] += (*V)["abij"];

//  Build Chi_klij intermediate
    Cklij["klij"] = (*V)["klij"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["cdij"];

//  Contract Chi_klij with T2 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["abkl"];

    if (V->abcd) {
      Tensor<> Cabcd(V->abcd);
  //  Build Chi_abcd intermediate
      Cabcd["abcd"] = (*V)["abcd"];
  //  Contract Chi_abcd with T2 Amplitudes
      Rabij["abij"] += Cabcd["abcd"] * (*T)["cdij"];
    } else {
  // Slicing:
      for (int b(0); b < options->nv; b += options->nw) {
        for (int a(b); a < options->nv; a += options->nw) {
          if (world->rank == 0) {
            std::cout << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          }
          Tensor<> Vxycd(V->getSlice(a, b));
          int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);
          int origin[] = {0, 0, 0, 0};
          int lens[] = {na, nb, options->no, options->no};
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
          Rxyij["xyij"] = Vxycd["xycd"] * (*T)["cdij"];

          int rBegin[] = {a, b, 0, 0};
          int rEnd[] = {a+na, b+nb, options->no, options->no};
          // R["abij"] += R["xyij"] at current x,y
          Rabij.slice(rBegin,rEnd,1.0, Rxyij,origin,lens,1.0);
          if (a>b) {
            // add the same slice at (b,a,j,i):
            rBegin[0] = b; rBegin[1] = a;
            rEnd[0] = b+nb; rEnd[1] = a+na;
            // note that na may be != nb
            lens[0] = nb; lens[1] = na;
            Tensor<> Ryxji(4,lens,syms,*world, "Ryxij", Vxycd.profile);
            Ryxji["yxji"] = Rxyij["xyij"];
            Rabij.slice(rBegin,rEnd,1.0, Ryxji,origin,lens,1.0);
          }
        }
      }
    }

    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fctr(&divide);
    T->abij->contract(1.0, Rabij, "abij", Dabij, "abij", 0.0, "abij", fctr);
  }
}


void Cc4s::iterateRccsd() {
  {
    int syms[] = {NS, NS, NS, NS};
    // Define Tensors
    //Allocate Tensors for T1 amplitude equations
    Tensor<> Rai = Tensor<>(T->ai);
    //Allocate Tensors for T2 amplitudes
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Tensor<> Rabij = Tensor<>(V->abij);
    Tensor<> Fba = Tensor<>(V->ab);
    Tensor<> Fji = Tensor<>(V->ij);
    Tensor<> Fai = Tensor<>(V->ai);
    //intermediates
    Tensor<> Cai = Tensor<>(T->ai);
    Tensor<> Lac = Tensor<>(V->ab);
    Tensor<> Kac = Tensor<>(V->ab);
    Tensor<> Lki = Tensor<>(V->ij);
    Tensor<> Kki = Tensor<>(V->ij);
    Tensor<> Cklij = Tensor<>(V->ijkl);
    Tensor<> Cabcd = Tensor<>(V->abcd);
    Tensor<> Cakic = Tensor<>(V->aijb);
    Tensor<> Cakci = Tensor<>(V->aibj);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);


//********************************************************************************
//  T2 amplitude equations
//********************************************************************************

    if (world->rank == 0) {
      std::cout << "Solving restricted T2 CCSD Amplitude Equations:" << std::endl;
    }



//Build Kac
    Kac["ac"] = Fba["ac"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["adkl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["adkl"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["ak"] * (*T)["dl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["ak"] * (*T)["dl"];

//Build Lac
    Lac["ac"] = Kac["ac"];
    Lac["ac"] -= Fai["ck"] * (*T)["ak"];
    Lac["ac"] += 2.0 * (*V)["akcd"] * (*T)["dk"];
    Lac["ac"] -= (*V)["akdc"] * (*T)["dk"];


//Build Kki
    Kki["ki"] = Fji["ki"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["cdil"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["cdil"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["ci"] * (*T)["dl"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["ci"] * (*T)["dl"];

//Build Lki
    Lki["ki"] = Kki["ki"];
    Lki["ki"] += Fai["ck"] * (*T)["ci"];
    Lki["ki"] += 2.0 * (*V)["klic"] * (*T)["cl"];
    Lki["ki"] -= (*V)["klci"] * (*T)["cl"];

//  Contract L_ac with T2 Amplitudes
    Rabij["abij"] += Lac["ac"] * (*T)["cbij"];

//  Contract L_ki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*T)["abkj"];

//  Contract Coulomb integrals with T2 amplitudes

    Rabij["abij"] += (*V)["abic"] * (*T)["cj"];

    Rabij["abij"] -= (*V)["kbic"] * (*T)["ak"] * (*T)["cj"];

    Rabij["abij"] -= (*V)["akij"] * (*T)["bk"];

    Rabij["abij"] += (*V)["akic"] * (*T)["cj"] * (*T)["bk"];

//Build C_akic
    Cakic["akic"] = (*V)["akic"];
    Cakic["akic"] -= (*V)["lkic"] * (*T)["al"];
    Cakic["akic"] += (*V)["akdc"] * (*T)["di"];
    Cakic["akic"] -= 0.5 * (*V)["lkdc"] * (*T)["dail"];
    Cakic["akic"] -= (*V)["lkdc"] * (*T)["di"] * (*T)["al"];
    Cakic["akic"] += (*V)["lkdc"] * (*T)["adil"];
    Cakic["akic"] -= 0.5 * (*V)["lkcd"] * (*T)["adil"];

//Build C_akci
    Cakci["akci"] = (*V)["akci"];
    Cakci["akci"] -= (*V)["lkci"] * (*T)["al"];
    Cakci["akci"] += (*V)["akcd"] * (*T)["di"];
    Cakci["akci"] -= 0.5 * (*V)["lkcd"] * (*T)["dail"];
    Cakci["akci"] -= (*V)["lkcd"] * (*T)["di"] * (*T)["al"];

//  Contract C_akic and C_akci intermediates with T2 amplitudes

    Rabij["abij"] += 2.0 * Cakic["akic"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakic["akic"] * (*T)["bckj"];

    Rabij["abij"] -= Cakci["akci"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakci["bkci"] * (*T)["ackj"];

// Symmetrize Rabij by applying permutation operator
    // to save memory we use Cakci as intermediate for the permutation operator 
    Cakci["aibj"] = Rabij["abij"];
    Rabij["abij"] += Cakci["bjai"]; 



//////////////////////////////////////////////////////////////////////
// Now add all terms to Rabij that do not need to be symmetrized with
// the permutation operator
//////////////////////////////////////////////////////////////////////

//  Rabij are the Tabij amplitudes for the next iteration and need to be build
    Rabij["abij"] += (*V)["abij"];

//  Build Chi_klij intermediate
    Cklij["klij"] = (*V)["klij"];
    Cklij["klij"] += (*V)["klic"] * (*T)["cj"];
    Cklij["klij"] += (*V)["klcj"] * (*T)["ci"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["cdij"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["ci"] * (*T)["dj"]; 

//  Contract Chi_klij with T2 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["abkl"];

//  Contract Chi_klij with T1 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["ak"] * (*T)["bl"];

//*****************************
//@Felix start slicing here
//*****************************
    if (V->abcd && V->aibc) {
  // copy from Vabcd including data
      Tensor<> Cabcd(V->abcd);
  //  Build Chi_abcd intermediate
      Cabcd["abcd"] -= (*V)["akcd"] * (*T)["bk"];
      Cabcd["abcd"] -= (*V)["kbcd"] * (*T)["ak"];

  //  Contract Chi_abcd with T2 Amplitudes
      Rabij["abij"] += Cabcd["abcd"] * (*T)["cdij"];
      Rabij["abij"] += Cabcd["abcd"] * (*T)["ci"] * (*T)["dj"];
    } else {
    }
//*****************************
//@Felix stop slicing here
//*****************************

    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fctr(&divide);
    T->abij->contract(1.0, Rabij, "abij", Dabij, "abij", 0.0, "abij", fctr);


//********************************************************************************
//  T1 amplitude equations
//********************************************************************************

    if (world->rank == 0) {
      std::cout << "Solving restricted T1 CCSD Amplitude Equations:" << std::endl;
    }

    Rai["ai"] = Fai["ai"];
    
    Rai["ai"] -= 2.0 * Fai["ck"] * (*T)["ak"] * (*T)["ci"];
    
    Rai["ai"] += Kac["ac"] * (*T)["ci"];

    Rai["ai"] -= Kki["ac"] * (*T)["ak"];
  
//  Calculate Kki
  
//    Rai["ai"] += Kki["ac"] * (*T)["ak"];
    
//   (*V)["acik"] * (*T)["cbkj"];


//    Rabij["abij"] += 2.0 * (*V)["acik"] * (*T)["cbkj"];
//    Cabij["abij"] =  2.0 * (*V)["cbkj"] * (*T)["acik"];
//    Rabij["abij"] += Cabij["abij"];
    //(*V)["cbkj"]*(*T)["acjk"];
//    Rabij["abij"] += 2.0 * Cabij["acik"] * (*T)["cbkj"];


  }
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

