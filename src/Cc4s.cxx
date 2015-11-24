/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

// TODO: change "" includes to <> includes
#include "Cc4s.hpp"
#include "util/Log.hpp"
#include "TextFtodReader.hpp"
#include "BinaryFtodReader.hpp"
#include "Exception.hpp"
#include "FtodRankDecomposition.hpp"
#include "CrossEntropyFtodRankDecomposition.hpp"
#include "util/CubicPolynomialRootFinder.hpp"
#include "util/ComplexPolynomialRootFinder.hpp"
#include "util/MathFunctions.hpp"
#include "util/IterativePseudoInverter.hpp"
#include <ctf.hpp>
#include <iostream>
#include <fstream>

using namespace CTF;

Cc4s::Cc4s(): flopCounter() {
}

Cc4s::~Cc4s() {
}

void Cc4s::run() {
  printBanner();

  // Read from disk
//  TextFtodReader textFtodReader;
  BinaryFtodReader binaryFtodReader;
//  textFtodReader.read();
  binaryFtodReader.read();
//  binaryFtodReader.write();

  IterativePseudoInverter<double>::test(world);
  IterativePseudoInverter<complex>::test(world);

  // experimental:
  std::vector<Argument const *> arguments;
  TensorData chiRData("chiR", *chiReal->gpq);
  InputArgument chiR("chiR", &chiRData);
  arguments.push_back(&chiR);
  TensorData chiIData("chiI", *chiImag->gpq);
  InputArgument chiI("chiI", &chiIData);
  arguments.push_back(&chiI);
  IntegerData rankData("rank", options->rank);
  InputArgument rank("rank", &rankData);
  arguments.push_back(&rank);
  RealData epsilonData("epsilon", options->accuracy);
  InputArgument epsilon("epsilon", &epsilonData);
  arguments.push_back(&epsilon);
//  CrossEntropyFtodRankDecomposition ftodRankDecomposition(arguments);
  FtodRankDecomposition ftodRankDecomposition(arguments);
//  util::CubicPolynomialRootFinder::test();
//  util::ComplexPolynomialRootFinder::test();
//  return;
  ftodRankDecomposition.run();
  (*chiReal->gpq)["Gqr"] = (*ftodRankDecomposition.chi0R)["Gqr"];
  (*chiImag->gpq)["Gqr"] = (*ftodRankDecomposition.chi0I)["Gqr"];

  // calculate Coulomb integrals from Fourier transformed overlap densities
  Cc4s::V->fetch();
  // write V(1,1,1,1) for testing
  int64_t readIndices[] = { 0 };
  double readValues[] = { 0.0 };
  Cc4s::V->ijkl->read(1l, readIndices, readValues);
  if (Cc4s::world->rank == 0) {
    std::cout << "V(1,1,1,1) = " << readValues[0] << std::endl;
  }

  // allocate and calculate the intial amplitudes
  Cc4s::T = new Amplitudes(Cc4s::V);

  Scalar<> energy(*world);
  double e, dire, exce, norm;
  (*T)["abij"] = 0.0;
  for (int i(0); i < options->niter; ++i) {
    double d = MPI_Wtime();
//    iterateMp2();
    iterateRccd();
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
  printStatistics();
}

void Cc4s::printBanner() {
  if (Cc4s::world->rank == 0) {
    std::cout << "====== Coupled Cluster for Solids ======" << std::endl;
    std::cout << "version " << CC4S_VERSION << " " << CC4S_DATE << std::endl;
    std::cout << "built " << __DATE__ << " " << __TIME__ <<
      " with c++ " << __cplusplus << std::endl;
  }
}

void Cc4s::printStatistics() {
  int64_t flops = flopCounter.count();
  std::string pid, comm, state, ppid, pgrp, session, ttyNr,
    tpgid, flags, minflt, cminflt, majflt, cmajflt,
    utime, stime, cutime, cstime, priority, nice,
    O, itrealvalue, starttime;
  int64_t vsize, rss;
  // assuming LINUX 
  std::ifstream statStream("/proc/self/stat", std::ios_base::in);
  statStream >> pid >> comm >> state >> ppid >> pgrp >> session >> ttyNr
    >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
    >> utime >> stime >> cutime >> cstime >> priority >> nice
    >> O >> itrealvalue >> starttime >> vsize >> rss;
  statStream.close();
  // in case x86-64 is configured to use 2MB pages
  int64_t pageSize = sysconf(_SC_PAGE_SIZE);
  if (world->rank == 0) {
    std::cout << "performance statistics:" << std::endl;
    std::cout << "  on root: " << flops / 1.e9 << " GFLOPS" << std::endl;
    std::cout << "    physical memory: " <<
      rss * pageSize / 1e9 << " GB" << std::endl;
    std::cout << "    virtual  memory: " <<
      vsize / 1e9 << " GB" << std::endl;
  }

  flops = flopCounter.count(world->comm);
  int64_t globalVSize, globalRss;
  MPI_Reduce(&vsize, &globalVSize, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  MPI_Reduce(&rss, &globalRss, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  if (world->rank == 0) {
    std::cout << "  total:   " << flops / 1.e9 << " GFLOPS" << std::endl;
    std::cout << "    physical memory: " <<
      globalRss * pageSize / 1e9 << " GB" << std::endl;
    std::cout << "    virtual  memory: " <<
      globalVSize / 1e9 << " GB" << std::endl;
  }
  // TODO: timing
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

    Bivar_Function<> fDivide(&MathFunctions::divide<>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
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
      for (int b(0); b < chiReal->nv; b += options->nw) {
        for (int a(b); a < chiReal->nv; a += options->nw) {
          if (world->rank == 0) {
            std::cout << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          }
          Tensor<> Vxycd(V->getSlice(a, b));
          int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);
          int origin[] = {0, 0, 0, 0};
          int lens[] = {na, nb, chiReal->no, chiReal->no};
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
          Rxyij["xyij"] = Vxycd["xycd"] * (*T)["cdij"];

          int rBegin[] = {a, b, 0, 0};
          int rEnd[] = {a+na, b+nb, chiReal->no, chiReal->no};
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

    Bivar_Function<> fDivide(&MathFunctions::divide<>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
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

    Bivar_Function<> fDivide(&MathFunctions::divide<>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);


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

    Bivar_Function<> fDivide(&MathFunctions::divide<>);
    T->abij->contract(1.0, *V->abij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
}

void Cc4s::iterateCcsd() {
  int no = chiReal->no;
  int nv = chiReal->nv;
  Tensor<> T21 = Tensor<>(T->abij);
  // NOTE: ctf double counts if lhs tensor is AS
  T21["abij"] += 0.5 * (*T)["ai"] * (*T)["bj"];
  Tensor<> tZabij = Tensor<>(V->abij);

  if (!V->abcd) {
    for (int b(0); b < chiReal->nv; b += chiReal->no) {
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
        int tzEnd[] = {a+na, b+nb, chiReal->no, chiReal->no};
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

    Bivar_Function<> fDivide(&MathFunctions::divide<>);
    T->abij->contract(1.0, tZabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
} 


World *Cc4s::world;
Options *Cc4s::options;
Chi *Cc4s::chiReal, *Cc4s::chiImag;
CoulombIntegrals *Cc4s::V;
Amplitudes *Cc4s::T;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  try {
    Cc4s::world = new World(argumentCount, arguments);
    Cc4s::options = new Options(argumentCount, arguments);
    Log::logLevel = Cc4s::world->rank == 0 ? Cc4s::options->logLevel : -1;
    Cc4s cc4s;
    cc4s.run();
  } catch (DetailedException *cause) {
    std::cout << std::endl << cause->getMessage() << std::endl;
  }

  MPI_Finalize();
  return 0;
}

