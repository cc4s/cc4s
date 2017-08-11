#include <algorithms/Mp2EquationOfMotion.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EquationOfMotion);

// TODO: Study the requirements to treat the problem with real numbers or
// complex

Mp2EquationOfMotion::Mp2EquationOfMotion(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2EquationOfMotion::~Mp2EquationOfMotion() {
}

void Mp2EquationOfMotion::run() {
  typedef CTF::Tensor<> T;

  // Get orbital energies
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int totalDimension(1 + Nv * No + No * No * Nv * Nv);
  LOG(1, "MP2_EOM") << "Nv " << Nv << std::endl;
  LOG(1, "MP2_EOM") << "No " << No << std::endl;
  LOG(1, "MP2_EOM") << "Problem dimension " << totalDimension << std::endl;

  // Get couloumb integrals (these shoul not be antisymmetrized)
  T *Vijkl(getTensorArgument<double, T>("HHHHCoulombIntegrals"));
  T *Vabcd(getTensorArgument<double, T>("PPPPCoulombIntegrals")); 

  T *Vabij(getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  // T *Vijab(getTensorArgument<double, T>("HHPPCoulombIntegrals")); // swap PPHH (done)

  T *Vijka(getTensorArgument<double, T>("HHHPCoulombIntegrals"));
  // T *Viajk(getTensorArgument<double, T>("HPHHCoulombIntegrals")); // swap HHHP (done)

  T *Vaibj(getTensorArgument<double, T>("PHPHCoulombIntegrals")); // not in eqs
  //T *Viajb(getTensorArgument<double, T>("HPHPCoulombIntegrals")); // swap PHPH (done)

  T *Vabci(getTensorArgument<double, T>("PPPHCoulombIntegrals")); // not in eqs
  //T *Viabc(getTensorArgument<double, T>("HPPPCoulombIntegrals")); // swap PPPH (done)
  //T *Vabic(getTensorArgument<double, T>("PPHPCoulombIntegrals")); // swap PPPH (done)

  int syms[] = {NS, NS, NS, NS};

  LOG(1, "MP2_EOM") << "Antisymmetrizing Vpqrs " << std::endl;

  //  Vijab
  int oovv[] = { No, No, Nv, Nv };
  T *Vijab(
    new T(4, oovv, syms, *Cc4s::world, "Vijab")
  );
  (*Vijab)["ijab"] =  (*Vabij)["abij"] - (*Vabij)["abji"];

  //  Viajk
  int ovoo[] = { No, Nv, No, No };
  T *Viajk(
    new T(4, ovoo, syms, *Cc4s::world, "Viajk")
  );
  (*Viajk)["iajk"] =  (*Vijka)["ijka"]  - (*Vijka)["ikja"];

  // Viajb
  int ovov[] = { No, Nv, No, Nv };
  T *Viajb(
    new T(4, ovov, syms, *Cc4s::world, "Viajb")
  );
  (*Viajb)["iajb"] =  (*Vaibj)["aibj"] - (*Vaibj)["aijb"];

  // Viabc
  int ovvv[] = { No, Nv, Nv, Nv };
  T *Viabc(
    new T(4, ovvv, syms, *Cc4s::world, "Viabc")
  );
  (*Viabc)["iabc"] =  (*Vabci)["abci"] - (*Vabci)["acbi"];

  // Vabic
  int vvov[] = { Nv, Nv, No, Nv };
  T *Vabic(
    new T(4, vvov, syms, *Cc4s::world, "Vabic")
  );
  (*Vabic)["abic"] =  (*Vabci)["abci"] - (*Vabci)["abic"];

  // Antisymmetrize integrals that are read in
  (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
  (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
  (*Vabij)["abij"] -= (*Vabij)["abji"];
  (*Vijka)["ijka"] -= (*Vijka)["ijak"];
  (*Vaibj)["aibj"] -= (*Vaibj)["aijb"];
  (*Vabci)["abci"] -= (*Vabci)["abic"];

  T Tabij(false, Vabij);
  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  LOG(1, "MP2_EOM") <<
    "Creating doubles amplitudes" << totalDimension << std::endl;
  CTF::Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*Vabij),"abij", Tabij,"abij", 0.0,"abij", fDivide);

  CTF::Scalar<> energy(0.0);
  double energy_val(0.0);

  LOG(2, "MP2_EOM") << "Calculating MP2 energy" << std::endl;
  energy[""] = ( 0.25 ) * Tabij["abij"] * (*Vabij)["abij"];
  energy_val = energy.get_val();
  LOG(1, "MP2_EOM") << " Mp2 energy = " << energy_val << std::endl;


  // Create L and R and intermediates
  int oneBodySyms[] = {NS, NS};
  int oneBodyLensL[] = {No, Nv};
  T *Lia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "Lia") );
  T Lijab(false, Vijab);
  T *LHia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "LHia") );
  T LHijab(false, Vijab);
  int oneBodyLensR[] = {Nv, No};
  T *Rai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "Rai") );
  T Rabij(false, Vabij);
  T *HRai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "Rai") );
  T HRabij(false, Vabij);

  // kinetic terms
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  T *Fab( new T(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab") );
  int kineticLensOccupied[] = {No, No};
  T *Fij( new T(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij") );

  (*Fab)["aa"] = (*epsa)["a"];
  (*Fij)["ii"] = (*epsi)["i"];

  //The totalDimension should be totalDimension, but the zero-particle part is
  //zero, so we restrict the hamiltonian to the space spanned by the singles
  //and doubles excitations
  int hLens[] = {totalDimension-1, totalDimension-1};
  int hSyms[] = {NS, NS};
  T *Hpq( new CTF::Tensor<>(2, hLens, hSyms, *Cc4s::world, "Hpq") );

  int64_t *hIndices;
  double *hValues;
  int oneBodyLength((*Lia).lens[0] * (*Lia).lens[1]);
  int twoBodyLength(
      Rabij.lens[0] * Rabij.lens[1] * Rabij.lens[2] *  Rabij.lens[3]
  );

  LOG(1,"Right EOM Vectors Lengths") << " Lengths= "
    << oneBodyLength << " " << twoBodyLength << std::endl;


  for (int64_t j = 0 ; j < totalDimension-1; j++) {
    getCanonicalPerturbationBasis(*Rai, Rabij, j);

    // Contruct HR (one body part)
    (*HRai)["bi"]  = 0.0;
    (*HRai)["bi"] += ( - 1.0 ) * (*Fij)["ki"] * (*Rai)["bk"];
    (*HRai)["bi"] += ( + 1.0 ) * (*Fab)["bc"] * (*Rai)["ci"];
    (*HRai)["bi"] += ( - 1.0 ) * (*Viajb)["kbid"] * (*Rai)["dk"];
    (*HRai)["bi"] += ( + 0.5 ) * (*Vijka)["klie"] * Rabij["ebkl"];
    (*HRai)["bi"] += ( + 0.5 ) * (*Viabc)["kbde"] * Rabij["deki"];
    (*HRai)["bi"] += ( + 1.0 ) * Tabij["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
    (*HRai)["bi"] += ( - 0.5 ) * Tabij["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
    (*HRai)["bi"] += ( - 0.5 ) * Tabij["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];

    // Contruct HR (two body part)
    HRabij["cdij"]  = 0.0;
    HRabij["cdij"] += ( - 1.0 ) * (*Viajk)["mdij"] * (*Rai)["cm"];
    HRabij["cdij"] += ( + 1.0 ) * (*Viajk)["mcij"] * (*Rai)["dm"];
    HRabij["cdij"] += ( + 1.0 ) * (*Vabic)["cdie"] * (*Rai)["ej"];
    HRabij["cdij"] += ( - 1.0 ) * (*Vabic)["cdje"] * (*Rai)["ei"];
    HRabij["cdij"] += ( - 1.0 ) * (*Fij)["mi"] * Rabij["cdmj"];
    HRabij["cdij"] += ( + 1.0 ) * (*Fij)["mj"] * Rabij["cdmi"];
    HRabij["cdij"] += ( - 1.0 ) * (*Fab)["de"] * Rabij["ecij"];
    HRabij["cdij"] += ( + 1.0 ) * (*Fab)["ce"] * Rabij["edij"];
    HRabij["cdij"] += ( + 0.5 ) * (*Vijkl)["mnij"] * Rabij["cdmn"];
    HRabij["cdij"] += ( + 1.0 ) * (*Viajb)["mdif"] * Rabij["fcmj"];
    HRabij["cdij"] += ( - 1.0 ) * (*Viajb)["mcif"] * Rabij["fdmj"];
    HRabij["cdij"] += ( - 1.0 ) * (*Viajb)["mdjf"] * Rabij["fcmi"];
    HRabij["cdij"] += ( + 1.0 ) * (*Viajb)["mcjf"] * Rabij["fdmi"];
    HRabij["cdij"] += ( + 0.5 ) * (*Vabcd)["cdef"] * Rabij["efij"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["cdmi"] * (*Vijka)["mnjg"] * (*Rai)["gn"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["edij"] * (*Vijab)["noeh"] * Rabij["hcno"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["ecij"] * (*Vijab)["noeh"] * Rabij["hdno"];
    HRabij["cdij"] += ( + 0.25 ) * Tabij["efij"] * (*Vijab)["opef"] * Rabij["cdop"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["cdmi"] * (*Vijab)["mngh"] * Rabij["ghnj"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["cdmj"] * (*Vijab)["mngh"] * Rabij["ghni"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["edni"] * (*Vijab)["noeh"] * Rabij["hcoj"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ecni"] * (*Vijab)["noeh"] * Rabij["hdoj"];
    HRabij["cdij"] += ( + 1.0 ) * Tabij["ednj"] * (*Vijab)["noeh"] * Rabij["hcoi"];
    HRabij["cdij"] += ( - 1.0 ) * Tabij["ecnj"] * (*Vijab)["noeh"] * Rabij["hdoi"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["efoi"] * (*Vijab)["opef"] * Rabij["cdpj"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["efoj"] * (*Vijab)["opef"] * Rabij["cdpi"];
    HRabij["cdij"] += ( + 0.25 ) * Tabij["cdmn"] * (*Vijab)["mngh"] * Rabij["ghij"];
    HRabij["cdij"] += ( + 0.5 ) * Tabij["edno"] * (*Vijab)["noeh"] * Rabij["hcij"];
    HRabij["cdij"] += ( - 0.5 ) * Tabij["ecno"] * (*Vijab)["noeh"] * Rabij["hdij"];

    writeEOMVectors(*HRai, 2);
    writeEOMVectors(HRabij, 4);
    bool new_line(false);
    if (HRai->wrld->rank == 0) {
      new_line = true;
    }
    if (new_line) {
      std::cout  << std::endl;
    }

  }

  LOG(1,"Right EOM END") << ' ' << std::endl;

  allocatedTensorArgument("SimlarityTransformedHamiltonianSD", Hpq);

  LOG(1,"Left EOM Vectors Lengths BEGIN") << " Lengths= "
    << oneBodyLength << " " << twoBodyLength << std::endl;

  for (int64_t i = 0 ; i < totalDimension-1 ; i++) {

     getCanonicalPerturbationBasis(*Lia, Lijab, i);

    // Contruct LH (one body part)
    (*LHia)["ja"]  = 0.0;
    (*LHia)["ja"] += ( - 1.0 ) * (*Fij)["jk"] * (*Lia)["ka"];
    (*LHia)["ja"] += ( + 1.0 ) * (*Fab)["ca"] * (*Lia)["jc"];
    (*LHia)["ja"] += ( - 1.0 ) * (*Viajb)["jcla"] * (*Lia)["lc"];
    (*LHia)["ja"] += ( - 0.5 ) * (*Viajk)["jclm"] * Lijab["mlca"];
    (*LHia)["ja"] += ( - 0.5 ) * (*Vabic)["cdma"] * Lijab["mjdc"];
    (*LHia)["ja"] += ( + 1.0 ) * Tabij["cdmn"] * (*Vijab)["njda"] * (*Lia)["mc"];
    (*LHia)["ja"] += ( + 0.5 ) * Tabij["cdmn"] * (*Vijab)["njcd"] * (*Lia)["ma"];
    (*LHia)["ja"] += ( + 0.5 ) * Tabij["cdmn"] * (*Vijab)["mnda"] * (*Lia)["jc"];
    (*LHia)["ja"] += ( - 0.5 ) * Tabij["cdmn"] * (*Vijka)["njoa"] * Lijab["omdc"];
    (*LHia)["ja"] += ( - 1.0 ) * Tabij["cdmn"] * (*Vijka)["njod"] * Lijab["omca"];
    (*LHia)["ja"] += ( - 0.25 ) * Tabij["cdmn"] * (*Vijka)["mnoa"] * Lijab["ojdc"];
    (*LHia)["ja"] += ( - 0.5 ) * Tabij["cdmn"] * (*Vijka)["mnod"] * Lijab["ojca"];
    (*LHia)["ja"] += ( - 0.5 ) * Tabij["cdmn"] * (*Viabc)["jgda"] * Lijab["nmgc"];
    (*LHia)["ja"] += ( - 0.25 ) * Tabij["cdmn"] * (*Viabc)["jgcd"] * Lijab["nmga"];
    (*LHia)["ja"] += ( - 1.0 ) * Tabij["cdmn"] * (*Viabc)["ngda"] * Lijab["mjgc"];
    (*LHia)["ja"] += ( - 0.5 ) * Tabij["cdmn"] * (*Viabc)["ngcd"] * Lijab["mjga"];

    // Contruct LH (two body part)
    LHijab["klab"]  = 0.0;
    LHijab["klab"] += ( - 1.0 ) * (*Vijka)["klmb"] * (*Lia)["ma"];
    LHijab["klab"] += ( + 1.0 ) * (*Vijka)["klma"] * (*Lia)["mb"];
    LHijab["klab"] += ( + 1.0 ) * (*Viabc)["keab"] * (*Lia)["le"];
    LHijab["klab"] += ( - 1.0 ) * (*Viabc)["leab"] * (*Lia)["ke"];
    LHijab["klab"] += ( - 1.0 ) * (*Fij)["km"] * Lijab["mlab"];
    LHijab["klab"] += ( + 1.0 ) * (*Fij)["lm"] * Lijab["mkab"];
    LHijab["klab"] += ( - 1.0 ) * (*Fab)["eb"] * Lijab["klea"];
    LHijab["klab"] += ( + 1.0 ) * (*Fab)["ea"] * Lijab["kleb"];
    LHijab["klab"] += ( - 0.5 ) * (*Vijkl)["klmn"] * Lijab["nmab"];
    LHijab["klab"] += ( + 1.0 ) * (*Viajb)["kenb"] * Lijab["nlea"];
    LHijab["klab"] += ( - 1.0 ) * (*Viajb)["kena"] * Lijab["nleb"];
    LHijab["klab"] += ( - 1.0 ) * (*Viajb)["lenb"] * Lijab["nkea"];
    LHijab["klab"] += ( + 1.0 ) * (*Viajb)["lena"] * Lijab["nkeb"];
    LHijab["klab"] += ( - 0.5 ) * (*Vabcd)["efab"] * Lijab["klfe"];
    LHijab["klab"] += ( + 0.5 ) * Tabij["efop"] * (*Vijab)["klfb"] * Lijab["poea"];
    LHijab["klab"] += ( - 0.5 ) * Tabij["efop"] * (*Vijab)["klfa"] * Lijab["poeb"];
    LHijab["klab"] += ( - 0.25 ) * Tabij["efop"] * (*Vijab)["klef"] * Lijab["poab"];
    LHijab["klab"] += ( - 0.5 ) * Tabij["efop"] * (*Vijab)["pkab"] * Lijab["olfe"];
    LHijab["klab"] += ( + 0.5 ) * Tabij["efop"] * (*Vijab)["plab"] * Lijab["okfe"];
    LHijab["klab"] += ( - 1.0 ) * Tabij["efop"] * (*Vijab)["pkfb"] * Lijab["olea"];
    LHijab["klab"] += ( + 1.0 ) * Tabij["efop"] * (*Vijab)["pkfa"] * Lijab["oleb"];
    LHijab["klab"] += ( + 1.0 ) * Tabij["efop"] * (*Vijab)["plfb"] * Lijab["okea"];
    LHijab["klab"] += ( - 1.0 ) * Tabij["efop"] * (*Vijab)["plfa"] * Lijab["okeb"];
    LHijab["klab"] += ( + 0.5 ) * Tabij["efop"] * (*Vijab)["pkef"] * Lijab["olab"];
    LHijab["klab"] += ( - 0.5 ) * Tabij["efop"] * (*Vijab)["plef"] * Lijab["okab"];
    LHijab["klab"] += ( - 0.25 ) * Tabij["efop"] * (*Vijab)["opab"] * Lijab["klfe"];
    LHijab["klab"] += ( - 0.5 ) * Tabij["efop"] * (*Vijab)["opfb"] * Lijab["klea"];
    LHijab["klab"] += ( + 0.5 ) * Tabij["efop"] * (*Vijab)["opfa"] * Lijab["kleb"];

    writeEOMVectors(*LHia, 2);
    writeEOMVectors(LHijab, 4);
    bool new_line(false);
    if (LHia->wrld->rank == 0) {
      new_line = true;
    }
    if (new_line) {
      std::cout  << std::endl;
    }

  }

  LOG(1,"Left EOM END") << ' ' << std::endl;

}

template <typename F>
void Mp2EquationOfMotion::writeEOMVectors(
    CTF::Tensor<F> &T,
    unsigned int rank
  ) {

  int nBodyLength(1);
  int arrayCount;
  int64_t *indices;
  F *values;

  for (int i = 0; i < rank; ++i) {
    nBodyLength *= T.lens[i];
  }

  if (T.wrld->rank == 0) {
    arrayCount = nBodyLength;
  } else {
    arrayCount = 0;
  }

  indices = new int64_t[arrayCount];
  values = new F[arrayCount];
  for (int64_t i = 0 ; i < arrayCount; i++) {
     indices[i] = i;
  }

  T.read(arrayCount, indices, values);

  for (int64_t i = 0 ; i < arrayCount; i++) {
     std::cout  << " " << values[i];
  }
}

// instantiate
template
void Mp2EquationOfMotion::writeEOMVectors(
  CTF::Tensor<double> &T,
  unsigned int rank
);


template <typename F>
void Mp2EquationOfMotion::getCanonicalPerturbationBasis(
    CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i) {
  int oneBodyLength(Tai.lens[0] * Tai.lens[1]);
  int twoBodyLength(
      Tabij.lens[0] * Tabij.lens[1] * Tabij.lens[2] *  Tabij.lens[3]
  );

  Tabij["abij"] = 0;
  Tai["ai"] = 0;
  int64_t *indices;
  F *values;
  int arrayCount;

  if (Tabij.wrld->rank == 0) {
    arrayCount = 1;
    values = (F*) malloc(arrayCount);
    indices = (int64_t*) malloc(arrayCount);
    indices[0] = i + 1 <= oneBodyLength ? i : i - oneBodyLength;
    values[0] = 1.0;
  } else {
    arrayCount = 0;
    values = (F*) malloc(arrayCount);
    indices = (int64_t*) malloc(arrayCount);
  }

  if (i+1 <= oneBodyLength) { // One body regime
    Tai.write(arrayCount, indices, values);
  } else { // Two body regime
    Tabij.write(arrayCount, indices, values);
  }

}

// instantiate
template
void Mp2EquationOfMotion::getCanonicalPerturbationBasis(
    CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, int64_t i);


