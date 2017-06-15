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
  T *Vabij(getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  T *Vijkl(getTensorArgument<double, T>("HHHHCoulombIntegrals"));
  T *Vijka(getTensorArgument<double, T>("HHHPCoulombIntegrals"));
  T *Vijab(getTensorArgument<double, T>("HHPPCoulombIntegrals"));
  T *Viajk(getTensorArgument<double, T>("HPHHCoulombIntegrals"));
  T *Viajb(getTensorArgument<double, T>("HPHPCoulombIntegrals"));
  T *Viabc(getTensorArgument<double, T>("HPPPCoulombIntegrals"));
  T *Vabic(getTensorArgument<double, T>("PPHPCoulombIntegrals"));
  T *Vabcd(getTensorArgument<double, T>("PPPPCoulombIntegrals"));

  LOG(1, "MP2_EOM") << "Antisymmetrizing Vabij " << std::endl;
  (*Vabij)["abij"] -= (*Vabij)["abji"];

  T Tabij(false, Vabij);
  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  LOG(1, "MP2_EOM") <<
    "Creating doubles amplitudes" << totalDimension << std::endl;
  CTF::Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*Vabij),"abij", Tabij,"abij", 0.0,"abij", fDivide);

  // Create L and R
  int oneBodyLens[] = {Nv, No};
  int oneBodySyms[] = {NS, NS};
  T *Lai( new T(2, oneBodyLens, oneBodySyms, *Cc4s::world, "Lai") );
  T Labij(false, Vabij);
  T *Rai( new T(2, oneBodyLens, oneBodySyms, *Cc4s::world, "Rai") );
  T Rabij(false, Vabij);

  // kinetic terms
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  T Fab( new T(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab") );
  int kineticLensOccupied[] = {No, No};
  T Fij( new T(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij") );

  Fab["aa"] = (*epsa)["a"];
  Fij["ii"] = (*epsi)["i"];

  //The totalDimension should be totalDimension, but the zero-particle part is
  //zero, so we restrict the hamiltonian to the space spanned by the singles
  //and doubles excitations
  int hLens[] = {totalDimension-1, totalDimension-1};
  int hSyms[] = {NS, NS};
  T *Hpq( new CTF::Tensor<>(2, hLens, hSyms, *Cc4s::world, "Hpq") );

  int64_t hIndices[1];
  double hValues[1];

  CTF::Scalar<> energy(0.0);
  double energy_val(0.0);


  for (int64_t i = 0 ; i < totalDimension-1 ; i++) {
    getCanonicalPerturbationBasis(*Lai, Labij, i);
    Labij["aaij"] = 0.0;
    Labij["abii"] = 0.0;
    for (int64_t j = 0 ; j < totalDimension-1; j++) {
      getCanonicalPerturbationBasis(*Rai, Rabij, j);
      Rabij["abii"] = 0.0;
      Rabij["aaij"] = 0.0;

      energy[""]  = ( - 1.0 ) * (*Lai)["ib"] * Fij["ki"] * (*Rai)["bk"];
      energy[""] += ( + 1.0 ) * (*Lai)["ib"] * Fab["bc"] * (*Rai)["ci"];
      energy[""] += ( - 1.0 ) * (*Lai)["ib"] * (*Viajb)["kbid"] * (*Rai)["dk"];
      energy[""] += ( + 1.0 ) * (*Lai)["ib"] * Tabij["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
      energy[""] += ( - 0.5 ) * (*Lai)["ib"] * Tabij["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
      energy[""] += ( - 0.5 ) * (*Lai)["ib"] * Tabij["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];
      energy[""] += ( + 0.5 ) * (*Lai)["ib"] * (*Vijka)["klie"] * Rabij["ebkl"];
      energy[""] += ( + 0.5 ) * (*Lai)["ib"] * (*Viabc)["kbde"] * Rabij["deki"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * (*Viajk)["mdij"] * (*Rai)["cm"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * (*Viajk)["mcij"] * (*Rai)["dm"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * (*Vabic)["cdie"] * (*Rai)["ej"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * (*Vabic)["cdje"] * (*Rai)["ei"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Fij["mi"] * Rabij["cdmj"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Fij["mj"] * Rabij["cdmi"];
      energy[""] += ( - 1.0 ) * Labij["ijdc"] * Fab["de"] * Rabij["ecij"];
      energy[""] += ( + 1.0 ) * Labij["ijdc"] * Fab["ce"] * Rabij["edij"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * (*Viajb)["mdif"] * Rabij["fcmj"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * (*Viajb)["mcif"] * Rabij["fdmj"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * (*Viajb)["mdjf"] * Rabij["fcmi"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * (*Viajb)["mcjf"] * Rabij["fdmi"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["edij"] * (*Vijab)["noeh"] * Rabij["hcno"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["ecij"] * (*Vijab)["noeh"] * Rabij["hdno"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["cdmi"] * (*Vijab)["mngh"] * Rabij["ghnj"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["cdmj"] * (*Vijab)["mngh"] * Rabij["ghni"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["edni"] * (*Vijab)["noeh"] * Rabij["hcoj"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ecni"] * (*Vijab)["noeh"] * Rabij["hdoj"];
      energy[""] += ( + 1.0 ) * Labij["ijcd"] * Tabij["ednj"] * (*Vijab)["noeh"] * Rabij["hcoi"];
      energy[""] += ( - 1.0 ) * Labij["ijcd"] * Tabij["ecnj"] * (*Vijab)["noeh"] * Rabij["hdoi"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["efoi"] * (*Vijab)["opef"] * Rabij["cdpj"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["efoj"] * (*Vijab)["opef"] * Rabij["cdpi"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * Tabij["edno"] * (*Vijab)["noeh"] * Rabij["hcij"];
      energy[""] += ( - 0.5 ) * Labij["ijcd"] * Tabij["ecno"] * (*Vijab)["noeh"] * Rabij["hdij"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * (*Vijkl)["mnij"] * Rabij["cdmn"];
      energy[""] += ( + 0.5 ) * Labij["ijcd"] * (*Vabcd)["cdef"] * Rabij["efij"];
      energy[""] += ( + 0.25 ) * Labij["ijcd"] * Tabij["efij"] * (*Vijab)["opef"] * Rabij["cdop"];
      energy[""] += ( + 0.25 ) * Labij["ijcd"] * Tabij["cdmn"] * (*Vijab)["mngh"] * Rabij["ghij"];

      energy_val = energy.get_val();
      hValues[0] = energy_val;
      hIndices[0] = i + j * (totalDimension - 1);
      (*Hpq).write(1, hIndices, hValues);
      LOG(1, "MP2_EOM") << "< " << i << " |H| " <<  j << " >"  << " = " << energy_val << std::endl;
    }
  }

  allocatedTensorArgument("SimlarityTransformedHamiltonianSD", Hpq);

}

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
    indices[0] = i;
    values[0] = 1.0;
  } else {
    arrayCount = 0;
    values = (F*) malloc(arrayCount);
    indices = (int64_t*) malloc(arrayCount);
  }

  if (i+1 <= oneBodyLength) { // One body regime
    Tai.write(arrayCount, indices, values);
  } else { // Two body regime
    indices[0] = i - oneBodyLength;
    Tabij.write(arrayCount, indices, values);
  }
  //Tai.print();
  //Tabij.print();

}

// instantiate
template
void Mp2EquationOfMotion::getCanonicalPerturbationBasis(
    CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, int64_t i);
