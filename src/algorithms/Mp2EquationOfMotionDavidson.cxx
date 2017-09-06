#include <algorithms/Mp2EquationOfMotionDavidson.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

#include <math/EigenSystemDavidson.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EquationOfMotionDavidson);

Mp2EquationOfMotionDavidson::Mp2EquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
Mp2EquationOfMotionDavidson::~Mp2EquationOfMotionDavidson() {}

void Mp2EquationOfMotionDavidson::run() {
  typedef CTF::Tensor<> T;

  // Get copy of couloumb integrals
  T Vijkl(*getTensorArgument<double, T>("HHHHCoulombIntegrals"));
  T Vabcd(*getTensorArgument<double, T>("PPPPCoulombIntegrals")); 

  T Vabij(*getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  // T *Vijab(getTensorArgument<double, T>("HHPPCoulombIntegrals")); // swap PPHH (done)

  T Vijka(*getTensorArgument<double, T>("HHHPCoulombIntegrals"));
  // T *Viajk(getTensorArgument<double, T>("HPHHCoulombIntegrals")); // swap HHHP (done)

  T Vaibj(*getTensorArgument<double, T>("PHPHCoulombIntegrals")); // not in eqs
  //T *Viajb(getTensorArgument<double, T>("HPHPCoulombIntegrals")); // swap PHPH (done)

  T Vabci(*getTensorArgument<double, T>("PPPHCoulombIntegrals")); // not in eqs
  //T *Viabc(getTensorArgument<double, T>("HPPPCoulombIntegrals")); // swap PPPH (done)
  //T *Vabic(getTensorArgument<double, T>("PPHPCoulombIntegrals")); // swap PPPH (done)


  // Get orbital energies
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int syms[] = {NS, NS, NS, NS};

  //  Tai
  int vo[] = { Nv, No };
  T Tai(2, vo, syms, *Cc4s::world, "Tai");
  T Tabij(false, Vabij);

  int totalDimension(1 + Nv * No + No * No * Nv * Nv);
  LOG(1, "MP2_EOM_DAVIDSON") << "Nv " << Nv << std::endl;
  LOG(1, "MP2_EOM_DAVIDSON") << "No " << No << std::endl;
  LOG(1, "MP2_EOM_DAVIDSON") << "Problem dimension " << totalDimension << std::endl;


  LOG(1, "MP2_EOM_DAVIDSON") << "Antisymmetrizing Vpqrs " << std::endl;

  //  Vijab
  int oovv[] = { No, No, Nv, Nv };
  T Vijab(4, oovv, syms, *Cc4s::world, "Vijab");
  Vijab["ijab"] =  Vabij["abij"] - Vabij["abji"];

  //  Viajk
  int ovoo[] = { No, Nv, No, No };
  T Viajk(4, ovoo, syms, *Cc4s::world, "Viajk");
  Viajk["iajk"] =  Vijka["ijka"]  - Vijka["ikja"];

  // Viajb
  int ovov[] = { No, Nv, No, Nv };
  T Viajb(4, ovov, syms, *Cc4s::world, "Viajb");
  Viajb["iajb"] =  Vaibj["aibj"] - Vabij["abji"];

  // Viabc
  int ovvv[] = { No, Nv, Nv, Nv };
  T Viabc(4, ovvv, syms, *Cc4s::world, "Viabc");
  Viabc["iabc"] =  Vabci["abci"] - Vabci["acbi"];

  // Vabic
  int vvov[] = { Nv, Nv, No, Nv };
  T Vabic(4, vvov, syms, *Cc4s::world, "Vabic");
  Vabic["abic"] =  Vabci["abci"] - Vabci["baci"];

  // Antisymmetrize integrals
  Vijkl["ijkl"] -= Vijkl["ijlk"];
  Vabcd["abcd"] -= Vabcd["abdc"];
  Vijka["ijka"] -= Vijka["jika"];
  Vaibj["aibj"] -= Vabij["baij"];
  Vabci["abci"] -= Vabci["baci"];
  Vabij["abij"] -= Vabij["abji"];

  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  LOG(1, "MP2_EOM_DAVIDSON") <<
    "Creating doubles amplitudes" << totalDimension << std::endl;
  CTF::Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, Vabij,"abij", Tabij,"abij", 0.0,"abij", fDivide);

  CTF::Scalar<> energy(0.0);
  double energy_val(0.0);

  LOG(2, "MP2_EOM_DAVIDSON") << "Calculating MP2 energy" << std::endl;
  energy[""] = ( 0.25 ) * Tabij["abij"] * Vabij["abij"];
  energy_val = energy.get_val();
  LOG(1, "MP2_EOM_DAVIDSON") << " Mp2 energy = " << energy_val << std::endl;


  // Create L and R and intermediates
  //int oneBodySyms[] = {NS, NS};
  //int oneBodyLensL[] = {No, Nv};
  //T *Lia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "Lia") );
  //T Lijab(false, Vijab);
  //T *LHia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "LHia") );
  //T LHijab(false, Vijab);
  //int oneBodyLensR[] = {Nv, No};
  //T *Rai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "Rai") );
  //T Rabij(false, Vabij);
  //T *HRai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "Rai") );
  //T HRabij(false, Vabij);

  // kinetic terms
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  T Fab(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab");
  int kineticLensOccupied[] = {No, No};
  T Fij(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij");
  Fab["aa"] = (*epsa)["a"];
  Fij["ii"] = (*epsi)["i"];

  Mp2PreConditioner<double, FockVector<double>> P(
    Fij, Fab, Tai, Tabij, Vabcd, Viajb, Vijab, Vijkl
  );
}

template <typename F>
void Mp2EquationOfMotionDavidson::getCanonicalPerturbationBasis(
    CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i) {
  int oneBodyLength(Tai.lens[-1] * Tai.lens[1]);
  //int twoBodyLength(
      //Tabij.lens[0] * Tabij.lens[1] * Tabij.lens[2] *  Tabij.lens[3]
  //);

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
void Mp2EquationOfMotionDavidson::getCanonicalPerturbationBasis(
    CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, int64_t i);

//template <typename F, typename V=FockVector<F>>
//std::vector<V> Mp2PreConditioner<F, V>::getInitialBasis(int eigenVectorsCount){
  //CTF::Scalar<F> nullies;
  //FockVector<F> fockA({{nullies}});
  //return fockA;
//}

template <typename F>
std::vector<FockVector<F>>
Mp2PreConditioner<F, FockVector<F>>::getInitialBasis(int eigenVectorsCount)
{
  CTF::Scalar<F> nullies;
  FockVector<F> fockA({{nullies}});
  return fockA;
}

template
std::vector<FockVector<double>>
Mp2PreConditioner<double, FockVector<double>>::getInitialBasis(int eigenVectorsCount);
