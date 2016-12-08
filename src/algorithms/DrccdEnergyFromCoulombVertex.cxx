#include <algorithms/DrccdEnergyFromCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombVertex);

DrccdEnergyFromCoulombVertex::DrccdEnergyFromCoulombVertex(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

DrccdEnergyFromCoulombVertex::~DrccdEnergyFromCoulombVertex() {
}

void DrccdEnergyFromCoulombVertex::iterate(int i) {
  // Read the DRCCD amplitudes Tabij
  Tensor<> *Tabij(&TabijMixer->getNext());

  // Construct intermediate Amplitudes
  Tensor<> Rabij(false, *Tabij);

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
                 abbreviation.begin(), ::toupper);

  LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;

    if (i == 0) {
      // For first iteration compute only the MP2 amplitudes 
      // Since Tabij = 0, Vabij is the only non-zero term

      // Read Vabij
      Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

      Rabij["abij"] += (*Vabij)["abij"];
    } 
    else {
      // For the rest iterations compute the DRCCD amplitudes

      // Read coulomb vertex GammaGai
      Tensor<complex> *GammaGai(getTensorArgument<complex>
                                ("ParticleHoleCoulombVertex"));

      // Allocate real and imag part of GammaGai
      Tensor<> realGammaGai(3, GammaGai->lens, GammaGai->sym, 
                            *GammaGai->wrld, "RealGammaGai");
      Tensor<> imagGammaGai(3, GammaGai->lens, GammaGai->sym, 
                            *GammaGai->wrld, "ImagGammaGai");

      // Split into real and imaginary parts
      fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

      // Construct left and right dressed Coulomb vertices GammaGai
      Tensor<> leftRealGammaGai (realGammaGai);
      Tensor<> leftImagGammaGai (imagGammaGai);
      Tensor<> rightRealGammaGai(realGammaGai);
      Tensor<> rightImagGammaGai(imagGammaGai);
      leftRealGammaGai ["Gai"] += 2.0 * realGammaGai["Gbk"] * (*Tabij)["abik"];
      leftImagGammaGai ["Gai"] += 2.0 * imagGammaGai["Gbk"] * (*Tabij)["abik"];
      rightRealGammaGai["Gai"] += 2.0 * realGammaGai["Gbk"] * (*Tabij)["baki"];
      rightImagGammaGai["Gai"] += 2.0 * imagGammaGai["Gbk"] * (*Tabij)["baki"];

      // Construct T2 amplitudes
      Rabij["abij"]  = leftRealGammaGai["Gai"] * rightRealGammaGai["Gbj"];
      Rabij["abij"] += leftImagGammaGai["Gai"] * rightImagGammaGai["Gbj"];
    }

    // Calculate the amplitudes from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // And append them to the mixer
    TabijMixer->append(Rabij);
}


void DrccdEnergyFromCoulombVertex::dryIterate() {
  {

    // TODO: the Mixer should provide a DryTensor in the future
    // Read the DRCCD amplitudes Tabij
    DryTensor<> *Tabij(getTensorArgument<double, DryTensor<double>>
                     ("DrccdDoublesAmplitudes"));

    // Read the Coulomb Vertex GammaGai
    DryTensor<complex> *GammaGai(getTensorArgument<complex, 
                                 DryTensor<complex>>
                                 ("ParticleHoleCoulombVertex"));

    // Compute the No,Nv,NG,Np
    int NG(GammaGai->lens[0]);
    int No(GammaGai->lens[2]);
    int Nv(GammaGai->lens[1]);
    int syms[] = { NS, NS, NS, };

    // Allocate and realGammaGai and imagGammaGai
    int GaiLens[] = {NG,Nv,No};

    DryTensor<> realGammaGai(3, GaiLens, syms);
    DryTensor<> imagGammaGai(3, GaiLens, syms);
  
    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Tabij);

    // TODO: implment dryDoublesAmplitudesFromResiduum
    // at the moment, assume usage of Dabij
    DryTensor<> Dabij(*Tabij);
  }
}

