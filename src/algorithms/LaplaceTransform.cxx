#include <algorithms/LaplaceTransform.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/LapackMatrix.hpp>
#include <util/Log.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(LaplaceTransform);

LaplaceTransform::LaplaceTransform(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

LaplaceTransform::~LaplaceTransform() {
}

void LaplaceTransform::run() {
  // get imaginary time and frequency grids on all nodes
  auto tn( getTensorArgument<real>("ImaginaryTimePoints") );
  std::vector<real> taus(tn->lens[0]);
  tn->read_all(taus.data());
  auto twn( getTensorArgument<real>("ImaginaryTimeWeights") );
  std::vector<real> tauWeights(twn->lens[0]);
  twn->read_all(tauWeights.data());
  auto nun( getTensorArgument<real>("EvenImaginaryFrequencyPoints") );
  std::vector<real> nus(nun->lens[0]);
  nun->read_all(nus.data());
  auto nuwn( getTensorArgument<real>("EvenImaginaryFrequencyWeights") );
  std::vector<real> nuWeights(nuwn->lens[0]);
  nuwn->read_all(nuWeights.data());

  // contour shift on real axis, such that s = eta + i*nu
  real eta( getRealArgument("EnergyShift", 1.0) );

  // forward and inverse cosine and sin transform weights
  auto ctfCTvn( getTensorArgument<real>("CosineTransform") );
  auto ctfSTvn( getTensorArgument<real>("SineTransform") );
  auto ctfICTnv( getTensorArgument<real>("InverseCosineTransform") );
  auto ctfISTnv( getTensorArgument<real>("InverseSineTransform") );

  // for testing
  CTF::Tensor<complex> cCTvn(2, ctfCTvn->lens, ctfCTvn->sym);
  toComplexTensor(*ctfCTvn, cCTvn);
  CTF::Tensor<complex> cSTvn(2, ctfSTvn->lens, ctfSTvn->sym);
  toComplexTensor(*ctfSTvn, cSTvn);
  CTF::Tensor<complex> cICTnv(2, ctfICTnv->lens, ctfICTnv->sym);
  toComplexTensor(*ctfICTnv, cICTnv);
  CTF::Tensor<complex> cISTnv(2, ctfISTnv->lens, ctfISTnv->sym);
  toComplexTensor(*ctfISTnv, cISTnv);

  LapackMatrix<real> CTvn(*ctfCTvn);
  LapackMatrix<real> STvn(*ctfSTvn);
  LapackMatrix<real> ICTnv(*ctfICTnv);
  LapackMatrix<real> ISTnv(*ctfISTnv);

  // cos(nu_v*tau_n) and sin(nu_v*tau_n), could also be done in CTF
  LapackMatrix<real> Cvn(CTvn);
  LapackMatrix<real> Svn(CTvn);
  for (size_t n(0); n < taus.size(); ++n) {
    for (size_t v(0); v < nus.size(); ++v) {
      Cvn(v,n) = std::cos(nus[v]*taus[n]);
      Svn(v,n) = std::sin(nus[v]*taus[n]);
    }
  }

  // convert to complex CTF tensors
  auto ctfCvn( *ctfCTvn );
  Cvn.write(ctfCvn);
  auto ctfSvn( *ctfSTvn );
  Svn.write(ctfSvn);
  CTF::Tensor<complex> cCvn(2, ctfCvn.lens, ctfCvn.sym);
  toComplexTensor(ctfCvn, cCvn);
  CTF::Tensor<complex> cSvn(2, ctfSvn.lens, ctfSvn.sym);
  toComplexTensor(ctfSvn, cSvn);

  auto No( getTensorArgument<real>("ThermalHoleEigenEnergies")->lens[0] );
  auto Nv( getTensorArgument<real>("ThermalParticleEigenEnergies")->lens[0] );
  auto Nt( tn->lens[0] );
  CTF::Tensor<real> Dai(2, std::vector<int>({Nv,No}).data());
  fetchDelta(Dai);
  auto Tain( new CTF::Tensor<complex>(3, std::vector<int>({Nv,No,Nt}).data()) );
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real Delta, complex &T) { T = Delta; }
    )
  ) (
    Dai["ai"], (*Tain)["ain"]
  );
  auto TRaiv( new CTF::Tensor<complex>(*Tain) );
  auto TIaiv( new CTF::Tensor<complex>(*Tain) );
  // propagator in imaginary time
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real tau, complex &T) { T = std::exp(-std::real(T)*tau); }
    )
  ) (
    (*tn)["n"], (*Tain)["ain"]
  );
  // propagator in imaginary frequency, real part
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real nu, complex &T) {
        real Delta_( std::real(T)+eta );
        T = Delta_ / (Delta_*Delta_ + nu*nu);
      }
    )
  ) (
    (*nun)["v"], (*TRaiv)["aiv"]
  );
  // propagator in imaginary frequency, imaginary part
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real nu, complex &T) {
        real Delta_( std::real(T)+eta );
        T = -nu / (Delta_*Delta_ + nu*nu);
      }
    )
  ) (
    (*nun)["v"], (*TIaiv)["aiv"]
  );
  // numericall transform from time to frequency
  auto T_ain( *Tain );
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(-eta*tau); }
    )
  ) (
    (*tn)["n"], T_ain["ain"]
  );
  auto NTRaiv( new CTF::Tensor<complex>(false, *TRaiv) );
  (*NTRaiv)["aiv"] = (+1.0) * T_ain["ain"] * cCTvn["vn"] * cCvn["vn"];
  auto NTIaiv( new CTF::Tensor<complex>(false, *TIaiv) );
  (*NTIaiv)["aiv"] = (-1.0) * T_ain["ain"] * cSTvn["vn"] * cSvn["vn"];
  {
    CTF::Tensor<complex> diffTaiv(*NTRaiv);
    diffTaiv["aiv"] += (-1.0) * (*TRaiv)["aiv"];
    real forwardError(frobeniusNorm(diffTaiv));
    LOG(0, "Laplace")
      << "forward Laplace transform error real=" << forwardError << std::endl;
  }
  {
    CTF::Tensor<complex> diffTaiv(*NTIaiv);
    diffTaiv["aiv"] += (-1.0) * (*TIaiv)["aiv"];
    real forwardError(frobeniusNorm(diffTaiv));
    LOG(0, "Laplace")
      << "forward Laplace transform error imag=" << forwardError << std::endl;
  }
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagatorsReal", TRaiv
  );
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagatorsImag", TIaiv
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryFrequencyPropagatorsReal", NTRaiv
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryFrequencyPropagatorsImag", NTIaiv
  );

  // numericall transform from frequency to time
  auto NRTain( new CTF::Tensor<complex>(false, *Tain) );
  (*NRTain)["ain"] = (+4.0) * (*TRaiv)["aiv"] * cICTnv["nv"] * cCvn["vn"];
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(+eta*tau); }
    )
  ) (
    (*tn)["n"], (*NRTain)["ain"]
  );
  auto NITain( new CTF::Tensor<complex>(false, *Tain) );
  (*NITain)["ain"] = (-4.0) * (*TIaiv)["aiv"] * cISTnv["nv"] * cSvn["vn"];
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(+eta*tau); }
    )
  ) (
    (*tn)["n"], (*NITain)["ain"]
  );
  {
    CTF::Tensor<complex> diffTain(*NRTain);
    diffTain["ain"] += (-1.0) * (*Tain)["ain"];
    real inverseError(frobeniusNorm(diffTain));
    LOG(0, "Laplace")
      << "inverse Laplace transform error real=" << inverseError << std::endl;
  }
  {
    CTF::Tensor<complex> diffTain(*NITain);
    diffTain["ain"] += (-1.0) * (*Tain)["ain"];
    real inverseError(frobeniusNorm(diffTain));
    LOG(0, "Laplace")
      << "inverse Laplace transform error imag=" << inverseError << std::endl;
  }
  allocatedTensorArgument<complex>(
    "ExactImaginaryTimePropagators", Tain
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryTimePropagatorsReal", NRTain
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryTimePropagatorsImag", NITain
  );
}

std::string LaplaceTransform::getAmplitudeIndices(Tensor<> &T) {
  char indices[T.order+1];
  const int excitationLevel(T.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    indices[i] = static_cast<char>('a'+i);
    indices[i+excitationLevel] = static_cast<char>('i'+i);
  }
  indices[T.order] = 0;
  return indices;
}

void LaplaceTransform::fetchDelta(Tensor<> &Delta) {
  Tensor<real> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<real> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  std::string indices(getAmplitudeIndices(Delta));
  real factor(0.0);
  const int excitationLevel(Delta.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    char aIndex[] = {static_cast<char>('a'+i), 0};
    Delta.sum(+1.0, *epsa,aIndex, factor,indices.c_str());
    factor = 1.0;
    char iIndex[] = {static_cast<char>('i'+i), 0};
    Delta.sum(-1.0, *epsi,iIndex, 1.0,indices.c_str());
  }
}

