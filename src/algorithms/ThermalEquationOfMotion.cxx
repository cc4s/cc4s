#include <algorithms/ThermalEquationOfMotion.hpp>

#include <tcc/Tcc.hpp>
#include <math/EigenSystemDavidson.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(ThermalEquationOfMotion);

ThermalEquationOfMotion::ThermalEquationOfMotion(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
ThermalEquationOfMotion::~ThermalEquationOfMotion() {}

void ThermalEquationOfMotion::run() {
  typedef CTF::Tensor<> T;

  T *Vabij(getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));
  CTF::Scalar<> E0(*Cc4s::world);
  E0[""] += 1.0;

  ThermalHamiltonian<double> H(
    &E0, epsi, epsa, Vabij
  );
  ThermalHamiltonianPreConditioner<double> P(
    E0, *epsi, *epsa, *Vabij
  );
  std::vector<FockVector<double>> basis( P.getInitialBasis(16) );
  allocatedTensorArgument(
    "VacuumHamiltonianDiagonal",
    new CTF::Tensor<>(P.getDiagonalH().componentTensors[0])
  );
  allocatedTensorArgument(
    "DoublesHamiltonianDiagonal",
    new CTF::Tensor<>(P.getDiagonalH().componentTensors[1])
  );
  allocatedTensorArgument(
    "VacuumBasis",
    new CTF::Tensor<>(basis[8].componentTensors[0])
  );
  allocatedTensorArgument(
    "DoublesBasis",
    new CTF::Tensor<>(basis[8].componentTensors[1])
  );
  // Davidson solver
  EigenSystemDavidson<FockVector<double>> eigenSystem(H, 16, P, 1E-14, 16*16);

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  for (auto &ev: eigenValues) {
    LOG(0, "FT_EOM_DAVIDSON") << "Eigenvalue=" << ev << std::endl;
  }
}



template <typename F>
ThermalHamiltonian<F>::ThermalHamiltonian(
  CTF::Tensor<F> *E0_,
  CTF::Tensor<F> *epsi_,
  CTF::Tensor<F> *epsa_,
  CTF::Tensor<F> *Vabij_
):
  E0(E0_),
  epsi(epsi_),
  epsa(epsa_),
  Vabij(Vabij_)
{
}

template <typename F>
FockVector<F> ThermalHamiltonian<F>::leftApply(
  FockVector<F> &L
) {
  FockVector<F> LH(L);
  // TODO: Implement left apply
  LH *= F(0);
  return LH;
}

template <typename F>
FockVector<F> ThermalHamiltonian<F>::rightApply(
  FockVector<F> &R
) {
  FockVector<F> HR(R);
  // get pointers to the component tensors
  CTF::Tensor<F> *R0( &R.componentTensors[0] );
  CTF::Tensor<F> *Rabij( &R.componentTensors[1] );
  CTF::Tensor<F> *HR0( &HR.componentTensors[0] );
  CTF::Tensor<F> *HRabij( &HR.componentTensors[1] );


  // HR (vacuum part)
  // from R
  (*HR0)[""]  = (*E0)[""] * (*R0)[""];
  // from Rabij
  (*HR0)[""] += 0.5 * (*Vabij)["abij"] * (*Rabij)["abij"];
  (*HR0)[""] -= 0.5 * (*Vabij)["abij"] * (*Rabij)["abji"];

  // HRabij (two body part)
  // from R
  (*HRabij)["abij"] = 0.5 * (*Vabij)["abij"] * (*R0)[""];
  (*HRabij)["abij"] = 0.5 * (*Vabij)["abji"] * (*R0)[""];

  // from Rabij
  (*HRabij)["abij"] += (*epsa)["a"] * (*Rabij)["abij"];
  (*HRabij)["abij"] += (*epsa)["b"] * (*Rabij)["abij"];
  (*HRabij)["abij"] -= (*epsi)["i"] * (*Rabij)["abij"];
  (*HRabij)["abij"] -= (*epsi)["j"] * (*Rabij)["abij"];

  return HR;
}

// instantiate class
template
class ThermalHamiltonian<double>;


template <typename F>
ThermalHamiltonianPreConditioner<F>::ThermalHamiltonianPreConditioner(
  CTF::Tensor<F> &E0,
  CTF::Tensor<F> &epsi,
  CTF::Tensor<F> &epsa,
  CTF::Tensor<F> &Vabij
): diagonalH({{E0, Vabij}, {"", "abij"}}) {
  // pointers to vacuum and doubles tensors of diagonal part
  CTF::Tensor<F> *D( &diagonalH.componentTensors[0] );
  CTF::Tensor<F> *Dabij( &diagonalH.componentTensors[1] );

  // calculate diagonal elements of H
  (*D)[""] = E0[""];

  (*Dabij)["abij"] =  epsa["a"];
  (*Dabij)["abij"] += epsa["b"];
  (*Dabij)["abij"] -= epsi["i"];
  (*Dabij)["abij"] -= epsi["j"];
}

template <typename F>
class EomDiagonalValueComparator;

template <>
class EomDiagonalValueComparator<double> {
public:
  bool operator ()(
    const std::pair<int, double> &a,
    const std::pair<int, double> &b
  ) {
    double A(
      std::abs(a.second) < 1E-13 ?
        std::numeric_limits<double>::infinity() : a.second
    );
    double B(
      std::abs(b.second) < 1E-13 ?
        std::numeric_limits<double>::infinity() : b.second
    );
    double diff(B-A);
    double magnitude(std::abs(A)+std::abs(B));
    if (std::real(diff) > +1E-13*magnitude) return true;
    if (std::real(diff) < -1E-13*magnitude) return false;
    return a.first < b.first;
  }
};


template <typename F>
std::vector<FockVector<F>> ThermalHamiltonianPreConditioner<F>::getInitialBasis(
  const int eigenVectorsCount
) {
  LOG(0, "FT_EOM_DAVIDSON") << "Get initial basis" << std::endl;
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<int64_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    EomDiagonalValueComparator<double>()
  );

  // gather all K lowest elements of each processor at root
  //   convert into homogeneous arrays for MPI gather
  std::vector<int64_t> localLowestElementIndices(eigenVectorsCount);
  std::vector<F> localLowestElementValues(eigenVectorsCount);
  for (int i(0); i < eigenVectorsCount; ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
  int lowestElementsCount(
    communicator.getRank() == 0 ?
      eigenVectorsCount * communicator.getProcesses() : 0
  );
  std::vector<int64_t> lowestElementIndices(lowestElementsCount);
  std::vector<F> lowestElementValues(lowestElementsCount);
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<int64_t, F>> lowestElements(lowestElementsCount);
  for (int i(0); i < lowestElementsCount; ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    EomDiagonalValueComparator<double>()
  );
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;
  //for (int b(0); b < eigenVectorsCount; ++b) {
  int bb(0);
  int b(0);
  while (bb < eigenVectorsCount) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<int64_t,F>> elements;
    if (communicator.getRank() == 0) {
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    b++;
    std::cout << "b" << b << std::endl;
    if (basisElement.dot(basisElement) < 1e-7) continue;
    bb++;
    basis.push_back(basisElement);
    std::cout << "bb" << bb << std::endl;
  }
  return basis;
}



template <typename F>
FockVector<F> ThermalHamiltonianPreConditioner<F>::getCorrection(
  const complex lambda, FockVector<F> &residuum
) {
  FockVector<F> w(diagonalH);
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_): lambda(lambda_) {
    }
    F operator ()(const F residuumElement, const F diagonalElement) {
      return std::abs(lambda - diagonalElement) < 1E-4 ?
        0.0 : residuumElement / (lambda - diagonalElement);
    }
  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  FockVector<F> correction(diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.componentTensors.size(); ++c) {
    const char *indices( correction.componentIndices[c].c_str() );
    correction.componentTensors[c].contract(
      1.0,
      residuum.componentTensors[c],indices,
      diagonalH.componentTensors[c],indices,
      0.0,indices,
      CTF::Bivar_Function<F>(diagonalCorrection)
    );
  }
  return correction;
}


// instantiate class
template
class ThermalHamiltonianPreConditioner<double>;

