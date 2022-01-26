/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithms/CoupledCluster.hpp>
#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>
#include <mixers/Mixer.hpp>
#include <MathFunctions.hpp>
#include <SharedPointer.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Options.hpp>
#include <Cc4s.hpp>
#include <Timer.hpp>


#include <array>
#include <initializer_list>
#include <iomanip>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CoupledCluster)

Ptr<MapNode> CoupledCluster::run(const Ptr<MapNode> &arguments_){
  this->arguments = arguments_;
  // multiplex calls to template methods
  Ptr<MapNode> result;
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    (result = run<Real<>,TE>()) || (result = run<Complex<>,TE>());
  } else {
    using TE = DefaultTensorEngine;
    (result = run<Real<>,TE>()) || (result = run<Complex<>,TE>());
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}


template <typename F, typename TE>
Ptr<MapNode> CoupledCluster::run() {
  using namespace std;

  auto coulombIntegrals(arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals"));
  if (!coulombIntegrals) return nullptr;

  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));

  Natural<> i(0);
  Ptr<TensorSet<F,TE>> amplitudes;
  if (arguments->isGiven("initialAmplitudes")) {
    amplitudes = arguments->getPtr<TensorSet<F,TE>>(
      "initialAmplitudes"
    );
    OUT() << "Using given initial amplitudes " << amplitudes << endl;
  }

  energy = New<MapNode>(SOURCE_LOCATION);

  // create a method handler, by default Ccsd
  auto methodName( arguments->getValue<std::string>("method", "Ccsd") );
  Ptr<CoupledClusterMethod<F,TE>> method(
    CoupledClusterMethodFactory<F,TE>::create(methodName, arguments)
  );
  ASSERT_LOCATION(
    method, std::string("Unknown method: '") + methodName + "'",
    arguments->get("method")->sourceLocation
  );
  OUT() << "Using method "
    << methodName << ". " << method->describeOptions() << endl;

  // create a mixer, by default use the linear one
  auto mixerArguments(arguments->getMap("mixer"));
  auto mixerType(mixerArguments->getValue<std::string>("type", "DiisMixer"));
  Ptr<Mixer<F,TE>> mixer(MixerFactory<F,TE>::create(mixerType, mixerArguments));
  ASSERT_LOCATION(
    mixer, std::string("Unknown mixer type: '") + mixerType + "'",
    mixerArguments->get("type")->sourceLocation
  );
  OUT() << "Using mixer "
    << mixerType << ". " << mixer->describeOptions() << endl;

  // number of iterations for determining the amplitudes
  auto maxIterationsCount(
    arguments->getValue<size_t>("maxIterations", DEFAULT_MAX_ITERATIONS)
  );

  auto amplitudesConvergence(
    arguments->getValue<Real<>>(
      "amplitudesConvergence", DEFAULT_AMPLITUDES_CONVERGENCE
    )
  );
  auto energyConvergence(
    arguments->getValue<Real<>>("energyConvergence", DEFAULT_ENERGY_CONVERGENCE)
  );

  OUT() << "Maximum number of iterations: " << maxIterationsCount << endl;
  OUT() << "Unless reaching energy convergence dE: " << energyConvergence << endl;
  OUT() << "and amplitudes convergence dR: " << amplitudesConvergence << endl;
  F e(0), previousE(0);
  Real<> residuumNorm;
  OUT()
    << "Iter         Energy         dE           dR         time   GF/s/core"
    << endl;

  bool isSecondOrder;
  for (; i < maxIterationsCount; ++i) {
    LOG() << "iteration: " << i+1 << endl;
    Time time;
    Natural<128> operations;
    {
      Timer timer(&time);
      OperationsCounter operationsCounter(&operations);
      auto residuum( method->getResiduum(amplitudes) );
      residuumToAmplitudes(residuum, amplitudes);
      auto amplitudesChange( New<TensorSet<F,TE>>(*residuum) );
      if (amplitudes) {
        *amplitudesChange -= *amplitudes;
        isSecondOrder = false;
      } else {
        isSecondOrder = true;
      }
      mixer->append(residuum, amplitudesChange);
      // get mixer's best guess for amplitudes
      amplitudes = mixer->get();
      residuumNorm = mixer->getResiduumNorm();
      e = getEnergy(amplitudes);
    }

    OUT()
      << setw(4) << i+1 << " "
      << scientific
      << setw(16) << setprecision(8) << real(e) << " "
      << setw(12) << setprecision(4) << real(e - previousE) << " "
      << setw(12) << setprecision(4) << real(residuumNorm) << " "
      << fixed
      << setw(8) << setprecision(1) << time.getFractionalSeconds() << " "
      << setw(6) << setprecision(1)
      << operations / 1e9 / time.getFractionalSeconds()
        / Cc4s::world->getProcesses()
      << endl;
    if (isSecondOrder) {
      energy->setValue("secondOrder", real(e));
    }
    if (
      !Cc4s::options->dryRun &&
      abs(e-previousE) < energyConvergence &&
      residuumNorm < amplitudesConvergence
    ) {
      break;
    }
    previousE = e;
  }

  if (maxIterationsCount == 0) {
    OUT() << "computing energy from given amplitudes" << endl;
  } else if (i == maxIterationsCount) {
    WARNING_LOCATION(arguments->sourceLocation) <<
      "energy or amplitudes convergence not reached." << endl;
  }

  e = getEnergy(amplitudes, true);
  bool convergenceReached = i < maxIterationsCount;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  result->setValue("convergenceReached", convergenceReached);
  result->setPtr("amplitudes", amplitudes);
  return result;
}

template <typename F, typename TE>
F CoupledCluster::getEnergy(
  const Ptr<TensorSet<F,TE>> &amplitudes,
  const bool isFinalReport
) {
  auto methodName( arguments->getValue<std::string>("method", "Ccsd") );

  // get the Coulomb integrals to compute the energy
  auto coulombIntegrals(arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals"));
  auto Vijab(coulombIntegrals->get("hhpp"));
  // TODO: get from size of spin properies
  Real<> degeneracy(2);

  auto Tai( amplitudes->get("ph") );
  auto Tabij( amplitudes->get("pphh") );
  // TODO: determine units in tcc
  Tai->inspect()->getUnit() = Tabij->inspect()->getUnit() = 1.0;
  F e;
  std::streamsize ss(std::cout.precision());
  {
    auto direct( Tcc<TE>::template tensor<F>("D") );
    auto exchange( Tcc<TE>::template tensor<F>("X") );
    COMPILE(
      (*direct)[""] <<=
        0.5*degeneracy*degeneracy * (*Tabij)["abij"] * (*Vijab)["ijab"],
      (*direct)[""] +=
        0.5*degeneracy*degeneracy * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijab"],
      (*exchange)[""] <<=
        -0.5*degeneracy * (*Tabij)["abij"] * (*Vijab)["ijba"],
      (*exchange)[""]
        += -0.5*degeneracy * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijba"]
    )->execute();
    F D(direct->read());
    F X(exchange->read());
    e = D+X;
    if (isFinalReport){
      OUT() << std::endl;
      OUT() << methodName << " correlation energy:          "
        << std::setprecision(10) << real(e) << std::endl;
      if (energy->get("secondOrder")) {
        OUT() << "2nd-order correlation energy:     " << std::setprecision(10)
          << energy->getValue<Real<>>("secondOrder") << std::endl;
      }
    }
    energy->setValue("correlation", real(e));
    energy->setValue("direct", real(D));
    energy->setValue("exchange", real(X));
    energy->setValue("unit", Vijab->inspect()->getUnit());
  }
  std::cout << std::setprecision(ss);

  return e;
}

template <typename F, typename TE>
void CoupledCluster::residuumToAmplitudes(
  const Ptr<TensorSet<F,TE>> &residuum,
  const Ptr<TensorSet<F,TE>> &amplitudes
) {
  auto levelShift(
    arguments->getValue<Real<>>("levelShift", DEFAULT_LEVEL_SHIFT)
  );

  if (amplitudes && levelShift != 0.0) {
    // apply level shifting on right hand side, if given and amplitudes present
    *residuum -= F(levelShift) * *amplitudes;
  }

  for (auto key: residuum->getKeys()) {
    auto R( residuum->get(key) );
    auto indices( residuum->generateIndices(key) );
    auto D( calculateEnergyDifferences<F,TE>(R->inspect()->getLens(),indices) );

    // divide by -Delta to get new estimate for T
    COMPILE(
      (*D)[indices] <<= map<F>(
        [levelShift](F delta) { return F(-1) / (delta + levelShift); },
        (*D)[indices]
      ),
      (*R)[indices] <<= (*R)[indices] * (*D)[indices]
    )->execute();
  }
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
) {
  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));
  auto Fepsh(Tcc<TE>::template tensor<F>(epsh->inspect()->getLens(), "Fepsh"));
  auto Fepsp(Tcc<TE>::template tensor<F>(epsp->inspect()->getLens(), "Fepsp"));
  // convert to type F (either complex or double)
  auto fromReal( [](Real<> eps) {return F(eps);} );
  COMPILE(
    (*Fepsp)["a"] <<= map<F>(fromReal, (*epsp)["a"]),
    (*Fepsh)["i"] <<= map<F>(fromReal, (*epsh)["i"])
  )->execute();

  auto D(
    Tcc<TE>::template tensor<F>(lens, std::string("D") + indices)
  );

  // create energy difference tensor
  int excitationLevel(indices.length()/2);
  for (int p(0); p < excitationLevel; ++p) {
    COMPILE(
      (*D)[indices] += (*Fepsp)[indices.substr(p,1)],
      (*D)[indices] -= (*Fepsh)[indices.substr(excitationLevel+p,1)]
    )->execute();
  }

  return D;
}

constexpr Real<64> CoupledCluster::DEFAULT_ENERGY_CONVERGENCE;
constexpr Real<64> CoupledCluster::DEFAULT_AMPLITUDES_CONVERGENCE;
constexpr Real<64> CoupledCluster::DEFAULT_LEVEL_SHIFT;

