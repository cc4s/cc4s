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
#include <math/MathFunctions.hpp>
#include <mixers/Mixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Options.hpp>
#include <Cc4s.hpp>
#include <util/Timer.hpp>


#include <array>
#include <initializer_list>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CoupledCluster)

Ptr<MapNode> CoupledCluster::run(const Ptr<MapNode> &arguments){
  this->arguments = arguments;
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto scalarType(coulombIntegrals->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return run<Real<>,TE>();
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>();
    }
  } else {
    using TE = DefaultTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return run<Real<>,TE>();
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>();
    }
  }
  ASSERT_LOCATION(
    false, "unsupported orbitals type '" + scalarType + "'",
    coulombIntegrals->get("scalarType")->sourceLocation
  );
}


template <typename F, typename TE>
Ptr<MapNode> CoupledCluster::run() {
  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<TensorExpression<Real<>,TE>>("h"));
  auto epsa(energySlices->getPtr<TensorExpression<Real<>,TE>>("p"));

  auto No(epsi->inspect()->getLen(0));
  auto Nv(epsa->inspect()->getLen(0));

  Natural<> i(0); bool restart(false);
  Ptr<TensorUnion<F,TE>> amplitudes;
  if (arguments->get("initialAmplitudes")) {
    amplitudes = arguments->getPtr<TensorUnion<F,TE>>(
      "initialAmplitudes"
    );
    OUT() << "Using given initial amplitudes " << amplitudes << std::endl;
    restart = true;
  } else {
    amplitudes = createAmplitudes<F,TE>(
      {{Nv,No}, {Nv,Nv,No,No}}, {"ai", "abij"}
    );
  }

  // TODO: conversion to eigen untis
  energy = New<MapNode>(SOURCE_LOCATION);
  energy->setValue<Real<>>("unit", eigenEnergies->getValue<Real<>>("unit"));

  // create a method handler, by default SinglesDoubles
  auto methodArguments(arguments->getMap("method"));
  auto methodType(
    methodArguments->getValue<std::string>("type", "SinglesDoubles")
  );
  Ptr<CoupledClusterMethod<F,TE>> method(
    CoupledClusterMethodFactory<F,TE>::create(methodType, arguments)
  );
  ASSERT_LOCATION(
    method, std::string("Unknown method type: '") + methodType + "'",
    methodArguments->get("type")->sourceLocation
  );
  OUT() << "Using method " << methodType << std::endl;

  // create a mixer, by default use the linear one
  auto mixerArguments(arguments->getMap("mixer"));
  auto mixerType(mixerArguments->getValue<std::string>("type", "LinearMixer"));
  Ptr<Mixer<F,TE>> mixer(MixerFactory<F,TE>::create(mixerType, mixerArguments));
  ASSERT_LOCATION(
    mixer, std::string("Unknown mixer type: '") + mixerType + "'",
    mixerArguments->get("type")->sourceLocation
  );
  // TODO: mixer should give state info
  if (mixerType.compare("DiisMixer") == 0 ){
    auto maxResidua(mixerArguments->getValue<std::string>("maxResidua", "4"));
    OUT() << "Using mixer " << mixerType << ", with maxResiua " << maxResidua << std::endl;
  }
  else if (mixerType.compare("LinearMixer") == 0){
    auto ratio(mixerArguments->getValue<Real<>>("ratio", 1.0));
    OUT() << "Using mixer " << mixerType << ", with ratio " << ratio << std::endl;
  }
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

  OUT() << "Maximum number of iterations: " << maxIterationsCount << std::endl;
  OUT() << "Unless reaching energy convergence dE: " << energyConvergence << std::endl;
  OUT() << "Or amplitudes convergence dR: " << amplitudesConvergence << std::endl;
  F e(0), previousE(0);
  char outstring[80];
  sprintf(outstring,"%4s %16s %11s %15s %6s\n",
          "Iter", "Energy  ", "dE   ", "dR      ", "time");
  OUT() << outstring;

  for (; i < maxIterationsCount; ++i) {
    auto startTime(Time::getCurrentRealTime());
    LOG() << "iteration: " << i+1 << std::endl;
    auto estimatedAmplitudes( method->getResiduum(i, restart, amplitudes) );
    estimateAmplitudesFromResiduum(estimatedAmplitudes, amplitudes);
    auto amplitudesChange( New<TensorUnion<F,TE>>(*estimatedAmplitudes) );
    *amplitudesChange -= *amplitudes;
    mixer->append(estimatedAmplitudes, amplitudesChange);
    // get mixer's best guess for amplitudes
    amplitudes = mixer->get();
    auto residuumNorm( mixer->getResiduumNorm());
    e = getEnergy(amplitudes);

    auto iterTime(Time::getCurrentRealTime() - startTime);

    sprintf(outstring,"%4ld %16.8f %12.4e %12.4e %8.1f\n",
            i+1, real(e), real(e - previousE), real(residuumNorm),
            iterTime.getFractionalSeconds());
    OUT() << outstring;
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
    OUT() << "computing energy from given amplitudes" << std::endl;
  } else if (i == maxIterationsCount) {
    WARNING_LOCATION(arguments->sourceLocation) <<
      "energy or amplitudes convergence not reached." << std::endl;
  }

  e = getEnergy(amplitudes, true);
  bool convergenceReached = i < maxIterationsCount;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  result->setValue<bool>("convergenceReached", convergenceReached);
  result->setPtr("amplitudes", amplitudes);
  return result;
}


template <typename F, typename TE>
F CoupledCluster::getEnergy(
  const Ptr<TensorUnion<F,TE>> &amplitudes,
  const bool finalReport
) {
  // get the Coulomb integrals to compute the energy
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vijab(coulombSlices->getPtr<TensorExpression<F,TE>>("hhpp"));
  auto orbitalType(
    coulombIntegrals->getMap(
      "indices"
    )->getMap("orbital")->getValue<std::string>("type")
  );
  Real<> spins;
  if (orbitalType == "spatial") {
    spins = 2;
  } else if (orbitalType == "spin") {
    spins = 1;
  } else {
    ASSERT_LOCATION(
      false, "unsupported orbital type '" + orbitalType + "'",
      coulombIntegrals->getMap(
        "indices"
      )->getMap("orbital")->get("type")->sourceLocation
    );
  }

  // singles amplitudes are optional
  auto Tai( amplitudes->get(0) );
  auto Tabij( amplitudes->get(1) );
  F e;
  std::streamsize ss(std::cout.precision());
  // TODO: antisymmetrized
/*
  if (antisymmetrized) {
    energy[""] += ( + 0.25  ) * (*Tabij)["abkl"] * (*Vijab)["klab"];
    energy[""] += ( + 0.5  ) * (*Tai)["aj"] * (*Tai)["cl"] * (*Vijab)["jlac"];
    e = energy.get_val();
    // FIXME: imaginary part ignored
    EMIT() << YAML::Key << "energy" << YAML::Value
      << YAML::BeginMap
      << YAML::Key << "value" << YAML::Value << real(e)
      << YAML::EndMap;
  } else
*/
  {
    auto direct( Tcc<TE>::template tensor<F>("D") );
    auto exchange( Tcc<TE>::template tensor<F>("X") );
    COMPILE(
      (*direct)[""] <<=
        0.5*spins*spins * (*Tabij)["abij"] * (*Vijab)["ijab"],
      (*direct)[""] +=
        0.5*spins*spins * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijab"],
      (*exchange)[""] <<=
        -0.5*spins * (*Tabij)["abij"] * (*Vijab)["ijba"],
      (*exchange)[""]
        += -0.5*spins * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijba"]
    )->execute();
    F D(direct->read());
    F X(exchange->read());
    e = D+X;
    if (finalReport){
      OUT() << std::endl;
      OUT() << "Total Energy: " << std::setprecision(10) << e << std::endl;
      OUT() << "Direct: "       << std::setprecision(10) << D << std::endl;
      OUT() << "Exchange: "     << std::setprecision(10) << X << std::endl;
    }
    energy->setValue<Real<>>("value", real<F>(e));
    energy->setValue<Real<>>("direct", real<F>(D));
    energy->setValue<Real<>>("exchange", real<F>(X));
    // TODO: energy units
  }



  std::cout << std::setprecision(ss);

  return e;
}

template <typename F, typename TE>
Ptr<TensorUnion<F,TE>> CoupledCluster::createAmplitudes(
  std::initializer_list<std::initializer_list<size_t>> amplitudeLens,
  std::initializer_list<std::string> amplitudeIndices
) {
  std::vector<Ptr<Tensor<F,TE>>> amplitudeTensors;
  for (auto lens: amplitudeLens) {
    amplitudeTensors.push_back(Tcc<TE>::template tensor<F>(lens, "T"));
  }
  return New<TensorUnion<F,TE>>(
    amplitudeTensors.begin(), amplitudeTensors.end(),
    amplitudeIndices.begin(), amplitudeIndices.end()
  );
}

template <typename F, typename TE>
void CoupledCluster::estimateAmplitudesFromResiduum(
  const Ptr<TensorUnion<F,TE>> &residuum,
  const Ptr<TensorUnion<F,TE>> &amplitudes
) {
  auto levelShift(
    arguments->getValue<Real<>>("levelShift", DEFAULT_LEVEL_SHIFT)
  );

  // apply level shifting on right hand side
  *residuum -= F(levelShift) * *amplitudes;

  for (unsigned int i(0); i < residuum->componentTensors.size(); ++i) {
    auto R( residuum->get(i) );
    auto indices( residuum->getIndices(i) );
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

// instantiate
template
void CoupledCluster::estimateAmplitudesFromResiduum(
  const Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>> &residuum,
  const Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>> &amplitudes
);
template
void CoupledCluster::estimateAmplitudesFromResiduum(
  const Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> &residuum,
  const Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> &amplitudes
);
template
void CoupledCluster::estimateAmplitudesFromResiduum(
  const Ptr<TensorUnion<Real<>, DefaultTensorEngine>> &residuum,
  const Ptr<TensorUnion<Real<>, DefaultTensorEngine>> &amplitudes
);
template
void CoupledCluster::estimateAmplitudesFromResiduum(
  const Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> &residuum,
  const Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> &amplitudes
);


template <typename F, typename TE>
Ptr<Tensor<F,TE>> CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
) {
  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<TensorExpression<Real<>,TE>>("h"));
  auto epsa(energySlices->getPtr<TensorExpression<Real<>,TE>>("p"));
  auto Fepsi(Tcc<TE>::template tensor<F>(epsi->inspect()->getLens(), "Fepsi"));
  auto Fepsa(Tcc<TE>::template tensor<F>(epsa->inspect()->getLens(), "Fepsa"));
  // convert to type F (either complex or double)
  COMPILE(
    (*Fepsa)["a"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"]),
    (*Fepsi)["i"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"])
  )->execute();

  auto D(
    Tcc<TE>::template tensor<F>(lens, std::string("D") + indices)
  );

  // create energy difference tensor
  int excitationLevel(indices.length()/2);
  for (int p(0); p < excitationLevel; ++p) {
    COMPILE(
      (*D)[indices] += (*Fepsa)[indices.substr(p,1)],
      (*D)[indices] -= (*Fepsi)[indices.substr(excitationLevel+p,1)]
    )->execute();
  }

  return D;
}

// instantiate
template
Ptr<Tensor<Real<>, DefaultDryTensorEngine>>
CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Complex<>, DefaultDryTensorEngine>>
CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Real<>, DefaultTensorEngine>>
CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Complex<>, DefaultTensorEngine>>
CoupledCluster::calculateEnergyDifferences(
  const std::vector<size_t> &lens, const std::string &indices
);


constexpr Real<64> CoupledCluster::DEFAULT_ENERGY_CONVERGENCE;
constexpr Real<64> CoupledCluster::DEFAULT_AMPLITUDES_CONVERGENCE;
constexpr Real<64> CoupledCluster::DEFAULT_LEVEL_SHIFT;

