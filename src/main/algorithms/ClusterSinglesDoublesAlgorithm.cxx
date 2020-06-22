#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <mixers/Mixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Options.hpp>
#include <Cc4s.hpp>

#include <array>
#include <initializer_list>

using namespace cc4s;

Ptr<MapNode> ClusterSinglesDoublesAlgorithm::run(const Ptr<MapNode> &arguments){
  this->arguments = arguments;
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto orbitals(coulombIntegrals->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (orbitals == "real64") {
    if (Cc4s::options->dryRun) {
      return run<Real<>,DryTensorEngine>();
    } else {
      return run<Real<>,DefaultTensorEngine>();
    }
  } else if (orbitals == "complex64") {
    if (Cc4s::options->dryRun) {
      return run<Complex<>,DryTensorEngine>();
    } else {
      return run<Complex<>,DefaultTensorEngine>();
    }
  } else {
    ASSERT_LOCATION(
      false, "unsupported orbitals type '" + orbitals + "'",
      coulombIntegrals->get("scalarType")->sourceLocation
    );
  }
}

template <typename F, typename TE>
Ptr<MapNode> ClusterSinglesDoublesAlgorithm::run() {
  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("h"));
  auto epsa(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("p"));

  auto No(epsi->getResult()->lens[0]);
  auto Nv(epsa->getResult()->lens[0]);

  Ptr<const FockVector<F,TE>> amplitudes(
    createAmplitudes<F,TE>(
      {"Singles", "Doubles"}, {{Nv,No}, {Nv,Nv,No,No}}, {"ai", "abij"}
    )
  );

  // TODO: conversion to eigen untis
  energy = New<MapNode>(SOURCE_LOCATION);
  energy->setValue<Real<>>("unit", eigenEnergies->getValue<Real<>>("unit"));

  // create a mixer, by default use the linear one
  auto mixerArguments(arguments->getMap("mixer"));
  auto mixerType(mixerArguments->getValue<std::string>("type", "LinearMixer"));
  Ptr<Mixer<F,TE>> mixer(MixerFactory<F,TE>::create(mixerType, mixerArguments));
  ASSERT_LOCATION(
    mixer, std::string("Unknown mixer type: '") + mixerType + "'",
    mixerArguments->get("type")->sourceLocation
  );

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


  F e(0), previousE(0);
  size_t i(0);
  for (; i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    // call the getResiduum of the actual algorithm,
    // which will be specified by inheriting classes
    auto estimatedAmplitudes( getResiduum(i, amplitudes) );
    estimateAmplitudesFromResiduum(estimatedAmplitudes, amplitudes);
    auto amplitudesChange( New<FockVector<F,TE>>(*estimatedAmplitudes) );
    *amplitudesChange -= *amplitudes;
    mixer->append(estimatedAmplitudes, amplitudesChange);
    // get mixer's best guess for amplitudes
    amplitudes = mixer->get();
    e = getEnergy(amplitudes);
    if (
      abs((e-previousE)/e) < energyConvergence &&
      abs(
        amplitudesChange->dot(*amplitudesChange) / amplitudes->dot(*amplitudes)
      ) < amplitudesConvergence * amplitudesConvergence
    ) {
      break;
    }
    previousE = e;
  }

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "computing energy from given amplitudes" << std::endl;
    e = getEnergy(amplitudes);
  } else if (i == maxIterationsCount) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "WARNING: energy or amplitudes convergence not reached." << std::endl;
  }
  bool convergenceReached = i < maxIterationsCount;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  result->setValue<bool>("convergenceReached", convergenceReached);
  result->get("amplitudes") = storeAmplitudes(arguments, amplitudes);
  return result;
}


template <typename F, typename TE>
F ClusterSinglesDoublesAlgorithm::getEnergy(
  const Ptr<const FockVector<F,TE>> &amplitudes
) {
  // get the Coulomb integrals to compute the energy
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vijab(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhpp"));
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
  std::streamsize ss = std::cout.precision();
  // TODO: antisymmetrized
/*
  if (antisymmetrized) {
    energy[""] += ( + 0.25  ) * (*Tabij)["abkl"] * (*Vijab)["klab"];
    energy[""] += ( + 0.5  ) * (*Tai)["aj"] * (*Tai)["cl"] * (*Vijab)["jlac"];
    e = energy.get_val();
    // FIXME: imaginary part ignored
    EMIT() << YAML::Key << "energy" << YAML::Value
      << YAML::BeginMap
      << YAML::Key << "value" << YAML::Value << std::real(e)
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
    F S(0.25*D - 0.5*X);
    F T(0.75*D + 1.5*X);
    e = D+X;
    LOG(1, getCapitalizedAbbreviation()) << std::setprecision(10) <<
      "dir= " << D << std::endl;
    LOG(1, getCapitalizedAbbreviation()) << std::setprecision(10) <<
      "exc= " << X << std::endl;
    LOG(1, getCapitalizedAbbreviation()) << std::setprecision(10) <<
      "sing= " << S << std::endl;
    LOG(1, getCapitalizedAbbreviation()) << std::setprecision(10) <<
      "trip= " << T << std::endl;

    energy->setValue<Real<>>("value", real<F>(e));
    energy->setValue<Real<>>("direct", real<F>(D));
    energy->setValue<Real<>>("exchange", real<F>(X));
    energy->setValue<Real<>>("singlet", real<F>(S));
    energy->setValue<Real<>>("triplet", real<F>(T));
  }

  LOG(0, getCapitalizedAbbreviation()) << std::setprecision(10) <<
    "energy= " << e << std::setprecision(ss) << std::endl;

  return e;
}

template <typename F, typename TE>
Ptr<FockVector<F,TE>> ClusterSinglesDoublesAlgorithm::createAmplitudes(
  std::initializer_list<std::string> amplitudeNames,
  std::initializer_list<std::initializer_list<size_t>> amplitudeLens,
  std::initializer_list<std::string> amplitudeIndices
) {
  std::vector<Ptr<Tensor<F,TE>>> amplitudeTensors;
  auto lensIterator( amplitudeLens.begin() );
  for (auto name: amplitudeNames) {
    std::stringstream initialDataName;
    initialDataName << "initial" << name << "Amplitudes";
// TODO: support initial amplitudes
/*
    if (isArgumentGiven(initialDataName.str())) {
      // use given amplitudes as initial amplitudes
      amplitudeTensors.push_back(
        NEW(CTF::Tensor<F>, *getTensorArgument<F>( initialDataName.str() ))
      );
    } else
*/
    {
      // otherwise, use zeros as initial amplitudes
      std::vector<size_t> lens(*lensIterator);
      amplitudeTensors.push_back(Tcc<TE>::template tensor<F>(lens, "T"));
    }
    ++lensIterator;
  }
  return New<FockVector<F,TE>>(
    amplitudeTensors.begin(), amplitudeTensors.end(),
    amplitudeIndices.begin(), amplitudeIndices.end()
  );
}


template <typename F, typename TE>
Ptr<MapNode> ClusterSinglesDoublesAlgorithm::storeAmplitudes(
  const Ptr<MapNode> &arguments,
  const Ptr<const FockVector<F,TE>> &amplitudes
) {
  // TODO: how should Fock-vectors be stored?
//  storeAmplitudes(amplitudes, {"Singles", "Doubles"});
  auto result(New<MapNode>(SOURCE_LOCATION));
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  result->get("scalarType") = coulombIntegrals->get("scalarType");
  result->setValue<Real<>>("unit", 1.0);
  result->get("indices") = coulombIntegrals->get("indices");
  auto components(New<MapNode>(SOURCE_LOCATION));
  components->get(0) = storeAmplitudesComponent(amplitudes->get(0));
  components->get(1) = storeAmplitudesComponent(amplitudes->get(1));
  result->get("components") = components;
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> ClusterSinglesDoublesAlgorithm::storeAmplitudesComponent(
  const Ptr<Tensor<F,TE>> &component
) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  auto dimensions(New<MapNode>(SOURCE_LOCATION));
  for (size_t d(0); d < component->lens.size(); ++d) {
    auto dimension(New<MapNode>(SOURCE_LOCATION));
    dimension->setValue<size_t>("length", component->lens[d]);
    dimension->setValue<std::string>("type", "orbital");
    dimensions->get(d) = dimension;
  }
  result->get("dimensions") = dimensions;
  result->setValue<Ptr<Tensor<F,TE>>>("data", component);
  return result;
}

template <typename F, typename TE>
void ClusterSinglesDoublesAlgorithm::estimateAmplitudesFromResiduum(
  const Ptr<FockVector<F,TE>> &residuum,
  const Ptr<const FockVector<F,TE>> &amplitudes
) {
  auto levelShift(
    arguments->getValue<Real<>>("levelShift", DEFAULT_LEVEL_SHIFT)
  );

  // apply level shifting on right hand side
  *residuum -= F(levelShift) * *amplitudes;

  for (unsigned int i(0); i < residuum->componentTensors.size(); ++i) {
    auto R( residuum->get(i) );
    auto indices( residuum->getIndices(i) );
    auto D( calculateExcitationEnergies<F,TE>(R->lens, indices) );

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
void ClusterSinglesDoublesAlgorithm::estimateAmplitudesFromResiduum(
  const Ptr<FockVector<Real<>, DryTensorEngine>> &residuum,
  const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
);
template
void ClusterSinglesDoublesAlgorithm::estimateAmplitudesFromResiduum(
  const Ptr<FockVector<Complex<>, DryTensorEngine>> &residuum,
  const Ptr<const FockVector<Complex<>, DryTensorEngine>> &amplitudes
);
template
void ClusterSinglesDoublesAlgorithm::estimateAmplitudesFromResiduum(
  const Ptr<FockVector<Real<>, DefaultTensorEngine>> &residuum,
  const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
);
template
void ClusterSinglesDoublesAlgorithm::estimateAmplitudesFromResiduum(
  const Ptr<FockVector<Complex<>, DefaultTensorEngine>> &residuum,
  const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
);


template <typename F, typename TE>
Ptr<Tensor<F,TE>> ClusterSinglesDoublesAlgorithm::calculateExcitationEnergies(
  const std::vector<size_t> &lens, const std::string &indices
) {
  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("h"));
  auto epsa(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("p"));
  auto Fepsi(Tcc<TE>::template tensor<F>(epsi->getResult()->lens, "Fepsi"));
  auto Fepsa(Tcc<TE>::template tensor<F>(epsa->getResult()->lens, "Fepsa"));
  // convert to type F (either complex or double)
  COMPILE(
    (*Fepsa)["a"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"]),
    (*Fepsi)["i"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"])
  )->execute();

  auto D(
    Tcc<TE>::template tensor<F>(lens, std::string("D") + indices)
  );

  // create excitation energy tensor
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
Ptr<Tensor<Real<>, DryTensorEngine>>
ClusterSinglesDoublesAlgorithm::calculateExcitationEnergies(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Complex<>, DryTensorEngine>>
ClusterSinglesDoublesAlgorithm::calculateExcitationEnergies(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Real<>, DefaultTensorEngine>>
ClusterSinglesDoublesAlgorithm::calculateExcitationEnergies(
  const std::vector<size_t> &lens, const std::string &indices
);
template
Ptr<Tensor<Complex<>, DefaultTensorEngine>>
ClusterSinglesDoublesAlgorithm::calculateExcitationEnergies(
  const std::vector<size_t> &lens, const std::string &indices
);



/*
Tensor<double> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(
  const Ptr<const FockVector<double>> &amplitudes,
  int a, int b, int integralsSliceSize
) {
  // Read the amplitudes Tai
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");

  // Read the Coulomb vertex GammaGqr
  auto GammaGqr( getTensorArgument<complex>("CoulombVertex"));
  GammaGqr->set_name("GammaGqr");

  // Compute No,Nv,NG,Np
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  int NG(GammaGqr->lens[0]);
  int Np = No + Nv;

  // Allocate and compute GammaGab,GammaGai from GammaGqr
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  int GabStart[] = {0 ,No,No};
  int GabEnd[]   = {NG,Np,Np};
  auto GammaGai(GammaGqr->slice(GaiStart,GaiEnd));
  auto GammaGab(GammaGqr->slice(GabStart,GabEnd));

  // Split GammaGab,GammaGai into real and imaginary parts
  Tensor<> realGammaGai(3, GammaGai.lens, GammaGai.sym,
                        *GammaGai.wrld, "RealGammaGai");
  Tensor<> imagGammaGai(3, GammaGai.lens, GammaGai.sym,
                        *GammaGai.wrld, "ImagGammaGai");
  fromComplexTensor(GammaGai, realGammaGai, imagGammaGai);

  Tensor<> realGammaGab(3, GammaGab.lens, GammaGab.sym,
                        *GammaGab.wrld, "RealGammaGab");
  Tensor<> imagGammaGab(3, GammaGab.lens, GammaGab.sym,
                        *GammaGab.wrld, "ImagGammaGab");
  fromComplexTensor(GammaGab, realGammaGab, imagGammaGab);

  // Construct dressed Coulomb vertex GammaGab
  realGammaGab["Gab"] += (-1.0) * realGammaGai["Gbk"] * (*Tai)["ak"];
  imagGammaGab["Gab"] += (-1.0) * imagGammaGai["Gbk"] * (*Tai)["ak"];
  toComplexTensor(realGammaGab, imagGammaGab, GammaGab);

  // Slice the respective parts from the dressed Coulomb vertex GammaGab
  int leftGammaStart[] = { 0, a, 0 };
  int leftGammaEnd[] = { NG, std::min(a+integralsSliceSize, Nv), Nv };
  int rightGammaStart[] = { 0, b, 0 };
  int rightGammaEnd[] = { NG, std::min(b+integralsSliceSize, Nv), Nv };

  auto leftGamma(GammaGab.slice(leftGammaStart, leftGammaEnd));
  auto rightGamma(GammaGab.slice(rightGammaStart, rightGammaEnd));

  // Split into real and imaginary parts
  Tensor<> realLeftGamma(
    3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld, "realLeftGamma"
  );
  Tensor<> imagLeftGamma(
    3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld, "imagLeftGamma"
  );
  fromComplexTensor(leftGamma, realLeftGamma, imagLeftGamma);
  Tensor<> realRightGamma(
    3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld, "realRightGamma"
  );
  Tensor<> imagRightGamma(
    3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld, "imagRightGamma"
  );
  fromComplexTensor(rightGamma, realRightGamma, imagRightGamma);

  // Allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  auto Vxycd(new Tensor<>(4, lens, syms, *GammaGqr->wrld, "Vxycd"));

  // Contract left and right slices of the dressed Coulomb vertices
  (*Vxycd)["xycd"]  = realLeftGamma["Gxc"] * realRightGamma["Gyd"];
  (*Vxycd)["xycd"] += imagLeftGamma["Gxc"] * imagRightGamma["Gyd"];
  return Vxycd;
}

Tensor<cc4s::complex> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(
  const Ptr<const FockVector<cc4s::complex>> &amplitudes,
  int a, int b, int integralsSliceSize
) {
  // Read the amplitudes Tai
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");

  // Read the Coulomb vertex GammaGqr
  auto GammaGqr( getTensorArgument<cc4s::complex>("CoulombVertex"));
  GammaGqr->set_name("GammaGqr");

  // Compute No,Nv,NG,Np
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  int NG(GammaGqr->lens[0]);
  int Np(No+Nv);

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int GiaStart[] = {0, iStart,aStart};
  int GiaEnd[]   = {NG,iEnd,  aEnd};
  int GaiStart[] = {0, aStart,iStart};
  int GaiEnd[]   = {NG,aEnd,  iEnd};
  int GabStart[] = {0, aStart,aStart};
  int GabEnd[]   = {NG,aEnd,  aEnd};
  auto GammaGia( new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd)) );
  auto GammaGai( new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd)) );
  auto GammaGab( new Tensor<complex>(GammaGqr->slice(GabStart, GabEnd)) );

  Univar_Function<complex> fConj(conj<complex>);

  Tensor<complex> conjTransposeGammaGia(false, *GammaGia);
  conjTransposeGammaGia.sum(1.0,*GammaGai,"Gai", 0.0,"Gia", fConj);
  Tensor<complex> conjTransposeGammaGab(false, *GammaGab);
  conjTransposeGammaGab.sum(1.0,*GammaGab,"Gba", 0.0,"Gab", fConj);

  // Construct dressed Coulomb vertex GammaGab
  Tensor<complex> DressedGammaGab(*GammaGab);
  DressedGammaGab.set_name("DressedGammaGab");
  DressedGammaGab["Gab"] += (-1.0) * (*GammaGia)["Gkb"] * (*Tai)["ak"];

  Tensor<complex> conjTransposeDressedGammaGab(conjTransposeGammaGab);
  conjTransposeDressedGammaGab.set_name("conjTransposeDressedGammaGab");
  conjTransposeDressedGammaGab["Gab"] += (-1.0) * conjTransposeGammaGia["Gkb"] * (*Tai)["ak"];

  // Slice the respective parts from the dressed Coulomb vertex GammaGab
  int leftGammaStart[] = { 0, a, 0 };
  int leftGammaEnd[] = { NG, std::min(a+integralsSliceSize, Nv), Nv };
  int rightGammaStart[] = { 0, b, 0 };
  int rightGammaEnd[] = { NG, std::min(b+integralsSliceSize, Nv), Nv };

  auto leftGamma(conjTransposeDressedGammaGab.slice(leftGammaStart, leftGammaEnd));
  auto rightGamma(DressedGammaGab.slice(rightGammaStart, rightGammaEnd));

  // Allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  auto Vxycd(new Tensor<complex>(4, lens, syms, *GammaGqr->wrld, "Vxycd"));

  // Contract left and right slices of the dressed Coulomb vertices
  (*Vxycd)["xycd"]  = leftGamma["Gxc"] * rightGamma["Gyd"];
  delete GammaGia;
  delete GammaGai;
  delete GammaGab;
  return Vxycd;
}


Tensor<double> *
  ClusterSinglesDoublesAlgorithm::sliceAmplitudesFromCoupledCoulombFactors
(
  const Ptr<const FockVector<double>> &amplitudes,
  int a, int b, int factorsSliceSize
) {
  auto PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  auto LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  auto epsi(getTensorArgument("HoleEigenEnergies"));
  auto epsa(getTensorArgument("ParticleEigenEnergies"));

  // Read the doubles amplitudes Tabij
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");
  auto Tabij( amplitudes->get(1) );
  Tabij->set_name("Tabij");

  // Intermediate tensor Iabij=T2+T1*T1
  auto Iabij(*Tabij);
  Iabij.set_name("Iabij");
  Iabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int Rx(std::min(factorsSliceSize, NR-a));
  int Ry(std::min(factorsSliceSize, NR-b));
  int Rvoo[] = { Rx, Nv, No, No };
  int RRoo[] = { Rx, Ry, No, No };
  int RR[] = { Rx, Ry };
  int syms[] = { NS, NS, NS, NS };

  Tensor<complex> VRS(2, RR, syms, *PirR->wrld, "VRS");

  Tensor<> realXRaij(4, Rvoo, syms, *PirR->wrld, "RealXRaij");
  Tensor<> imagXRaij(4, Rvoo, syms, *PirR->wrld, "ImagXRaij");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Slice the respective parts from PiaR
  int leftPiStart[]  = { 0 ,                                a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };
  int rightPiStart[] = { 0 ,                                b };
  int rightPiEnd[]   = { Nv, std::min(b+factorsSliceSize, NR) };

  auto leftPiaR (PiaR.slice(leftPiStart  ,  leftPiEnd));
  leftPiaR.set_name("leftPiaR");
  auto rightPiaR(PiaR.slice(rightPiStart , rightPiEnd));
  rightPiaR.set_name("rightPiaR");

  // Split left and right PiaR into real and imaginary parts
  Tensor<> realLeftPiaR(
    2, leftPiaR.lens, leftPiaR.sym, *leftPiaR.wrld, "RealLeftPiaR"
  );
  Tensor<> imagLeftPiaR(
    2, leftPiaR.lens, leftPiaR.sym, *leftPiaR.wrld, "ImagRightPiaR"
  );
  fromComplexTensor(leftPiaR, realLeftPiaR, imagLeftPiaR);

  Tensor<> realRightPiaR(
    2, rightPiaR.lens, rightPiaR.sym, *rightPiaR.wrld, "RealLeftPiaR"
  );
  Tensor<> imagRightPiaR(
    2, rightPiaR.lens, rightPiaR.sym, *rightPiaR.wrld, "ImagRightPiaR"
  );
  fromComplexTensor(leftPiaR, realLeftPiaR, imagLeftPiaR);

  // Slice the respective parts from LambdaGR
  int leftLambdaStart[]  = { 0  ,                            a };
  int leftLambdaEnd[]    = { NG , std::min(a+factorsSliceSize, NR) };
  auto leftLambdaGR(LambdaGR->slice(leftLambdaStart, leftLambdaEnd));
  leftLambdaGR.set_name("leftLambdaGR");

  int rightLambdaStart[]  = { 0  ,                            b };
  int rightLambdaEnd[]    = { NG , std::min(b+factorsSliceSize, NR) };
  auto rightLambdaGR(LambdaGR->slice(rightLambdaStart, rightLambdaEnd));
  rightLambdaGR.set_name("rightLambdaGR");

  // TODO: specify how the vertex should be computed
  // assuming GammaGqr = PirR*PirR*LambdaGR (first Pi not conjugated)
  realXRaij["Rdij"] = (+1.0) * Iabij["cdij"] * realLeftPiaR["cR"];
  imagXRaij["Rdij"] = (-1.0) * Iabij["cdij"] * imagLeftPiaR["cR"];
  Tensor<complex> XRaij(4, Rvoo, syms, *PirR->wrld, "XRaij");
  toComplexTensor(realXRaij, imagXRaij, XRaij);

  Tensor<complex> XRSij(4, RRoo, syms, *PirR->wrld, "XRSij");
  XRSij["RSij"] = XRaij["Rdij"] * rightPiaR["dS"];

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);
  VRS["RS"] = conjLeftLambdaGR["GR"] * rightLambdaGR["GS"];

  XRSij["RSij"] = XRSij["RSij"]  * VRS["RS"];

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  auto PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");

  // Split PiiR into real and imaginary parts
  Tensor<> realPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "RealPiiR");
  Tensor<> imagPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "ImagPiiR");
  fromComplexTensor(PiiR, realPiiR, imagPiiR);

  // Initialize dressedPiaR
  auto dressedPiaR(PiaR);
  dressedPiaR.set_name("dressedPiaR");

  // Split dressedPiaR into real and imaginary parts
  Tensor<> realDressedPiaR(
    2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "RealDressedPiaR"
  );
  Tensor<> imagDressedPiaR(
    2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "ImagDressedPiaR"
  );
  fromComplexTensor(dressedPiaR, realDressedPiaR, imagDressedPiaR);

  // Construct dressedPiaR
  realDressedPiaR["aR"] += (-1.0) * realPiiR["kR"] * (*Tai)["ak"];
  imagDressedPiaR["aR"] += (-1.0) * imagPiiR["kR"] * (*Tai)["ak"];
  toComplexTensor(realDressedPiaR, imagDressedPiaR, dressedPiaR);

  // Slice the respective parts from dressedPiaR
  auto dressedLeftPiaR (dressedPiaR.slice(leftPiStart  ,  leftPiEnd));
  dressedLeftPiaR.set_name("dressedLeftPiaR");
  auto dressedRightPiaR(dressedPiaR.slice(rightPiStart , rightPiEnd));
  dressedRightPiaR.set_name("dressedRightPiaR");

  // Split dressed left PiaR into real and imaginary parts
  Tensor<> dressedRealLeftPiaR(
    2, dressedLeftPiaR.lens, dressedLeftPiaR.sym, *dressedLeftPiaR.wrld,
    "dressedRealLeftPiaR"
  );
  Tensor<> dressedImagLeftPiaR(
    2, dressedLeftPiaR.lens, dressedLeftPiaR.sym, *dressedLeftPiaR.wrld,
    "dressedImagLeftPiaR"
  );
  fromComplexTensor(dressedLeftPiaR, dressedRealLeftPiaR, dressedImagLeftPiaR);

  XRaij["Rbij"] = XRSij["RSij"]  * dressedRightPiaR["bS"];

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  auto Fabij(new Tensor<>(4, vvoo, syms, *PirR->wrld, "Fabij"));

  // compute sliced amplitudes
  fromComplexTensor(XRaij, realXRaij, imagXRaij);
  (*Fabij)["abij"]  = realXRaij["Rbij"]  * dressedRealLeftPiaR["aR"];
  (*Fabij)["abij"] += imagXRaij["Rbij"]  * dressedImagLeftPiaR["aR"];

  // return sliced amplitudes
  return Fabij;
}

Tensor<cc4s::complex> *
  ClusterSinglesDoublesAlgorithm::sliceAmplitudesFromCoupledCoulombFactors
(
  const Ptr<const FockVector<complex>> &amplitudes,
  int a, int b, int factorsSliceSize
) {
  auto PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  auto PiqR(getTensorArgument<complex>("OutgoingFactorOrbitals"));
  PiqR->set_name("PiqR");
  auto LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  auto epsi(getTensorArgument("HoleEigenEnergies"));
  auto epsa(getTensorArgument("ParticleEigenEnergies"));

  // Read the doubles amplitudes Tabij
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");
  auto Tabij( amplitudes->get(1) );
  Tabij->set_name("Tabij");

  // Intermediate tensor Iabij=T2+T1*T1
  auto Iabij(*Tabij);
  Iabij.set_name("Iabij");
  Iabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int Rx(std::min(factorsSliceSize, NR-a));
  int Ry(std::min(factorsSliceSize, NR-b));
  int Rvoo[] = { Rx, Nv, No, No };
  int RRoo[] = { Rx, Ry, No, No };
  int RR[] = { Rx, Ry };
  int syms[] = { NS, NS, NS, NS };

  Univar_Function<complex> fConj(&cc4s::conj<complex>);

  Tensor<complex> VRS(2, RR, syms, *PirR->wrld, "VRS");

  Tensor<complex> XRaij(4, Rvoo, syms, *PirR->wrld, "XRaij");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  auto PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");
  auto PicR(PiqR->slice(aRStart,aREnd));
  PicR.set_name("PicR");

  Tensor<complex> conjPiaR(false, PiaR);
  conjPiaR.set_name("ConjPiaR");
  conjPiaR.sum(1.0, PiaR,"aR", 0.0,"aR", fConj);
  Tensor<complex> conjPicR(false, PicR);
  conjPicR.set_name("ConjPicR");
  conjPicR.sum(1.0, PicR,"aR", 0.0,"aR", fConj);

  // Slice the respective parts from PiaR
  int leftPiStart[]  = { 0 ,                                a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };
  int rightPiStart[] = { 0 ,                                b };
  int rightPiEnd[]   = { Nv, std::min(b+factorsSliceSize, NR) };

  auto leftPiaR(conjPicR.slice(leftPiStart  ,  leftPiEnd));
  leftPiaR.set_name("leftPiaR");
  auto rightPiaR(PiaR.slice(rightPiStart , rightPiEnd));
  rightPiaR.set_name("rightPiaR");

  // Slice the respective parts from LambdaGR
  int leftLambdaStart[]  = { 0  ,                                a };
  int leftLambdaEnd[]    = { NG , std::min(a+factorsSliceSize, NR) };
  auto leftLambdaGR (LambdaGR->slice(leftLambdaStart , leftLambdaEnd));
  leftLambdaGR.set_name("leftLambdaGR");

  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);

  int rightLambdaStart[]  = { 0  ,                                b };
  int rightLambdaEnd[]    = { NG , std::min(b+factorsSliceSize, NR) };
  auto rightLambdaGR (LambdaGR->slice(rightLambdaStart , rightLambdaEnd));
  rightLambdaGR.set_name("rightLambdaGR");

  // TODO: specify how the vertex should be computed
  // assuming GammaGqr = (PiqR*)*(PirR)*(LambdaGR) (first Pi conjugated)
  XRaij["Rdij"] = (+1.0) * Iabij["cdij"] * leftPiaR["cR"];

  Tensor<complex> XRSij(4, RRoo, syms, *PirR->wrld, "XRSij");
  XRSij["RSij"] = XRaij["Rdij"] * rightPiaR["dS"];

  VRS["RS"] = conjLeftLambdaGR["GR"] * rightLambdaGR["GS"];

  XRSij["RSij"] = XRSij["RSij"]  * VRS["RS"];

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  auto PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");
  auto PijR(PiqR->slice(iRStart,iREnd));
  PijR.set_name("PijR");
  Tensor<complex> conjPiiR(false, PiiR);
  conjPiiR.set_name("ConjPiiR");
  conjPiiR.sum(1.0, PiiR,"iR", 0.0,"iR", fConj);
  Tensor<complex> conjPijR(false, PijR);
  conjPijR.set_name("ConjPijR");
  conjPijR.sum(1.0, PijR,"iR", 0.0,"iR", fConj);

  // Construct dressedPiaR
  auto dressedPiaR(PicR);
  dressedPiaR.set_name("dressedPiaR");
  dressedPiaR["aR"] += (-1.0) * PijR["kR"] * (*Tai)["ak"];

  auto conjDressedPiaR(conjPiaR);
  conjDressedPiaR.set_name("conjDressedPiaR");
  conjDressedPiaR["aR"] += (-1.0) * conjPiiR["kR"] * (*Tai)["ak"];

  // Slice the respective parts from dressedPiaR
  auto dressedLeftPiaR (conjDressedPiaR.slice(leftPiStart  ,  leftPiEnd));
  dressedLeftPiaR.set_name("dressedLeftPiaR");

  auto dressedRightPiaR(dressedPiaR.slice(rightPiStart , rightPiEnd));
  dressedRightPiaR.set_name("dressedrightPiaR");

  XRaij["Rbij"] = XRSij["RSij"]  * dressedRightPiaR["bS"];

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  auto Fabij(new Tensor<complex>(4, vvoo, syms, *PirR->wrld, "Fabij"));

  // compute sliced amplitudes
  (*Fabij)["abij"]  = XRaij["Rbij"]  * dressedLeftPiaR["aR"];

  // return sliced amplitudes
  return Fabij;
}


template <typename F>
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<F> &Rxyij, int a, int b, Tensor<F> &Rabij
) {
  int Nx(Rxyij.lens[0]);
  int Ny(Rxyij.lens[1]);
  int No(Rxyij.lens[2]);
  int dstStart[] = { a, b, 0, 0 };
  int dstEnd[] = { a+Nx, b+Ny, No, No };
  int srcStart[] = { 0, 0, 0, 0 };
  int srcEnd[] = { Nx, Ny, No, No };
  // R["abij"] += R["xyij"] at current x,y
  Rabij.slice(dstStart,dstEnd,1.0, Rxyij,srcStart,srcEnd,1.0);
  if (a>b) {
    // Add the same slice at (b,a,j,i):
    dstStart[0] = b; dstStart[1] = a;
    dstEnd[0] = b+Ny; dstEnd[1] = a+Nx;
    srcEnd[0] = Ny; srcEnd[1] = Nx;
    // Swap xy and ij simultaneously
    Tensor<F> Ryxji(4, srcEnd, Rxyij.sym, *Rxyij.wrld, "Ryxji");
    Ryxji["yxji"] = Rxyij["xyij"];
    // Add Ryxij to Rabij
    Rabij.slice(dstStart,dstEnd,1.0, Ryxji,srcStart,srcEnd,1.0);
  }
}

// instantiate:
template
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<double> &Rxyij, int a, int b, Tensor<double> &Rabij
);
template
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<complex> &Rxyij, int a, int b, Tensor<complex> &Rabij
);
*/

std::string ClusterSinglesDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string capitalizedAbbreviation(getAbbreviation());
  std::transform(
    capitalizedAbbreviation.begin(), capitalizedAbbreviation.end(),
    capitalizedAbbreviation.begin(), ::toupper
  );
  return capitalizedAbbreviation;
}


std::string ClusterSinglesDoublesAlgorithm::getDataName(
  const std::string &type, const std::string &data
) {
  std::stringstream dataName;
  dataName << getAbbreviation() << type << data;
  return dataName.str();
}

constexpr double ClusterSinglesDoublesAlgorithm::DEFAULT_ENERGY_CONVERGENCE;
constexpr double ClusterSinglesDoublesAlgorithm::DEFAULT_AMPLITUDES_CONVERGENCE;
constexpr double ClusterSinglesDoublesAlgorithm::DEFAULT_LEVEL_SHIFT;

