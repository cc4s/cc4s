/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <CrossEntropyFtodRankDecomposition.hpp>
#include <Exception.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <iostream>
#include <limits>


using namespace CTF;

RankDecomposition::RankDecomposition(
  int rank, int np, int nG, World *world
):
  X(rank, np, NS, *world),
  GamR(rank, nG, NS, *world), GamI(rank, nG, NS, *world),
  residuum(0)
{
}

CrossEntropyFtodRankDecomposition::CrossEntropyFtodRankDecomposition(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
}

CrossEntropyFtodRankDecomposition::~CrossEntropyFtodRankDecomposition() {
}

void CrossEntropyFtodRankDecomposition::run() {
  rank = getIntegerArgument("rank");
  chiR = const_cast<Tensor<> *>(getTensorArgument("chiR"));
  chiI = const_cast<Tensor<> *>(getTensorArgument("chiI"));
  double epsilon = getRealArgument("epsilon");
  LOG(3) << "rank=" << rank << std::endl;
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  LOG(3) << "nG=" << nG << ", np=" << np << std::endl;
  // allocate decompisiton tensors
  X = new Matrix<>(int(rank), np, NS, *chiR->wrld, "XRp", chiR->profile);
  gamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamRRG", chiR->profile);
  gamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamIRG", chiR->profile);

  // allocate intermediate tensors:
  // XX
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  XX = new Tensor<>(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);
  // tensor rank approximation
  chi0R = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0RGqr", chiR->profile);
  chi0I = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0IGqr", chiR->profile);
  // residuum
  RR = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RRGqr", chiR->profile);
  RI = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RIGqr", chiR->profile);

  random.seed(chiR->wrld->rank);

  samplesCount = 1000;
  // 0.1 acceptance probability
  estimatorsCount = samplesCount/10;
  // allocate estimators
  estimators = new RankDecomposition *[estimatorsCount];
  for (int i(0); i < estimatorsCount; ++i) {
    estimators[i] = new RankDecomposition(rank, np, nG, chiR->wrld);
  }
  estimator = new RankDecomposition(rank, np, nG, chiR->wrld);
  mu = new RankDecomposition(rank, np, nG, chiR->wrld);
  sigma = new RankDecomposition(rank, np, nG, chiR->wrld);
  sigma->X["Rq"] += 1.0;
  sigma->GamR["RG"] += 100.0;
  sigma->GamI["RG"] += 100.0;

  for (int i(0); i < estimatorsCount; ++i) {
    setRandom(*estimators[i], *mu, *sigma);
    calculateChi0(*estimators[i]);
    calculateResiduum();
    estimators[i]->residuum = R;
  }
  double worstResiduum(findWorstEstimator());

  while (worstResiduum > epsilon) {
    for (int i(0); i < samplesCount; ++i) {
      LOG(4) << "  R=";
      setRandom(*estimator, *mu, *sigma);
      calculateChi0(*estimator);
      calculateResiduum();
      LOG(4) << R << std::endl;
      if (R < worstResiduum) {
        LOG(3) << "  better R=" << R << std::endl;
        // found better estimator
        estimator->residuum = R;
        RankDecomposition *previousWorst(estimators[worstEstimator]);
        estimators[worstEstimator] = estimator;
        estimator = previousWorst;
        worstResiduum = findWorstEstimator();
      }
    }
    calculateMu();
    LOG(3) << "  |sigma|^2=";
    double variance(calculateSigma());
    LOG(3) << variance << std::endl;
  }
}

// TODO: maybe use a queue for sorted access
double CrossEntropyFtodRankDecomposition::findWorstEstimator() {
  double worstResiduum(0);
  for (int i(0); i < estimatorsCount; ++i) {
    if (estimators[i]->residuum > worstResiduum) {
      worstEstimator = i;
      worstResiduum = R;
    }
  }
  return worstResiduum;
}

void CrossEntropyFtodRankDecomposition::calculateChi0(
  RankDecomposition &d
) {
  (*XX)["Rqr"] = d.X["Rq"] * d.X["Rr"];
  (*chi0R)["Gqr"] = (*XX)["Rqr"] * d.GamR["RG"];
  (*chi0I)["Gqr"] = (*XX)["Rqr"] * d.GamI["RG"];
}

void CrossEntropyFtodRankDecomposition::calculateResiduum() {
  // calculate the residuum
  (*RR)["Gqr"] = (*chi0R)["Gqr"] - (*chiR)["Gqr"];
  (*RI)["Gqr"] = (*chi0I)["Gqr"] - (*chiI)["Gqr"];
  Scalar<> s(*chiR->wrld);
  s[""] =  (*RR)["Gqr"] * (*RR)["Gqr"];
  s[""] += (*RI)["Gqr"] * (*RI)["Gqr"];
  R = s.get_val();
}

void CrossEntropyFtodRankDecomposition::setRandom(
  RankDecomposition &d, RankDecomposition &mu, RankDecomposition &sigma
) {
  setRandom(d.X, mu.X, sigma.X);
  setRandom(d.GamR, mu.GamR, sigma.GamR);
  setRandom(d.GamI, mu.GamI, sigma.GamI);
}

void CrossEntropyFtodRankDecomposition::setRandom(
  Tensor<> &t, Tensor<> &mu, Tensor<> &sigma
) {
  int64_t indicesCount, *indices;
  double *values;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  t.read_local(&indicesCount, &indices, &values);
  for (int64_t i(0); i < indicesCount; ++i) {
    values[i] = normalDistribution(random);
  }
  t.write(indicesCount, indices, values);
  free(indices); free(values);
  t["ab"] *= sigma["ab"];
  t["ab"] += mu["ab"];
}

// TODO: implement in RankDecomposition class
void CrossEntropyFtodRankDecomposition::calculateMu() {
  mu->X["Rq"] = estimators[0]->X["Rq"];
  mu->GamR["RG"] = estimators[0]->GamR["RG"];
  mu->GamI["RG"] = estimators[0]->GamI["RG"];

  for (int i(1); i < estimatorsCount; ++i) {
    mu->X["Rq"] += estimators[i]->X["Rq"];
    mu->GamR["RG"] += estimators[i]->GamR["RG"];
    mu->GamI["RG"] += estimators[i]->GamI["RG"];
  }
  // TODO: use division function for non-power-2 estimatorsCount
  mu->X["Rq"] *= 1.0 / estimatorsCount;
  mu->GamR["RG"] *= 1.0 / estimatorsCount;
  mu->GamI["RG"] *= 1.0 / estimatorsCount;
}

double CrossEntropyFtodRankDecomposition::calculateSigma() {
  sigma->X["Rq"] *= 0.0;
  sigma->GamR["RG"] *= 0.0;
  sigma->GamI["RG"] = 0.0;

  for (int i(0); i < estimatorsCount; ++i) {
    estimator->X["Rq"] = estimators[i]->X["Rq"];
    estimator->X["Rq"] -= mu->X["Rq"];
    estimator->X["Rq"] *= estimator->X["Rq"];
    estimator->GamR["RG"] = estimators[i]->GamR["RG"];
    estimator->GamR["RG"] -= mu->GamR["RG"];
    estimator->GamR["RG"] *= estimator->GamR["RG"];
    estimator->GamI["RG"] = estimators[i]->GamI["RG"];
    estimator->GamI["RG"] -= mu->GamI["RG"];
    estimator->GamI["RG"] *= estimator->GamI["RG"];
    sigma->X["Rq"] += estimator->X["Rq"];
    sigma->GamR["RG"] += estimator->GamR["RG"];
    sigma->GamI["RG"] += estimator->GamI["RG"];
  }
  // TODO: use division function for non-power-2 estimatorsCount
  Univar_Function<> fSqrt(&MathFunctions::sqrt<>);
  Scalar<> s(*chiR->wrld);
  sigma->X.sum(
    1.0/(estimatorsCount-1), sigma->X, "Rq",
    0.0, "Rq", fSqrt
  );
  s[""]  = sigma->X["Rq"] * sigma->X["Rq"];
  sigma->GamR.sum(
    1.0/(estimatorsCount-1), sigma->GamR, "RG",
    0.0, "RG", fSqrt
  );
  s[""] += sigma->GamR["RG"] * sigma->GamR["RG"];
  sigma->GamI.sum(
    1.0/(estimatorsCount-1), sigma->GamI, "RG",
    0.0, "RG", fSqrt
  );
  s[""] += sigma->GamI["RG"] * sigma->GamI["RG"];
  return s.get_val();
}

