#include <algorithms/Mp2ImaginaryFrequencyGrid.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>

#include <ctf.hpp>
#include <algorithm>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Mp2ImaginaryFrequencyGrid);

Mp2ImaginaryFrequencyGrid::Mp2ImaginaryFrequencyGrid(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2ImaginaryFrequencyGrid::~Mp2ImaginaryFrequencyGrid() {
}

void Mp2ImaginaryFrequencyGrid::run() {
  auto epsi(getTensorArgument<>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  Dai = NEW(CTF::Tensor<double>, 2, std::vector<int>({Nv,No}).data());
  (*Dai)["ai"]  = (*epsa)["a"];
  (*Dai)["ai"] -= (*epsi)["i"];

  if (isArgumentGiven("ImaginaryFrequencyPoints")) {
    auto nun(getTensorArgument<>("ImaginaryFrequencyPoints"));
    auto wn(getTensorArgument<>("ImaginaryFrequencyWeights"));
    grid = IntegrationGrid(nun->lens[0]);
    nun->read_all(grid.points.data(), true);
    wn->read_all(grid.weights.data(), true);
  } else {
    size_t N(getIntegerArgument("imaginaryFrequencies", 6));
    grid = IntegrationGrid(N);
    std::vector<double> deltas(No*Nv);
    Dai->read_all(deltas.data(), true);
    std::sort(deltas.begin(), deltas.end());

    // scale all Deltas to the range [1,R]
    scale = deltas[0];
    LOG(1, "Mp2Grid") << "scale=" << scale << std::endl;
    for (size_t i(0); i < deltas.size(); ++i) {
//      deltas[i] /= scale;
    }
//    (*Dai)["ai"] *= 1.0 / scale;

    // first half: use equidistant grid for head ~ 1/Delta^2 in [0,Delta_0]
    size_t n(0);
    for (; n < N/2; ++n) {
      grid.points[n] = scale * (n+0.5) / (N/2);
    }
    // second half: use 1/nu^(1/3) grid for tail ~ Delta^2/nu^4 in [Delta_0,inf]
    // where Delta_0 = min(eps_a-eps_i)
    for (; n < N; ++n) {
      grid.points[n] = scale * std::pow((N-n+0.5) / (N-N/2), -1.0/3);
    }

    // initialize weights
    double lastNu(0);
    for (size_t n(0); n < N; ++n) {
      grid.weights[n] = grid.points[n] - lastNu;
      lastNu = grid.points[n];
    }
  }
  Eai = NEW(CTF::Tensor<double>, false, *Dai);

  optimize(getIntegerArgument("stepCount", 1024));

  size_t N(grid.points.size());
  std::vector<int64_t> indices(N);
  for (size_t n(0); n < indices.size(); ++n) { indices[n] = n; }

  auto nun( new CTF::Tensor<double>(1, std::vector<int>{int(N)}.data()) );
  auto wn( new CTF::Tensor<double>(1, std::vector<int>{int(N)}.data()) );
  int64_t count(nun->wrld->rank == 0 ? N : 0);
  nun->write(count, indices.data(), grid.points.data());
  wn->write(count, indices.data(), grid.weights.data());

  allocatedTensorArgument<>("ImaginaryFrequencyPoints", nun);
  allocatedTensorArgument<>("ImaginaryFrequencyWeights", wn);
}

void Mp2ImaginaryFrequencyGrid::optimize(const int stepCount) {
  LOG(1, "Mp2Grid") << "optimizing grid for " << grid.points.size() <<
    " points" << std::endl;
  double E(getError(grid));
  NEW_FILE("E.dat") << 0 << " " << E << std::endl;
  NEW_FILE("nu.dat");
  NEW_FILE("w.dat");
  writeGrid(0);
  IntegrationGrid lastDelta(-getGradient(grid));
  IntegrationGrid lastDirection(lastDelta);
  applyConstraints(lastDelta);
  applyConstraints(lastDirection);
  for (int m(0); m < stepCount; ++m) {
    IntegrationGrid Delta(-getGradient(grid));
    if (Delta.dot(Delta) < 1e-16) break;
    applyConstraints(Delta);
    double beta(
      std::max(0.0, Delta.dot(Delta-lastDelta) / lastDelta.dot(lastDelta))
    );
    IntegrationGrid direction(Delta + beta*lastDirection);
    E = lineSearch(grid, direction);
    FILE("E.dat") << m+1 << " " << E << std::endl;
    writeGrid(m+1);
    if (m % 100 == 0) {
      LOG(1, "Mp2Grid") <<
        "iteration " << m+1 << ": error=" << E << ", beta=" << beta << std::endl;
    }
    lastDirection = direction;
    lastDelta = Delta;
  }
  LOG(0, "Mp2Grid") << "error=" << E << std::endl;

  // scale the grid back from the range [1,R] to the original range
  // [min(eps_a-eps_i),max(eps_a-eps-i)]
  for (size_t n(0); n < grid.points.size(); ++n) {
//    grid.points[n] *= scale;
//    grid.weights[n] *= scale;
  }
}

void Mp2ImaginaryFrequencyGrid::applyConstraints(IntegrationGrid &direction) {
  for (int n(0); n < getIntegerArgument("firstPoint", 0); ++n) {
    direction.points[n] = 0.0;
    direction.weights[n] = 0.0;
  }
}

double Mp2ImaginaryFrequencyGrid::lineSearch(
  IntegrationGrid &grid, const IntegrationGrid &direction
) {
  // starting vectors
  double alpha0(0.0);
  double alpha1(1.0);
  double E0(getError(grid));
  double E1(getError(grid + direction));
  // increase alpha until error also increases
  while (E1 <= E0) {
    alpha1 *= 2.0;
    E1 = getError(grid + alpha1*direction);
  };
  // interval bisection to machine precision
  double alpham(0.5*(alpha0+alpha1));
  double Em(getError(grid + alpham*direction));
  while (alpham-alpha0 > 1e-16 && alpha1-alpham > 1e-16) {
    if (E0 < E1) {
      alpha1 = alpham;
      E1 = Em;
    } else {
      alpha0 = alpham;
      E0 = Em;
    }
    alpham = 0.5*(alpha0+alpha1);
    Em = getError(grid + alpham*direction);
  };
  grid += alpham*direction;
  return Em;
}

double Mp2ImaginaryFrequencyGrid::gradientLineSearch(
  IntegrationGrid &grid, const IntegrationGrid &direction
) {
  // starting vectors
  double alpha0(0.0);
  double alpha1(1.0);
  IntegrationGrid newGrid(grid + direction);
  double E(getError(newGrid));
  IntegrationGrid grad(getGradient(newGrid));
  // increase alpha until gradient also increases
  while (direction.dot(grad) < 0) {
    alpha1 *= 2.0;
    newGrid = grid + alpha1*direction;
    E = getError(newGrid);
    grad = getGradient(newGrid);
//    LOG(1, "Mp2Grid") << "increasing to " << alpha1  << std::endl;
  };
  // interval bisection to machine precision
  double alpham(0.5*(alpha0+alpha1));
  newGrid = grid + alpham*direction;
  E = getError(newGrid);
  grad = getGradient(newGrid);
  while (alpham-alpha0 > 1e-16 && alpha1-alpham > 1e-16) {
    if (direction.dot(grad) < 0) {
      alpha0 = alpham;
    } else {
      alpha1 = alpham;
    }
    alpham = 0.5*(alpha0+alpha1);
    newGrid = grid + alpham*direction;
    E = getError(newGrid);
    grad = getGradient(newGrid);
//    LOG(1, "Mp2Grid") << "bisecting to " << alpham << std::endl;
  };
  grid = newGrid;
  return E;
}

void Mp2ImaginaryFrequencyGrid::testGradient(const double stepSize) {
  double E0(getError(grid));
  IntegrationGrid grad(getGradient(grid));
  IntegrationGrid numGrad(grid);
  for (size_t n(0); n < grid.points.size(); ++n) {
    grid.weights[n] += stepSize;
    double E1(getError(grid));
    grid.weights[n] -= stepSize;
    numGrad.weights[n] = (E1-E0) / stepSize;
    LOG(1, "Mp2Grid") << "dE/dw[" << n << "] = " <<
      grad.weights[n] << " ~ " << numGrad.weights[n] << std::endl;
  }
  for (size_t n(0); n < grid.points.size(); ++n) {
    grid.points[n] += stepSize;
    double E1(getError(grid));
    grid.points[n] -= stepSize;
    numGrad.points[n] = (E1-E0) / stepSize;
    LOG(1, "Mp2Grid") << "dE/dnu[" << n << "] = " <<
      grad.points[n] << " ~ " << numGrad.points[n] << std::endl;
  }
}

double Mp2ImaginaryFrequencyGrid::getError(
  const IntegrationGrid &grid
) {
  // error of numerical quadrature with current grid wn[n] & points[n] from
  // analytic value 1/(eps_a-eps_i) for each ai
  CTF::Transform<double, double>(
    std::function<void(double, double &)>(
      [&grid](double delta, double &error) {
        error = 0.0;
        // use quadrature to integrate Mp2 propagator
        for (size_t n(0); n < grid.points.size(); ++n) {
          error += grid.weights[n] * propagator(delta, grid.points[n]);
        }
        error /= Pi();
        // subtract analytic result = 1/delta
        error -= 1.0 / delta;
      }
    )
  ) (
    (*Dai)["ai"], (*Eai)["ai"]
  );

  // get frobenius norm of errors, which is to be optimized
  CTF::Scalar<double> E;
  E[""] = 0.5 * (*Eai)["ai"] * (*Eai)["ai"];
  return E.get_val();
}

// expects a call of getError with the same grid first
IntegrationGrid Mp2ImaginaryFrequencyGrid::getGradient(
  const IntegrationGrid &grid
) {
  IntegrationGrid gradGrid(grid);
  for (size_t n(0); n < grid.points.size(); ++n) {
    // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
    // integration weight weights[n] for each ai
    auto DWai(*Dai);
    CTF::Transform<double, double>(
      std::function<void(double, double &)>(
        [&grid, n](double error, double &delta) {
          double diff(error * propagator(delta, grid.points[n]) / Pi());
          delta = diff;
        }
      )
    ) (
      (*Eai)["ai"], DWai["ai"]
    );
    // sum contributions from all energy gaps eps_a-eps_i
    CTF::Scalar<double> DW;
    DW[""] = DWai["ai"];
    gradGrid.weights[n] = DW.get_val();

    // derivative of 0.5*error(eps_a-eps_i)^2 with respect to n-th
    // integration point points[n] for each ai
    auto DNuai(*Dai);
    CTF::Transform<double, double>(
      std::function<void(double, double &)>(
        [&grid, n](double error, double &delta) {
          double diff(
            error *
            (-2) * grid.weights[n] * propagator(delta, grid.points[n], 3) *
            2*grid.points[n]
          );
          delta = diff / Pi();
        }
      )
    ) (
      (*Eai)["ai"], DNuai["ai"]
    );
    // sum contributions from all energy gaps eps_a-eps_i
    CTF::Scalar<double> DNu;
    DNu[""] = DNuai["ai"];
    gradGrid.points[n] = DNu.get_val();
  }
  return gradGrid;
}

void Mp2ImaginaryFrequencyGrid::writeGrid(const int m) {
  FILE("nu.dat") << m;
  FILE("w.dat") << m;
  for (size_t n(0); n < grid.points.size(); ++n) {
    FILE("nu.dat") << " " << grid.points[n];
    FILE("w.dat") << " " << grid.weights[n];
  }
  FILE("nu.dat") << std::endl;
  FILE("w.dat") << std::endl;
}

