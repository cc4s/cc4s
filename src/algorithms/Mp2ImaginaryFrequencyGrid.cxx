#include <algorithms/Mp2ImaginaryFrequencyGrid.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

#include <ctf.hpp>
#include <algorithm>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(Mp2ImaginaryFrequencyGrid);

Mp2ImaginaryFrequencyGrid::Mp2ImaginaryFrequencyGrid(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2ImaginaryFrequencyGrid::~Mp2ImaginaryFrequencyGrid() {
}

void Mp2ImaginaryFrequencyGrid::run() {
  Mp2ImaginaryFrequencyGridOptimizer optimizer(
    getIntegerArgument("imaginaryFrequencies", 6),
    *getTensorArgument<>("HoleEigenEnergies"),
    *getTensorArgument<>("ParticleEigenEnergies")
  );
  optimizer.optimize(getIntegerArgument("stepCount", 1024));
}

Mp2ImaginaryFrequencyGridOptimizer::Mp2ImaginaryFrequencyGridOptimizer(
  const int N, CTF::Tensor<double> epsi, CTF::Tensor<double> epsa
):
  grid(N)
{
  int No(epsi.lens[0]);
  int Nv(epsa.lens[0]);
  Dai = new CTF::Tensor<double>(2, std::vector<int>({Nv,No}).data());
  (*Dai)["ai"]  = epsa["a"];
  (*Dai)["ai"] -= epsi["i"];

  // take N quantiles as initial estimates for grid points
  std::vector<double> deltas(No*Nv);
  Dai->read_all(deltas.data(), true);
  std::sort(deltas.begin(), deltas.end());
  double lastNu(0);
  for (size_t n(0); n < grid.points.size(); ++n) {
    grid.points[n] = deltas[No*Nv*n/N + 0*No*Nv/(2*N)];
    grid.weights[n] = grid.points[n] - lastNu;
    lastNu = grid.points[n];
  }
  Eai = new CTF::Tensor<double>(false, *Dai);
}

void Mp2ImaginaryFrequencyGridOptimizer::optimize(const int stepCount) {
  LOG(1, "RPA") << "optimizing grid" << std::endl;
  double E(getError(grid));
  IntegrationGrid lastDelta(-getGradient(grid));
  IntegrationGrid lastDirection(lastDelta);
  E = lineSearch(grid, lastDirection);
  NEW_FILE("E.dat");
  NEW_FILE("nu.dat");
  NEW_FILE("w.dat");
  for (int m(0); m < stepCount; ++m) {
    IntegrationGrid Delta(-getGradient(grid));
    double beta(
      std::max(0.0, Delta.dot(Delta-lastDelta) / lastDelta.dot(lastDelta))
    );
    IntegrationGrid direction(Delta + beta*lastDirection);
    E = lineSearch(grid, direction);
    FILE("E.dat") << m << " " << E << std::endl;
    FILE("nu.dat") << m;
    FILE("w.dat") << m;
    for (size_t n(0); n < grid.points.size(); ++n) {
      FILE("nu.dat") << " " << grid.points[n];
      FILE("w.dat") << " " << grid.weights[n];
    }
    FILE("nu.dat") << std::endl;
    FILE("w.dat") << std::endl;
    LOG(1, "RPA") << "error=" << E << ", beta=" << beta << std::endl;
    lastDirection = direction;
    lastDelta = Delta;
  }
}

double Mp2ImaginaryFrequencyGridOptimizer::lineSearch(
  IntegrationGrid &grid, const IntegrationGrid &direction
) {
  // starting vectors
  double E0(getError(grid));
  double E1(getError(grid + direction));
  double Em;
  double alpham(0.5);
  double alpha1(1.0);
  // increase alpha until error also increases
  while (E1 <= E0) {
    alpham = alpha1;
    alpha1 *= 2.0;
    Em = E1;
    E1 = getError(grid + alpha1*direction);
  } while (E1 <= E0);
  double alpha0(0.0);
  // interval bisection to machine precision
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
  }
  grid += alpham*direction;
  return Em;
}

void Mp2ImaginaryFrequencyGridOptimizer::testGradient(const double stepSize) {
  double E0(getError(grid));
  IntegrationGrid grad(getGradient(grid));
  IntegrationGrid numGrad(grid);
  for (size_t n(0); n < grid.points.size(); ++n) {
    grid.weights[n] += stepSize;
    double E1(getError(grid));
    grid.weights[n] -= stepSize;
    numGrad.weights[n] = (E1-E0) / stepSize;
    LOG(1, "RPA") << "dE/dw[" << n << "] = " <<
      grad.weights[n] << " ~ " << numGrad.weights[n] << std::endl;
  }
  for (size_t n(0); n < grid.points.size(); ++n) {
    grid.points[n] += stepSize;
    double E1(getError(grid));
    grid.points[n] -= stepSize;
    numGrad.points[n] = (E1-E0) / stepSize;
    LOG(1, "RPA") << "dE/dnu[" << n << "] = " <<
      grad.points[n] << " ~ " << numGrad.points[n] << std::endl;
  }
}

double Mp2ImaginaryFrequencyGridOptimizer::getError(
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
        error -= 1 / delta;
      }
    )
  ) (
    (*Dai)["ai"], (*Eai)["ai"]
  );

  // get frobenius norm of errors, which is to be optimized
  double E(Eai->norm2());
  return 0.5*E*E;
}

// expects a call of getError with the same grid first
IntegrationGrid Mp2ImaginaryFrequencyGridOptimizer::getGradient(
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

