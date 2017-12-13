/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MP2_IMAGINARY_FREQUENCY_GRID_DEFINED 
#define MP2_IMAGINARY_FREQUENCY_GRID_DEFINED

#include <algorithms/Algorithm.hpp>

#include <ctf.hpp>

namespace cc4s {
  class Mp2ImaginaryFrequencyGrid: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2ImaginaryFrequencyGrid);
    Mp2ImaginaryFrequencyGrid(
      std::vector<Argument> const &argumentList
    );
    virtual ~Mp2ImaginaryFrequencyGrid();

    virtual void run();
  };

  /**
   * \brief Provides vector space structure for integration grids
   **/
  class IntegrationGrid {
  public:
    IntegrationGrid(const int N): points(N), weights(N) { }
    IntegrationGrid(const IntegrationGrid &g):
      points(g.points), weights(g.weights)
    { }
    IntegrationGrid &operator =(const IntegrationGrid &g) {
      points = g.points; weights = g.weights;
      return *this;
    }
    IntegrationGrid &operator +=(const IntegrationGrid &g) {
      for (size_t n(0); n < points.size(); ++n) {
        points[n] += g.points[n];
        weights[n] += g.weights[n];
      }
      return *this;
    }
    IntegrationGrid &operator -=(const IntegrationGrid &g) {
      for (size_t n(0); n < points.size(); ++n) {
        points[n] -= g.points[n];
        weights[n] -= g.weights[n];
      }
      return *this;
    }
    IntegrationGrid &operator *=(double s) {
      for (size_t n(0); n < points.size(); ++n) {
        points[n] *= s;
        weights[n] *= s;
      }
      return *this;
    }
    IntegrationGrid &operator /=(double s) {
      for (size_t n(0); n < points.size(); ++n) {
        points[n] /= s;
        weights[n] /= s;
      }
      return *this;
    }
    IntegrationGrid &operator -() {
      for (size_t n(0); n < points.size(); ++n) {
        points[n] = -points[n];
        weights[n] = -weights[n];
      }
      return *this;
    }
    double dot(const IntegrationGrid &g) {
      double result(0.0);
      for (size_t n(0); n < points.size(); ++n) {
        result += points[n] * g.points[n];
        result += weights[n] * g.weights[n];
      }
      return result;
    }
    std::vector<double> points, weights;
  };

  inline IntegrationGrid operator +(
    const IntegrationGrid &a, const IntegrationGrid &b
  ) {
    IntegrationGrid result(a);
    return result += b;
  }

  inline IntegrationGrid operator -(
    const IntegrationGrid &a, const IntegrationGrid &b
  ) {
    IntegrationGrid result(a);
    return result -= b;
  }

  inline IntegrationGrid operator *(
    const double s, const IntegrationGrid &g
  ) {
    IntegrationGrid result(g);
    return result *= s;
  }

  inline IntegrationGrid operator /(
    const IntegrationGrid &g, const double s
  ) {
    IntegrationGrid result(g);
    return result /= s;
  }

  class Mp2ImaginaryFrequencyGridOptimizer {
  public:
    Mp2ImaginaryFrequencyGridOptimizer(
      const int N, CTF::Tensor<double> epsi, CTF::Tensor<double> epsa
    );
    void optimize(const int stepCount);

    double lineSearch(IntegrationGrid &grid, const IntegrationGrid &direction);
    void testGradient(const double stepSize);
    double getError(const IntegrationGrid &grid);
    // expects a call of getError with the same grid first
    IntegrationGrid getGradient(const IntegrationGrid &grid);

    inline static double propagator(
      const double delta, const double nu, const int e = 2
    ) {
      return 4*delta*delta / std::pow(delta*delta + nu*nu, e);
    }

    IntegrationGrid grid;

    // Dai = eps_a-eps_i for each ai
    CTF::Tensor<double> *Dai;

    // error of numerical quadrature of current grid wn[n] & nus[n] from
    // analytic value 1/(eps_a-eps_i) for each ai.
    // The current grid is the argument of last getError(grid).
    CTF::Tensor<double> *Eai;
  };
}

#endif

