#ifndef FOCK_VECTOR_DEFINED
#define FOCK_VECTOR_DEFINED

#include <math/ComplexTensor.hpp>

#include <vector>
#include <string>
#include <initializer_list>

#include <ctf.hpp>

namespace cc4s {
  template <typename F=double>
  class FockVector {
  public:
    std::vector<CTF::Tensor<F>> componentTensors;
    std::vector<std::string> componentIndices;

    FockVector(
      FockVector<F> const &a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices)
    {
    }
  
    FockVector(
      std::initializer_list<CTF::Tensor<F>> componentTensors_
      std::initializer_list<std::string> componentIndices_
    ):
      componentTensors(componentTensors_),
      componentIndices(componentIndices_)
    {
    }

    FockVector<F> &operator += (FockVector<F> &a) {
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        // TODO: check whether given vector a is compatible w.r.t.
        // number and shape components
        char const *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator *= (F const s) {
      // TODO: implement scalar multiplication
    }

    F dot(FockVector<F> &a) {
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        // TODO: check whether given vector a is compatible w.r.t.
        // number and shape components
        char const *indices(componentIndices[i].c_str());
        CTF::Tensor<F> conjThis(componentTensors[i]);
        cc4s::conjugate(conjThis);
        // add to result
        result[""] += conjThis[indices] * a.componentTensors[i][indices];
      }
      return result.get_val();
    }

  };

  template <typename F>
  FockVector<F> inline operator +(FockVector<F> &a, FockVector<F> &b) {
    FockVector<F> result(a);
    result += b;
    return result;
  }
}

#endif

/*
  example usage:

  CTF::Scalar<> nulliesA;
  CTF::Tensor<> singlesA(2, vo, ...);
  CTF::Tensor<> doublesA(4, vvoo, ...);
  
  FockVector<double> fockA({{nulliesA, singlesA, doublesA}}, {{"", "ai", "abij"}});
  FockVector<double> fockB({{nulliesB, singlesB, doublesB}}, {{"", "ai", "abij"}});

  fockA += fockB;
  fockA.dot(fockB);
*/
