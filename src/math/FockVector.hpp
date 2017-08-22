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
    typedef F FieldType;

    std::vector<CTF::Tensor<F>> componentTensors;
    std::vector<std::string> componentIndices;

    FockVector(
      const FockVector<F> &a
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

    FockVector<F> &operator += (const FockVector<F> &a) {
      checkCompatabilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator -= (const FockVector<F> &a) {
      checkCompatabilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += F(-1) * a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator *= (coonst F s) {
      checkCompatabilityTo(a);
      CTF::Scalar<F> scalar(s);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] *= scalar[""];
      }
      return *this;
    }

    F dot(const FockVector<F> &a) const {
      checkCompatabilityTo(a);
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        CTF::Bivar_Function<F> fDot(&cc4s::dot<F>);
        // add to result
        result.contract(
          1.0, componentTensors[i], indices, a.componentTensors[i], indices,
          1.0, "", fDot
        );
      }
      return result.get_val();
    }
  protected:
    void checkCompatabilityTo(const FockVector<F> &a) const {
      if (
        componentTensors.size() != a.componentTensors.size() ||
        componentIndices.size() != a.componentIndices.size() 
      ) {
        throw EXCEPTION("Number of component tensors does no match");
      }
      // TODO: check shapes.
    }
  };

  template <typename F>
  FockVector<F> inline operator +(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result += b;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator -(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result -= b;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator *(const FockVector<F> &a, const F &s) {
    FockVector<F> result(a);
    result *= s;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator *(const F &s, const FockVector<F> &a) {
    FockVector<F> result(a);
    result *= s;
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
