#include <math/ComplexTensor.hpp>

#include <util/Exception.hpp>
#include <complex>

#define AssertCompatibleTensorShapes(c,r,i) \
{ \
  Assert( \
    (c).order==(r).order && (c).order==(i).order, \
    "Incompatible tensor orders" \
  ); \
  for (int k(0); k < c.order; ++k) \
    Assert( \
      (c).lens[k]==(r).lens[k] && (c).lens[k]==(i).lens[k], \
      "Incompatible tensor shapes" \
    ); \
}

#define AssertCompatibleTensorShape(c,r) \
{ \
  Assert( \
    (c).order==(r).order, \
    "Incompatible tensor orders" \
  ); \
  for (int k(0); k < c.order; ++k) \
    Assert( \
      (c).lens[k]==(r).lens[k], \
      "Incompatible tensor shapes" \
    ); \
}

void cc4s::fromComplexTensor(
  CTF::Tensor<complex> &C,
  CTF::Tensor<double> &R, CTF::Tensor<double> &I
) {
  AssertCompatibleTensorShapes(C,R,I);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  fromComplexTensor(C, R);
  I[inds] = CTF::Function<complex, double>(
    std::function<double(complex)>([](complex c){ return c.imag(); })
  ) (
    C[inds]
  );
}

void cc4s::fromComplexTensor(
  CTF::Tensor<double> &C,
  CTF::Tensor<double> &R, CTF::Tensor<double> &I
) {
  AssertCompatibleTensorShapes(C,R,I);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  R[inds] = C[inds];
  I[inds] = 0.0;
}

void cc4s::fromComplexTensor(
  CTF::Tensor<complex> &C,
  CTF::Tensor<double> &R
) {
  AssertCompatibleTensorShape(C,R);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  R[inds] = CTF::Function<complex, double>(
    std::function<double(complex)>([](complex c){ return c.real(); })
  ) (
    C[inds]
  );
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> &R, CTF::Tensor<double> &I,
  CTF::Tensor<complex> &C
) {
  AssertCompatibleTensorShapes(C,R,I);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  toComplexTensor(R, C);
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double i, complex & c) {
//#ifdef INTEL_COMPILER
//        c.imag() = i;
//#else
        c.imag(i);
//#endif
      }
    )
  ) (
    I[inds], C[inds]
  );
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> &R,
  CTF::Tensor<complex> &C
) {
  AssertCompatibleTensorShape(C,R);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double r, complex & c) {
//#ifdef INTEL_COMPILER
//        c.real() = r;
//#else
        c.real(r);
        c.imag(0);
//#endif
      }
    )
  ) (
    R[inds], C[inds]
  );
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> &R,
  CTF::Tensor<double> &C
) {
  AssertCompatibleTensorShape(C,R);
  char inds[C.order];
  for (int i(0); i < C.order; ++i) { inds[i] = 'a'+i; }
  C[inds] = R[inds];
}

void cc4s::conjugate(
  CTF::Tensor<complex> &C
) {
  char inds[C.order];
  for (int i=0; i<C.order; i++){ inds[i] = 'a'+i; }
  CTF::Transform<complex>(
    std::function<void(complex &)> (
      [](complex &c){ c = std::conj(c); }
    )
  ) (
    C[inds]
  );
}

void cc4s::conjugate(
  CTF::Tensor<double> &C
) {
  // ;-)
}

