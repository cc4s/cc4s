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
  CTF::Tensor<complex> &c,
  CTF::Tensor<double> &r, CTF::Tensor<double> &i
) {
  AssertCompatibleTensorShapes(c,r,i);
  char inds[c.order];
  for (int i=0; i<c.order; i++){ inds[i] = 'a'+i; }
  r[inds] = CTF::Function<complex, double>([](complex c){ return c.real(); }) (
    c[inds]
  );
  i[inds] = CTF::Function<complex, double>([](complex c){ return c.imag(); }) (
    c[inds]
  );
}

void cc4s::fromComplexTensor(
  CTF::Tensor<complex> &c,
  CTF::Tensor<double> &r
) {
  AssertCompatibleTensorShapes(c,r,r);
  char inds[c.order];
  for (int i=0; i<c.order; i++){ inds[i] = 'a'+i; }
  r[inds] = CTF::Function<complex, double>([](complex c){ return c.real(); }) (
    c[inds]
  );
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> &r, CTF::Tensor<double> &i,
  CTF::Tensor<complex> &c
) {
  AssertCompatibleTensorShapes(c,r,i);
  char inds[c.order];
  for (int i=0; i<c.order; i++){ inds[i] = 'a'+i; }
  toComplexTensor(r, c);
  CTF::Transform<double, complex>([](double d, complex & c){ 
    c.imag(d); 
  })(i[inds], c[inds]);
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> &r,
  CTF::Tensor<complex> &c
) {
  AssertCompatibleTensorShapes(c,r,r);
    char inds[c.order];
  for (int i=0; i<c.order; i++){ inds[i] = 'a'+i; }
  CTF::Transform<double, complex>([](double d, complex & c){ 
    c.real(d); 
  })(r[inds], c[inds]);
}

