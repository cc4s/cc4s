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
  CTF::Tensor<complex> const &c,
  CTF::Tensor<double> &r, CTF::Tensor<double> &i
) {
  AssertCompatibleTensorShapes(c,r,i);
  int64_t indicesCount, *indices;
  complex *values;
  c.read_local(&indicesCount, &indices, &values);
  double *reals(new double[indicesCount]), *imags(new double[indicesCount]);
  for (int64_t i(0); i < indicesCount; ++i) {
    reals[i] = std::real(values[i]);
    imags[i] = std::imag(values[i]);
  }
  free(values);
  r.write(indicesCount, indices, reals);
  i.write(indicesCount, indices, imags);
  free(indices);
  delete[] reals; delete[] imags;
}

void cc4s::fromComplexTensor(
  CTF::Tensor<complex> const &c,
  CTF::Tensor<double> &r
) {
  AssertCompatibleTensorShape(c,r);
  int64_t indicesCount, *indices;
  complex *values;
  c.read_local(&indicesCount, &indices, &values);
  double *reals(new double[indicesCount]);
  for (int64_t i(0); i < indicesCount; ++i) {
    reals[i] = std::real(values[i]);
  }
  free(values);
  r.write(indicesCount, indices, reals);
  free(indices);
  delete[] reals;
}


void cc4s::toComplexTensor(
  CTF::Tensor<double> const &r, CTF::Tensor<double> &i,
  CTF::Tensor<complex> &c
) {
  AssertCompatibleTensorShapes(c,r,i);
  int64_t indicesCount, *indices;
  double *reals;
  r.read_local(&indicesCount, &indices, &reals);
  double *imags(new double[indicesCount]);
  i.read(indicesCount, indices, imags);
  complex *values(new complex[indicesCount]);
  for (int64_t i(0); i < indicesCount; ++i) {
#ifdef INTEL_COMPILER
    values[i].real() = reals[i];
    values[i].imag() = imags[i];
#else
    values[i].real(reals[i]);
    values[i].imag(imags[i]);
#endif
  }
  free(reals); delete[] imags;
  c.write(indicesCount, indices, values);
  free(indices);
  delete[] values;
}

void cc4s::toComplexTensor(
  CTF::Tensor<double> const &r, CTF::Tensor<double> const &i,
  CTF::Tensor<complex> &c
) {
  AssertCompatibleTensorShapes(c,r,i);
  int64_t indicesCount, *indices;
  double *components;
  r.read_local(&indicesCount, &indices, &components);
  complex *values(new complex[indicesCount]);
  for (int64_t i(0); i < indicesCount; ++i) {
#ifdef INTEL_COMPILER
    values[i].real() = components[i];
#else
    values[i].real(components[i]);
#endif
  }
  free(components);
  c.write(indicesCount, indices, values);
  free(indices);
  delete[] values;

  i.read_local(&indicesCount, &indices, &components);
  values = new complex[indicesCount];
  // read complex values that already contain the real part in the
  // index order of the imaginary part. This may redistribute the complex
  // tensor.
  c.read(indicesCount, indices, values);
  for (int64_t i(0); i < indicesCount; ++i) {
#ifdef INTEL_COMPILER
    values[i].imag() = components[i];
#else
    values[i].imag(components[i]);
#endif
  }
  free(components);
  c.write(indicesCount, indices, values);
  free(indices);
  delete[] values;
}

