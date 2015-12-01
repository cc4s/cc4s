#include "Log.hpp"

#include <util/Complex.hpp>

int cc4s::Log::logLevel(0);
std::string cc4s::Log::indent("\t");

std::ostream &cc4s::Log::logStream(int const level) {
  for (int i(0); i < level; ++i) {
    std::cout << indent.c_str();
  }
  return std::cout;
}

template <typename F>
void cc4s::logMatrix(int level, CTF::Matrix<F> &m) {
  F *values(new F[m.lens[0]*m.lens[1]]);
  m.read_all(values);
  for (int i(0); i < m.lens[0]; ++i) {
    for (int j(0); j < m.lens[1]; ++j) {
      LOG(level) << " " << values[i+j*m.lens[0]];
    }
    LOG(level) << std::endl;
  }
}

// instantiate
template
void cc4s::logMatrix(int level, CTF::Matrix<double> &m);
template
void cc4s::logMatrix(int level, CTF::Matrix<complex> &m);

