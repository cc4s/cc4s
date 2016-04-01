#include <mixers/TrivialMixer.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(TrivialMixer);

template <typename F>
TrivialMixer<F>::TrivialMixer(): last(nullptr) {
}

template <typename F>
TrivialMixer<F>::~TrivialMixer() {
}

template <typename F>
void TrivialMixer<F>::append(Tensor<F> &A) {
  if (!last) {
    // create new, copying A
    last = new Tensor<F>(A);
  } else {
    // overwrite last with A
    std::string idx(Mixer<F>::indices(A));
    (*last)[idx.c_str()] = A[idx.c_str()];
  }
}

template <typename F>
Tensor<F> &TrivialMixer<F>::getNext() {
  return *last;
}

// instantiate
template class TrivialMixer<double>;
template class TrivialMixer<complex>;

