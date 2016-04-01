#include <mixers/TrivialMixer.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(TrivialMixer);

template <typename F>
TrivialMixer<F>::TrivialMixer(
  Algorithm *algorithm
):
  Mixer<F>(algorithm), last(nullptr)
{
  ratio = (algorithm->getRealArgument("mixingRatio", 1.0));
  LOG(1,"TrivialMixer") << "ratio=" << ratio << std::endl;
}

template <typename F>
TrivialMixer<F>::~TrivialMixer() {
  if (last) delete last;
}

template <typename F>
void TrivialMixer<F>::append(Tensor<F> &A) {
  if (!last) {
    // create new, copying A
    last = new Tensor<F>(A);
  } else {
    // overwrite last with A
    std::string idx(Mixer<F>::indices(A));
    // (*last)[] = ratio*A[idx.c_str()] + (1-ratio)*(*last);
    last->sum(ratio, A, idx.c_str(), 1-ratio, idx.c_str());
  }
}

template <typename F>
Tensor<F> &TrivialMixer<F>::getNext() {
  return *last;
}

// instantiate
template class TrivialMixer<double>;
template class TrivialMixer<complex>;

