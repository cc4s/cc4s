#include <mixers/LinearMixer.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(LinearMixer);

template <typename F>
LinearMixer<F>::LinearMixer(
  Algorithm *algorithm
):
  Mixer<F>(algorithm), last(nullptr)
{
  ratio = (algorithm->getRealArgument("mixingRatio", 1.0));
  LOG(1,"LinearMixer") << "ratio=" << ratio << std::endl;
}

template <typename F>
LinearMixer<F>::~LinearMixer() {
  if (last) delete last;
}

template <typename F>
void LinearMixer<F>::append(Tensor<F> &A) {
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
Tensor<F> &LinearMixer<F>::getNext() {
  return *last;
}

// instantiate
template class LinearMixer<double>;
template class LinearMixer<complex>;

