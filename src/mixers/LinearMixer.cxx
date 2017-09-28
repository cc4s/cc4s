#include <mixers/LinearMixer.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>

#include <memory>

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
}

template <typename F>
void LinearMixer<F>::append(FockVector<F> &A) {
  if (!last) {
    // create new, copying A
    last = std::make_shared<FockVector<F>>(A);
  } else {
    // overwrite last with A
    *last *= 1-ratio;
    *last += A;
  }
}

template <typename F>
FockVector<F> &LinearMixer<F>::getNext() {
    return *last;
}

// instantiate
template class LinearMixer<double>;
template class LinearMixer<complex>;

