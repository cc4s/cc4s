#include <mixers/LinearMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>

using namespace CTF;
using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(LinearMixer);

template <typename F>
LinearMixer<F>::LinearMixer(
  Algorithm *algorithm
):
  Mixer<F>(algorithm), last(nullptr), lastResiduum(nullptr)
{
  ratio = (algorithm->getRealArgument("mixingRatio", 1.0));
  LOG(1,"LinearMixer") << "ratio=" << ratio << std::endl;
}

template <typename F>
LinearMixer<F>::~LinearMixer() {
}

template <typename F>
void LinearMixer<F>::append(
  const PTR(FockVector<F>) &next, const PTR(FockVector<F>) &nextResiduum
) {
  if (last) {
    // mix accordingly
    *last *= 1-ratio;
    *next *= ratio;
    *next += *last;

    *lastResiduum *= 1-ratio;
    *nextResiduum *= ratio;
    *nextResiduum += *lastResiduum;
  }
  last = next;
  lastResiduum = nextResiduum;
}

template <typename F>
PTR(const FockVector<F>) LinearMixer<F>::get() {
    return last;
}

template <typename F>
PTR(const FockVector<F>) LinearMixer<F>::getResiduum() {
    return lastResiduum;
}

// instantiate
template class LinearMixer<double>;
template class LinearMixer<complex>;

