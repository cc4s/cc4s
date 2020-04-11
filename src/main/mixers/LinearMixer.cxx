#include <mixers/LinearMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>
#include <Data.hpp>

using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(LinearMixer);

template <typename F, typename TE>
LinearMixer<F,TE>::LinearMixer(
  Algorithm *algorithm
):
  Mixer<F,TE>(algorithm), last(nullptr), lastResiduum(nullptr)
{
  ratio = (algorithm->getRealArgument("mixingRatio", 1.0));
  LOG(1,"LinearMixer") << "ratio=" << std::real(ratio) << std::endl;
  EMIT() << YAML::Key << "mixer" << YAML::Value;
  EMIT() << YAML::BeginMap;
  EMIT() << YAML::Key << "type" << YAML::Value << "linear";
  EMIT() << YAML::Key << "ratio" << YAML::Value << std::real(ratio);
  EMIT() << YAML::EndMap;
}

template <typename F, typename TE>
LinearMixer<F,TE>::~LinearMixer() {
}

template <typename F, typename TE>
void LinearMixer<F,TE>::append(
  const PTR(ESC(FockVector<F,TE>)) &next,
  const PTR(ESC(FockVector<F,TE>)) &nextResiduum
) {
  if (last) {
    // mix accordingly
    *last *= F(1)-ratio;
    *next *= ratio;
    *next += *last;

    *lastResiduum *= F(1)-ratio;
    *nextResiduum *= ratio;
    *nextResiduum += *lastResiduum;
  }
  last = next;
  lastResiduum = nextResiduum;
}

template <typename F, typename TE>
PTR(ESC(const FockVector<F,TE>)) LinearMixer<F,TE>::get() {
    return last;
}

template <typename F, typename TE>
PTR(ESC(const FockVector<F,TE>)) LinearMixer<F,TE>::getResiduum() {
    return lastResiduum;
}

// instantiate
template class LinearMixer<Real<64>, DryTensorEngine>;
template class LinearMixer<Complex<64>, DryTensorEngine>;
template class LinearMixer<Real<64>, DefaultTensorEngine>;
template class LinearMixer<Complex<64>, DefaultTensorEngine>;

