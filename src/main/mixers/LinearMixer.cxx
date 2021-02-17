#include <mixers/LinearMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <algorithms/Algorithm.hpp>
#include <Data.hpp>

using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(LinearMixer)

template <typename F, typename TE>
LinearMixer<F,TE>::LinearMixer(
  const Ptr<MapNode> &arguments
):
  Mixer<F,TE>(arguments), last(nullptr), lastResiduum(nullptr)
{
  ratio = arguments->getValue<Real<>>("ratio", 1.0);
  LOG() << "ratio=" << real(ratio) << std::endl;
}

template <typename F, typename TE>
void LinearMixer<F,TE>::append(
  const Ptr<FockVector<F,TE>> &next,
  const Ptr<FockVector<F,TE>> &nextResiduum
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
  residuumNorm = sqrt(real(lastResiduum->dot(*lastResiduum)));
}

template <typename F, typename TE>
Ptr<const FockVector<F,TE>> LinearMixer<F,TE>::get() {
  return last;
}

template <typename F, typename TE>
double LinearMixer<F,TE>::getResiduumNorm() {
  return residuumNorm;
}

// instantiate
template class cc4s::LinearMixer<Real<64>, DefaultDryTensorEngine>;
template class cc4s::LinearMixer<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::LinearMixer<Real<64>, DefaultTensorEngine>;
template class cc4s::LinearMixer<Complex<64>, DefaultTensorEngine>;

