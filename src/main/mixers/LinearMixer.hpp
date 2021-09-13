#ifndef LINEAR_MIXER_DEFINED 
#define LINEAR_MIXER_DEFINED

#include <mixers/Mixer.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class LinearMixer: public Mixer<F,TE> {
  public:
    MIXER_REGISTRAR_DECLARATION(LinearMixer)
    LinearMixer(const Ptr<MapNode> &arguments);

    virtual void append(
      const Ptr<TensorUnion<F,TE>> &A, const Ptr<TensorUnion<F,TE>> &R
    );
    virtual Ptr<const TensorUnion<F,TE>> get();
    virtual Real<> getResiduumNorm();

    Ptr<TensorUnion<F,TE>> last;
    Ptr<TensorUnion<F,TE>> lastResiduum;

    F ratio;
    Real<> residuumNorm;
  };
}

#endif

