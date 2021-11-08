#ifndef TCC_MAP_DEFINED
#define TCC_MAP_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/MapOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

namespace cc4s {
  template <typename Target, typename Domain, typename TE>
  class Map: public IndexedTensorExpression<Target,TE> {
  public:
    /**
     * \brief Creates a map expression of a unary map f and one tensor
     * expressions source.
     **/
    static Ptr<Map<Target,Domain,TE>> create(
      const std::function<Target(const Domain)> &f,
      const Ptr<IndexedTensorExpression<Domain,TE>> &source
    ) {
      return New<Map<Target,Domain,TE>>(
        f, source,
        typename Expression<TE>::ProtectedToken()
      );
    }

    /**
     * \brief Creates a map expression of a unary map f and one tensor
     * expressions source.
     * Not indended for direct invocation. Use Map::create instead.
     **/
    Map(
      const std::function<Target(const Domain)> &f_,
      const Ptr<IndexedTensorExpression<Domain,TE>> &source_,
      const typename Expression<TE>::ProtectedToken &
    ): f(f_), source(source_) {
    }

    virtual ~Map() {
    }

    Ptr<Operation<TE>> compile(Scope &scope) override {
      auto sourceOperation(
        dynamicPtrCast<IndexedTensorOperation<Domain,TE>>(
          source->compile(scope)
        )
      );
      return MapOperation<Target,Domain,TE>::create(f, sourceOperation, scope);
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    void countIndices(Scope &scope) override {
      source->countIndices(scope);
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Map( " << "f" << ", " << std::string(*source) << " )";
      return stream.str();
    }

  protected:
    std::function<Target(const Domain)> f;
    Ptr<IndexedTensorExpression<Domain,TE>> source;
  };

  /**
   * \brief Creates a map expression of a unary map f and one tensor
   * expressions A.
   **/
  template <
    typename Target, typename RHS
  >
  inline
  Ptr<Map<Target,typename RHS::FieldType,typename RHS::TensorEngine>> map(
    const std::function<Target(typename RHS::FieldType)> &f,
    const Ptr<RHS> &A
  ) {
    return
    Map<Target,typename RHS::FieldType,typename RHS::TensorEngine>::create(
      f, A
    );
  }
}

#endif

