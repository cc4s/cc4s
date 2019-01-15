/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MAP_DEFINED
#define TCC_MAP_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/MapOperation.hpp>

#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

namespace tcc {
  template <typename F>
  class Map: public Expression<F> {
  public:
    /**
     * \brief Creates a map expression of a unary map f and one tensor
     * expressions source.
     **/
    static PTR(Map<F>) create(
      const std::function<F(const F)> &f,
      const PTR(Expression<F>) &source
    ) {
      auto map(
        NEW(Map<F>,
          f, source,
          typename Expression<F>::ProtectedToken()
        )
      );
      return map;
    }

    /**
     * \brief Creates a map expression of a unary map f and one tensor
     * expressions source.
     * Not indended for direct invocation. Use Map::create instead.
     **/
    Map(
      const std::function<F(const F)> &f_,
      const PTR(Expression<F>) &source_,
      const typename Expression<F>::ProtectedToken &
    ): f(f_), source(source_) {
    }

    virtual ~Map() {
    }

    virtual PTR(Operation<F>) compile(IndexCounts &indexCounts) {
      auto sourceOperation(source->compile(indexCounts));
      return MapOperation<F>::create(f, sourceOperation);
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      source->countIndices(indexCounts);
    }

  protected:
    std::function<F(const F)> f;
    PTR(Expression<F>) source;
  };

  /**
   * \brief Creates a map expression of a unary map f and one tensor
   * expressions A.
   **/
  template <typename F>
  inline PTR(Map<F>) map(
    const std::function<F(const F)> &f,
    const PTR(Expression<F>) &A
  ) {
    return Map<F>::create(f, A);
  }
}

#endif

