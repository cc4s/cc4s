/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/IndexedTensor.hpp>
#include <util/Log.hpp>
#include <util/StaticAssert.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace tcc {
  class Tcc;

  template <typename F>
  class Tensor: public std::enable_shared_from_this<Tensor<F>> {
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };

  public:
    Tensor(
      const std::vector<int> &lens_,
      const std::string &name_,
      Tcc *tcc_,
      const ProtectedToken &
    ): lens(lens_), name(name_), tcc(tcc_) {
    }

    void setName(const std::string &name_) {
      name = name_;
    }
    const std::string &getName() const {
      return name;
    }

    Tcc *getTcc() {
      return tcc;
    }

    int64_t getElementsCount() const {
      int64_t elementsCount(1);
      for (unsigned int i(0); i < lens.size(); ++i) {
        elementsCount *= lens[i];
      }
      return elementsCount;
    }

    /**
     * \brief Specify named indices of this tensor to be used in a
     * tensor expression. Indexed tensors are atomic types of tensor
     * expressions.
     **/
    std::shared_ptr<IndexedTensor<F>> operator[](std::string const &indices) {
      return IndexedTensor<F>::create(
        this->shared_from_this(), indices
      );
    }

    std::vector<int> lens;

    friend class Tcc;

  protected:
    std::string name;
    Tcc *tcc;
  };

  class Tcc {
  public:
    /**
     * \brief Creates an abstract tensor object which is elementary to
     * all tensor expression compilation and execution functionality.
     * Symmetryies are not supported at the moment.
     * Note that tensor objects should only be created by the Tcc object
     * which specifies the environment the tensor lives in.
     */
    template <typename F=double>
    std::shared_ptr<Tensor<F>> createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::make_shared<Tensor<F>>(
        lens, name, this, typename Tensor<F>::ProtectedToken()
      );
    }
/*
    void *operator new(const size_t size) {
      static_assert(
        cc4s::StaticAssert<Tcc>::False,
        "Tcc objects cannot be allocated dynamically"
      );
      return nullptr;
    }
*/
  };
}

#endif

