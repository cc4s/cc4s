/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/IndexedTensor.hpp>
#include <util/Log.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace tcc {
  template <typename F>
  class Tensor: public std::enable_shared_from_this<Tensor<F>> {
  public:
    /**
     * \brief Creates an abstract tensor object which is elementary to
     * all tensor expression compilation and execution functionality.
     * Symmetryies are not supported at the moment.
     */
    Tensor(
      const std::vector<int> lens_,
      const std::string &name_
    ): lens(lens_), name(name_) {
    }

    void setName(const std::string &name_) {
      name = name_;
    }
    const std::string &getName() const {
      return name;
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
      return std::make_shared<IndexedTensor<F>>(
        this->shared_from_this(), indices
      );
    }

    std::vector<int> lens;

  protected:
    std::string name;
  };
}

#endif

