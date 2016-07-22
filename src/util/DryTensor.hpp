/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_DEFINED
#define DRY_TENSOR_DEFINED

#include <util/Log.hpp>
#include <cstdint>
#include <vector>

namespace cc4s {
  class DryMemory {
  public:
    static void allocate(int64_t size) {
      currentTotalSize += size;
      if (currentTotalSize > maxTotalSize) maxTotalSize = currentTotalSize;
    }
    static void free(int64_t size) {
      currentTotalSize -= size;
    }
    static int64_t currentTotalSize, maxTotalSize;
  };

  template <typename F=double>
  class DryTensor {
  public:
    /**
     * \brief Creates a dry tensor for resource consumption estimation
     * without actually allocating its data.
     * Note that currently only the memory requirement of a DryTensor
     * is estimated. The number of FLOPS should also be estimated in the
     * future.
     * Symmetryies are also ignored at the moment.
     */
    DryTensor(
      int order_, int const *lens_, int const *syms_
    ): order(order_) {
      for (int i(0); i < order_; ++i) {
        lens.push_back(lens_[i]);
        syms.push_back(syms_[i]);
      }
      allocate();
    }
    DryTensor(DryTensor const &A): order(A.order), lens(A.lens), syms(A.syms) {
      allocate();
    }
    virtual ~DryTensor() {
      free();
    }
    virtual void use() {}

    int order;
    std::vector<int> lens, syms;

  protected:
    void allocate() {
      size = sizeof(F);
      for (int i(0); i < order; ++i) {
        size *= lens[i];
      }
      DryMemory::allocate(size);
    }
    void free() {
      DryMemory::free(size);
    }
    int64_t size;
  };

  template <typename F=double>
  class DryMatrix: public DryTensor<F> {
  public:
    DryMatrix(int rowCount, int columnCount, int sym);
  };

  template <typename F=double>
  class DryVector: public DryTensor<F> {
  public:
    DryVector(int elementsCount);
  };

  template <typename F=double>
  class DryScalar: public DryTensor<F> {
  public:
    DryScalar();
    DryScalar(F const value);
  };
}

#endif

