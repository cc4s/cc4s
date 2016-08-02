#include <algorithms/PartitionTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PartitionTensor);

PartitionTensor::PartitionTensor(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PartitionTensor::~PartitionTensor() {
}

/**
 * \brief Partition the given tensor along the specified coordinate.
 */
void PartitionTensor::run() {
  Tensor<> *A(getTensorArgument<>("A"));
  int dimension(getIntegerArgument("dimension", 0));
  if (dimension == A->order-1) {
    partitionLastDimension();
  } else {
    int *start(new int[A->order]);
    int *end(new int[A->order]);
    int *sliceStart(new int[A->order]);
    int *sliceLens(new int[A->order]);
    for (int dim(0); dim < A->order; ++dim) {
      start[dim] = 0;
      end[dim] = A->lens[dim];
      sliceStart[dim] = 0;
      sliceLens[dim] = A->lens[dim];
    }
    sliceLens[dimension] = 1;
    std::string prefix(getTextArgument("prefix", A->get_name()));
    for (int i(0); i < A->lens[dimension]; ++i) {
      LOG(1, "PartitionTensor") << "Slicing part " << i << std::endl;
      start[dimension] = i;
      end[dimension] = i+1;
      std::stringstream sliceName;
      sliceName << prefix << i;
      Tensor<> *ASlice(
        new Tensor<>(A->order, sliceLens, A->sym, *A->wrld, sliceName.str().c_str())
      );
      ASlice->slice(sliceStart, sliceLens, 0.0, *A, start, end, 1.0);
      allocatedTensorArgument(sliceName.str(), ASlice);
    }
  }
  if (getIntegerArgument("discardSource", 0) == 1) {
    LOG(0, "PartitionTensor") << "Discarding source tensor" << std::endl;
    delete A; // this also removes A from the data list
  }
}

void PartitionTensor::partitionLastDimension() {
  Tensor<> *A(getTensorArgument<>("A"));
  int *sliceLens(new int[A->order]);
  int64_t stride(1);
  for (int dim(0); dim < A->order-1; ++dim) {
    stride *= sliceLens[dim] = A->lens[dim];
  }
  sliceLens[A->order-1] = 1;
  int64_t offset(0);
  int64_t *indices(new int64_t[stride]);
  double *values(new double[stride]);
  std::string prefix(getTextArgument("prefix", A->get_name()));
  for (int i(0); i < A->lens[A->order-1]; ++i) {
    LOG(1, "PartitionTensor") << "Directrly slicing part " << i << std::endl;
    std::stringstream sliceName;
    sliceName << prefix << i;
    for (int64_t j(0); j < stride; ++j) {
      indices[j] = offset + j;
    }
    A->read(stride, indices, values);
    LOG(1, "PartitionTensor") << "Last value=" << values[stride-1] << std::endl;
    for (int64_t j(0); j < stride; ++j) {
      indices[j] = j;
    }
    Tensor<> *ASlice(
      new Tensor<>(A->order, sliceLens, A->sym, *A->wrld, sliceName.str().c_str())
    );
    ASlice->write(stride, indices, values);
    allocatedTensorArgument(sliceName.str(), ASlice);
    offset += stride;
  }
}

