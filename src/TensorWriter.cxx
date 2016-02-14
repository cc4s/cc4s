/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <TensorWriter.hpp>
#include <util/Log.hpp>
#include <fstream> 
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorWriter);

TensorWriter::TensorWriter(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorWriter::~TensorWriter() {
}

/**
 * \brief Writes the real tensor data given as Data argument to a file.
 */
void TensorWriter::run() {
  Tensor<> *A(getTensorArgument<>("Data"));
  char defaultIndexOrder[A->order + 1];
  for (int dim(0); dim < A->order; ++dim) defaultIndexOrder[dim] = 'i' + dim;
  defaultIndexOrder[A->order] = 0;
  std::string rowIndexOrder(getTextArgument("rowIndexOrder",defaultIndexOrder));
  std::string columnIndexOrder(getTextArgument("columnIndexOrder", ""));
  Assert(
    rowIndexOrder.length()+columnIndexOrder.length() == unsigned(A->order),
    "Number of indices in rowIndexOrder and columnIndexOrder must match tensor order"
  );
  char indexOrder[A->order+1];
  int lens[A->order];
  int syms[A->order];
  int64_t columnElementsCount(1);
  int columnOrder(columnIndexOrder.length());
  for (int dim(0); dim < columnOrder; ++dim) {
    indexOrder[dim] = columnIndexOrder[dim];
    lens[dim] = A->lens[columnIndexOrder[dim]-'i'];
    syms[dim] = NS;
    columnElementsCount *= lens[dim];
  }
  int64_t rowElementsCount(1);
  int rowOrder(rowIndexOrder.length());
  for (int dim(0); dim < rowOrder; ++dim) {
    indexOrder[columnOrder+dim] = rowIndexOrder[dim];
    lens[columnOrder+dim] = A->lens[rowIndexOrder[dim]-'i'];
    syms[columnOrder+dim] = NS;
    rowElementsCount *= lens[columnOrder+dim];
  }
  indexOrder[A->order] = 0;
  Tensor<> B(A->order, lens, syms, *A->wrld, "DataOrdered");
  // reorder indices for writing:
  B[indexOrder] = (*A)[defaultIndexOrder];
  int64_t valuesCount;
  double *values;
  // and unpack symmetries for writing
  // TODO: implement memory scalable version
  B.read_all(&valuesCount, &values, true);
  Assert(
    rowElementsCount*columnElementsCount == valuesCount,
    "Wrong number of elements read"
  );

  // only the root writes the file
  if (A->wrld->rank == 0) {
    // create file and write header
    std::string dataName(getArgumentData("Data")->getName());
    // by default the file is named after the written data
    std::ofstream file(getTextArgument("file", dataName + ".dat"));
    std::string delimiter(getTextArgument("delimiter", " "));
    file << dataName << delimiter << A->order;
    for (int i(0); i < A->order; ++i) {
      file << delimiter << A->lens[i];
    }
    file << std::endl;
    file << "\"" << rowIndexOrder << "\"" << delimiter
      << "\"" << columnIndexOrder << "\"" << std::endl;

    // write the actual data
    int64_t index(0);
    LOG(4) << "rows=" << rowElementsCount
      << ", columns=" << columnElementsCount << std::endl;
    for (int64_t row(0); row < rowElementsCount; ++row) {
      file << values[index++];
      for (int64_t column(1); column < columnElementsCount; ++column) {
        file << delimiter << values[index++];
      }
      file << std::endl;
    }
  }
  free(values);
}

