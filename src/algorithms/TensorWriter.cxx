/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorWriter.hpp>
#include <util/BinaryTensorFormat.hpp>
#include <util/Log.hpp>
#include <fstream> 
#include <iomanip>
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

void TensorWriter::run() {
  std::string mode(getTextArgument("mode", "text"));

  if (mode == "binary") writeBinary();
  else writeText();
}

void TensorWriter::writeBinary() {
  Tensor<> *A(getTensorArgument<>("Data"));
  std::string dataName(getArgumentData("Data")->getName());
  // by default the file is named after the written data
  std::string fileName(getTextArgument("file", dataName + ".bin"));

  MPI_File file;
  MPI_Status status;
  MPI_File_open(
    A->wrld->comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &file
  );
  MPI_File_set_size(file, 0);
  BinaryTensorHeader header(*A);
  MPI_File_write(file, &header, sizeof(header), MPI_BYTE, &status);
  // FIXME: status checking
  A->write_dense_to_file(file, sizeof(header));
  MPI_File_close(&file);
}

void TensorWriter::writeText() {
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
    std::ofstream file(getTextArgument("file", dataName + ".dat").c_str());
    std::string delimiter(getTextArgument("delimiter", " "));
    file << dataName << delimiter << A->order;
    for (int i(0); i < A->order; ++i) {
      file << delimiter << A->lens[i];
    }
    file << std::endl;
    file << rowIndexOrder << delimiter << columnIndexOrder << std::endl;

    // write the actual data
    int64_t index(0);
    LOG(1, "Writer") << "rows=" << rowElementsCount
      << ", columns=" << columnElementsCount << std::endl;
    file << std::setprecision(16);
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

