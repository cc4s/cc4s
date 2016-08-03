/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorReader.hpp>
#include <util/BinaryTensorFormat.hpp>
#include <util/LineNumberStream.hpp>
#include <util/Scanner.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <fstream>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorReader);

TensorReader::TensorReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorReader::~TensorReader() {
}

void TensorReader::run() {
  std::string mode(getTextArgument("mode", "text"));

  // make sure all processes start reading the file at the same time in case
  // it has been modified before
  MPI_Barrier(Cc4s::world->comm);

  if (mode == "binary") readBinary();
  else readText();
}

void TensorReader::readBinary() {
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

void TensorReader::readText() {
  int order;
  int *lens;
  int *syms;

  // open the file and read header on all processes
  std::string dataName(getArgumentData("Data")->getName());
  // by default the file is named after the written data
  std::string fileName(getTextArgument("file", dataName + ".dat").c_str());
  std::ifstream stream(fileName.c_str());
  Scanner scanner(&stream);
  std::string name(scanner.nextLine(' '));
  std::stringstream lineStream(scanner.nextLine());
  lineStream >> order;
  lens = new int[order];
  syms = new int[order];
  int64_t indexCount(1);
  for (int dim(0); dim < order; ++dim) {
    lineStream >> lens[dim];
    syms[dim] = NS;
    indexCount *= lens[dim];
  }
  std::string rowIndexOrder(scanner.nextLine(' '));
  std::string columnIndexOrder(scanner.nextLine());

  int *storedLens(new int[order]);
  int storedIndex(0);
  for (unsigned int dim(0); dim < columnIndexOrder.length(); ++dim) {
    // FIXME: i,j,k ... assumed in indexOrder strings
    storedLens[storedIndex] = lens[columnIndexOrder[dim] - 'i'];
    ++storedIndex;
  }
  for (unsigned int dim(0); dim < rowIndexOrder.length(); ++dim) {
    storedLens[storedIndex] = lens[rowIndexOrder[dim] - 'i'];
    ++storedIndex;
  }

  Tensor<> *B(new Tensor<> (order, storedLens, syms, *Cc4s::world, "B"));
  int64_t bufferSize(getIntegerArgument("bufferSize", 128l*1024*1024));
  int64_t localBufferSize(B->wrld->rank == 0 ? bufferSize : 0);
  int64_t *indices(new int64_t[localBufferSize]);
  double *values(new double[localBufferSize]);
  // read the values only on root
  int64_t index(0);
  LOG(1, "TensorReader") << "indexCount=" << indexCount << std::endl;
  while (index < indexCount) {
    int64_t elementsCount(std::min(bufferSize, indexCount-index));
    int64_t localElementsCount(B->wrld->rank == 0 ? elementsCount : 0);
    for (int64_t i(0); i < localElementsCount; ++i) {
      indices[i] = index+i;
      values[i] = scanner.nextReal();
    }
    // wait until all processes finished reading this buffer
    MPI_Barrier(Cc4s::world->comm);
    B->write(localElementsCount, indices, values);
    index += elementsCount;
  }
  delete[] indices;
  delete[] values;

  char indexOrder[order + 1];
  for (int dim(0); dim < order; ++dim) {
    indexOrder[dim] = 'i' + dim;
  }
  indexOrder[order] = 0;
  std::string storedIndexOrder(columnIndexOrder + rowIndexOrder);
  if (std::string(indexOrder) != storedIndexOrder) {
    Tensor<> *A(new Tensor<>(order, lens, syms, *B->wrld, name.c_str()));
    (*A)[indexOrder] = (*B)[storedIndexOrder.c_str()];
    allocatedTensorArgument("Data", A);
    delete B;
  } else {
    B->set_name(name.c_str());
    allocatedTensorArgument("Data", B);
  }
}

