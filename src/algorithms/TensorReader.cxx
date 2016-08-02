/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorReader.hpp>
#include <util/BinaryTensorFormat.hpp>
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

  // open the file and read header on all ranks
  std::string dataName(getArgumentData("Data")->getName());
  // by default the file is named after the written data
  std::ifstream file(getTextArgument("file", dataName + ".dat").c_str());
  // FIXME: space assumed as delimiter
  char line[1024];
  file.getline(line, sizeof(line), ' ');
  std::string name(line);
  file.getline(line, sizeof(line));
  std::stringstream lineStream(line);
  lineStream >> order;
  lens = new int[order];
  syms = new int[order];
  int64_t indexCount(1);
  for (int dim(0); dim < order; ++dim) {
    lineStream >> lens[dim];
    syms[dim] = NS;
    indexCount *= lens[dim];
  }
  file.getline(line, sizeof(line), ' ');
  std::string rowIndexOrder(line);
  file.getline(line, sizeof(line));
  std::string columnIndexOrder(line);

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

  Tensor<> B(order, storedLens, syms, *Cc4s::world, "B");
  if (B.wrld->rank > 0) indexCount = 0;
  int64_t *indices(new int64_t[indexCount]);
  double *values(new double[indexCount]);
  // read the values only on root
  for (int64_t i(0); i < indexCount; ++i) {
    // FIXME: this algorithm uses a lot of memory: better to read in blocks
    indices[i] = i;
    file >> values[i];
  }
  B.write(indexCount, indices, values);
  LOG(1, "TensorReader") << "Last value=" << values[indexCount-1] << std::endl;
  delete[] indices;
  delete[] values;

  Tensor<> *A(new Tensor<>(order, lens, syms, *B.wrld, name.c_str()));
  char indexOrder[A->order + 1];
  for (int dim(0); dim < A->order; ++dim) {
    indexOrder[dim] = 'i' + dim;
  }
  indexOrder[A->order] = 0;
  (*A)[indexOrder] = B[(columnIndexOrder + rowIndexOrder).c_str()];
  allocatedTensorArgument("Data", A);
}

