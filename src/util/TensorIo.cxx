/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <util/TensorIo.hpp>
#include <util/BinaryTensorFormat.hpp>
#include <util/Scanner.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <fstream>
#include <iomanip>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

template <typename F>
Tensor<F> *TensorIo::readBinary(std::string const &fileName) {
  // open the file
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Cc4s::world->comm, fileName.c_str(), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &file
    )
  );
  if (mpiError) throw new EXCEPTION("Failed to open file");

  int64_t offset(0);
  MPI_Status status;
  // reade header
  BinaryTensorHeader header;
  MPI_File_read_at(file, offset, &header, sizeof(header), MPI_BYTE, &status);
  offset += sizeof(header);
  if (strncmp(header.magic, header.MAGIC, sizeof(header.magic)) != 0)
    throw new EXCEPTION("Invalid file format");
  if (header.version > header.VERSION)
    throw new EXCEPTION("Incompatible file format version");

  // read dimension headers
  int lens[header.order];
  int syms[header.order];
  for (int dim(0); dim < header.order; ++dim) {
    BinaryTensorDimensionHeader dimensionHeader;
    MPI_File_read_at(
      file, offset, &dimensionHeader, sizeof(dimensionHeader), MPI_BYTE, &status
    );
    offset += sizeof(dimensionHeader);
    lens[dim] = dimensionHeader.length;
    syms[dim] = NS;
  }

  // allocate tensor
  Tensor<F> *A(
    new Tensor<F>(header.order, lens, syms, *Cc4s::world)
  );

  // read dense data
  A->read_dense_from_file(file, offset);

  // done
  MPI_File_close(&file);
  return A;
}

template <typename F>
Tensor<F> *TensorIo::readText(
  std::string const &fileName,
  std::string const &delimiter,
  int64_t const bufferSize
) {
  int order;
  int *lens;
  int *syms;

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

  Tensor<F> *B(new Tensor<F> (order, storedLens, syms, *Cc4s::world, "B"));
  int64_t localBufferSize(B->wrld->rank == 0 ? bufferSize : 0);
  int64_t *indices(new int64_t[localBufferSize]);
  F *values(new F[localBufferSize]);
  // read the values only on root
  int64_t index(0);
  LOG(1, "TensorReader") << "indexCount=" << indexCount << std::endl;
  NumberScanner<F> numberScanner(&scanner);
  while (index < indexCount) {
    int64_t elementsCount(std::min(bufferSize, indexCount-index));
    int64_t localElementsCount(B->wrld->rank == 0 ? elementsCount : 0);
    for (int64_t i(0); i < localElementsCount; ++i) {
      indices[i] = index+i;
      values[i] = numberScanner.nextNumber();
    }
    // wait until all processes finished reading this buffer
    MPI_Barrier(Cc4s::world->comm);
    LOG(1, "TensorReader") << "writing " << elementsCount << " values to tensor..." << std::endl;
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
    Tensor<F> *A(new Tensor<F>(order, lens, syms, *B->wrld, name.c_str()));
    (*A)[indexOrder] = (*B)[storedIndexOrder.c_str()];
    delete B;
    return A;
  } else {
    B->set_name(name.c_str());
    return B;
  }
}

template <typename F>
void TensorIo::writeBinary(std::string const &fileName, Tensor<F> &A) {
  MPI_File file;
  MPI_Status status;
  MPI_File_open(
    A.wrld->comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &file
  );
  int64_t offset(0);

  // truncate possibly existing file
  MPI_File_set_size(file, offset);
  // write header
  BinaryTensorHeader header(A);
  MPI_File_write_at(file, offset, &header, sizeof(header), MPI_BYTE, &status);
  offset += sizeof(header);
  // FIXME: status checking

  // write dimension header for each dimension
  for (int dim(0); dim < A.order; ++dim) {
    BinaryTensorDimensionHeader dimensionHeader(A.lens[dim], 'a'+dim);
    MPI_File_write_at(
      file, offset, &dimensionHeader, sizeof(dimensionHeader), MPI_BYTE, &status
    );
    offset += sizeof(dimensionHeader);
  }

  // write dense data
  A.write_dense_to_file(file, offset);

  // done
  MPI_File_close(&file);
}

template <typename F>
void TensorIo::writeText(
  std::string const &fileName, Tensor<F> &A,
  std::string const &givenRowIndexOrder, std::string const &columnIndexOrder,
  std::string const &delimiter
) {
  char defaultIndexOrder[A.order + 1];
  for (int dim(0); dim < A.order; ++dim) defaultIndexOrder[dim] = 'i' + dim;
  defaultIndexOrder[A.order] = 0;
  std::string rowIndexOrder(
    givenRowIndexOrder == "" && columnIndexOrder == "" ?
      defaultIndexOrder : givenRowIndexOrder
  );
  Assert(
    rowIndexOrder.length()+columnIndexOrder.length() == unsigned(A.order),
    "Number of indices in rowIndexOrder and columnIndexOrder must match tensor order"
  );
  char indexOrder[A.order+1];
  int lens[A.order];
  int syms[A.order];
  int64_t columnElementsCount(1);
  int columnOrder(columnIndexOrder.length());
  for (int dim(0); dim < columnOrder; ++dim) {
    indexOrder[dim] = columnIndexOrder[dim];
    lens[dim] = A.lens[columnIndexOrder[dim]-'i'];
    syms[dim] = NS;
    columnElementsCount *= lens[dim];
  }
  int64_t rowElementsCount(1);
  int rowOrder(rowIndexOrder.length());
  for (int dim(0); dim < rowOrder; ++dim) {
    indexOrder[columnOrder+dim] = rowIndexOrder[dim];
    lens[columnOrder+dim] = A.lens[rowIndexOrder[dim]-'i'];
    syms[columnOrder+dim] = NS;
    rowElementsCount *= lens[columnOrder+dim];
  }
  indexOrder[A.order] = 0;
  Tensor<F> B(A.order, lens, syms, *A.wrld, "DataOrdered");
  // reorder indices for writing:
  B[indexOrder] = A[defaultIndexOrder];
  int64_t valuesCount;
  F *values;
  // and unpack symmetries for writing
  // TODO: implement memory scalable version
  B.read_all(&valuesCount, &values, true);
  Assert(
    rowElementsCount*columnElementsCount == valuesCount,
    "Wrong number of elements read"
  );

  // only the root writes the file
  if (A.wrld->rank == 0) {
    std::ofstream file(fileName.c_str());
    file << A.get_name() << delimiter << A.order;
    for (int i(0); i < A.order; ++i) {
      file << delimiter << A.lens[i];
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

// instantiate
template
Tensor<double> *TensorIo::readBinary(std::string const &fileName);
template
Tensor<complex> *TensorIo::readBinary(std::string const &fileName);

template
Tensor<double> *TensorIo::readText(
  std::string const &fileName,
  std::string const &delimiter,
  int64_t const bufferSize
);
template
Tensor<complex> *TensorIo::readText(
  std::string const &fileName,
  std::string const &delimiter,
  int64_t const bufferSize
);

template void TensorIo::writeBinary<double>(
  std::string const &fileName, Tensor<double> &A
);
template void TensorIo::writeBinary<complex>(
  std::string const &fileName, Tensor<complex> &A
);

template void TensorIo::writeText<double>(
  std::string const &fileName, Tensor<double> &A,
  std::string const &rowIndexOrder, std::string const &columnIndexOrder,
  std::string const &delimiter
);
template void TensorIo::writeText<complex>(
  std::string const &fileName, Tensor<complex> &A,
  std::string const &rowIndexOrder, std::string const &columnIndexOrder,
  std::string const &delimiter
);

