/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <util/TensorIo.hpp>
#include <util/BinaryTensorFormat.hpp>
#include <util/Scanner.hpp>
#include <util/Log.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <Data.hpp>
#include <fstream>
#include <iomanip>

using namespace cc4s;

template <typename F, typename TE>
Ptr<Tensor<F,TE>> TensorIo::readBinary(const std::string &fileName) {
  // open the file
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Cc4s::world->getComm(), fileName.c_str(), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &file
    )
  );
  if (mpiError) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw New<Exception>(explanation.str(), SOURCE_LOCATION);
  }

  size_t offset(0);
  auto A( readBinaryHeader<F,TE>(file, offset) );

  // read dense data
  // FIXME: invoke tcc read from file calls
//  A->read_dense_from_file(file, offset);

  // done
  MPI_File_close(&file);
  return A;
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> TensorIo::readText(
  const std::string &fileName,
  const std::string &delimiter,
  const size_t bufferSize
) {
  std::ifstream stream(fileName.c_str());
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw New<Exception>(explanation.str(), SOURCE_LOCATION);
  }
  Scanner scanner(&stream);
  std::string name(scanner.nextLine(' '));
  std::stringstream lineStream(scanner.nextLine());
  unsigned int order;
  lineStream >> order;
  std::vector<size_t> lens(order);
  for (unsigned int dim(0); dim < order; ++dim) {
    lineStream >> lens[dim];
  }

  // dismiss index order string
  std::string columnIndexOrder(scanner.nextLine());

  // create tensor
  auto A( Tcc<TE>::template tensor<F>(lens, name) );

  // read the values only on root
  size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
  std::vector<int64_t> indices(localBufferSize);
  std::vector<F> values(localBufferSize);

  size_t index(0);
  LOG() << "indexCount=" << A->getElementsCount() << std::endl;
  NumberScanner<F> numberScanner(&scanner);
  while (index < A->getElementsCount()) {
    size_t elementsCount(std::min(bufferSize, A->getElementsCount()-index));
    size_t localElementsCount(Cc4s::world->getRank() == 0 ? elementsCount : 0);
    for (size_t i(0); i < localElementsCount; ++i) {
      indices[i] = index+i;
      values[i] = numberScanner.nextNumber();
    }
    // wait until all processes finished reading this buffer into the tensor
    Cc4s::world->barrier();
    LOG() << "writing " << elementsCount << " values to tensor..." << std::endl;
    // FIXME: invoke tcc io
    //A->write(localElementsCount, indices.data(), values.data());
    index += elementsCount;
  }

  return A;
}

template <typename F, typename TE>
void TensorIo::writeBinary(
  const std::string &fileName, const Ptr<Tensor<F,TE>> &A
) {
  MPI_File file;
  MPI_File_open(
    Cc4s::world->getComm(), fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
    MPI_INFO_NULL, &file
  );
  size_t offset(0);

  // truncate possibly existing file
  MPI_File_set_size(file, offset);

  // write dense data
  // FIXME: invoke tcc write to file calls
  // A.write_dense_to_file(file, offset);

  // done
  MPI_File_close(&file);
}

template <typename F, typename TE>
void TensorIo::writeText(
  const std::string &fileName,
  const Ptr<Tensor<F,TE>> &A,
  const std::string &givenRowIndexOrder, const std::string &columnIndexOrder,
  const std::string &delimiter
) {
  // TODO: implement memory scalable version
  size_t elementsCount(Cc4s::world->getRank() == 0 ? A->getElementsCount() : 0);
  std::vector<size_t> indices(elementsCount);
  for (size_t i(0); i < elementsCount; ++i) {
    indices[i] = i;
  }
  std::vector<F> values(elementsCount);
  // FIXME: invoke tcc read
  //A->read(elementsCount, indices.data(), values.data());

  // only the root writes the file
  if (elementsCount > 0) {
    std::ofstream file(fileName.c_str());
    file << A->getName() << delimiter << A->getLens().size();
    for (size_t len: A->getLens()) {
      file << delimiter << len;
    }
    file << '\n';
    file << "" << delimiter << "" << '\n';

    // write the actual data
    file << std::setprecision(16);
    for (size_t index(0); index < elementsCount; ++index) {
      file << values[index++];
      file << '\n';
    }
  }
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> TensorIo::readBinaryHeader(
  MPI_File &file, size_t &offset
) {
  MPI_Status status;
  // reade header
  BinaryTensorHeader header;
  MPI_File_read_at(file, offset, &header, sizeof(header), MPI_BYTE, &status);
  offset += sizeof(header);
  if (strncmp(header.magic, header.MAGIC, sizeof(header.magic)) != 0)
    throw New<Exception>("Invalid file format", SOURCE_LOCATION);
  if (header.version > header.VERSION)
    throw New<Exception>("Incompatible file format version", SOURCE_LOCATION);

  // read dimension headers
  std::vector<size_t> lens(header.order);
  for (unsigned int dim(0); dim < header.order; ++dim) {
    BinaryTensorDimensionHeader dimensionHeader;
    MPI_File_read_at(
      file, offset, &dimensionHeader, sizeof(dimensionHeader), MPI_BYTE, &status
    );
    offset += sizeof(dimensionHeader);
    lens[dim] = dimensionHeader.length;
  }

  // allocate tensor
  return Tcc<TE>::template tensor<F>(lens, "binary");
}



// instantiate
template
Ptr<Tensor<Real<64>,DefaultTensorEngine>>
TensorIo::readBinary<Real<64>,DefaultTensorEngine>(
  const std::string &fileName
);
template
Ptr<Tensor<Complex<64>,DefaultTensorEngine>>
TensorIo::readBinary<Complex<64>,DefaultTensorEngine>(
  const std::string &fileName
);
template
Ptr<Tensor<Real<64>,DryTensorEngine>>
TensorIo::readBinary<Real<64>,DryTensorEngine>(
  const std::string &fileName
);
template
Ptr<Tensor<Complex<64>,DryTensorEngine>>
TensorIo::readBinary<Complex<64>,DryTensorEngine>(
  const std::string &fileName
);
// TODO: 128 bit tensors

template
Ptr<Tensor<Real<64>,DefaultTensorEngine>>
TensorIo::readText<Real<64>,DefaultTensorEngine>(
  const std::string &fileName,
  const std::string &delimiter,
  const size_t bufferSize
);
template
Ptr<Tensor<Complex<64>,DefaultTensorEngine>>
TensorIo::readText<Complex<64>>(
  const std::string &fileName,
  const std::string &delimiter,
  const size_t bufferSize
);
template
Ptr<Tensor<Real<64>,DryTensorEngine>>
TensorIo::readText<Real<64>,DryTensorEngine>(
  const std::string &fileName,
  const std::string &delimiter,
  const size_t bufferSize
);
template
Ptr<Tensor<Complex<64>,DryTensorEngine>>
TensorIo::readText<Complex<64>,DryTensorEngine>(
  const std::string &fileName,
  const std::string &delimiter,
  const size_t bufferSize
);
// TODO: 128 bit tensors

template
void TensorIo::writeBinary<Real<64>,DefaultTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Real<64>,DefaultTensorEngine>> &A
);
template
void TensorIo::writeBinary<Complex<64>,DefaultTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Complex<64>,DefaultTensorEngine>> &A
);
template
void TensorIo::writeBinary<Real<64>,DryTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Real<64>,DryTensorEngine>> &A
);
template
void TensorIo::writeBinary<Complex<64>,DryTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Complex<64>,DryTensorEngine>> &A
);
// TODO: 128 bit tensors

template
void TensorIo::writeText<Real<64>,DefaultTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Real<64>,DefaultTensorEngine>> &A,
  const std::string &rowIndexOrder, const std::string &columnIndexOrder,
  const std::string &delimiter
);
template
void TensorIo::writeText<Complex<64>,DefaultTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Complex<64>,DefaultTensorEngine>> &A,
  const std::string &rowIndexOrder, const std::string &columnIndexOrder,
  const std::string &delimiter
);
template
void TensorIo::writeText<Real<64>,DryTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Real<64>,DryTensorEngine>> &A,
  const std::string &rowIndexOrder, const std::string &columnIndexOrder,
  const std::string &delimiter
);
template
void TensorIo::writeText<Complex<64>,DryTensorEngine>(
  const std::string &fileName,
  const Ptr<Tensor<Complex<64>,DryTensorEngine>> &A,
  const std::string &rowIndexOrder, const std::string &columnIndexOrder,
  const std::string &delimiter
);
// TODO: 128 bit tensors

