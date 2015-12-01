#include <BinaryFtodReader.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <fstream>

using namespace cc4s;

char const *BinaryFtodReader::Header::MAGIC = "cc4sFTOD";
char const *BinaryFtodReader::Chunk::REALS_MAGIC = "FTODreal";
char const *BinaryFtodReader::Chunk::IMAGS_MAGIC = "FTODimag";
char const *BinaryFtodReader::Chunk::REALSIA_MAGIC = "FTIAreal";
char const *BinaryFtodReader::Chunk::IMAGSIA_MAGIC = "FTIAimag";
char const *BinaryFtodReader::Chunk::EPSILONS_MAGIC = "FTODepsi";

BinaryFtodReader::BinaryFtodReader(bool stridedIo_): stridedIo(stridedIo_) {
}

BinaryFtodReader::~BinaryFtodReader() {
}

/**
 * \brief Reads the Fourier transformed overlap densities from disk.
 */
void BinaryFtodReader::read() {
  LOG(4) <<
    "using " << (stridedIo ? "strided" : "blocked") << " IO" << std::endl;
  LOG(0) <<
    "Reading Fourier transformed overlap densities from binary FTODDUMP...";
  std::ifstream file("FTODDUMP", std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open FTODDUMP file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format of FTODDUMP file");
  nG = header.nG;
  no = header.no;
  nv = header.nv;
  np = no+nv;
  // allocate chi and Coulomb integral tensors
  Cc4s::chiReal = new Chi(nG, no, nv);
  Cc4s::chiImag = new Chi(nG, no, nv);
  Cc4s::chiaiReal = new Chiai(nG, no, nv);
  Cc4s::chiaiImag = new Chiai(nG, no, nv);
  Cc4s::V = new CoulombIntegrals(Cc4s::chiReal, Cc4s::chiImag);

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALS_MAGIC, sizeof(chunk.magic)) == 0) {
      readChiChunk(file, Cc4s::chiReal);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {
      readChiChunk(file, Cc4s::chiImag);
    } else
    if (strncmp(chunk.magic, Chunk::REALSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(0) << "Found ia chunk";
      //readChiChunk(file, Cc4s::chiIAReal);
      readChiaiChunkBlocked(file, Cc4s::chiaiReal);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(0) << "Found ia chunk";
      //readChiChunk(file, Cc4s::chiIAImag);
      readChiaiChunkBlocked(file, Cc4s::chiaiImag);
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      readEpsChunk(file);
    }
  }
  file.close();
  LOG(0) << " OK" << std::endl;

  // TODO: factorize out of reading methods
  double realNorm = Cc4s::chiReal->gpq->norm2();
  double imagNorm = Cc4s::chiImag->gpq->norm2();
  double iNorm = Cc4s::V->i->norm2();
  double aNorm = Cc4s::V->a->norm2();
  LOG(4) <<
    "2-Norm of FTOD = (" << realNorm << "," << imagNorm << ")" << std::endl;
  LOG(4) <<
    "2-Norm of (eps_i,eps_a) = (" << iNorm << "," << aNorm << ")" << std::endl;
}

void BinaryFtodReader::readChiChunk(std::ifstream &file, Chi *chi) {
  if (stridedIo) {
    readChiChunkStrided(file, chi);
  } else {
    readChiChunkBlocked(file, chi);
  }
}

// TODO: use several write calls instead of one big to reduce int64 requirement
void BinaryFtodReader::readChiChunkBlocked(std::ifstream &file, Chi *chi) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  int64_t npPerNode(np / Cc4s::world->np);
  int64_t npLocal(
    Cc4s::world->rank+1 < Cc4s::world->np ?
      npPerNode : np - Cc4s::world->rank * npPerNode
  );
  int64_t npToSkipBefore(Cc4s::world->rank * npPerNode);
  int64_t npToSkipAfter(np - (npToSkipBefore + npLocal));
  double *values(new double[npLocal*np*nG]);
  int64_t *indices(new int64_t[npLocal*np*nG]);
  file.seekg(sizeof(double)*npToSkipBefore*np*nG, file.cur);
  file.read(reinterpret_cast<char *>(values), sizeof(double)*npLocal*np*nG);
  file.seekg(sizeof(double)*npToSkipAfter*np*nG, file.cur);
  for (int64_t i(0); i < npLocal*np*nG; ++i) {
    indices[i] = i + npToSkipBefore*np*nG;
  }
  chi->gpq->write(npLocal*np*nG, indices, values);
  delete[] values; delete[] indices;
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void BinaryFtodReader::readChiaiChunkBlocked(std::ifstream &file, Chiai *chiai) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  int64_t nvPerNode(nv / chiai->nv);
  int64_t nvLocal(
    Cc4s::world->rank+1 < chiai->nv ?
      nvPerNode : nv - Cc4s::world->rank * nvPerNode
  );
  int64_t nvToSkipBefore(Cc4s::world->rank * nvPerNode);
  int64_t nvToSkipAfter(nv - (nvToSkipBefore + nvLocal));
  double *values(new double[nvLocal*no*nG]);
  int64_t *indices(new int64_t[nvLocal*no*nG]);
  file.seekg(sizeof(double)*nvToSkipBefore*no*nG, file.cur);
  file.read(reinterpret_cast<char *>(values), sizeof(double)*nvLocal*no*nG);
  file.seekg(sizeof(double)*nvToSkipAfter*no*nG, file.cur);
  for (int64_t i(0); i < nvLocal*no*nG; ++i) {
    indices[i] = i + nvToSkipBefore*no*nG;
  }
  chiai->gai->write(nvLocal*no*nG, indices, values);
  delete[] values; delete[] indices;
}



// TODO: use intrinsic tensor distribution rather than enforce own distribution
// TODO: use several write calls instead of one big to reduce int64 requirement
void BinaryFtodReader::readChiChunkStrided(std::ifstream &file, Chi *chi) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  // distributed along g: round up g/nprocs 
  int64_t maxValuesCount((nG+Cc4s::world->np-1)/Cc4s::world->np * np*np);
  // each process can have at most maxValuesCount entires
  double *values(new double[maxValuesCount]);
  int64_t *indices(new int64_t[maxValuesCount]);
  int64_t valuesCount(0), index(0);
  for (int q(0); q < np; ++q) {
    for (int p(0); p < np; ++p) {
      for (int g(0); g < nG; ++g) {
        // TODO: use counter to avoid division
        // distributed along g: current rank is responsible, only
        if (g % Cc4s::world->np == Cc4s::world->rank) {
          file.read(
            reinterpret_cast<char *>(&values[valuesCount]), sizeof(double)
          );
          indices[valuesCount] = index;
          ++valuesCount;
        } else {
          // skip the data otherwise
          file.seekg(sizeof(double), file.cur);
        }
        ++index;
      }
    }
  }
  chi->gpq->write(valuesCount, indices, values);
  delete[] values; delete[] indices;
}

void BinaryFtodReader::readEpsChunk(std::ifstream &file) {
  // allocate local indices and values of eigenenergies
  double *iValues(new double[no]);
  double *aValues(new double[nv]);
  int64_t *iIndices(new int64_t[no]);
  int64_t *aIndices(new int64_t[nv]);

  if (Cc4s::world->rank == 0) {
    file.read(reinterpret_cast<char *>(iValues), no*sizeof(double));
    for (int i(0); i < no; ++i) iIndices[i] = i;
    file.read(reinterpret_cast<char *>(aValues), nv*sizeof(double));
    for (int a(0); a < nv; ++a) aIndices[a] = a;
  } else {
    // skip the data otherwise
    file.seekg(sizeof(double)*np, file.cur);
  }
  int64_t iValuesCount(Cc4s::world->rank == 0 ? no : 0);
  int64_t aValuesCount(Cc4s::world->rank == 0 ? nv : 0);
  Cc4s::V->i->write(iValuesCount, iIndices, iValues);
  Cc4s::V->a->write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

void BinaryFtodReader::write() {
  // collect the Ftod data on all nodes (required by read_all)
  // NOTE: not memory scalable
  no = Cc4s::chiReal->no;
  nv = Cc4s::chiReal->nv;
  nG = Cc4s::chiReal->nG;
  np = no+nv;
  double *reals = new double[nG*np*np];
  double *imags = new double[nG*np*np];
  double *epsilons = new double[np];
  Cc4s::chiReal->gpq->read_all(reals);
  Cc4s::chiImag->gpq->read_all(imags);
  Cc4s::V->i->read_all(&epsilons[0]);
  Cc4s::V->a->read_all(&epsilons[no]);

  // only the root writes to file
  if (Cc4s::world->rank == 0) {
    std::ofstream file("FTODDUMP", std::ios::binary);
    // write header
    Header header;
    strncpy(header.magic, Header::MAGIC, 8);
    header.nG = nG;
    header.no = no;
    header.nv = nv;
    header.nSpins = 1;
    header.kPoints = 1;
    header.reserved_ = 0;
    file.write(reinterpret_cast<const char *>(&header), sizeof(header));

    // write reals chunk
    Chunk chunk;
    strncpy(chunk.magic, Chunk::REALS_MAGIC, sizeof(header.magic));
    chunk.size = sizeof(chunk)+sizeof(double)*(nG*np*np);
    file.write(reinterpret_cast<const char *>(&chunk), sizeof(chunk));
    file.write(reinterpret_cast<const char *>(reals), chunk.size-sizeof(chunk));
    // write imags chunk
    strncpy(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic));
    file.write(reinterpret_cast<const char *>(&chunk), sizeof(chunk));
    file.write(reinterpret_cast<const char *>(imags), chunk.size-sizeof(chunk));
    // write epsilons chunk
    strncpy(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic));
    chunk.size = sizeof(chunk)+sizeof(double)*np;
    file.write(reinterpret_cast<const char *>(&chunk), sizeof(chunk));
    file.write(
      reinterpret_cast<const char *>(epsilons), chunk.size-sizeof(chunk)
    );
    file.close();
  }
  delete[] reals; delete[] imags; delete[] epsilons;
}

