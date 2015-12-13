#include <ParticleHoleCoulombVertexReader.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <fstream>

using namespace cc4s;
using namespace CTF;

char const *ParticleHoleCoulombVertexReader::Header::MAGIC = "cc4sFTOD";
char const *ParticleHoleCoulombVertexReader::Chunk::REALS_MAGIC = "FTODreal";
char const *ParticleHoleCoulombVertexReader::Chunk::IMAGS_MAGIC = "FTODimag";
char const *ParticleHoleCoulombVertexReader::Chunk::REALSIA_MAGIC = "FTIAreal";
char const *ParticleHoleCoulombVertexReader::Chunk::IMAGSIA_MAGIC = "FTIAimag";
char const *ParticleHoleCoulombVertexReader::Chunk::EPSILONS_MAGIC = "FTODepsi";

ParticleHoleCoulombVertexReader::ParticleHoleCoulombVertexReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
  
}

ParticleHoleCoulombVertexReader::~ParticleHoleCoulombVertexReader() {
}

/**
 * \brief Reads the Fourier transformed overlap densities from disk.
 */
void ParticleHoleCoulombVertexReader::run() {
  std::string fileName(getTextArgument("file"));
  LOG(0) <<
    "Reading particle hole Coulomb vertex from file " << fileName << " ...";
  std::ifstream file(fileName, std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open FTODDUMP file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format of FTODDUMP file");
  nG = header.nG;
  no = header.no;
  nv = header.nv;
  np = no + nv;

  // allocate output tensors
  int vertexLens[] = { nG, nv, no };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *iEps(new Vector<>(no, *Cc4s::world, "iEps"));
  Tensor<> *aEps(new Vector<>(nv, *Cc4s::world, "aEps"));
  Tensor<> *aiCoulombVertexReal(
    new Tensor<>(
      3, vertexLens, vertexSyms, *Cc4s::world, "SvRgai"
    )
  );
  Tensor<> *aiCoulombVertexImag(
    new Tensor<>(
      3, vertexLens, vertexSyms, *Cc4s::world, "SvIgai"
    )
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("iEps", iEps);
  allocatedTensorArgument("aEps", aEps);
  allocatedTensorArgument("aiCoulombVertexReal", aiCoulombVertexReal);
  allocatedTensorArgument("aiCoulombVertexImag", aiCoulombVertexImag);

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(4) << "Found ia chunk. ";
      //readChiChunk(file, Cc4s::chiIAReal);
      readChiAiChunkBlocked(file, aiCoulombVertexReal);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(4) << "Found ia chunk. ";
      //readChiChunk(file, Cc4s::chiIAImag);
      readChiAiChunkBlocked(file, aiCoulombVertexImag);
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(4) << "Found eps chunk.";
      readEpsChunk(file, iEps, aEps);
    }
  }
  file.close();
  LOG(0) << " OK" << std::endl;
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void ParticleHoleCoulombVertexReader::readChiAiChunkBlocked(
  std::ifstream &file, Tensor<> *chiAi
) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  // FIXME: continue here...
  int64_t nvPerNode(nv / Cc4s::world->np);
  int64_t nvLocal(
    Cc4s::world->rank+1 < Cc4s::world->np ?
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
  chiAi->write(nvLocal*no*nG, indices, values);
  delete[] values; delete[] indices;
}


void ParticleHoleCoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> *ieps, Tensor<> *aeps
) {
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
  ieps->write(iValuesCount, iIndices, iValues);
  aeps->write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

