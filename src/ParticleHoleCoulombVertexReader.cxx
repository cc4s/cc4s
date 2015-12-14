#include <ParticleHoleCoulombVertexReader.hpp>
#include <util/ComplexTensor.hpp>
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

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombVertexReader);

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
  if (!file.is_open()) throw new Exception("Failed to open file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  nG = header.nG;
  no = header.no;
  nv = header.nv;
  np = no + nv;

  // allocate output tensors
  int vertexLens[] = { nG, nv, no };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(no, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(nv, *Cc4s::world, "epsa"));
  Tensor<complex> *gammaGai(
    new Tensor<complex>(
      3, vertexLens, vertexSyms, *Cc4s::world, "GammaGai"
    )
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument("ParticleHoleCoulombVertex", gammaGai);

  // real and imaginary parts are read in seperately
  Tensor<> realGammaGai(
    3, vertexLens, vertexSyms, *Cc4s::world, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, vertexLens, vertexSyms, *Cc4s::world, "ImagGammaGai"
  );

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      readGammaGaiChunkBlocked(file, &realGammaGai);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      readGammaGaiChunkBlocked(file, &imagGammaGai);
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      readEpsChunk(file, epsi, epsa);
    }
  }
  file.close();
  // combine to complex tensor
  toComplexTensor(realGammaGai, imagGammaGai, *gammaGai);
  LOG(0) << " OK" << std::endl;
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void ParticleHoleCoulombVertexReader::readGammaGaiChunkBlocked(
  std::ifstream &file, Tensor<> *gammaGai
) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
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
  gammaGai->write(nvLocal*no*nG, indices, values);
  delete[] values; delete[] indices;
}


void ParticleHoleCoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> *epsi, Tensor<> *epsa
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
  epsi->write(iValuesCount, iIndices, iValues);
  epsa->write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

