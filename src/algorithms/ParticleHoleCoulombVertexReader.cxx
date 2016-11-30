#include <algorithms/ParticleHoleCoulombVertexReader.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
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

void ParticleHoleCoulombVertexReader::run() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading particle hole Coulomb vertex from file " << fileName << " ...";
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  NG = header.NG;
  No = header.No;
  Nv = header.Nv;
  Np = No + Nv;

  // allocate output tensors
  int vertexLens[] = { NG, Nv, No };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(No, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(Nv, *Cc4s::world, "epsa"));
  Tensor<complex> *GammaGai(
    new Tensor<complex>(
      3, vertexLens, vertexSyms, *Cc4s::world, "GammaGai"
    )
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex>("ParticleHoleCoulombVertex", GammaGai);

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
  toComplexTensor(realGammaGai, imagGammaGai, *GammaGai);
  LOG(0, "Reader") << " OK" << std::endl;

  // print the number of NG's, Nv's, No's, and Np's at the level of LOG(1)
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;
}

void ParticleHoleCoulombVertexReader::dryRun() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading particle hole Coulomb vertex from file " << fileName << std::endl;
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  file.close();
  NG = header.NG;
  No = header.No;
  Nv = header.Nv;
  Np = No + Nv;

  // allocate output tensors
  int vertexLens[] = { NG, Nv, No };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(No, SOURCE_LOCATION));
  DryTensor<> *epsa(new DryVector<>(Nv, SOURCE_LOCATION));
  DryTensor<complex> *GammaGai(
    new DryTensor<complex>(3, vertexLens, vertexSyms, SOURCE_LOCATION)
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex>("ParticleHoleCoulombVertex", GammaGai);

  // real and imaginary parts are read in seperately
  DryTensor<> realGammaGai(3, vertexLens, vertexSyms, SOURCE_LOCATION);
  DryTensor<> imagGammaGai(3, vertexLens, vertexSyms, SOURCE_LOCATION);
  dryReadGammaGaiChunkBlocked(&realGammaGai);
  dryReadGammaGaiChunkBlocked(&imagGammaGai);

  // print the number of NG's, Nv's, No's, and Np's at the level of LOG(1)
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void ParticleHoleCoulombVertexReader::readGammaGaiChunkBlocked(
  std::ifstream &file, Tensor<> *GammaGai
) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  int64_t NvPerNode(Nv / Cc4s::world->np);
  int64_t NvLocal(
    Cc4s::world->rank+1 < Cc4s::world->np ?
      NvPerNode : Nv - Cc4s::world->rank * NvPerNode
  );
  int64_t NvToSkipBefore(Cc4s::world->rank * NvPerNode);
  int64_t NvToSkipAfter(Nv - (NvToSkipBefore + NvLocal));
  double *values(new double[NvLocal*No*NG]);
  int64_t *indices(new int64_t[NvLocal*No*NG]);
  file.seekg(sizeof(double)*NvToSkipBefore*No*NG, file.cur);
  file.read(reinterpret_cast<char *>(values), sizeof(double)*NvLocal*No*NG);
  file.seekg(sizeof(double)*NvToSkipAfter*No*NG, file.cur);
  for (int64_t i(0); i < NvLocal*No*NG; ++i) {
    indices[i] = i + NvToSkipBefore*No*NG;
  }
  GammaGai->write(NvLocal*No*NG, indices, values);
  delete[] values; delete[] indices;
}

void ParticleHoleCoulombVertexReader::dryReadGammaGaiChunkBlocked(
  DryTensor<> *GammaGai
) {
  int64_t elements(Nv*No*NG);
  // values array
  DryMemory::allocate(sizeof(double)*elements, SOURCE_LOCATION);
  // indices array
  DryMemory::allocate(sizeof(int64_t)*elements, SOURCE_LOCATION);
  DryMemory::free(sizeof(int64_t)*elements);
  DryMemory::free(sizeof(double)*elements);
}

void ParticleHoleCoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> *epsi, Tensor<> *epsa
) {
  // allocate local indices and values of eigenenergies
  double *iValues(new double[No]);
  double *aValues(new double[Nv]);
  int64_t *iIndices(new int64_t[No]);
  int64_t *aIndices(new int64_t[Nv]);

  if (Cc4s::world->rank == 0) {
    file.read(reinterpret_cast<char *>(iValues), No*sizeof(double));
    for (int i(0); i < No; ++i) iIndices[i] = i;
    file.read(reinterpret_cast<char *>(aValues), Nv*sizeof(double));
    for (int a(0); a < Nv; ++a) aIndices[a] = a;
  } else {
    // skip the data otherwise
    file.seekg(sizeof(double)*Np, file.cur);
  }
  int64_t iValuesCount(Cc4s::world->rank == 0 ? No : 0);
  int64_t aValuesCount(Cc4s::world->rank == 0 ? Nv : 0);
  epsi->write(iValuesCount, iIndices, iValues);
  epsa->write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

