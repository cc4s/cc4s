#include <algorithms/CoulombVertexReader.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <fstream>

using namespace cc4s;
using namespace CTF;

char const *CoulombVertexReader::Header::MAGIC = "cc4sFTOD";
char const *CoulombVertexReader::Chunk::REALS_MAGIC = "FTODreal";
char const *CoulombVertexReader::Chunk::IMAGS_MAGIC = "FTODimag";
char const *CoulombVertexReader::Chunk::EPSILONS_MAGIC = "FTODepsi";

ALGORITHM_REGISTRAR_DEFINITION(CoulombVertexReader);

CoulombVertexReader::CoulombVertexReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
  
}

CoulombVertexReader::~CoulombVertexReader() {
}

void CoulombVertexReader::run() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "CoulombVertexReader") <<
    "Reading Coulomb vertex from file " << fileName << " ..." << std::endl;
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
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
  int vertexLens[] = { nG, np, np };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(no, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(nv, *Cc4s::world, "epsa"));
  Tensor<complex> *GammaGpq(
    new Tensor<complex>(
      3, vertexLens, vertexSyms, *Cc4s::world, "GammaGpq"
    )
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex>("CoulombVertex", GammaGpq);

  // real and imaginary parts are read in seperately
  Tensor<> realGammaGpq(
    3, vertexLens, vertexSyms, *Cc4s::world, "RealGammaGpq"
  );
  Tensor<> imagGammaGpq(
    3, vertexLens, vertexSyms, *Cc4s::world, "ImagGammaGpq"
  );

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALS_MAGIC, sizeof(chunk.magic)) == 0) {
      readGammaGpqChunkBlocked(file, realGammaGpq);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {
      readGammaGpqChunkBlocked(file, imagGammaGpq);
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      readEpsChunk(file, *epsi, *epsa);
    }
  }
  file.close();
  // combine to complex tensor
  toComplexTensor(realGammaGpq, imagGammaGpq, *GammaGpq);
  
  // Print Ok
  //LOG(0, "CoulombVertexReader") << " OK" << std::endl;
  
  // print nG, no, nv, np
  LOG(1, "CoulombVertexReader") << "nG = " << nG << std::endl;
  LOG(1, "CoulombVertexReader") << "no = " << no << std::endl;
  LOG(1, "CoulombVertexReader") << "nv = " << nv << std::endl;
  LOG(1, "CoulombVertexReader") << "np = " << np << std::endl;

  // Test print the norm of GammaGpq
  //double result(realGammaGpq.norm2());
  //LOG(4) << "|GammaGpq| = " << result << std::endl;
  //print the norm of epsi and epsa
  //result = epsi->norm2();
  //LOG(4) << "|epsi| = " << result << std::endl;
  //result = epsa->norm2();
  //LOG(4) << "|epsa| = " << result << std::endl;

  // Test print the norm of Vpqrs
  //int lens[] = { np, np, np, np };
  //int syms[] = { NS, NS, NS, NS };
  //Tensor<> vpqrs(4, lens, syms, *Cc4s::world, "Vpqrs");
  //vpqrs["prqs"]  = realGammaGpq["Gpq"] * realGammaGpq["Grs"];
  //vpqrs["prqs"] += imagGammaGpq["Gpq"] * imagGammaGpq["Grs"];
  //double error(vpqrs.norm2());
  //LOG(4) << "|Vpqrs| = " << error << std::endl;
}

void CoulombVertexReader::dryRun() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "CoulombVertexReader") <<
    "Reading Coulomb vertex from file " << fileName << std::endl;
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  file.close();
  nG = header.nG;
  no = header.no;
  nv = header.nv;
  np = no + nv;
  
  // print nG, no, nv, np
  LOG(1, "CoulombVertexReader") << "nG = " << nG << std::endl;
  LOG(1, "CoulombVertexReader") << "no = " << no << std::endl;
  LOG(1, "CoulombVertexReader") << "nv = " << nv << std::endl;
  LOG(1, "CoulombVertexReader") << "np = " << np << std::endl;

  // allocate output tensors
  int vertexLens[] = { nG, np, np };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(no));
  DryTensor<> *epsa(new DryVector<>(nv));
  DryTensor<complex> *GammaGpq(
    new DryTensor<complex>(3, vertexLens, vertexSyms)
  );
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<
    complex, DryTensor<complex>
  >("CoulombVertex", GammaGpq);

  // real and imaginary parts are read in seperately
  DryTensor<> realGammaGpq(3, vertexLens, vertexSyms);
  DryTensor<> imagGammaGpq(3, vertexLens, vertexSyms);
  //realGammaGpq.use();
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void CoulombVertexReader::readGammaGpqChunkBlocked(
  std::ifstream &file, Tensor<> &GammaGpq
) {
  //LOG(4) << "Reading " << GammaGpq.get_name() << std::endl;
  // TODO: separate distribution from reading
  // allocate local indices and values of the GammaGpq tensors
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
  GammaGpq.write(npLocal*np*nG, indices, values);
  //for (int ind=0; ind<npLocal*np*nG; ind++)
  //{
  //LOG(6) << values[ind] << std::endl;
  //}
  delete[] values; delete[] indices;
}


void CoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> &epsi, Tensor<> &epsa
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
  epsi.write(iValuesCount, iIndices, iValues);
  epsa.write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

