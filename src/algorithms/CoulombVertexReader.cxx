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
  LOG(0, "Reader") <<
    "Reading Coulomb vertex from file " << fileName << std::endl;
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
  
  // print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // allocate output tensors
  int vertexLens[] = { NG, Np, Np };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(No, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(Nv, *Cc4s::world, "epsa"));
  Tensor<complex> *GammaGpq(new Tensor<complex>
			    (3, vertexLens, vertexSyms, 
			     *Cc4s::world, "GammaGpq"));

  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex>("CoulombVertex", GammaGpq);

  // real and imaginary parts are read in seperately
  Tensor<> realGammaGpq(3, vertexLens, vertexSyms, 
			*Cc4s::world, "RealGammaGpq");
  Tensor<> imagGammaGpq(3, vertexLens, vertexSyms, 
			*Cc4s::world, "ImagGammaGpq");

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALS_MAGIC, sizeof(chunk.magic)) == 0) {

      // choose reading algorithm
      readGammaGpqChunkSequential(file, realGammaGpq);
      // TODO: find out why Blocked routine produces incorrect numbers with
      // different number of processors.
      //readGammaGpqChunkBlocked(file, realGammaGpq);

    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {

      // choose reading algorithm
      readGammaGpqChunkSequential(file, imagGammaGpq);
      // TODO: find out why Blocked routine produces incorrect numbers with
      // different number of processors.
      //readGammaGpqChunkBlocked(file, imagGammaGpq);

    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      readEpsChunk(file, *epsi, *epsa);
    }
  }
  file.close();

  // combine to complex tensor
  toComplexTensor(realGammaGpq, imagGammaGpq, *GammaGpq);
}

void CoulombVertexReader::dryRun() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading Coulomb vertex from file " << fileName << std::endl;
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
  
  // print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // allocate output tensors
  int vertexLens[] = { NG, Np, Np };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(No));
  DryTensor<> *epsa(new DryVector<>(Nv));
  DryTensor<complex> *GammaGpq(new DryTensor<complex>(3, vertexLens, 
						      vertexSyms));
  // enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex, DryTensor<complex>>("CoulombVertex", 
						       GammaGpq);

  // real and imaginary parts are read in seperately
  DryTensor<> realGammaGpq(3, vertexLens, vertexSyms);
  DryTensor<> imagGammaGpq(3, vertexLens, vertexSyms);
  //realGammaGpq.use();
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void CoulombVertexReader::readGammaGpqChunkBlocked(std::ifstream &file, 
						   Tensor<> &GammaGpq)
{
  LOG(1, "Reader") << "Reading " 
		   << GammaGpq.get_name() << std::endl;

  // TODO: separate distribution from reading
  // allocate local indices and values of the GammaGpq tensors
  int64_t NpPerNode(Np / Cc4s::world->np);
  int64_t NpLocal(Cc4s::world->rank+1 < Cc4s::world->np ? 
		  NpPerNode : Np - Cc4s::world->rank * NpPerNode);
  int64_t NpToSkipBefore(Cc4s::world->rank * NpPerNode);
  int64_t NpToSkipAfter(Np - (NpToSkipBefore + NpLocal));
  double *values(new double[NpLocal*Np*NG]);
  int64_t *indices(new int64_t[NpLocal*Np*NG]);
  file.seekg(sizeof(double)*NpToSkipBefore*Np*NG, file.cur);
  file.read(reinterpret_cast<char *>(values), sizeof(double)*NpLocal*Np*NG);
  file.seekg(sizeof(double)*NpToSkipAfter*Np*NG, file.cur);
  for (int64_t i(0); i < NpLocal*Np*NG; ++i) {
    indices[i] = i + NpToSkipBefore*Np*NG;
  }
  GammaGpq.write(NpLocal*Np*NG, indices, values);
  delete[] values; delete[] indices;
}


// Sequential reading routine with only the master process
void CoulombVertexReader::readGammaGpqChunkSequential(std::ifstream &file, 
						      Tensor<> &GammaGpq) 
{
  int64_t index(0);

  int readRank(getIntegerArgument
		("readRank",No));
  
  for (int p(1); p <= Np; p+=readRank) {
    LOG(1, "Reader") << "Reading " << GammaGpq.get_name()
				  << " at p=" << p << std::endl;
    double *values(new double[std::min(readRank,Np-p)*Np*NG]);
    int64_t *indices(new int64_t[std::min(readRank,Np-p)*Np*NG]);
    if (Cc4s::world->rank == 0) {
      file.read(reinterpret_cast<char *>(values), 
		std::min(readRank,Np-p)*Np*NG*sizeof(double));
      for (int i(0); i < std::min(readRank,Np-p)*Np*NG; ++i) {
        indices[i] = index;
        index++;
      }
    }
    else {
      // skip the data otherwise
      file.seekg(sizeof(double)*std::min(readRank,Np-p)*Np*NG, file.cur);
    }
    int64_t valuesCount(Cc4s::world->rank == 0 ? 
			std::min(readRank,Np-p)*Np*NG : 0);
    GammaGpq.write(valuesCount, indices, values);
    delete[] values; delete[] indices;
  }
}


void CoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> &epsi, Tensor<> &epsa
) {
  LOG(1, "Reader") << "Reading " << epsi.get_name() << ", " 
				<< epsa.get_name() <<  std::endl;

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
  epsi.write(iValuesCount, iIndices, iValues);
  epsa.write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}

