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
    "Reading Coulomb vertex from file " << fileName << std::endl;
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
  
  // print nG, no, nv, np
  LOG(1, "CoulombVertexReader") << "nG=" << nG << std::endl;
  LOG(1, "CoulombVertexReader") << "no=" << no << std::endl;
  LOG(1, "CoulombVertexReader") << "nv=" << nv << std::endl;
  LOG(1, "CoulombVertexReader") << "np=" << np << std::endl;

  // allocate output tensors
  int vertexLens[] = { nG, np, np };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(no, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(nv, *Cc4s::world, "epsa"));
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

      //readGammaGpqChunkBlocked(file, realGammaGpq);
      readGammaGpqChunkSequential(file, realGammaGpq);

    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {

      //readGammaGpqChunkBlocked(file, imagGammaGpq);
      readGammaGpqChunkSequential(file, imagGammaGpq);

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
  LOG(1, "CoulombVertexReader") << "nG=" << nG << std::endl;
  LOG(1, "CoulombVertexReader") << "no=" << no << std::endl;
  LOG(1, "CoulombVertexReader") << "nv=" << nv << std::endl;
  LOG(1, "CoulombVertexReader") << "np=" << np << std::endl;

  // allocate output tensors
  int vertexLens[] = { nG, np, np };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(no));
  DryTensor<> *epsa(new DryVector<>(nv));
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
void CoulombVertexReader::readGammaGpqChunkBlocked(
  std::ifstream &file, Tensor<> &GammaGpq
) {
  LOG(1, "CoulombVertexReader") << "Reading " 
				<< GammaGpq.get_name() << std::endl;

  // TODO: separate distribution from reading
  // allocate local indices and values of the GammaGpq tensors
  int64_t npPerNode(np / Cc4s::world->np);
  int64_t npLocal(Cc4s::world->rank+1 < Cc4s::world->np ? 
		  npPerNode : np - Cc4s::world->rank * npPerNode);
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
  delete[] values; delete[] indices;
}


// Sequential reading routine with only the master process
void CoulombVertexReader::readGammaGpqChunkSequential(std::ifstream &file, 
						      Tensor<> &GammaGpq) 
{
  int64_t index(0);

  int readRank(getIntegerArgument
		("readRank",no));
  
  for (int p(1); p <= np; p+=readRank) {
    LOG(1, "CoulombVertexReader") << "Reading " << GammaGpq.get_name()
				  << " at p=" << p << std::endl;
    double *values(new double[std::min(readRank,np-p)*np*nG]);
    int64_t *indices(new int64_t[std::min(readRank,np-p)*np*nG]);
    if (Cc4s::world->rank == 0) {
      file.read(reinterpret_cast<char *>(values), 
		std::min(readRank,np-p)*np*nG*sizeof(double));
      for (int i(0); i < std::min(readRank,np-p)*np*nG; ++i) {
        indices[i] = index;
        index++;
      }
    }
    else {
      // skip the data otherwise
      file.seekg(sizeof(double)*std::min(readRank,np-p)*np*nG, file.cur);
    }
    int64_t valuesCount(Cc4s::world->rank == 0 ? 
			std::min(readRank,np-p)*np*nG : 0);
    GammaGpq.write(valuesCount, indices, values);
    delete[] values; delete[] indices;
  }
}


void CoulombVertexReader::readEpsChunk(
  std::ifstream &file, Tensor<> &epsi, Tensor<> &epsa
) {
  LOG(1, "CoulombVertexReader") << "Reading " << epsi.get_name() << ", " 
				<< epsa.get_name() <<  std::endl;

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

