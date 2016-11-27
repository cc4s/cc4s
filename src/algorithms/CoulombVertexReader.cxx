#include <algorithms/CoulombVertexReader.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
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
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Cc4s::world->comm, fileName.c_str(), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &file
    )
  );
  if (mpiError) throw new Exception("Failed to open file");
  // Read header: obtain NG,No,Nv,Np
  Header header;
  MPI_Status status;
  MPI_File_read(file, &header, sizeof(header), MPI_BYTE, &status);
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  int NG(header.NG);
  int No(header.No);
  int Nv(header.Nv);
  int Np(No + Nv);
  
  // Print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // Allocate output tensors
  int vertexLens[] = { NG, Np, Np };
  int vertexSyms[] = { NS, NS, NS };
  Tensor<> *epsi(new Vector<>(No, *Cc4s::world, "epsi"));
  Tensor<> *epsa(new Vector<>(Nv, *Cc4s::world, "epsa"));
  Tensor<complex> *GammaGqr(
    new Tensor<complex>(3, vertexLens, vertexSyms, *Cc4s::world, "GammaGqr")
  );

  // Enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex>("CoulombVertex", GammaGqr);

  // Real and imaginary parts are read in seperately
  Tensor<> realGammaGqr(
    3, vertexLens, vertexSyms, *Cc4s::world, "RealGammaGqr"
  );
  Tensor<> imagGammaGqr(
    3, vertexLens, vertexSyms, *Cc4s::world, "ImagGammaGqr"
  );

  int64_t offset(sizeof(header));
  MPI_Offset fileSize;
  MPI_File_get_size(file, &fileSize);
  Chunk chunk;
  while (offset < fileSize) {
    MPI_File_read_at(file, offset, &chunk, sizeof(chunk), MPI_BYTE, &status);
    if (strncmp(chunk.magic, Chunk::REALS_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(1, "Reader") << "reading " << realGammaGqr.get_name() << std::endl;
      realGammaGqr.read_dense_from_file(file, offset+sizeof(chunk));
    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(1, "Reader") << "reading " << imagGammaGqr.get_name() << std::endl;
      imagGammaGqr.read_dense_from_file(file, offset+sizeof(chunk));
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(1, "Reader") << "reading " << epsi->get_name() << ", "
        << epsa->get_name() << std::endl;
      epsi->read_dense_from_file(file, offset+sizeof(chunk));
      epsa->read_dense_from_file(file, offset+sizeof(chunk)+No*sizeof(double));
    }
    offset += chunk.size;
  }
  MPI_File_close(&file);

  // Combine to complex tensor
  toComplexTensor(realGammaGqr, imagGammaGqr, *GammaGqr);
}

void CoulombVertexReader::dryRun() {
  std::string fileName(getTextArgument("file"));
  LOG(0, "Reader") <<
    "Reading Coulomb vertex from file " << fileName << std::endl;
  std::ifstream file(fileName.c_str(), std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open file");
  // Read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format");
  file.close();
  int NG(header.NG);
  int No(header.No);
  int Nv(header.Nv);
  int Np(No + Nv);
  
  // Print NG, No, Nv, Np
  LOG(1, "Reader") << "NG=" << NG << std::endl;
  LOG(1, "Reader") << "No=" << No << std::endl;
  LOG(1, "Reader") << "Nv=" << Nv << std::endl;
  LOG(1, "Reader") << "Np=" << Np << std::endl;

  // Allocate output tensors
  int vertexLens[] = { NG, Np, Np };
  int vertexSyms[] = { NS, NS, NS };
  DryTensor<> *epsi(new DryVector<>(No, SOURCE_LOCATION));
  DryTensor<> *epsa(new DryVector<>(Nv, SOURCE_LOCATION));
  DryTensor<complex> *GammaGqr(
    new DryTensor<complex>(3, vertexLens, vertexSyms, SOURCE_LOCATION)
  );
  // Enter the allocated data (and by that type the output data to tensors)
  allocatedTensorArgument("HoleEigenEnergies", epsi);
  allocatedTensorArgument("ParticleEigenEnergies", epsa);
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "CoulombVertex", GammaGqr
  );

  // Real and imaginary parts are read in seperately
  DryTensor<> realGammaGqr(3, vertexLens, vertexSyms, SOURCE_LOCATION);
  DryTensor<> imagGammaGqr(3, vertexLens, vertexSyms, SOURCE_LOCATION);
  //  realGammaGqr.use();
}

