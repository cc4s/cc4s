#include <TextFtodReader.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <Exception.hpp>
#include <ctf.hpp>
#include <fstream>

using namespace cc4s;

/**
 * \brief Reads the Fourier transformed overlap densities from disk.
 */
void TextFtodReader::read() {
  LOG(0) <<
    "Reading Fourier transformed overlap densities from text FTOD...";
  std::ifstream file("FTOD");
  if (!file.is_open()) throw new Exception("Failed to open FTOD file");
  std::string line;
  // read a comment line
  std::getline(file, line);
  // NOTE: currently unused
  int nG, no, nv, nSpins, nk;
  // read the size data
  std::getline(file, line);
  std::stringstream lineStream(line);
  lineStream >> no >> nv >> nG >> nSpins >> nk;

  // allocate chi and Coulomb integral tensors
  Cc4s::chiReal = new Chi(nG, no, nv);
  Cc4s::chiImag = new Chi(nG, no, nv);
  Cc4s::V = new CoulombIntegrals(Cc4s::chiReal, Cc4s::chiImag);

  // allocate local indices and values of the chi tensors
  int64_t np(no+nv);
  // distributed along g: round up g/nprocs 
  int64_t maxValuesCount((nG+Cc4s::world->np-1)/Cc4s::world->np * np*np);
  // each process can have at most maxValuesCount entires
  double *reals(new double[maxValuesCount]);
  double *imags(new double[maxValuesCount]);
  int64_t *indices(new int64_t[maxValuesCount]);
  int64_t valuesCount(0);
  // allocate local indices and values of eigenenergies
  double *iValues(new double[no]);
  double *aValues(new double[nv]);
  int64_t *iIndices(new int64_t[no]);
  int64_t *aIndices(new int64_t[nv]);

  // read another comment line
  std::getline(file, line);
  // start reading the file line by line
  while (std::getline(file, line)) {
    std::stringstream lineStream(line);
    double real, imag;
    int g, p, q, spin;
    // parse the line
    lineStream >> real >> imag >> g >> p >> q >> spin;
    if (g > 0) {
      // chi_q^p(g), spin is ignored
      // g, p and q are zero based
      --g; --p; --q;
      // distributed along g: current rank is responsible, only
      if (g % Cc4s::world->np == Cc4s::world->rank) {
        reals[valuesCount] = real;
        imags[valuesCount] = imag;
        indices[valuesCount] = g + nG*(p + np*q);
        ++valuesCount;
      }
    } else {
      // eigenenergy with eps_p, all other indices are to be ignored
      // they are written only on the root
      if (Cc4s::world->rank == 0) {
        --p;
        if (p < Cc4s::chiReal->no) {
          iValues[p] = real;
          iIndices[p] = p;
        } else {
          aValues[p-no] = real;
          aIndices[p-no] = p-no;
        }
      }
    }
  }
  Cc4s::chiReal->gpq->write(valuesCount, indices, reals);
  Cc4s::chiImag->gpq->write(valuesCount, indices, imags);
  int64_t iValuesCount(Cc4s::world->rank == 0 ? no : 0);
  int64_t aValuesCount(Cc4s::world->rank == 0 ? nv : 0);
  Cc4s::V->i->write(iValuesCount, iIndices, iValues);
  Cc4s::V->a->write(aValuesCount, aIndices, aValues);
  delete[] indices; delete[] reals; delete[] imags;
  delete[] iIndices; delete[] aIndices; delete[] iValues; delete[] aValues;
  file.close();
  LOG(0) << " OK" << std::endl;

  double realNorm = Cc4s::chiReal->gpq->norm2();
  double imagNorm = Cc4s::chiImag->gpq->norm2();
  double iNorm = Cc4s::V->i->norm2();
  double aNorm = Cc4s::V->a->norm2();
  LOG(4) <<
    "2-Norm of FTOD = (" << realNorm << "," << imagNorm << ")" << std::endl;
  LOG(4) <<
    "2-Norm of (eps_i,eps_a) = (" << iNorm << "," << aNorm << ")" << std::endl;
}

