/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_READER_DEFINED
#define COULOMB_VERTEX_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>
#include <cstdint>
#include <fstream>

namespace cc4s {
  /**
   * \brief Reads the Coulomb vertex \f$\Gamma_{pG}^q\f$ and the occupied and
   * virtual orbital energies \f$\varepsilon_i, \varepsilon_a\f$ from binary
   * data file, and stores them in the CTF Tensors GammaGpq, epsi, epsa.
   */
  class CoulombVertexReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexReader);
    CoulombVertexReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexReader();

    /**
     * \brief Reads the Full Coulomb Vertex GammaGpq from FTODDUMP file.
     */
    virtual void run();

    /** \brief The occupied orbital energies  */
    CTF::Tensor<> *epsi;
    /** \brief The virtual orbital energies  */
    CTF::Tensor<> *epsa;
    /** \brief The Coulomb Vertex GammaGpq  */
    CTF::Tensor<complex> *GammaGpq;

    /**
     * \brief Reads the size of the Full Coulomb Vertex GammaGpq from FTODDUMP
     * file.
     */
    virtual void dryRun();

  protected:
    /** \brief number of occupied orbitals */
    int no;
    /** \brief number of unoccupied orbitals */
    int nv;
    /** \brief number of occupied plus unoccupied orbitals */
    int np;
    /** \brief number of G-vectors */
    int nG;

    /**
     * \brief Reads a chunk from the Full Coulomb Vertex GammaGpq
     * \param[in] file name of file of Coulomb Vertex
     * \param[in,out] GammaGpq tensor to read the Coulomb vertex
     */
    void readGammaGpqChunkBlocked(std::ifstream &file, CTF::Tensor<> &GammaGpq);
    /**
     * \brief Reads a chunk from the Full Coulomb Vertex GammaGpq sequentially
     * using one processor only.
     * \param[in] file name of file of Coulomb Vertex
     * \param[in,out] GammaGpq tensor to read the Coulomb vertex
     */
    void readGammaGpqChunkSequential(std::ifstream &file, CTF::Tensor<> &GammaGpq);
    /**
     * \brief Reads a chunk from the orbital energies epsi, epsa
     * \param[in] file name of file of Coulomb Vertex
     * \param[in,out] epsi energies of occupied orbitals
     * \param[in,out] epsa energies of unoccupied orbitals
     */
    void readEpsChunk(
      std::ifstream &file,  CTF::Tensor<> &epsi, CTF::Tensor<> &epsa
    );

    class Header {
    public:
      char magic[8];
      int32_t no, nv, nG, nSpins, kPoints, reserved_;
      static char const *MAGIC;
    };
    class Chunk {
    public:
      char magic[8];
      int64_t size;
      static char const *REALS_MAGIC;
      static char const *IMAGS_MAGIC;
      static char const *REALSIA_MAGIC;
      static char const *IMAGSIA_MAGIC;
      static char const *EPSILONS_MAGIC;
    };
  };
}

#endif

