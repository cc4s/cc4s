/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_VERTEX_READER_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/DryTensor.hpp>
#include <ctf.hpp>
#include <cstdint>
#include <fstream>

namespace cc4s {
  /**
   * \brief Reads the particle-hole Coulomb vertex \f$\Gamma_{iG}^a\f$ and
   * the occupied and
   * virtual orbital energies \f$\varepsilon_i, \varepsilon_a\f$ from binary
   * data file, and stores them in the CTF Tensors GammaGai, epsi, epsa.
   */
  class ParticleHoleCoulombVertexReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ParticleHoleCoulombVertexReader);
    ParticleHoleCoulombVertexReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~ParticleHoleCoulombVertexReader();
    /**
     * \brief Reads the Full Coulomb Vertex \f$\Gamma_{iG}^a\f$ from FTODDUMP file.
     */
    virtual void run();
    /**
     * \brief Performs a dry run on reading the
     * binary Coulomb vertex file from disk.
     * Only the header information is read containing size information.
     */
    virtual void dryRun();

  protected:
    /** \brief number of occupied orbitals */
    int No;
    /** \brief number of unoccupied orbitals */
    int Nv;
    /** \brief number of occupied plus unoccupied orbitals */
    int Np;
    /** \brief number of G-vectors */
    int NG;

    /**
     * \brief Reads a chunk from the ParticleHole Coulomb Vertex GammaGai
     * \param[in] file name of file of Coulomb Vertex
     * \param[in,out] GammaGai tensor to read the ParticleHole Coulomb Vertex
     */
    void readGammaGaiChunkBlocked(std::ifstream &file, CTF::Tensor<> *GammaGai);
    /**
     * \brief Performs a dry run of reading the ParticleHole Coulomb Vertex.
     * \param[in,out] GammaGai DryTensor to read the ParticleHole Coulomb Vertex
     */
    void dryReadGammaGaiChunkBlocked(DryTensor<> *GammaGai);
    /**
     * \brief Reads a chunk from the orbital energies epsi, epsa
     * \param[in] file name of file of Coulomb Vertex
     * \param[in,out] epsi energies of occupied orbitals
     * \param[in,out] epsa energies of unoccupied orbitals
     */
    void readEpsChunk(
      std::ifstream &file,  CTF::Tensor<> *epsi, CTF::Tensor<> *epsa
    );

    class Header {
    public:
      char magic[8];
      int32_t No, Nv, NG, NSpins, kPoints, reserved_;
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

