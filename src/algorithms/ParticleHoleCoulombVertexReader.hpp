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

