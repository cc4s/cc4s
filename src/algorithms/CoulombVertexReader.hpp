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
   * data file, and stores them in the CTF Tensors GammaGqr, epsi, epsa.
   */
  class CoulombVertexReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexReader);
    CoulombVertexReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexReader();

    /**
     * \brief Reads the Full Coulomb Vertex \f$\Gamma_{pG}^q\f$ and the occupied and
     * virtual orbital energies \f$\varepsilon_i, \varepsilon_a\f$ from binary
     * data file, and stores them in the CTF Tensors GammaGqr, epsi, epsa.
     */
    virtual void run();

    /**
     * \brief Dry run for reading the Coulomb vertex from binary data file.
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
  protected:
    void handleUnrestricted();
    void unrestrictVertex();
    void unrestrictEigenEnergies(const std::string &name);
  };
}

#endif

