/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CC4S_DEFINED
#define CC4S_DEFINED

#include <util/SharedPointer.hpp>
#include <util/MpiCommunicator.hpp>
#include <Options.hpp>
#include <Data.hpp>

namespace cc4s {
  class Cc4s {
  public:
    void run();

    static bool isDebugged();

    // static properties, accessible from everywhere
    static Ptr<MpiCommunicator> world;
    static Ptr<Options> options;

  protected:
    void fetch(const Ptr<MapNode> &arguments);
    void store(const Ptr<MapNode> &result, const Ptr<MapNode> &variables);

    void printBanner(const Ptr<MapNode> &job);
    void printStatistics(const Ptr<MapNode> &job);
    void listHosts(const Ptr<MapNode> &job);

    Ptr<MapNode> storage;
  };
}

#endif

/*! \mainpage 
 * \section intro Introduction
 * 
 * Coupled Cluster For Solids (CC4S) is is a parallel quantum
 * chemistry package built on the Cyclops Tensor Framework which
 * provides high-performance structured tensor operations. CC4S is
 * primarily focused on iterative CC methods such as CCD, CCSD, and
 * CCSD(T) for periodic systems in a plane-wave basis. The code is
 * currently interfaced to the Vienna Ab-inition Simulation Package.
 * 
 * The software is available on
 * https://hulldamage.org/cc4s.org/doku.php and maybe obtained
 * via the command 
 *
 * git clone git@gitlab.com:grueneis/cc4s.git
 *
 * For information or access to the repository please contact <mailto:andreas.grueneis@tuwien.ac.at>
 *
 * CC4S depends on the CTF and requires MPI to be built as the main parallel execution and communication mechanism.
 *
 * \section algorithms Algorithms
 *
 *The main functionality is based on different algorithms.
 */
