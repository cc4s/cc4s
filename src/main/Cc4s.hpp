/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CC4S_DEFINED
#define CC4S_DEFINED

#include <util/SharedPointer.hpp>
#include <util/MpiCommunicator.hpp>
#include <math/Integer.hpp>
#include <Options.hpp>
#include <Node.hpp>

namespace cc4s {
  class Cc4s {
  public:
    void run(const Ptr<MapNode> &report);

    static bool isDebugged();

    // static properties, accessible from everywhere
    static Ptr<MpiCommunicator> world;
    static Ptr<Options> options;

    static Natural<128> getFloatingPointOperations();
    static void addFloatingPointOperations(const Natural<128> ops);

  protected:
    void runStep(const Natural<> i, const Ptr<MapNode> &step);
    void fetchSymbols(const Ptr<MapNode> &arguments);
    void storeSymbols(const Ptr<MapNode> &result,const Ptr<MapNode> &variables);

    void printBanner(const Ptr<MapNode> &job);
    void printStatistics(const Ptr<MapNode> &job);
    void listHosts(const Ptr<MapNode> &job);

    Ptr<MapNode> storage;
  };

  class OperationsCounter {
  public:
    OperationsCounter(Natural<128> *_counter): counter(_counter) {
      *counter = -Cc4s::getFloatingPointOperations();
    }
    ~OperationsCounter() {
      *counter += Cc4s::getFloatingPointOperations();
    }
  protected:
    Natural<128> *counter;
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
