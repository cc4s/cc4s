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

#define CATCH_CONFIG_RUNNER
#include <test/Test.hpp>
#include <Cc4s.hpp>
#include <Parser.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;

World *Cc4s::world;
Options *Cc4s::options;

void printBanner() {
  OUT() << std::endl
        << "       UNIT TESTING__      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  OUT() << "version=" << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  OUT() << "build date=" << __DATE__ << " " << __TIME__ << std::endl;
  OUT() << "compiler=" << COMPILER_VERSION << std::endl;
  OUT() << "number of processes=" << Cc4s::world->np << std::endl << std::endl;
}


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Cc4s::world = new World(argumentCount, arguments);
  Cc4s::options = new Options(argumentCount, arguments);
  Log::setRank(Cc4s::world->rank);

  printBanner();

  const int result = Catch::Session().run(argumentCount, arguments);

  MPI_Finalize();

  return result;
}
